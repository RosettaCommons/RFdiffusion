"""SO(3) diffusion methods."""
import numpy as np
import os
from functools import cached_property
import torch
from scipy.spatial.transform import Rotation
import scipy.linalg


### First define geometric operations on the SO3 manifold

# hat map from vector space R^3 to Lie algebra so(3)
def hat(v):
    hat_v = torch.zeros([v.shape[0], 3, 3])
    hat_v[:, 0, 1], hat_v[:, 0, 2], hat_v[:, 1, 2] = -v[:, 2], v[:, 1], -v[:, 0]
    return hat_v + -hat_v.transpose(2, 1)

def hat_batch(v):
    """Batch hat map: [..., 3] -> [..., 3, 3] (cross-product / skew-symmetric matrix)."""
    bshape = v.shape[:-1]
    h = torch.zeros(*bshape, 3, 3, device=v.device, dtype=v.dtype)
    h[..., 0, 1] = -v[..., 2]
    h[..., 0, 2] =  v[..., 1]
    h[..., 1, 0] =  v[..., 2]
    h[..., 1, 2] = -v[..., 0]
    h[..., 2, 0] = -v[..., 1]
    h[..., 2, 1] =  v[..., 0]
    return h

def Log_torch(R):
    """On-device rotation matrix -> rotation vector. R: [..., 3, 3] -> [..., 3].
    Stays on the original device — no scipy or CPU transfers.
    Numerically stable across the full [0, pi] range:
      - Uses ||skew|| = 2*sin(theta) for theta when cos(theta) < 0 (avoids trace
        instability near pi where float32 R loses precision in the trace but skew
        elements remain accurate).
      - Falls back to R+I decomposition only when sin(theta) is sub-epsilon
        (skew elements are below float32 resolution, i.e. theta very close to pi).
    """
    orig_dtype = R.dtype
    R64 = R.to(torch.float64)
    trace = R64[..., 0, 0] + R64[..., 1, 1] + R64[..., 2, 2]
    cos_theta = torch.clamp((trace - 1.0) / 2.0, -1.0, 1.0)

    # Skew-symmetric part: (R - R^T)_vee = 2*sin(theta)*n_vec (sign-correct for all theta)
    skew = torch.stack([
        R64[..., 2, 1] - R64[..., 1, 2],
        R64[..., 0, 2] - R64[..., 2, 0],
        R64[..., 1, 0] - R64[..., 0, 1],
    ], dim=-1)
    skew_norm = torch.norm(skew, dim=-1)  # = 2*|sin(theta)|
    axis = skew / torch.clamp(skew_norm, min=1e-12)[..., None]

    # Theta: acos(cos_theta) for small angles; pi - asin(skew_norm/2) near pi.
    # The asin estimate uses the skew magnitude directly, avoiding trace instability.
    theta_trace = torch.acos(cos_theta)
    theta_asin  = skew.new_full(skew_norm.shape, np.pi) - torch.asin(torch.clamp(skew_norm / 2.0, 0.0, 1.0))
    theta = torch.where(cos_theta < 0.0, theta_asin, theta_trace)
    rotvec_std = theta[..., None] * axis

    # Near-pi fallback: when sin(theta) < float32 noise floor in R, skew -> 0 but
    # R + I = 2*outer(n,n) is still readable.  Use R+I decomposition for axis.
    diag = torch.stack([R64[..., 0, 0] + 1.0, R64[..., 1, 1] + 1.0, R64[..., 2, 2] + 1.0], dim=-1)
    ax_mags = torch.sqrt(torch.clamp(diag / 2.0, min=0.0))
    ref = torch.argmax(diag, dim=-1, keepdim=True)
    ref_row = torch.gather(R64, -2, ref.unsqueeze(-1).expand(*ref.shape[:-1], 1, 3)).squeeze(-2)
    signs = torch.sign(ref_row + 1e-30)
    ref_mask = torch.zeros_like(signs).scatter_(-1, ref, 1.0)
    signs = signs * (1.0 - ref_mask) + ref_mask
    ax_pi = ax_mags * signs
    ax_pi = ax_pi / torch.norm(ax_pi, dim=-1, keepdim=True).clamp(min=1e-15)
    rotvec_nearpi = ax_pi * theta[..., None]

    # Branch thresholds (based on cos_theta, stable for float32 R inputs):
    #   near_zero: cos ≈ 1  → identity rotation
    #   near_pi:   cos ≈ -1 → skew magnitude below float32 noise (~3.5e-4)
    near_zero = cos_theta[..., None] > (1.0 - 1e-10)
    near_pi   = cos_theta[..., None] < -(1.0 - 6.25e-8)

    rotvec = torch.where(near_zero, torch.zeros_like(rotvec_std),
             torch.where(near_pi,   rotvec_nearpi, rotvec_std))
    return rotvec.to(orig_dtype)

def Exp_torch(v):
    """On-device rotation vector -> rotation matrix. v: [..., 3] -> [..., 3, 3].
    Rodrigues formula. Stays on the original device/dtype."""
    theta = torch.norm(v, dim=-1)
    theta_safe = torch.clamp(theta, min=1e-7)
    axis = v / theta_safe[..., None]
    K = hat_batch(axis)
    I = torch.eye(3, device=v.device, dtype=v.dtype).expand(*v.shape[:-1], 3, 3)
    sin_t = torch.sin(theta)[..., None, None]
    cos_t = torch.cos(theta)[..., None, None]
    R = I + sin_t * K + (1.0 - cos_t) * (K @ K)
    return torch.where(theta[..., None, None] < 1e-7, I, R)

# Logarithmic map from SO(3) to R^3 (i.e. rotation vector) — legacy CPU version
def Log(R): return torch.tensor(Rotation.from_matrix(R.numpy()).as_rotvec())

# logarithmic map from SO(3) to so(3), this is the matrix logarithm
def log(R): return hat(Log(R))

# Exponential map from vector space of so(3) to SO(3) — legacy CPU version
def Exp(A): return torch.tensor(Rotation.from_rotvec(A.numpy()).as_matrix())

# Angle of rotation SO(3) to R^+
def Omega(R): return np.linalg.norm(log(R), axis=[-2, -1])/np.sqrt(2.)

L_default = 2000
def f_igso3(omega, t, L=L_default):
    """Truncated sum of IGSO(3) distribution.

    This function approximates the power series in equation 5 of
    "DENOISING DIFFUSION PROBABILISTIC MODELS ON SO(3) FOR ROTATIONAL
    ALIGNMENT"
    Leach et al. 2022

    This expression diverges from the expression in Leach in that here, sigma =
    sqrt(2) * eps, if eps_leach were the scale parameter of the IGSO(3).

    With this reparameterization, IGSO(3) agrees with the Brownian motion on
    SO(3) with t=sigma^2 when defined for the canonical inner product on SO3,
    <u, v>_SO3 = Trace(u v^T)/2

    Args:
        omega: i.e. the angle of rotation associated with rotation matrix
        t: variance parameter of IGSO(3), maps onto time in Brownian motion
        L: Truncation level
    """
    ls = torch.arange(L)[None]  # of shape [1, L]
    return ((2*ls + 1) * torch.exp(-ls*(ls+1)*t/2) *
             torch.sin(omega[:, None]*(ls+1/2)) / torch.sin(omega[:, None]/2)).sum(dim=-1)

def d_logf_d_omega(omega, t, L=L_default):
    omega = torch.tensor(omega, requires_grad=True)
    log_f = torch.log(f_igso3(omega, t, L))
    return torch.autograd.grad(log_f.sum(), omega)[0].numpy()

# IGSO3 density with respect to the volume form on SO(3)
def igso3_density(Rt, t, L=L_default):
    return f_igso3(torch.tensor(Omega(Rt)), t, L).numpy()

def igso3_density_angle(omega, t, L=L_default): 
    return f_igso3(torch.tensor(omega), t, L).numpy()*(1-np.cos(omega))/np.pi

# grad_R log IGSO3(R; I_3, t)
def igso3_score(R, t, L=L_default):
    omega = Omega(R)
    unit_vector = np.einsum('Nij,Njk->Nik', R, log(R))/omega[:, None, None]
    return unit_vector * d_logf_d_omega(omega, t, L)[:, None, None]

def calculate_igso3(*, num_sigma, num_omega, min_sigma, max_sigma):
    """calculate_igso3 pre-computes numerical approximations to the IGSO3 cdfs
    and score norms and expected squared score norms.

    Args:
        num_sigma: number of different sigmas for which to compute igso3
            quantities.
        num_omega: number of point in the discretization in the angle of
            rotation.
        min_sigma, max_sigma: the upper and lower ranges for the angle of
            rotation on which to consider the IGSO3 distribution.  This cannot
            be too low or it will create numerical instability.
    """
    # Discretize omegas for calculating CDFs. Skip omega=0.
    discrete_omega = np.linspace(0, np.pi, num_omega+1)[1:]

    # Exponential noise schedule.  This choice is closely tied to the
    # scalings used when simulating the reverse time SDE. For each step n,
    # discrete_sigma[n] = min_eps^(1-n/num_eps) * max_eps^(n/num_eps)
    discrete_sigma = 10 ** np.linspace(np.log10(min_sigma), np.log10(max_sigma), num_sigma + 1)[1:]

    # Compute the pdf and cdf values for the marginal distribution of the angle
    # of rotation (which is needed for sampling)
    pdf_vals = np.asarray(
        [igso3_density_angle(discrete_omega, sigma**2) for sigma in discrete_sigma])
    cdf_vals = np.asarray(
        [pdf.cumsum() / num_omega * np.pi for pdf in pdf_vals])

    # Compute the norms of the scores.  This are used to scale the rotation axis when
    # computing the score as a vector.
    score_norm = np.asarray(
        [d_logf_d_omega(discrete_omega, sigma**2) for sigma in discrete_sigma])

    # Compute the standard deviation of the score norm for each sigma
    exp_score_norms = np.sqrt(
        np.sum(
            score_norm**2 * pdf_vals, axis=1) / np.sum(
                pdf_vals, axis=1))
    return {
        'cdf': cdf_vals,
        'score_norm': score_norm,
        'exp_score_norms': exp_score_norms,
        'discrete_omega': discrete_omega,
        'discrete_sigma': discrete_sigma,
    }
