"""
Benchmark script for PR #454 performance improvements.

Measures:
  1. Flash Attention (F.scaled_dot_product_attention) vs manual einsum attention
  2. Torch-native SO(3) ops vs scipy (isolated + pipelined context)
  3. Log_torch accuracy including near-theta=pi

Run with: python scripts/benchmark_pr454.py
"""
import time
import math
import torch
import torch.nn.functional as F
import numpy as np
from scipy.spatial.transform import Rotation

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DTYPE  = torch.float32
N_WARMUP = 50
N_TRIALS = 500

def timer(fn, warmup=N_WARMUP, trials=N_TRIALS):
    for _ in range(warmup):
        fn()
    if DEVICE == "cuda":
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(trials):
        fn()
    if DEVICE == "cuda":
        torch.cuda.synchronize()
    return (time.perf_counter() - t0) / trials * 1000  # ms

# ─── SO(3) helpers ────────────────────────────────────────────────────────────

def hat_batch(v):
    bshape = v.shape[:-1]
    h = torch.zeros(*bshape, 3, 3, device=v.device, dtype=v.dtype)
    h[..., 0, 1] = -v[..., 2];  h[..., 0, 2] =  v[..., 1]
    h[..., 1, 0] =  v[..., 2];  h[..., 1, 2] = -v[..., 0]
    h[..., 2, 0] = -v[..., 1];  h[..., 2, 1] =  v[..., 0]
    return h

def Log_torch(R):
    """Three-branch Log: near-identity / standard / near-pi."""
    trace = R[..., 0, 0] + R[..., 1, 1] + R[..., 2, 2]
    theta = torch.acos(torch.clamp((trace - 1.0) / 2.0, -1.0, 1.0))
    skew = torch.stack([
        R[..., 2, 1] - R[..., 1, 2],
        R[..., 0, 2] - R[..., 2, 0],
        R[..., 1, 0] - R[..., 0, 1],
    ], dim=-1)
    rotvec_std = (theta / (2.0 * torch.clamp(torch.sin(theta), min=1e-7)))[..., None] * skew
    # Near-pi branch: R + I = 2 * outer(n, n)
    diag = torch.stack([R[..., 0, 0]+1.0, R[..., 1, 1]+1.0, R[..., 2, 2]+1.0], dim=-1)
    ax_mags = torch.sqrt(torch.clamp(diag / 2.0, min=0.0))
    ref = torch.argmax(diag, dim=-1, keepdim=True)
    ref_row = torch.gather(R, -2, ref.unsqueeze(-1).expand(*ref.shape[:-1], 1, 3)).squeeze(-2)
    signs = torch.sign(ref_row + 1e-10)
    ref_mask = torch.zeros_like(signs).scatter_(-1, ref, 1.0)
    signs = signs * (1.0 - ref_mask) + ref_mask
    ax_pi = ax_mags * signs
    ax_pi = ax_pi / torch.norm(ax_pi, dim=-1, keepdim=True).clamp(min=1e-7)
    rotvec_pi = ax_pi * theta[..., None]
    near_zero = theta[..., None] < 1e-6
    near_pi   = theta[..., None] > (math.pi - 1e-3)
    return torch.where(near_zero, torch.zeros_like(rotvec_std),
           torch.where(near_pi, rotvec_pi, rotvec_std))

def Exp_torch(v):
    theta = torch.norm(v, dim=-1)
    axis = v / torch.clamp(theta, min=1e-7)[..., None]
    K = hat_batch(axis)
    I = torch.eye(3, device=v.device, dtype=v.dtype).expand(*v.shape[:-1], 3, 3)
    R = I + torch.sin(theta)[..., None, None] * K + (1.0 - torch.cos(theta))[..., None, None] * (K @ K)
    return torch.where(theta[..., None, None] < 1e-7, I, R)

def scipy_Log(R_gpu):
    return torch.from_numpy(Rotation.from_matrix(R_gpu.cpu().numpy()).as_rotvec()).to(R_gpu.device)

def scipy_Exp(v_gpu):
    return torch.from_numpy(Rotation.from_rotvec(v_gpu.cpu().numpy()).as_matrix()).float().to(v_gpu.device)

# ─── Benchmark 1: Attention ────────────────────────────────────────────────────

def bench_attention():
    print("\n" + "="*66)
    print("Flash Attention  (F.scaled_dot_product_attention vs einsum)")
    print("="*66)

    configs = [
        (1, 4,  64, 32, "L=64  (short)"),
        (1, 4, 200, 32, "L=200 (typical)"),
        (1, 4, 500, 32, "L=500 (long)"),
        (4, 4, 200, 32, "L=200 batch=4"),
    ]

    for label, (B, h, L, d) in [(c[-1], c[:-1]) for c in configs]:
        q = torch.randn(B, h, L, d, device=DEVICE, dtype=DTYPE)
        k = torch.randn(B, h, L, d, device=DEVICE, dtype=DTYPE)
        v = torch.randn(B, h, L, d, device=DEVICE, dtype=DTYPE)
        bias = torch.randn(B, h, L, L, device=DEVICE, dtype=DTYPE)
        scale = d ** -0.5

        def old_plain():
            a = torch.einsum('bhid,bhjd->bhij', q * scale, k)
            a = F.softmax(a, dim=-1)
            return torch.einsum('bhij,bhjd->bhid', a, v)

        def new_plain():
            return F.scaled_dot_product_attention(q, k, v)

        def old_bias():
            a = torch.einsum('bhid,bhjd->bhij', q * scale, k) + bias
            a = F.softmax(a, dim=-1)
            return torch.einsum('bhij,bhjd->bhid', a, v)

        def new_bias():
            return F.scaled_dot_product_attention(q, k, v, attn_mask=bias)

        t_op = timer(old_plain); t_np = timer(new_plain)
        t_ob = timer(old_bias);  t_nb = timer(new_bias)
        print(f"  {label:<18}  plain {t_op:.3f}ms -> {t_np:.3f}ms  ({t_op/t_np:.1f}x)   "
              f"biased {t_ob:.3f}ms -> {t_nb:.3f}ms  ({t_ob/t_nb:.1f}x)")

# ─── Benchmark 2: SO(3) ───────────────────────────────────────────────────────

def bench_so3():
    print("\n" + "="*66)
    print("SO(3) Log+Exp: scipy CPU-roundtrip vs torch-native (on GPU)")
    print("Note: isolated timings; GPU pipeline benefit not captured here")
    print("="*66)

    for N in [50, 200, 500]:
        vn = np.random.randn(N, 3).astype(np.float32)
        vn = vn / np.linalg.norm(vn, axis=1, keepdims=True) * np.random.uniform(0.01, 3.0, (N, 1))
        Rn = Rotation.from_rotvec(vn.astype(np.float64)).as_matrix().astype(np.float32)
        R_gpu = torch.from_numpy(Rn).to(DEVICE)
        v_gpu = torch.from_numpy(vn).to(DEVICE)

        t_sc = timer(lambda: (scipy_Log(R_gpu), scipy_Exp(v_gpu)), warmup=10, trials=200)
        t_th = timer(lambda: (Log_torch(R_gpu), Exp_torch(v_gpu)), warmup=10, trials=200)
        note = "(*)" if t_th > t_sc else ""
        print(f"  N={N:<5}  scipy {t_sc:.3f}ms  torch {t_th:.3f}ms  {t_sc/t_th:.1f}x{note}")

    print("  (*) For small N, GPU kernel launch overhead exceeds scipy cost in isolation.")
    print("      The torch path eliminates the CPU sync point that stalls the GPU pipeline")
    print("      between denoising steps, and keeps the full trajectory on-device.")

# ─── Benchmark 3: GPU pipeline stall ──────────────────────────────────────────

def bench_pipeline_stall():
    print("\n" + "="*66)
    print("GPU pipeline stall from .cpu() transfer (100-step denoising loop)")
    print("="*66)
    N = 200
    Rn = Rotation.from_rotvec(np.random.randn(N, 3)).as_matrix().astype(np.float32)
    R_gpu = torch.from_numpy(Rn).to(DEVICE)
    v_gpu = torch.randn(N, 3, device=DEVICE)
    steps = 100

    if DEVICE == "cuda":
        torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(steps):
        dummy = torch.randn(N, 3, device=DEVICE) * 0.01 + v_gpu
        scipy_Log(R_gpu + 0)   # forces CPU sync each step
        scipy_Exp(dummy)
    torch.cuda.synchronize()
    t_scipy = (time.perf_counter() - t0) * 1000

    torch.cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(steps):
        dummy = torch.randn(N, 3, device=DEVICE) * 0.01 + v_gpu
        Log_torch(R_gpu + 0)   # stays on GPU
        Exp_torch(dummy)
    torch.cuda.synchronize()
    t_torch = (time.perf_counter() - t0) * 1000

    print(f"  {steps}-step loop (N={N} residues):")
    print(f"    scipy  (CPU sync each step): {t_scipy:.1f}ms")
    print(f"    torch  (no sync):            {t_torch:.1f}ms")
    print(f"    speedup: {t_scipy/t_torch:.1f}x")

# ─── Accuracy ─────────────────────────────────────────────────────────────────

def check_accuracy():
    """Round-trip test: R -> Log_torch -> Exp_torch -> R.
    A perfect Log would recover R exactly; errors here are the combined
    Log+Exp error against the ground-truth rotation.
    """
    import sys; sys.path.insert(0, '.')
    from rfdiffusion.igso3 import Log_torch as Log_igso3, Exp_torch as Exp_igso3

    print("\n" + "="*66)
    print("Accuracy: R -> Log_torch -> Exp_torch -> R  (round-trip |dR|)")
    print("="*66)

    np.random.seed(0)
    N = 100000
    mags = np.random.uniform(0.01, math.pi, N)
    axes = np.random.randn(N, 3); axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    vn = (axes * mags[:, None]).astype(np.float32)
    Rn = Rotation.from_rotvec(vn.astype(np.float64)).as_matrix().astype(np.float32)
    R_gpu = torch.from_numpy(Rn).to(DEVICE)
    R_rec = Exp_igso3(Log_igso3(R_gpu))
    errs = (R_rec - R_gpu).abs().amax(dim=(-1, -2)).cpu().numpy()
    print(f"  Log_torch + Exp_torch, full [0,pi]: max|dR|={errs.max():.2e}  mean={errs.mean():.2e}")
    for lo, hi, label in [(0,2,'0..2'), (2,3,'2..3'), (3,math.pi-0.01,'3..pi-0.01'), (math.pi-0.01,math.pi,'pi-0.01..pi')]:
        m = (mags>=lo)&(mags<=hi)
        if m.sum(): print(f"    [{label:<12}] N={m.sum():6d}  max={errs[m].max():.2e}  mean={errs[m].mean():.2e}")

    log_s = torch.from_numpy(Rotation.from_matrix(Rn.astype(np.float64)).as_rotvec().astype(np.float32)).to(DEVICE)
    R_sc  = Exp_igso3(log_s)
    errs_s = (R_sc - R_gpu).abs().amax(dim=(-1,-2)).cpu().numpy()
    print(f"  scipy Log + Exp_torch baseline:    max|dR|={errs_s.max():.2e}  mean={errs_s.mean():.2e}")

# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\nDevice : {DEVICE}")
    if DEVICE == "cuda":
        print(f"GPU    : {torch.cuda.get_device_name(0)}")
    print(f"PyTorch: {torch.__version__}")

    check_accuracy()
    bench_attention()
    bench_so3()
    if DEVICE == "cuda":
        bench_pipeline_stall()

    print("\nDone.")
