import numpy as np
import os
from omegaconf import DictConfig
import torch
import torch.nn.functional as nn
from rfdiffusion.diffusion import get_beta_schedule
from scipy.spatial.transform import Rotation as scipy_R
from rfdiffusion.util import rigid_from_3_points
from rfdiffusion.util_module import ComputeAllAtomCoords
from rfdiffusion import util
import random
import logging
from rfdiffusion.inference import model_runners
import glob

###########################################################
#### Functions which can be called outside of Denoiser ####
###########################################################


def get_next_frames(xt, px0, t, diffuser, so3_type, diffusion_mask, noise_scale=1.0):
    """
    get_next_frames gets updated frames using IGSO(3) + score_based reverse diffusion.


    based on self.so3_type use score based update.

    Generate frames at t-1
    Rather than generating random rotations (as occurs during forward process), calculate rotation between xt and px0

    Args:
        xt: noised coordinates of shape [L, 14, 3]
        px0: prediction of coordinates at t=0, of shape [L, 14, 3]
        t: integer time step
        diffuser: Diffuser object for reverse igSO3 sampling
        so3_type: The type of SO3 noising being used ('igso3')
        diffusion_mask: of shape [L] of type bool, True means not to be
            updated (e.g. mask is true for motif residues)
        noise_scale: scale factor for the noise added (IGSO3 only)

    Returns:
        backbone coordinates for step x_t-1 of shape [L, 3, 3]
    """
    N_0 = px0[None, :, 0, :]
    Ca_0 = px0[None, :, 1, :]
    C_0 = px0[None, :, 2, :]

    R_0, Ca_0 = rigid_from_3_points(N_0, Ca_0, C_0)

    N_t = xt[None, :, 0, :]
    Ca_t = xt[None, :, 1, :]
    C_t = xt[None, :, 2, :]

    R_t, Ca_t = rigid_from_3_points(N_t, Ca_t, C_t)

    # this must be to normalize them or something
    R_0 = scipy_R.from_matrix(R_0.squeeze().numpy()).as_matrix()
    R_t = scipy_R.from_matrix(R_t.squeeze().numpy()).as_matrix()

    L = R_t.shape[0]
    all_rot_transitions = np.broadcast_to(np.identity(3), (L, 3, 3)).copy()
    # Sample next frame for each residue
    if so3_type == "igso3":
        # don't do calculations on masked positions since they end up as identity matrix
        all_rot_transitions[
            ~diffusion_mask
        ] = diffuser.so3_diffuser.reverse_sample_vectorized(
            R_t[~diffusion_mask],
            R_0[~diffusion_mask],
            t,
            noise_level=noise_scale,
            mask=None,
            return_perturb=True,
        )
    else:
        assert False, "so3 diffusion type %s not implemented" % so3_type

    all_rot_transitions = all_rot_transitions[:, None, :, :]

    # Apply the interpolated rotation matrices to the coordinates
    next_crds = (
        np.einsum(
            "lrij,laj->lrai",
            all_rot_transitions,
            xt[:, :3, :] - Ca_t.squeeze()[:, None, ...].numpy(),
        )
        + Ca_t.squeeze()[:, None, None, ...].numpy()
    )

    # (L,3,3) set of backbone coordinates with slight rotation
    return next_crds.squeeze(1)


def get_mu_xt_x0(xt, px0, t, beta_schedule, alphabar_schedule, eps=1e-6):
    """
    Given xt, predicted x0 and the timestep t, give mu of x(t-1)
    Assumes t is 0 indexed
    """
    # sigma is predefined from beta. Often referred to as beta tilde t
    t_idx = t - 1
    sigma = (
        (1 - alphabar_schedule[t_idx - 1]) / (1 - alphabar_schedule[t_idx])
    ) * beta_schedule[t_idx]

    xt_ca = xt[:, 1, :]
    px0_ca = px0[:, 1, :]

    a = (
        (torch.sqrt(alphabar_schedule[t_idx - 1] + eps) * beta_schedule[t_idx])
        / (1 - alphabar_schedule[t_idx])
    ) * px0_ca
    b = (
        (
            torch.sqrt(1 - beta_schedule[t_idx] + eps)
            * (1 - alphabar_schedule[t_idx - 1])
        )
        / (1 - alphabar_schedule[t_idx])
    ) * xt_ca

    mu = a + b

    return mu, sigma


def get_next_ca(
    xt,
    px0,
    t,
    diffusion_mask,
    crd_scale,
    beta_schedule,
    alphabar_schedule,
    noise_scale=1.0,
):
    """
    Given full atom x0 prediction (xyz coordinates), diffuse to x(t-1)

    Parameters:

        xt (L, 14/27, 3) set of coordinates

        px0 (L, 14/27, 3) set of coordinates

        t: time step. Note this is zero-index current time step, so are generating t-1

        logits_aa (L x 20 ) amino acid probabilities at each position

        seq_schedule (L): Tensor of bools, True is unmasked, False is masked. For this specific t

        diffusion_mask (torch.tensor, required): Tensor of bools, True means NOT diffused at this residue, False means diffused

        noise_scale: scale factor for the noise being added

    """
    get_allatom = ComputeAllAtomCoords().to(device=xt.device)
    L = len(xt)

    # bring to origin after global alignment (when don't have a motif) or replace input motif and bring to origin, and then scale
    px0 = px0 * crd_scale
    xt = xt * crd_scale

    # get mu(xt, x0)
    mu, sigma = get_mu_xt_x0(
        xt, px0, t, beta_schedule=beta_schedule, alphabar_schedule=alphabar_schedule
    )

    sampled_crds = torch.normal(mu, torch.sqrt(sigma * noise_scale))
    delta = sampled_crds - xt[:, 1, :]  # check sign of this is correct

    if not diffusion_mask is None:
        # Don't move motif
        delta[diffusion_mask, ...] = 0

    out_crds = xt + delta[:, None, :]

    return out_crds / crd_scale, delta / crd_scale


def get_noise_schedule(T, noiseT, noise1, schedule_type):
    """
    Function to create a schedule that varies the scale of noise given to the model over time

    Parameters:

        T: The total number of timesteps in the denoising trajectory

        noiseT: The inital (t=T) noise scale

        noise1: The final (t=1) noise scale

        schedule_type: The type of function to use to interpolate between noiseT and noise1

    Returns:

        noise_schedule: A function which maps timestep to noise scale

    """

    noise_schedules = {
        "constant": lambda t: noiseT,
        "linear": lambda t: ((t - 1) / (T - 1)) * (noiseT - noise1) + noise1,
    }

    assert (
        schedule_type in noise_schedules
    ), f"noise_schedule must be one of {noise_schedules.keys()}. Received noise_schedule={schedule_type}. Exiting."

    return noise_schedules[schedule_type]


class Denoise:
    """
    Class for getting x(t-1) from predicted x0 and x(t)
    Strategy:
        Ca coordinates: Rediffuse to x(t-1) from predicted x0
        Frames: Approximate update from rotation score
        Torsions: 1/t of the way to the x0 prediction

    """

    def __init__(
        self,
        T,
        L,
        diffuser,
        b_0=0.001,
        b_T=0.1,
        min_b=1.0,
        max_b=12.5,
        min_sigma=0.05,
        max_sigma=1.5,
        noise_level=0.5,
        schedule_type="linear",
        so3_schedule_type="linear",
        schedule_kwargs={},
        so3_type="igso3",
        noise_scale_ca=1.0,
        final_noise_scale_ca=1,
        ca_noise_schedule_type="constant",
        noise_scale_frame=0.5,
        final_noise_scale_frame=0.5,
        frame_noise_schedule_type="constant",
        crd_scale=1 / 15,
        potential_manager=None,
        partial_T=None,
    ):
        """

        Parameters:
            noise_level: scaling on the noise added (set to 0 to use no noise,
                to 1 to have full noise)

        """
        self.T = T
        self.L = L
        self.diffuser = diffuser
        self.b_0 = b_0
        self.b_T = b_T
        self.noise_level = noise_level
        self.schedule_type = schedule_type
        self.so3_type = so3_type
        self.crd_scale = crd_scale
        self.noise_scale_ca = noise_scale_ca
        self.final_noise_scale_ca = final_noise_scale_ca
        self.ca_noise_schedule_type = ca_noise_schedule_type
        self.noise_scale_frame = noise_scale_frame
        self.final_noise_scale_frame = final_noise_scale_frame
        self.frame_noise_schedule_type = frame_noise_schedule_type
        self.potential_manager = potential_manager
        self._log = logging.getLogger(__name__)

        self.schedule, self.alpha_schedule, self.alphabar_schedule = get_beta_schedule(
            self.T, self.b_0, self.b_T, self.schedule_type, inference=True
        )

        self.noise_schedule_ca = get_noise_schedule(
            self.T,
            self.noise_scale_ca,
            self.final_noise_scale_ca,
            self.ca_noise_schedule_type,
        )
        self.noise_schedule_frame = get_noise_schedule(
            self.T,
            self.noise_scale_frame,
            self.final_noise_scale_frame,
            self.frame_noise_schedule_type,
        )

    @property
    def idx2steps(self):
        return self.decode_scheduler.idx2steps.numpy()

    def align_to_xt_motif(self, px0, xT, diffusion_mask, eps=1e-6):
        """
        Need to align px0 to motif in xT. This is to permit the swapping of residue positions in the px0 motif for the true coordinates.
        First, get rotation matrix from px0 to xT for the motif residues.
        Second, rotate px0 (whole structure) by that rotation matrix
        Third, centre at origin
        """

        def rmsd(V, W, eps=0):
            # First sum down atoms, then sum down xyz
            N = V.shape[-2]
            return np.sqrt(np.sum((V - W) * (V - W), axis=(-2, -1)) / N + eps)

        assert (
            xT.shape[1] == px0.shape[1]
        ), f"xT has shape {xT.shape} and px0 has shape {px0.shape}"

        L, n_atom, _ = xT.shape  # A is number of atoms
        atom_mask = ~torch.isnan(px0)
        # convert to numpy arrays
        px0 = px0.cpu().detach().numpy()
        xT = xT.cpu().detach().numpy()
        diffusion_mask = diffusion_mask.cpu().detach().numpy()

        # 1 centre motifs at origin and get rotation matrix
        px0_motif = px0[diffusion_mask, :3].reshape(-1, 3)
        xT_motif = xT[diffusion_mask, :3].reshape(-1, 3)
        px0_motif_mean = np.copy(px0_motif.mean(0))  # need later
        xT_motif_mean = np.copy(xT_motif.mean(0))

        # center at origin
        px0_motif = px0_motif - px0_motif_mean
        xT_motif = xT_motif - xT_motif_mean

        # A = px0_motif
        # B = xT_motif
        A = xT_motif
        B = px0_motif

        C = np.matmul(A.T, B)

        # compute optimal rotation matrix using SVD
        U, S, Vt = np.linalg.svd(C)

        # ensure right handed coordinate system
        d = np.eye(3)
        d[-1, -1] = np.sign(np.linalg.det(Vt.T @ U.T))

        # construct rotation matrix
        R = Vt.T @ d @ U.T

        # get rotated coords
        rB = B @ R

        # calculate rmsd
        rms = rmsd(A, rB)
        self._log.info(f"Sampled motif RMSD: {rms:.2f}")

        # 2 rotate whole px0 by rotation matrix
        atom_mask = atom_mask.cpu()
        px0[~atom_mask] = 0  # convert nans to 0
        px0 = px0.reshape(-1, 3) - px0_motif_mean
        px0_ = px0 @ R

        # 3 put in same global position as xT
        px0_ = px0_ + xT_motif_mean
        px0_ = px0_.reshape([L, n_atom, 3])
        px0_[~atom_mask] = float("nan")
        return torch.Tensor(px0_)

    def get_potential_gradients(self, xyz, diffusion_mask):
        """
        This could be moved into potential manager if desired - NRB

        Function to take a structure (x) and get per-atom gradients used to guide diffusion update

        Inputs:

            xyz (torch.tensor, required): [L,27,3] Coordinates at which the gradient will be computed

        Outputs:

            Ca_grads (torch.tensor): [L,3] The gradient at each Ca atom
        """

        if self.potential_manager == None or self.potential_manager.is_empty():
            return torch.zeros(xyz.shape[0], 3)

        use_Cb = False

        # seq.requires_grad = True
        xyz.requires_grad = True

        if not xyz.grad is None:
            xyz.grad.zero_()

        current_potential = self.potential_manager.compute_all_potentials(xyz)
        current_potential.backward()

        # Since we are not moving frames, Cb grads are same as Ca grads
        # Need access to calculated Cb coordinates to be able to get Cb grads though
        Ca_grads = xyz.grad[:, 1, :]

        if not diffusion_mask == None:
            Ca_grads[diffusion_mask, :] = 0

        # check for NaN's
        if torch.isnan(Ca_grads).any():
            print("WARNING: NaN in potential gradients, replacing with zero grad.")
            Ca_grads[:] = 0

        return Ca_grads

    def get_next_pose(
        self,
        xt,
        px0,
        t,
        diffusion_mask,
        fix_motif=True,
        align_motif=True,
        include_motif_sidechains=True,
    ):
        """
        Wrapper function to take px0, xt and t, and to produce xt-1
        First, aligns px0 to xt
        Then gets coordinates, frames and torsion angles

        Parameters:

            xt (torch.tensor, required): Current coordinates at timestep t

            px0 (torch.tensor, required): Prediction of x0

            t (int, required): timestep t

            diffusion_mask (torch.tensor, required): Mask for structure diffusion

            fix_motif (bool): Fix the motif structure

            align_motif (bool): Align the model's prediction of the motif to the input motif

            include_motif_sidechains (bool): Provide sidechains of the fixed motif to the model
        """

        get_allatom = ComputeAllAtomCoords().to(device=xt.device)
        L, n_atom = xt.shape[:2]
        assert (xt.shape[1] == 14) or (xt.shape[1] == 27)
        assert (px0.shape[1] == 14) or (px0.shape[1] == 27)

        ###############################
        ### Align pX0 onto Xt motif ###
        ###############################

        if align_motif and diffusion_mask.any():
            px0 = self.align_to_xt_motif(px0, xt, diffusion_mask)
        # xT_motif_aligned = self.align_to_xt_motif(px0, xt, diffusion_mask)

        px0 = px0.to(xt.device)
        # Now done with diffusion mask. if fix motif is False, just set diffusion mask to be all True, and all coordinates can diffuse
        if not fix_motif:
            diffusion_mask[:] = False

        # get the next set of CA coordinates
        noise_scale_ca = self.noise_schedule_ca(t)
        _, ca_deltas = get_next_ca(
            xt,
            px0,
            t,
            diffusion_mask,
            crd_scale=self.crd_scale,
            beta_schedule=self.schedule,
            alphabar_schedule=self.alphabar_schedule,
            noise_scale=noise_scale_ca,
        )

        # get the next set of backbone frames (coordinates)
        noise_scale_frame = self.noise_schedule_frame(t)
        frames_next = get_next_frames(
            xt,
            px0,
            t,
            diffuser=self.diffuser,
            so3_type=self.so3_type,
            diffusion_mask=diffusion_mask,
            noise_scale=noise_scale_frame,
        )

        # Apply gradient step from guiding potentials
        # This can be moved to below where the full atom representation is calculated to allow for potentials involving sidechains

        grad_ca = self.get_potential_gradients(
            xt.clone(), diffusion_mask=diffusion_mask
        )

        ca_deltas += self.potential_manager.get_guide_scale(t) * grad_ca

        # add the delta to the new frames
        frames_next = torch.from_numpy(frames_next) + ca_deltas[:, None, :]  # translate

        fullatom_next = torch.full_like(xt, float("nan")).unsqueeze(0)
        fullatom_next[:, :, :3] = frames_next[None]
        # This is never used so just make it a fudged tensor - NRB
        torsions_next = torch.zeros(1, 1)

        if include_motif_sidechains:
            fullatom_next[:, diffusion_mask, :14] = xt[None, diffusion_mask]

        return fullatom_next.squeeze()[:, :14, :], px0


def sampler_selector(conf: DictConfig):
    if conf.scaffoldguided.scaffoldguided:
        sampler = model_runners.ScaffoldedSampler(conf)
    else:
        if conf.inference.model_runner == "default":
            sampler = model_runners.Sampler(conf)
        elif conf.inference.model_runner == "SelfConditioning":
            sampler = model_runners.SelfConditioning(conf)
        elif conf.inference.model_runner == "ScaffoldedSampler":
            sampler = model_runners.ScaffoldedSampler(conf)
        else:
            raise ValueError(f"Unrecognized sampler {conf.model_runner}")
    return sampler


def parse_pdb(filename, **kwargs):
    """extract xyz coords for all heavy atoms"""
    with open(filename,"r") as f:
        lines=f.readlines()
    return parse_pdb_lines(lines, **kwargs)


def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res, pdb_idx = [],[]
    for l in lines:
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            res.append((l[22:26], l[17:20]))
            # chain letter, res num
            pdb_idx.append((l[21:22].strip(), int(l[22:26].strip())))
    seq = [util.aa2num[r[1]] if r[1] in util.aa2num.keys() else 20 for r in res]
    pdb_idx = [
        (l[21:22].strip(), int(l[22:26].strip()))
        for l in lines
        if l[:4] == "ATOM" and l[12:16].strip() == "CA"
    ]  # chain letter, res num

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 14, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = (
            l[21:22],
            int(l[22:26]),
            " " + l[12:16].strip().ljust(3),
            l[17:20],
        )
        if (chain,resNo) in pdb_idx:
            idx = pdb_idx.index((chain, resNo))
            # for i_atm, tgtatm in enumerate(util.aa2long[util.aa2num[aa]]):
            for i_atm, tgtatm in enumerate(
                util.aa2long[util.aa2num[aa]][:14]
                ):
                if (
                    tgtatm is not None and tgtatm.strip() == atom.strip()
                    ):  # ignore whitespace
                    xyz[idx, i_atm, :] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    break

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[..., 0]))
    xyz[np.isnan(xyz[..., 0])] = 0.0

    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i, idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)

    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]

    seq = np.array(seq)[i_unique]

    out = {
        "xyz": xyz,  # cartesian coordinates, [Lx14]
        "mask": mask,  # mask showing which atoms are present in the PDB file, [Lx14]
        "idx": np.array(
            [i[1] for i in pdb_idx]
        ),  # residue numbers in the PDB file, [L]
        "seq": np.array(seq),  # amino acid sequence, [L]
        "pdb_idx": pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
    }

    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6] == "HETATM" and not (ignore_het_h and l[77] == "H"):
                info_het.append(
                    dict(
                        idx=int(l[7:11]),
                        atom_id=l[12:16],
                        atom_type=l[77],
                        name=l[16:20],
                    )
                )
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out["xyz_het"] = np.array(xyz_het)
        out["info_het"] = info_het

    return out


def process_target(pdb_path, parse_hetatom=False, center=True):
    # Read target pdb and extract features.
    target_struct = parse_pdb(pdb_path, parse_hetatom=parse_hetatom)

    # Zero-center positions
    ca_center = target_struct["xyz"][:, :1, :].mean(axis=0, keepdims=True)
    if not center:
        ca_center = 0
    xyz = torch.from_numpy(target_struct["xyz"] - ca_center)
    seq_orig = torch.from_numpy(target_struct["seq"])
    atom_mask = torch.from_numpy(target_struct["mask"])
    seq_len = len(xyz)

    # Make 27 atom representation
    xyz_27 = torch.full((seq_len, 27, 3), np.nan).float()
    xyz_27[:, :14, :] = xyz[:, :14, :]

    mask_27 = torch.full((seq_len, 27), False)
    mask_27[:, :14] = atom_mask
    out = {
        "xyz_27": xyz_27,
        "mask_27": mask_27,
        "seq": seq_orig,
        "pdb_idx": target_struct["pdb_idx"],
    }
    if parse_hetatom:
        out["xyz_het"] = target_struct["xyz_het"]
        out["info_het"] = target_struct["info_het"]
    return out


def get_idx0_hotspots(mappings, ppi_conf, binderlen):
    """
    Take pdb-indexed hotspot resudes and the length of the binder, and makes the 0-indexed tensor of hotspots
    """

    hotspot_idx = None
    if binderlen > 0:
        if ppi_conf.hotspot_res is not None:
            assert all(
                [i[0].isalpha() for i in ppi_conf.hotspot_res]
            ), "Hotspot residues need to be provided in pdb-indexed form. E.g. A100,A103"
            hotspots = [(i[0], int(i[1:])) for i in ppi_conf.hotspot_res]
            hotspot_idx = []
            for i, res in enumerate(mappings["receptor_con_ref_pdb_idx"]):
                if res in hotspots:
                    hotspot_idx.append(mappings["receptor_con_hal_idx0"][i])
    return hotspot_idx


class BlockAdjacency:
    """
    Class for handling PPI design inference with ss/block_adj inputs.
    Basic idea is to provide a list of scaffolds, and to output ss and adjacency
    matrices based off of these, while sampling additional lengths.
    Inputs:
        - scaffold_list: list of scaffolds (e.g. ['2kl8','1cif']). Can also be a .txt file.
        - scaffold dir: directory where scaffold ss and adj are precalculated
        - sampled_insertion: how many additional residues do you want to add to each loop segment? Randomly sampled 0-this number (or within given range)
        - sampled_N: randomly sample up to this number of additional residues at N-term
        - sampled_C: randomly sample up to this number of additional residues at C-term
        - ss_mask: how many residues do you want to mask at either end of a ss (H or E) block. Fixed value
        - num_designs: how many designs are you wanting to generate? Currently only used for bookkeeping
        - systematic: do you want to systematically work through the list of scaffolds, or randomly sample (default)
        - num_designs_per_input: Not really implemented yet. Maybe not necessary
    Outputs:
        - L: new length of chain to be diffused
        - ss: all loops and insertions, and ends of ss blocks (up to ss_mask) set to mask token (3). Onehot encoded. (L,4)
        - adj: block adjacency with equivalent masking as ss (L,L)
    """

    def __init__(self, conf, num_designs):
        """
        Parameters:
          inputs:
             conf.scaffold_list as conf
             conf.inference.num_designs for sanity checking
        """
       
        self.conf=conf 
        # either list or path to .txt file with list of scaffolds
        if self.conf.scaffoldguided.scaffold_list is not None:
            if type(self.conf.scaffoldguided.scaffold_list) == list:
                self.scaffold_list = scaffold_list
            elif self.conf.scaffoldguided.scaffold_list[-4:] == ".txt":
                # txt file with list of ids
                list_from_file = []
                with open(self.conf.scaffoldguided.scaffold_list, "r") as f:
                    for line in f:
                        list_from_file.append(line.strip())
                self.scaffold_list = list_from_file
            else:
                raise NotImplementedError
        else:
            self.scaffold_list = [
                os.path.split(i)[1][:-6]
                for i in glob.glob(f"{self.conf.scaffoldguided.scaffold_dir}/*_ss.pt")
            ]
            self.scaffold_list.sort()

        # path to directory with scaffolds, ss files and block_adjacency files
        self.scaffold_dir = self.conf.scaffoldguided.scaffold_dir

        # maximum sampled insertion in each loop segment
        if "-" in str(self.conf.scaffoldguided.sampled_insertion):
            self.sampled_insertion = [
                int(str(self.conf.scaffoldguided.sampled_insertion).split("-")[0]),
                int(str(self.conf.scaffoldguided.sampled_insertion).split("-")[1]),
            ]
        else:
            self.sampled_insertion = [0, int(self.conf.scaffoldguided.sampled_insertion)]

        # maximum sampled insertion at N- and C-terminus
        if "-" in str(self.conf.scaffoldguided.sampled_N):
            self.sampled_N = [
                int(str(self.conf.scaffoldguided.sampled_N).split("-")[0]),
                int(str(self.conf.scaffoldguided.sampled_N).split("-")[1]),
            ]
        else:
            self.sampled_N = [0, int(self.conf.scaffoldguided.sampled_N)]
        if "-" in str(self.conf.scaffoldguided.sampled_C):
            self.sampled_C = [
                int(str(self.conf.scaffoldguided.sampled_C).split("-")[0]),
                int(str(self.conf.scaffoldguided.sampled_C).split("-")[1]),
            ]
        else:
            self.sampled_C = [0, int(self.conf.scaffoldguided.sampled_C)]

        # number of residues to mask ss identity of in H/E regions (from junction)
        # e.g. if ss_mask = 2, L,L,L,H,H,H,H,H,H,H,L,L,E,E,E,E,E,E,L,L,L,L,L,L would become\
        # M,M,M,M,M,H,H,H,M,M,M,M,M,M,E,E,M,M,M,M,M,M,M,M where M is mask
        self.ss_mask = self.conf.scaffoldguided.ss_mask

        # whether or not to work systematically through the list
        self.systematic = self.conf.scaffoldguided.systematic

        self.num_designs = num_designs

        if len(self.scaffold_list) > self.num_designs:
            print(
                "WARNING: Scaffold set is bigger than num_designs, so not every scaffold type will be sampled"
            )

        # for tracking number of designs
        self.num_completed = 0
        if self.systematic:
            self.item_n = 0

        # whether to mask loops or not
        if not self.conf.scaffoldguided.mask_loops:
            assert self.conf.scaffoldguided.sampled_N == 0, "can't add length if not masking loops"
            assert self.conf.scaffoldguided.sampled_C == 0, "can't add lemgth if not masking loops"
            assert self.conf.scaffoldguided.sampled_insertion == 0, "can't add length if not masking loops"
            self.mask_loops = False
        else:
            self.mask_loops = True

    def get_ss_adj(self, item):
        """
        Given at item, get the ss tensor and block adjacency matrix for that item
        """
        ss = torch.load(os.path.join(self.scaffold_dir, f'{item.split(".")[0]}_ss.pt'))
        adj = torch.load(
            os.path.join(self.scaffold_dir, f'{item.split(".")[0]}_adj.pt')
        )

        return ss, adj

    def mask_to_segments(self, mask):
        """
        Takes a mask of True (loop) and False (non-loop), and outputs list of tuples (loop or not, length of element)
        """
        segments = []
        begin = -1
        end = -1
        for i in range(mask.shape[0]):
            # Starting edge case
            if i == 0:
                begin = 0
                continue

            if not mask[i] == mask[i - 1]:
                end = i
                if mask[i - 1].item() is True:
                    segments.append(("loop", end - begin))
                else:
                    segments.append(("ss", end - begin))
                begin = i

        # Ending edge case: last segment is length one
        if not end == mask.shape[0]:
            if mask[i].item() is True:
                segments.append(("loop", mask.shape[0] - begin))
            else:
                segments.append(("ss", mask.shape[0] - begin))
        return segments

    def expand_mask(self, mask, segments):
        """
        Function to generate a new mask with dilated loops and N and C terminal additions
        """
        N_add = random.randint(self.sampled_N[0], self.sampled_N[1])
        C_add = random.randint(self.sampled_C[0], self.sampled_C[1])

        output = N_add * [False]
        for ss, length in segments:
            if ss == "ss":
                output.extend(length * [True])
            else:
                # randomly sample insertion length
                ins = random.randint(
                    self.sampled_insertion[0], self.sampled_insertion[1]
                )
                output.extend((length + ins) * [False])
        output.extend(C_add * [False])
        assert torch.sum(torch.tensor(output)) == torch.sum(~mask)
        return torch.tensor(output)

    def expand_ss(self, ss, adj, mask, expanded_mask):
        """
        Given an expanded mask, populate a new ss and adj based on this
        """
        ss_out = torch.ones(expanded_mask.shape[0]) * 3  # set to mask token
        adj_out = torch.full((expanded_mask.shape[0], expanded_mask.shape[0]), 0.0)
        ss_out[expanded_mask] = ss[~mask]
        expanded_mask_2d = torch.full(adj_out.shape, True)
        # mask out loops/insertions, which is ~expanded_mask
        expanded_mask_2d[~expanded_mask, :] = False
        expanded_mask_2d[:, ~expanded_mask] = False

        mask_2d = torch.full(adj.shape, True)
        # mask out loops. This mask is True=loop
        mask_2d[mask, :] = False
        mask_2d[:, mask] = False
        adj_out[expanded_mask_2d] = adj[mask_2d]
        adj_out = adj_out.reshape((expanded_mask.shape[0], expanded_mask.shape[0]))

        return ss_out, adj_out

    def mask_ss_adj(self, ss, adj, expanded_mask):
        """
        Given an expanded ss and adj, mask some number of residues at either end of non-loop ss
        """
        original_mask = torch.clone(expanded_mask)
        if self.ss_mask > 0:
            for i in range(1, self.ss_mask + 1):
                expanded_mask[i:] *= original_mask[:-i]
                expanded_mask[:-i] *= original_mask[i:]

        if self.mask_loops:
            ss[~expanded_mask] = 3
            adj[~expanded_mask, :] = 0
            adj[:, ~expanded_mask] = 0

        # mask adjacency
        adj[~expanded_mask] = 2
        adj[:, ~expanded_mask] = 2

        return ss, adj

    def get_scaffold(self):
        """
        Wrapper method for pulling an item from the list, and preparing ss and block adj features
        """
        
        # Handle determinism. Useful for integration tests
        if self.conf.inference.deterministic:
            torch.manual_seed(self.num_completed)
            np.random.seed(self.num_completed)
            random.seed(self.num_completed)
  
        if self.systematic:
            # reset if num designs > num_scaffolds
            if self.item_n >= len(self.scaffold_list):
                self.item_n = 0
            item = self.scaffold_list[self.item_n]
            self.item_n += 1
        else:
            item = random.choice(self.scaffold_list)
        print("Scaffold constrained based on file: ", item)
        # load files
        ss, adj = self.get_ss_adj(item)
        adj_orig = torch.clone(adj)
        # separate into segments (loop or not)
        mask = torch.where(ss == 2, 1, 0).bool()
        segments = self.mask_to_segments(mask)

        # insert into loops to generate new mask
        expanded_mask = self.expand_mask(mask, segments)

        # expand ss and adj
        ss, adj = self.expand_ss(ss, adj, mask, expanded_mask)

        # finally, mask some proportion of the ss at either end of the non-loop ss blocks
        ss, adj = self.mask_ss_adj(ss, adj, expanded_mask)

        # and then update num_completed
        self.num_completed += 1

        return ss.shape[0], torch.nn.functional.one_hot(ss.long(), num_classes=4), adj


class Target:
    """
    Class to handle targets (fixed chains).
    Inputs:
        - path to pdb file
        - hotspot residues, in the form B10,B12,B60 etc
        - whether or not to crop, and with which method
    Outputs:
        - Dictionary of xyz coordinates, indices, pdb_indices, pdb mask
    """

    def __init__(self, conf: DictConfig, hotspots=None):
        self.pdb = parse_pdb(conf.target_path)

        if hotspots is not None:
            self.hotspots = hotspots
        else:
            self.hotspots = []
        self.pdb["hotspots"] = np.array(
            [
                True if f"{i[0]}{i[1]}" in self.hotspots else False
                for i in self.pdb["pdb_idx"]
            ]
        )

        if conf.contig_crop:
            self.contig_crop(conf.contig_crop)

    def parse_contig(self, contig_crop):
        """
        Takes contig input and parses
        """
        contig_list = []
        for contig in contig_crop[0].split(" "):
            subcon = []
            for crop in contig.split("/"):
                if crop[0].isalpha():
                    subcon.extend(
                        [
                            (crop[0], p)
                            for p in np.arange(
                                int(crop.split("-")[0][1:]), int(crop.split("-")[1]) + 1
                            )
                        ]
                    )
            contig_list.append(subcon)
        return contig_list

    def contig_crop(self, contig_crop, residue_offset=200) -> None:
        """
        Method to take a contig string referring to the receptor and output a pdb dictionary with just this crop
        NB there are two ways to provide inputs:
            - 1) e.g. B1-30,0 B50-60,0. This will add a residue offset between each chunk
            - 2) e.g. B1-30,B50-60,B80-100. This will keep the original indexing of the pdb file.
        Can handle the target being on multiple chains
        """

        # add residue offset between chains if multiple chains in receptor file
        for idx, val in enumerate(self.pdb["pdb_idx"]):
            if idx != 0 and val != self.pdb["pdb_idx"][idx - 1]:
                self.pdb["idx"][idx:] += residue_offset + idx

        # convert contig to mask
        contig_list = self.parse_contig(contig_crop)

        # add residue offset to different parts of contig_list
        for contig in contig_list[1:]:
            start = int(contig[0][1])
            self.pdb["idx"][start:] += residue_offset
        # flatten list
        contig_list = [i for j in contig_list for i in j]
        mask = np.array(
            [True if i in contig_list else False for i in self.pdb["pdb_idx"]]
        )

        # sanity check
        assert np.sum(self.pdb["hotspots"]) == np.sum(
            self.pdb["hotspots"][mask]
        ), "Supplied hotspot residues are missing from the target contig!"
        # crop pdb
        for key, val in self.pdb.items():
            try:
                self.pdb[key] = val[mask]
            except:
                self.pdb[key] = [i for idx, i in enumerate(val) if mask[idx]]
        self.pdb["crop_mask"] = mask

    def get_target(self):
        return self.pdb
