import scipy.sparse
from rfdiffusion.chemical import *
from rfdiffusion.scoring import *


def generate_Cbeta(N, Ca, C):
    # recreate Cb given N,Ca,C
    b = Ca - N
    c = C - Ca
    a = torch.cross(b, c, dim=-1)
    # These are the values used during training
    Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca
    # fd: below matches sidechain generator (=Rosetta params)
    # Cb = -0.57910144 * a + 0.5689693 * b - 0.5441217 * c + Ca

    return Cb


def th_ang_v(ab, bc, eps: float = 1e-8):
    def th_norm(x, eps: float = 1e-8):
        return x.square().sum(-1, keepdim=True).add(eps).sqrt()

    def th_N(x, alpha: float = 0):
        return x / th_norm(x).add(alpha)

    ab, bc = th_N(ab), th_N(bc)
    cos_angle = torch.clamp((ab * bc).sum(-1), -1, 1)
    sin_angle = torch.sqrt(1 - cos_angle.square() + eps)
    dih = torch.stack((cos_angle, sin_angle), -1)
    return dih


def th_dih_v(ab, bc, cd):
    def th_cross(a, b):
        a, b = torch.broadcast_tensors(a, b)
        return torch.cross(a, b, dim=-1)

    def th_norm(x, eps: float = 1e-8):
        return x.square().sum(-1, keepdim=True).add(eps).sqrt()

    def th_N(x, alpha: float = 0):
        return x / th_norm(x).add(alpha)

    ab, bc, cd = th_N(ab), th_N(bc), th_N(cd)
    n1 = th_N(th_cross(ab, bc))
    n2 = th_N(th_cross(bc, cd))
    sin_angle = (th_cross(n1, bc) * n2).sum(-1)
    cos_angle = (n1 * n2).sum(-1)
    dih = torch.stack((cos_angle, sin_angle), -1)
    return dih


def th_dih(a, b, c, d):
    return th_dih_v(a - b, b - c, c - d)


# More complicated version splits error in CA-N and CA-C (giving more accurate CB position)
# It returns the rigid transformation from local frame to global frame
def rigid_from_3_points(N, Ca, C, non_ideal=False, eps=1e-8):
    # N, Ca, C - [B,L, 3]
    # R - [B,L, 3, 3], det(R)=1, inv(R) = R.T, R is a rotation matrix
    B, L = N.shape[:2]

    v1 = C - Ca
    v2 = N - Ca
    e1 = v1 / (torch.norm(v1, dim=-1, keepdim=True) + eps)
    u2 = v2 - (torch.einsum("bli, bli -> bl", e1, v2)[..., None] * e1)
    e2 = u2 / (torch.norm(u2, dim=-1, keepdim=True) + eps)
    e3 = torch.cross(e1, e2, dim=-1)
    R = torch.cat(
        [e1[..., None], e2[..., None], e3[..., None]], axis=-1
    )  # [B,L,3,3] - rotation matrix

    if non_ideal:
        v2 = v2 / (torch.norm(v2, dim=-1, keepdim=True) + eps)
        cosref = torch.sum(e1 * v2, dim=-1)  # cosine of current N-CA-C bond angle
        costgt = cos_ideal_NCAC.item()
        cos2del = torch.clamp(
            cosref * costgt
            + torch.sqrt((1 - cosref * cosref) * (1 - costgt * costgt) + eps),
            min=-1.0,
            max=1.0,
        )
        cosdel = torch.sqrt(0.5 * (1 + cos2del) + eps)
        sindel = torch.sign(costgt - cosref) * torch.sqrt(1 - 0.5 * (1 + cos2del) + eps)
        Rp = torch.eye(3, device=N.device).repeat(B, L, 1, 1)
        Rp[:, :, 0, 0] = cosdel
        Rp[:, :, 0, 1] = -sindel
        Rp[:, :, 1, 0] = sindel
        Rp[:, :, 1, 1] = cosdel

        R = torch.einsum("blij,bljk->blik", R, Rp)

    return R, Ca


def get_tor_mask(seq, torsion_indices, mask_in=None):
    B, L = seq.shape[:2]
    tors_mask = torch.ones((B, L, 10), dtype=torch.bool, device=seq.device)
    tors_mask[..., 3:7] = torsion_indices[seq, :, -1] > 0
    tors_mask[:, 0, 1] = False
    tors_mask[:, -1, 0] = False

    # mask for additional angles
    tors_mask[:, :, 7] = seq != aa2num["GLY"]
    tors_mask[:, :, 8] = seq != aa2num["GLY"]
    tors_mask[:, :, 9] = torch.logical_and(seq != aa2num["GLY"], seq != aa2num["ALA"])
    tors_mask[:, :, 9] = torch.logical_and(tors_mask[:, :, 9], seq != aa2num["UNK"])
    tors_mask[:, :, 9] = torch.logical_and(tors_mask[:, :, 9], seq != aa2num["MAS"])

    if mask_in != None:
        # mask for missing atoms
        # chis
        ti0 = torch.gather(mask_in, 2, torsion_indices[seq, :, 0])
        ti1 = torch.gather(mask_in, 2, torsion_indices[seq, :, 1])
        ti2 = torch.gather(mask_in, 2, torsion_indices[seq, :, 2])
        ti3 = torch.gather(mask_in, 2, torsion_indices[seq, :, 3])
        is_valid = torch.stack((ti0, ti1, ti2, ti3), dim=-2).all(dim=-1)
        tors_mask[..., 3:7] = torch.logical_and(tors_mask[..., 3:7], is_valid)
        tors_mask[:, :, 7] = torch.logical_and(
            tors_mask[:, :, 7], mask_in[:, :, 4]
        )  # CB exist?
        tors_mask[:, :, 8] = torch.logical_and(
            tors_mask[:, :, 8], mask_in[:, :, 4]
        )  # CB exist?
        tors_mask[:, :, 9] = torch.logical_and(
            tors_mask[:, :, 9], mask_in[:, :, 5]
        )  # XG exist?

    return tors_mask


def get_torsions(
    xyz_in, seq, torsion_indices, torsion_can_flip, ref_angles, mask_in=None
):
    B, L = xyz_in.shape[:2]

    tors_mask = get_tor_mask(seq, torsion_indices, mask_in)

    # torsions to restrain to 0 or 180degree
    tors_planar = torch.zeros((B, L, 10), dtype=torch.bool, device=xyz_in.device)
    tors_planar[:, :, 5] = seq == aa2num["TYR"]  # TYR chi 3 should be planar

    # idealize given xyz coordinates before computing torsion angles
    xyz = xyz_in.clone()
    Rs, Ts = rigid_from_3_points(xyz[..., 0, :], xyz[..., 1, :], xyz[..., 2, :])
    Nideal = torch.tensor([-0.5272, 1.3593, 0.000], device=xyz_in.device)
    Cideal = torch.tensor([1.5233, 0.000, 0.000], device=xyz_in.device)
    xyz[..., 0, :] = torch.einsum("brij,j->bri", Rs, Nideal) + Ts
    xyz[..., 2, :] = torch.einsum("brij,j->bri", Rs, Cideal) + Ts

    torsions = torch.zeros((B, L, 10, 2), device=xyz.device)
    # avoid undefined angles for H generation
    torsions[:, 0, 1, 0] = 1.0
    torsions[:, -1, 0, 0] = 1.0

    # omega
    torsions[:, :-1, 0, :] = th_dih(
        xyz[:, :-1, 1, :], xyz[:, :-1, 2, :], xyz[:, 1:, 0, :], xyz[:, 1:, 1, :]
    )
    # phi
    torsions[:, 1:, 1, :] = th_dih(
        xyz[:, :-1, 2, :], xyz[:, 1:, 0, :], xyz[:, 1:, 1, :], xyz[:, 1:, 2, :]
    )
    # psi
    torsions[:, :, 2, :] = -1 * th_dih(
        xyz[:, :, 0, :], xyz[:, :, 1, :], xyz[:, :, 2, :], xyz[:, :, 3, :]
    )

    # chis
    ti0 = torch.gather(xyz, 2, torsion_indices[seq, :, 0, None].repeat(1, 1, 1, 3))
    ti1 = torch.gather(xyz, 2, torsion_indices[seq, :, 1, None].repeat(1, 1, 1, 3))
    ti2 = torch.gather(xyz, 2, torsion_indices[seq, :, 2, None].repeat(1, 1, 1, 3))
    ti3 = torch.gather(xyz, 2, torsion_indices[seq, :, 3, None].repeat(1, 1, 1, 3))
    torsions[:, :, 3:7, :] = th_dih(ti0, ti1, ti2, ti3)

    # CB bend
    NC = 0.5 * (xyz[:, :, 0, :3] + xyz[:, :, 2, :3])
    CA = xyz[:, :, 1, :3]
    CB = xyz[:, :, 4, :3]
    t = th_ang_v(CB - CA, NC - CA)
    t0 = ref_angles[seq][..., 0, :]
    torsions[:, :, 7, :] = torch.stack(
        (torch.sum(t * t0, dim=-1), t[..., 0] * t0[..., 1] - t[..., 1] * t0[..., 0]),
        dim=-1,
    )

    # CB twist
    NCCA = NC - CA
    NCp = xyz[:, :, 2, :3] - xyz[:, :, 0, :3]
    NCpp = (
        NCp
        - torch.sum(NCp * NCCA, dim=-1, keepdim=True)
        / torch.sum(NCCA * NCCA, dim=-1, keepdim=True)
        * NCCA
    )
    t = th_ang_v(CB - CA, NCpp)
    t0 = ref_angles[seq][..., 1, :]
    torsions[:, :, 8, :] = torch.stack(
        (torch.sum(t * t0, dim=-1), t[..., 0] * t0[..., 1] - t[..., 1] * t0[..., 0]),
        dim=-1,
    )

    # CG bend
    CG = xyz[:, :, 5, :3]
    t = th_ang_v(CG - CB, CA - CB)
    t0 = ref_angles[seq][..., 2, :]
    torsions[:, :, 9, :] = torch.stack(
        (torch.sum(t * t0, dim=-1), t[..., 0] * t0[..., 1] - t[..., 1] * t0[..., 0]),
        dim=-1,
    )

    mask0 = torch.isnan(torsions[..., 0]).nonzero()
    mask1 = torch.isnan(torsions[..., 1]).nonzero()
    torsions[mask0[:, 0], mask0[:, 1], mask0[:, 2], 0] = 1.0
    torsions[mask1[:, 0], mask1[:, 1], mask1[:, 2], 1] = 0.0

    # alt chis
    torsions_alt = torsions.clone()
    torsions_alt[torsion_can_flip[seq, :]] *= -1

    return torsions, torsions_alt, tors_mask, tors_planar


def get_tips(xyz, seq):
    B, L = xyz.shape[:2]

    xyz_tips = torch.gather(
        xyz, 2, tip_indices.to(xyz.device)[seq][:, :, None, None].expand(-1, -1, -1, 3)
    ).reshape(B, L, 3)
    mask = ~(torch.isnan(xyz_tips[:, :, 0]))
    if torch.isnan(xyz_tips).any():  # replace NaN tip atom with virtual Cb atom
        # three anchor atoms
        N = xyz[:, :, 0]
        Ca = xyz[:, :, 1]
        C = xyz[:, :, 2]

        # recreate Cb given N,Ca,C
        b = Ca - N
        c = C - Ca
        a = torch.cross(b, c, dim=-1)
        Cb = -0.58273431 * a + 0.56802827 * b - 0.54067466 * c + Ca

        xyz_tips = torch.where(torch.isnan(xyz_tips), Cb, xyz_tips)
    return xyz_tips, mask


# process ideal frames
def make_frame(X, Y):
    Xn = X / torch.linalg.norm(X)
    Y = Y - torch.dot(Y, Xn) * Xn
    Yn = Y / torch.linalg.norm(Y)
    Z = torch.cross(Xn, Yn)
    Zn = Z / torch.linalg.norm(Z)

    return torch.stack((Xn, Yn, Zn), dim=-1)


def cross_product_matrix(u):
    B, L = u.shape[:2]
    matrix = torch.zeros((B, L, 3, 3), device=u.device)
    matrix[:, :, 0, 1] = -u[..., 2]
    matrix[:, :, 0, 2] = u[..., 1]
    matrix[:, :, 1, 0] = u[..., 2]
    matrix[:, :, 1, 2] = -u[..., 0]
    matrix[:, :, 2, 0] = -u[..., 1]
    matrix[:, :, 2, 1] = u[..., 0]
    return matrix


# writepdb
def writepdb(
    filename, atoms, seq, binderlen=None, idx_pdb=None, bfacts=None, chain_idx=None
):
    f = open(filename, "w")
    ctr = 1
    scpu = seq.cpu().squeeze()
    atomscpu = atoms.cpu().squeeze()
    if bfacts is None:
        bfacts = torch.zeros(atomscpu.shape[0])
    if idx_pdb is None:
        idx_pdb = 1 + torch.arange(atomscpu.shape[0])

    Bfacts = torch.clamp(bfacts.cpu(), 0, 1)
    for i, s in enumerate(scpu):
        if chain_idx is None:
            if binderlen is not None:
                if i < binderlen:
                    chain = "A"
                else:
                    chain = "B"
            elif binderlen is None:
                chain = "A"
        else:
            chain = chain_idx[i]
        if len(atomscpu.shape) == 2:
            f.write(
                "%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                % (
                    "ATOM",
                    ctr,
                    " CA ",
                    num2aa[s],
                    chain,
                    idx_pdb[i],
                    atomscpu[i, 0],
                    atomscpu[i, 1],
                    atomscpu[i, 2],
                    1.0,
                    Bfacts[i],
                )
            )
            ctr += 1
        elif atomscpu.shape[1] == 3:
            for j, atm_j in enumerate([" N  ", " CA ", " C  "]):
                f.write(
                    "%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                    % (
                        "ATOM",
                        ctr,
                        atm_j,
                        num2aa[s],
                        chain,
                        idx_pdb[i],
                        atomscpu[i, j, 0],
                        atomscpu[i, j, 1],
                        atomscpu[i, j, 2],
                        1.0,
                        Bfacts[i],
                    )
                )
                ctr += 1
        elif atomscpu.shape[1] == 4:
            for j, atm_j in enumerate([" N  ", " CA ", " C  ", " O  "]):
                f.write(
                    "%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                    % (
                        "ATOM",
                        ctr,
                        atm_j,
                        num2aa[s],
                        chain,
                        idx_pdb[i],
                        atomscpu[i, j, 0],
                        atomscpu[i, j, 1],
                        atomscpu[i, j, 2],
                        1.0,
                        Bfacts[i],
                    )
                )
                ctr += 1

        else:
            natoms = atomscpu.shape[1]
            if natoms != 14 and natoms != 27:
                print("bad size!", atoms.shape)
                assert False
            atms = aa2long[s]
            # his prot hack
            if (
                s == 8
                and torch.linalg.norm(atomscpu[i, 9, :] - atomscpu[i, 5, :]) < 1.7
            ):
                atms = (
                    " N  ",
                    " CA ",
                    " C  ",
                    " O  ",
                    " CB ",
                    " CG ",
                    " NE2",
                    " CD2",
                    " CE1",
                    " ND1",
                    None,
                    None,
                    None,
                    None,
                    " H  ",
                    " HA ",
                    "1HB ",
                    "2HB ",
                    " HD2",
                    " HE1",
                    " HD1",
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                )  # his_d

            for j, atm_j in enumerate(atms):
                if (
                    j < natoms and atm_j is not None
                ):  # and not torch.isnan(atomscpu[i,j,:]).any()):
                    f.write(
                        "%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                        % (
                            "ATOM",
                            ctr,
                            atm_j,
                            num2aa[s],
                            chain,
                            idx_pdb[i],
                            atomscpu[i, j, 0],
                            atomscpu[i, j, 1],
                            atomscpu[i, j, 2],
                            1.0,
                            Bfacts[i],
                        )
                    )
                    ctr += 1


# resolve tip atom indices
tip_indices = torch.full((22,), 0)
for i in range(22):
    tip_atm = aa2tip[i]
    atm_long = aa2long[i]
    tip_indices[i] = atm_long.index(tip_atm)

# resolve torsion indices
torsion_indices = torch.full((22, 4, 4), 0)
torsion_can_flip = torch.full((22, 10), False, dtype=torch.bool)
for i in range(22):
    i_l, i_a = aa2long[i], aa2longalt[i]
    for j in range(4):
        if torsions[i][j] is None:
            continue
        for k in range(4):
            a = torsions[i][j][k]
            torsion_indices[i, j, k] = i_l.index(a)
            if i_l.index(a) != i_a.index(a):
                torsion_can_flip[i, 3 + j] = True  ##bb tors never flip
# HIS is a special case
torsion_can_flip[8, 4] = False

# build the mapping from atoms in the full rep (Nx27) to the "alternate" rep
allatom_mask = torch.zeros((22, 27), dtype=torch.bool)
long2alt = torch.zeros((22, 27), dtype=torch.long)
for i in range(22):
    i_l, i_lalt = aa2long[i], aa2longalt[i]
    for j, a in enumerate(i_l):
        if a is None:
            long2alt[i, j] = j
        else:
            long2alt[i, j] = i_lalt.index(a)
            allatom_mask[i, j] = True

# bond graph traversal
num_bonds = torch.zeros((22, 27, 27), dtype=torch.long)
for i in range(22):
    num_bonds_i = np.zeros((27, 27))
    for bnamei, bnamej in aabonds[i]:
        bi, bj = aa2long[i].index(bnamei), aa2long[i].index(bnamej)
        num_bonds_i[bi, bj] = 1
    num_bonds_i = scipy.sparse.csgraph.shortest_path(num_bonds_i, directed=False)
    num_bonds_i[num_bonds_i >= 4] = 4
    num_bonds[i, ...] = torch.tensor(num_bonds_i)


# LJ/LK scoring parameters
ljlk_parameters = torch.zeros((22, 27, 5), dtype=torch.float)
lj_correction_parameters = torch.zeros(
    (22, 27, 4), dtype=bool
)  # donor/acceptor/hpol/disulf
for i in range(22):
    for j, a in enumerate(aa2type[i]):
        if a is not None:
            ljlk_parameters[i, j, :] = torch.tensor(type2ljlk[a])
            lj_correction_parameters[i, j, 0] = (type2hb[a] == HbAtom.DO) + (
                type2hb[a] == HbAtom.DA
            )
            lj_correction_parameters[i, j, 1] = (type2hb[a] == HbAtom.AC) + (
                type2hb[a] == HbAtom.DA
            )
            lj_correction_parameters[i, j, 2] = type2hb[a] == HbAtom.HP
            lj_correction_parameters[i, j, 3] = a == "SH1" or a == "HS"


# hbond scoring parameters
def donorHs(D, bonds, atoms):
    dHs = []
    for i, j in bonds:
        if i == D:
            idx_j = atoms.index(j)
            if idx_j >= 14:  # if atom j is a hydrogen
                dHs.append(idx_j)
        if j == D:
            idx_i = atoms.index(i)
            if idx_i >= 14:  # if atom j is a hydrogen
                dHs.append(idx_i)
    assert len(dHs) > 0
    return dHs


def acceptorBB0(A, hyb, bonds, atoms):
    if hyb == HbHybType.SP2:
        for i, j in bonds:
            if i == A:
                B = atoms.index(j)
                if B < 14:
                    break
            if j == A:
                B = atoms.index(i)
                if B < 14:
                    break
        for i, j in bonds:
            if i == atoms[B]:
                B0 = atoms.index(j)
                if B0 < 14:
                    break
            if j == atoms[B]:
                B0 = atoms.index(i)
                if B0 < 14:
                    break
    elif hyb == HbHybType.SP3 or hyb == HbHybType.RING:
        for i, j in bonds:
            if i == A:
                B = atoms.index(j)
                if B < 14:
                    break
            if j == A:
                B = atoms.index(i)
                if B < 14:
                    break
        for i, j in bonds:
            if i == A and j != atoms[B]:
                B0 = atoms.index(j)
                break
            if j == A and i != atoms[B]:
                B0 = atoms.index(i)
                break

    return B, B0


hbtypes = torch.full(
    (22, 27, 3), -1, dtype=torch.long
)  # (donortype, acceptortype, acchybtype)
hbbaseatoms = torch.full(
    (22, 27, 2), -1, dtype=torch.long
)  # (B,B0) for acc; (D,-1) for don
hbpolys = torch.zeros(
    (HbDonType.NTYPES, HbAccType.NTYPES, 3, 15)
)  # weight,xmin,xmax,ymin,ymax,c9,...,c0

for i in range(22):
    for j, a in enumerate(aa2type[i]):
        if a in type2dontype:
            j_hs = donorHs(aa2long[i][j], aabonds[i], aa2long[i])
            for j_h in j_hs:
                hbtypes[i, j_h, 0] = type2dontype[a]
                hbbaseatoms[i, j_h, 0] = j
        if a in type2acctype:
            j_b, j_b0 = acceptorBB0(
                aa2long[i][j], type2hybtype[a], aabonds[i], aa2long[i]
            )
            hbtypes[i, j, 1] = type2acctype[a]
            hbtypes[i, j, 2] = type2hybtype[a]
            hbbaseatoms[i, j, 0] = j_b
            hbbaseatoms[i, j, 1] = j_b0

for i in range(HbDonType.NTYPES):
    for j in range(HbAccType.NTYPES):
        weight = dontype2wt[i] * acctype2wt[j]

        pdist, pbah, pahd = hbtypepair2poly[(i, j)]
        xrange, yrange, coeffs = hbpolytype2coeffs[pdist]
        hbpolys[i, j, 0, 0] = weight
        hbpolys[i, j, 0, 1:3] = torch.tensor(xrange)
        hbpolys[i, j, 0, 3:5] = torch.tensor(yrange)
        hbpolys[i, j, 0, 5:] = torch.tensor(coeffs)
        xrange, yrange, coeffs = hbpolytype2coeffs[pahd]
        hbpolys[i, j, 1, 0] = weight
        hbpolys[i, j, 1, 1:3] = torch.tensor(xrange)
        hbpolys[i, j, 1, 3:5] = torch.tensor(yrange)
        hbpolys[i, j, 1, 5:] = torch.tensor(coeffs)
        xrange, yrange, coeffs = hbpolytype2coeffs[pbah]
        hbpolys[i, j, 2, 0] = weight
        hbpolys[i, j, 2, 1:3] = torch.tensor(xrange)
        hbpolys[i, j, 2, 3:5] = torch.tensor(yrange)
        hbpolys[i, j, 2, 5:] = torch.tensor(coeffs)

# kinematic parameters
base_indices = torch.full((22, 27), 0, dtype=torch.long)
xyzs_in_base_frame = torch.ones((22, 27, 4))
RTs_by_torsion = torch.eye(4).repeat(22, 7, 1, 1)
reference_angles = torch.ones((22, 3, 2))

for i in range(22):
    i_l = aa2long[i]
    for name, base, coords in ideal_coords[i]:
        idx = i_l.index(name)
        base_indices[i, idx] = base
        xyzs_in_base_frame[i, idx, :3] = torch.tensor(coords)

    # omega frame
    RTs_by_torsion[i, 0, :3, :3] = torch.eye(3)
    RTs_by_torsion[i, 0, :3, 3] = torch.zeros(3)

    # phi frame
    RTs_by_torsion[i, 1, :3, :3] = make_frame(
        xyzs_in_base_frame[i, 0, :3] - xyzs_in_base_frame[i, 1, :3],
        torch.tensor([1.0, 0.0, 0.0]),
    )
    RTs_by_torsion[i, 1, :3, 3] = xyzs_in_base_frame[i, 0, :3]

    # psi frame
    RTs_by_torsion[i, 2, :3, :3] = make_frame(
        xyzs_in_base_frame[i, 2, :3] - xyzs_in_base_frame[i, 1, :3],
        xyzs_in_base_frame[i, 1, :3] - xyzs_in_base_frame[i, 0, :3],
    )
    RTs_by_torsion[i, 2, :3, 3] = xyzs_in_base_frame[i, 2, :3]

    # chi1 frame
    if torsions[i][0] is not None:
        a0, a1, a2 = torsion_indices[i, 0, 0:3]
        RTs_by_torsion[i, 3, :3, :3] = make_frame(
            xyzs_in_base_frame[i, a2, :3] - xyzs_in_base_frame[i, a1, :3],
            xyzs_in_base_frame[i, a0, :3] - xyzs_in_base_frame[i, a1, :3],
        )
        RTs_by_torsion[i, 3, :3, 3] = xyzs_in_base_frame[i, a2, :3]

    # chi2~4 frame
    for j in range(1, 4):
        if torsions[i][j] is not None:
            a2 = torsion_indices[i, j, 2]
            if (i == 18 and j == 2) or (
                i == 8 and j == 2
            ):  # TYR CZ-OH & HIS CE1-HE1 a special case
                a0, a1 = torsion_indices[i, j, 0:2]
                RTs_by_torsion[i, 3 + j, :3, :3] = make_frame(
                    xyzs_in_base_frame[i, a2, :3] - xyzs_in_base_frame[i, a1, :3],
                    xyzs_in_base_frame[i, a0, :3] - xyzs_in_base_frame[i, a1, :3],
                )
            else:
                RTs_by_torsion[i, 3 + j, :3, :3] = make_frame(
                    xyzs_in_base_frame[i, a2, :3],
                    torch.tensor([-1.0, 0.0, 0.0]),
                )
            RTs_by_torsion[i, 3 + j, :3, 3] = xyzs_in_base_frame[i, a2, :3]

    # CB/CG angles
    NCr = 0.5 * (xyzs_in_base_frame[i, 0, :3] + xyzs_in_base_frame[i, 2, :3])
    CAr = xyzs_in_base_frame[i, 1, :3]
    CBr = xyzs_in_base_frame[i, 4, :3]
    CGr = xyzs_in_base_frame[i, 5, :3]
    reference_angles[i, 0, :] = th_ang_v(CBr - CAr, NCr - CAr)
    NCp = xyzs_in_base_frame[i, 2, :3] - xyzs_in_base_frame[i, 0, :3]
    NCpp = NCp - torch.dot(NCp, NCr) / torch.dot(NCr, NCr) * NCr
    reference_angles[i, 1, :] = th_ang_v(CBr - CAr, NCpp)
    reference_angles[i, 2, :] = th_ang_v(CGr, torch.tensor([-1.0, 0.0, 0.0]))

N_BACKBONE_ATOMS = 3
N_HEAVY = 14


def writepdb_multi(
    filename,
    atoms_stack,
    bfacts,
    seq_stack,
    backbone_only=False,
    chain_ids=None,
    use_hydrogens=True,
):
    """
    Function for writing multiple structural states of the same sequence into a single
    pdb file.
    """

    f = open(filename, "w")

    if seq_stack.ndim != 2:
        T = atoms_stack.shape[0]
        seq_stack = torch.tile(seq_stack, (T, 1))
    seq_stack = seq_stack.cpu()
    for atoms, scpu in zip(atoms_stack, seq_stack):
        ctr = 1
        atomscpu = atoms.cpu()
        Bfacts = torch.clamp(bfacts.cpu(), 0, 1)
        for i, s in enumerate(scpu):
            atms = aa2long[s]
            for j, atm_j in enumerate(atms):
                if backbone_only and j >= N_BACKBONE_ATOMS:
                    break
                if not use_hydrogens and j >= N_HEAVY:
                    break
                if (atm_j is None) or (torch.all(torch.isnan(atomscpu[i, j]))):
                    continue
                chain_id = "A"
                if chain_ids is not None:
                    chain_id = chain_ids[i]
                f.write(
                    "%-6s%5s %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                    % (
                        "ATOM",
                        ctr,
                        atm_j,
                        num2aa[s],
                        chain_id,
                        i + 1,
                        atomscpu[i, j, 0],
                        atomscpu[i, j, 1],
                        atomscpu[i, j, 2],
                        1.0,
                        Bfacts[i],
                    )
                )
                ctr += 1

        f.write("ENDMDL\n")

def calc_rmsd(xyz1, xyz2, eps=1e-6):
    """
    Calculates RMSD between two sets of atoms (L, 3)
    """
    # center to CA centroid
    xyz1 = xyz1 - xyz1.mean(0)
    xyz2 = xyz2 - xyz2.mean(0)

    # Computation of the covariance matrix
    C = xyz2.T @ xyz1

    # Compute otimal rotation matrix using SVD
    V, S, W = np.linalg.svd(C)

    # get sign to ensure right-handedness
    d = np.ones([3,3])
    d[:,-1] = np.sign(np.linalg.det(V)*np.linalg.det(W))

    # Rotation matrix U
    U = (d*V) @ W

    # Rotate xyz2
    xyz2_ = xyz2 @ U
    L = xyz2_.shape[0]
    rmsd = np.sqrt(np.sum((xyz2_-xyz1)*(xyz2_-xyz1), axis=(0,1)) / L + eps)

    return rmsd, U
