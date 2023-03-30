#!/usr/bin/env python
import os,sys,glob,torch,random
import numpy as np
import argparse
try:
    import pyrosetta
    pyrosetta.init()
    APPROX = False
except:
    print("WARNING: pyRosetta not found, will use an approximate SSE calculation")
    APPROX = True

def main():
    args=get_args()
    assert args.input_pdb or args.pdb_dir is not None, 'Need to provide either an input pdb (--input_pdb) or a path to pdbs (--pdb_dir)'
    assert not (args.input_pdb is not None and args.pdb_dir is not None), 'Need to provide either --input_pdb or --pdb_dir, not both'

    os.makedirs(args.out_dir, exist_ok=True)
    if args.pdb_dir is not None:
        pdbs=glob.glob(f'{args.pdb_dir}/*pdb')
    else:
        pdbs=[args.input_pdb]
    for pdb in pdbs:
        name=os.path.split(pdb)[1][:-4]
        secstruc_dict=extract_secstruc(pdb)
        xyz,_,_ = parse_pdb_torch(pdb)
        ss, idx = ss_to_tensor(secstruc_dict)
        block_adj = construct_block_adj_matrix(torch.FloatTensor(ss), torch.tensor(xyz)).float()
        ss_tens, mask = mask_ss(ss, idx, max_mask=0)
        ss_argmax = torch.argmax(ss_tens[:,:4], dim=1).float()
        torch.save(ss_argmax, os.path.join(args.out_dir, f'{name}_ss.pt'))
        torch.save(block_adj, os.path.join(args.out_dir, f'{name}_adj.pt'))

def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb_dir",required=False, help="path to directory of pdbs. Either pass this or the path to a specific pdb (--input_pdb)", default=None)
    parser.add_argument("--input_pdb", required=False, help="path to input pdb. Either provide this of path to directory of pdbs (--pdb_dir)", default=None)
    parser.add_argument("--out_dir",dest="out_dir", required=True, help='need to specify an output path')
    args = parser.parse_args()
    return args


def extract_secstruc(fn):
    pdb=parse_pdb(fn)
    idx = pdb['idx']
    if APPROX:
        aa_sequence = pdb["seq"]
        secstruct = get_sse(pdb["xyz"][:,1])
    else:
        dssp = pyrosetta.rosetta.core.scoring.dssp
        pose = pyrosetta.io.pose_from_pdb(fn)
        dssp.Dssp(pose).insert_ss_into_pose(pose, True)
        aa_sequence = pose.sequence()
        secstruct = pose.secstruct()
    secstruc_dict = {'sequence':[i for i in aa_sequence],
                     'idx':[int(i) for i in idx],
                     'ss':[i for i in secstruct]}
    return secstruc_dict

def ss_to_tensor(ss):
    """
    Function to convert ss files to indexed tensors
    0 = Helix
    1 = Strand
    2 = Loop
    3 = Mask/unknown
    4 = idx for pdb
    """
    ss_conv = {'H':0,'E':1,'L':2}
    idx = np.array(ss['idx'])
    ss_int = np.array([int(ss_conv[i]) for i in ss['ss']])
    return ss_int, idx

def mask_ss(ss, idx, min_mask = 0, max_mask = 1.0):
    mask_prop = random.uniform(min_mask, max_mask)
    transitions = np.where(ss[:-1] - ss[1:] != 0)[0] #gets last index of each block of ss
    stuck_counter = 0
    while len(ss[ss == 3])/len(ss) < mask_prop or stuck_counter > 100:
        width = random.randint(1,9)
        start = random.choice(transitions)
        offset = random.randint(-8,1)
        try:

            ss[start+offset:start+offset+width] = 3
        except:
            stuck_counter += 1
            pass
    ss = torch.tensor(ss)
    ss = torch.nn.functional.one_hot(ss, num_classes=4)
    ss = torch.cat((ss, torch.tensor(idx)[...,None]), dim=-1)
#     mask = torch.where(torch.argmax(ss[:,:-1], dim=-1) == 3, False, True)
    mask=torch.tensor(np.where(np.argmax(ss[:,:-1].numpy(), axis=-1) == 3))
    return ss, mask

def generate_Cbeta(N,Ca,C):
    # recreate Cb given N,Ca,C
    b = Ca - N 
    c = C - Ca
    a = torch.cross(b, c, dim=-1)
    #Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca
    # fd: below matches sidechain generator (=Rosetta params)
    Cb = -0.57910144*a + 0.5689693*b - 0.5441217*c + Ca

    return Cb

def get_pair_dist(a, b): 
    """calculate pair distances between two sets of points
    
    Parameters
    ----------
    a,b : pytorch tensors of shape [batch,nres,3]
          store Cartesian coordinates of two sets of atoms
    Returns
    -------
    dist : pytorch tensor of shape [batch,nres,nres]
           stores paitwise distances between atoms in a and b
    """

    dist = torch.cdist(a, b, p=2)
    return dist


def construct_block_adj_matrix( sstruct, xyz, cutoff=6, include_loops=False ):
    ''' 
    Given a sstruct specification and backbone coordinates, build a block adjacency matrix.

    Input:
    
        sstruct (torch.FloatTensor): (L) length tensor with numeric encoding of sstruct at each position

        xyz (torch.FloatTensor): (L,3,3) tensor of Cartesian coordinates of backbone N,Ca,C atoms

        cutoff (float): The Cb distance cutoff under which residue pairs are considered adjacent
                        By eye, Nate thinks 6A is a good Cb distance cutoff

    Output:

        block_adj (torch.FloatTensor): (L,L) boolean matrix where adjacent secondary structure contacts are 1
    '''

    L = xyz.shape[0]
    
    # three anchor atoms
    N  = xyz[:,0]
    Ca = xyz[:,1]
    C  = xyz[:,2]
    
    # recreate Cb given N,Ca,C
    Cb = generate_Cbeta(N,Ca,C)
    
    # May need a batch dimension - NRB
    dist = get_pair_dist(Cb,Cb) # [L,L]
    dist[torch.isnan(dist)] = 999.9

    dist += 999.9*torch.eye(L,device=xyz.device)
    # Now we have dist matrix and sstruct specification, turn this into a block adjacency matrix
    # There is probably a way to do this in closed-form with a beautiful einsum but I am going to do the loop approach
    
    # First: Construct a list of segments and the index at which they begin and end
    in_segment = True
    segments = []

    begin = -1
    end = -1

    for i in range(sstruct.shape[0]):
        # Starting edge case
        if i == 0:
            begin = 0 
            continue

        if not sstruct[i] == sstruct[i-1]:
            end = i 
            segments.append( (sstruct[i-1], begin, end) )

            begin = i

    # Ending edge case: last segment is length one
    if not end == sstruct.shape[0]:
        segments.append( (sstruct[-1], begin, sstruct.shape[0]) )


    block_adj = torch.zeros_like(dist)
    for i in range(len(segments)):
        curr_segment = segments[i]

        if curr_segment[0] == 2 and not include_loops: continue

        begin_i = curr_segment[1]
        end_i = curr_segment[2]
        for j in range(i+1, len(segments)):
            j_segment = segments[j]

            if j_segment[0] == 2 and not include_loops: continue

            begin_j = j_segment[1]
            end_j = j_segment[2]

            if torch.any( dist[begin_i:end_i, begin_j:end_j] < cutoff ):
                # Matrix is symmetic
                block_adj[begin_i:end_i, begin_j:end_j] = torch.ones(end_i - begin_i, end_j - begin_j)
                block_adj[begin_j:end_j, begin_i:end_i] = torch.ones(end_j - begin_j, end_i - begin_i)
    return block_adj

def parse_pdb_torch(filename):
    lines = open(filename,'r').readlines()
    return parse_pdb_lines_torch(lines)

#'''
def parse_pdb_lines_torch(lines):

    # indices of residues observed in the structure
    pdb_idx = []
    for l in lines:
      if l[:4]=="ATOM" and l[12:16].strip()=="CA":
        idx = ( l[21:22].strip(), int(l[22:26].strip()) )
        if idx not in pdb_idx:
          pdb_idx.append(idx)
 
    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(pdb_idx), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm == atom:
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break
    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0

    return xyz,mask,np.array(pdb_idx)

def parse_pdb(filename, **kwargs):
    '''extract xyz coords for all heavy atoms'''
    lines = open(filename,'r').readlines()
    return parse_pdb_lines(lines, **kwargs)

def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res = [(l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    seq = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res]
    pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num
    
    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm is not None and tgtatm.strip() == atom.strip(): # ignore whitespace
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break
        
    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0 
    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i,idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)
    
    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]
    seq = np.array(seq)[i_unique]

    out = {'xyz':xyz, # cartesian coordinates, [Lx14]
            'mask':mask, # mask showing which atoms are present in the PDB file, [Lx14]
            'idx':np.array([i[1] for i in pdb_idx]), # residue numbers in the PDB file, [L]
            'seq':np.array(seq), # amino acid sequence, [L]
            'pdb_idx': pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
           }
    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6]=='HETATM' and not (ignore_het_h and l[77]=='H'):
                info_het.append(dict(
                    idx=int(l[7:11]),
                    atom_id=l[12:16],
                    atom_type=l[77],
                    name=l[16:20]
                ))
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out['xyz_het'] = np.array(xyz_het)
        out['info_het'] = info_het

    return out

num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    ]   

aa2num= {x:i for i,x in enumerate(num2aa)}
# full sc atom representation (Nx14)
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE1"," HE2",  None,  None,  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HD2"," HE1"," HE2"," HZ ",  None,  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HE2"," HD2"," HH ",  None,  None,  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), # val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

def get_sse(ca_coord):
  '''
  calculates the SSE of a peptide chain based on the P-SEA algorithm (Labesse 1997)
  code borrowed from biokite: https://github.com/biokit/biokit
  '''
  def vector_dot(v1,v2): return (v1*v2).sum(-1)
  def norm_vector(v): return v / np.linalg.norm(v, axis=-1, keepdims=True)
  def displacement(atoms1, atoms2):
    v1 = np.asarray(atoms1)
    v2 = np.asarray(atoms2)
    if len(v1.shape) <= len(v2.shape):
      diff = v2 - v1
    else:
      diff = -(v1 - v2)
    return diff
  def distance(atoms1, atoms2):
    diff = displacement(atoms1, atoms2)
    return np.sqrt(vector_dot(diff, diff))

  def angle(atoms1, atoms2, atoms3):
    v1 = norm_vector(displacement(atoms1, atoms2))
    v2 = norm_vector(displacement(atoms3, atoms2))
    return np.arccos(vector_dot(v1,v2))

  def dihedral(atoms1, atoms2, atoms3, atoms4):
    v1 = norm_vector(displacement(atoms1, atoms2))
    v2 = norm_vector(displacement(atoms2, atoms3))
    v3 = norm_vector(displacement(atoms3, atoms4))
    
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    
    # Calculation using atan2, to ensure the correct sign of the angle 
    x = vector_dot(n1,n2)
    y = vector_dot(np.cross(n1,n2), v2)
    return np.arctan2(y,x)

  _radians_to_angle = 2*np.pi/360

  _r_helix = ((89-12)*_radians_to_angle, (89+12)*_radians_to_angle)
  _a_helix = ((50-20)*_radians_to_angle, (50+20)*_radians_to_angle)
  _d2_helix = ((5.5-0.5), (5.5+0.5))
  _d3_helix = ((5.3-0.5), (5.3+0.5))
  _d4_helix = ((6.4-0.6), (6.4+0.6))

  _r_strand = ((124-14)*_radians_to_angle, (124+14)*_radians_to_angle)
  _a_strand = ((-180)*_radians_to_angle, (-125)*_radians_to_angle,
              (145)*_radians_to_angle, (180)*_radians_to_angle)
  _d2_strand = ((6.7-0.6), (6.7+0.6))
  _d3_strand = ((9.9-0.9), (9.9+0.9))
  _d4_strand = ((12.4-1.1), (12.4+1.1))

  # Filter all CA atoms in the relevant chain.

  d2i_coord = np.full(( len(ca_coord), 2, 3 ), np.nan)
  d3i_coord = np.full(( len(ca_coord), 2, 3 ), np.nan)
  d4i_coord = np.full(( len(ca_coord), 2, 3 ), np.nan)
  ri_coord = np.full(( len(ca_coord), 3, 3 ), np.nan)
  ai_coord = np.full(( len(ca_coord), 4, 3 ), np.nan)
  
  # The distances and angles are not defined for the entire interval,
  # therefore the indices do not have the full range
  # Values that are not defined are NaN
  for i in range(1, len(ca_coord)-1): d2i_coord[i] = (ca_coord[i-1], ca_coord[i+1])
  for i in range(1, len(ca_coord)-2): d3i_coord[i] = (ca_coord[i-1], ca_coord[i+2])
  for i in range(1, len(ca_coord)-3): d4i_coord[i] = (ca_coord[i-1], ca_coord[i+3])
  for i in range(1, len(ca_coord)-1): ri_coord[i]  = (ca_coord[i-1], ca_coord[i], ca_coord[i+1])
  for i in range(1, len(ca_coord)-2): ai_coord[i]  = (ca_coord[i-1], ca_coord[i], ca_coord[i+1], ca_coord[i+2])
  
  d2i = distance(d2i_coord[:,0], d2i_coord[:,1])
  d3i = distance(d3i_coord[:,0], d3i_coord[:,1])
  d4i = distance(d4i_coord[:,0], d4i_coord[:,1])
  ri = angle(ri_coord[:,0], ri_coord[:,1], ri_coord[:,2])
  ai = dihedral(ai_coord[:,0], ai_coord[:,1], ai_coord[:,2], ai_coord[:,3])
  
  sse = ["L"] * len(ca_coord)
  
  # Annotate helices
  # Find CA that meet criteria for potential helices
  is_pot_helix = np.zeros(len(sse), dtype=bool)
  for i in range(len(sse)):
    if (
            d3i[i] >= _d3_helix[0] and d3i[i] <= _d3_helix[1]
        and d4i[i] >= _d4_helix[0] and d4i[i] <= _d4_helix[1]
        ) or (
            ri[i] >= _r_helix[0] and ri[i] <= _r_helix[1]
        and ai[i] >= _a_helix[0] and ai[i] <= _a_helix[1]
        ):
          is_pot_helix[i] = True
  # Real helices are 5 consecutive helix elements
  is_helix = np.zeros(len(sse), dtype=bool)
  counter = 0
  for i in range(len(sse)):
    if is_pot_helix[i]:
      counter += 1
    else:
      if counter >= 5:
        is_helix[i-counter : i] = True
      counter = 0
  # Extend the helices by one at each end if CA meets extension criteria
  i = 0
  while i < len(sse):
    if is_helix[i]:
      sse[i] = "H"
      if (
          d3i[i-1] >= _d3_helix[0] and d3i[i-1] <= _d3_helix[1]
          ) or (
          ri[i-1] >= _r_helix[0] and ri[i-1] <= _r_helix[1]
          ):
            sse[i-1] = "H"
      sse[i] = "H"
      if (
          d3i[i+1] >= _d3_helix[0] and d3i[i+1] <= _d3_helix[1]
          ) or (
          ri[i+1] >= _r_helix[0] and ri[i+1] <= _r_helix[1]
          ):
            sse[i+1] = "H"
    i += 1
  
  # Annotate sheets
  # Find CA that meet criteria for potential strands
  is_pot_strand = np.zeros(len(sse), dtype=bool)
  for i in range(len(sse)):
    if (    d2i[i] >= _d2_strand[0] and d2i[i] <= _d2_strand[1]
        and d3i[i] >= _d3_strand[0] and d3i[i] <= _d3_strand[1]
        and d4i[i] >= _d4_strand[0] and d4i[i] <= _d4_strand[1]
        ) or (
          ri[i] >= _r_strand[0] and ri[i] <= _r_strand[1]
        and (   (ai[i] >= _a_strand[0] and ai[i] <= _a_strand[1])
              or (ai[i] >= _a_strand[2] and ai[i] <= _a_strand[3]))
        ):
          is_pot_strand[i] = True
  # Real strands are 5 consecutive strand elements,
  # or shorter fragments of at least 3 consecutive strand residues,
  # if they are in hydrogen bond proximity to 5 other residues
  pot_strand_coord = ca_coord[is_pot_strand]
  is_strand = np.zeros(len(sse), dtype=bool)
  counter = 0
  contacts = 0
  for i in range(len(sse)):
    if is_pot_strand[i]:
      counter += 1
      coord = ca_coord[i]
      for strand_coord in ca_coord:
        dist = distance(coord, strand_coord)
        if dist >= 4.2 and dist <= 5.2:
          contacts += 1
    else:
      if counter >= 4:
        is_strand[i-counter : i] = True
      elif counter == 3 and contacts >= 5:
        is_strand[i-counter : i] = True
      counter = 0
      contacts = 0
  # Extend the strands by one at each end if CA meets extension criteria
  i = 0
  while i < len(sse):
    if is_strand[i]:
      sse[i] = "E"
      if d3i[i-1] >= _d3_strand[0] and d3i[i-1] <= _d3_strand[1]:
        sse[i-1] = "E"
      sse[i] = "E"
      if d3i[i+1] >= _d3_strand[0] and d3i[i+1] <= _d3_strand[1]:
        sse[i+1] = "E"
    i += 1
  return sse

if __name__ == "__main__":
    main()
