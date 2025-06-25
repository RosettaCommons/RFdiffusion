import torch
import numpy as np 
from rfdiffusion.util import generate_Cbeta

class Potential:
    '''
        Interface class that defines the functions a potential must implement
    '''

    def compute(self, xyz):
        '''
            Given the current structure of the model prediction, return the current
            potential as a PyTorch tensor with a single entry

            Args:
                xyz (torch.tensor, size: [L,27,3]: The current coordinates of the sample
            
            Returns:
                potential (torch.tensor, size: [1]): A potential whose value will be MAXIMIZED
                                                     by taking a step along it's gradient
        '''
        raise NotImplementedError('Potential compute function was not overwritten')

class monomer_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging monomer compactness

        Written by DJ and refactored into a class by NRB
    '''

    def __init__(self, weight=1, min_dist=15):

        self.weight   = weight
        self.min_dist = min_dist

    def compute(self, xyz):
        Ca = xyz[:,1] # [L,3]

        centroid = torch.mean(Ca, dim=0, keepdim=True) # [1,3]

        dgram = torch.cdist(Ca[None,...].contiguous(), centroid[None,...].contiguous(), p=2) # [1,L,1,3]

        dgram = torch.maximum(self.min_dist * torch.ones_like(dgram.squeeze(0)), dgram.squeeze(0)) # [L,1,3]

        rad_of_gyration = torch.sqrt( torch.sum(torch.square(dgram)) / Ca.shape[0] ) # [1]

        return -1 * self.weight * rad_of_gyration

class binder_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging binder compactness

        Author: NRB
    '''

    def __init__(self, binderlen, weight=1, min_dist=15):

        self.binderlen = binderlen
        self.min_dist  = min_dist
        self.weight    = weight

    def compute(self, xyz):
        
        # Only look at binder residues
        Ca = xyz[:self.binderlen,1] # [Lb,3]

        centroid = torch.mean(Ca, dim=0, keepdim=True) # [1,3]

        # cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), centroid[None,...].contiguous(), p=2) # [1,Lb,1,3]

        dgram = torch.maximum(self.min_dist * torch.ones_like(dgram.squeeze(0)), dgram.squeeze(0)) # [Lb,1,3]

        rad_of_gyration = torch.sqrt( torch.sum(torch.square(dgram)) / Ca.shape[0] ) # [1]

        return -1 * self.weight * rad_of_gyration


class dimer_ROG(Potential):
    '''
        Radius of Gyration potential for encouraging compactness of both monomers when designing dimers

        Author: PV
    '''

    def __init__(self, binderlen, weight=1, min_dist=15):

        self.binderlen = binderlen
        self.min_dist  = min_dist
        self.weight    = weight

    def compute(self, xyz):

        # Only look at monomer 1 residues
        Ca_m1 = xyz[:self.binderlen,1] # [Lb,3]
        
        # Only look at monomer 2 residues
        Ca_m2 = xyz[self.binderlen:,1] # [Lb,3]

        centroid_m1 = torch.mean(Ca_m1, dim=0, keepdim=True) # [1,3]
        centroid_m2 = torch.mean(Ca_m1, dim=0, keepdim=True) # [1,3]

        # cdist needs a batch dimension - NRB
        #This calculates RoG for Monomer 1
        dgram_m1 = torch.cdist(Ca_m1[None,...].contiguous(), centroid_m1[None,...].contiguous(), p=2) # [1,Lb,1,3]
        dgram_m1 = torch.maximum(self.min_dist * torch.ones_like(dgram_m1.squeeze(0)), dgram_m1.squeeze(0)) # [Lb,1,3]
        rad_of_gyration_m1 = torch.sqrt( torch.sum(torch.square(dgram_m1)) / Ca_m1.shape[0] ) # [1]

        # cdist needs a batch dimension - NRB
        #This calculates RoG for Monomer 2
        dgram_m2 = torch.cdist(Ca_m2[None,...].contiguous(), centroid_m2[None,...].contiguous(), p=2) # [1,Lb,1,3]
        dgram_m2 = torch.maximum(self.min_dist * torch.ones_like(dgram_m2.squeeze(0)), dgram_m2.squeeze(0)) # [Lb,1,3]
        rad_of_gyration_m2 = torch.sqrt( torch.sum(torch.square(dgram_m2)) / Ca_m2.shape[0] ) # [1]

        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return -1 * self.weight * (rad_of_gyration_m1 + rad_of_gyration_m2)/2

class binder_ncontacts(Potential):
    '''
        Differentiable way to maximise number of contacts within a protein
        
        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

    '''

    def __init__(self, binderlen, weight=1, r_0=8, d_0=4):

        self.binderlen = binderlen
        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, xyz):

        # Only look at binder Ca residues
        Ca = xyz[:self.binderlen,1] # [Lb,3]
        
        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        binder_ncontacts = (1 - numerator) / (1 - denominator)
        
        print("BINDER CONTACTS:", binder_ncontacts.sum())
        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * binder_ncontacts.sum()

class interface_ncontacts(Potential):

    '''
        Differentiable way to maximise number of contacts between binder and target
        
        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html

        Author: PV
    '''


    def __init__(self, binderlen, weight=1, r_0=8, d_0=6):

        self.binderlen = binderlen
        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0

    def compute(self, xyz):

        # Extract binder Ca residues
        Ca_b = xyz[:self.binderlen,1] # [Lb,3]

        # Extract target Ca residues
        Ca_t = xyz[self.binderlen:,1] # [Lt,3]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca_b[None,...].contiguous(), Ca_t[None,...].contiguous(), p=2) # [1,Lb,Lt]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)
        interface_ncontacts = (1 - numerator) / (1 - denominator)
        #Potential is the sum of values in the tensor
        interface_ncontacts = interface_ncontacts.sum()

        print("INTERFACE CONTACTS:", interface_ncontacts.sum())

        return self.weight * interface_ncontacts


class monomer_contacts(Potential):
    '''
        Differentiable way to maximise number of contacts within a protein

        Motivation is given here: https://www.plumed.org/doc-v2.7/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html
        Author: PV

        NOTE: This function sometimes produces NaN's -- added check in reverse diffusion for nan grads
    '''

    def __init__(self, weight=1, r_0=8, d_0=2, eps=1e-6):

        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0
        self.eps       = eps

    def compute(self, xyz):

        Ca = xyz[:,1] # [L,3]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), Ca[None,...].contiguous(), p=2) # [1,Lb,Lb]
        divide_by_r_0 = (dgram - self.d_0) / self.r_0
        numerator = torch.pow(divide_by_r_0,6)
        denominator = torch.pow(divide_by_r_0,12)

        ncontacts = (1 - numerator) / ((1 - denominator))


        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * ncontacts.sum()


class loop_contacts(Potential):

    def __init__(self, weight=0.01, ideal_dist_Ca=7, res1=27, res2=39):
        
        self.weight = weight
        self.ideal_dist_Ca = ideal_dist_Ca
        self.res1 = res1
        self.res2 = res2

    def compute(self, xyz):
        
        Ca_1 = xyz[self.res1,1]
        Ca_2 = xyz[self.res2,1]
        dgram = torch.cdist(Ca_1[None,None],Ca_2[None,None], p=2)
        dgram.squeeze(0).squeeze(0)
        dist = (dgram - self.ideal_dist_Ca)**2
        print('Ca dist:', dist)

        return -self.weight * dist


class olig_contacts(Potential):
    """
    Applies PV's num contacts potential within/between chains in symmetric oligomers 

    Author: DJ 
    """

    def __init__(self, 
                 contact_matrix, 
                 weight_intra=1, 
                 weight_inter=1,
                 r_0=8, d_0=2):
        """
        Parameters:
            chain_lengths (list, required): List of chain lengths, length is (Nchains)

            contact_matrix (torch.tensor/np.array, required): 
                square matrix of shape (Nchains,Nchains) whose (i,j) enry represents 
                attractive (1), repulsive (-1), or non-existent (0) contact potentials 
                between chains in the complex

            weight (int/float, optional): Scaling/weighting factor
        """
        self.contact_matrix = contact_matrix
        self.weight_intra = weight_intra 
        self.weight_inter = weight_inter 
        self.r_0 = r_0
        self.d_0 = d_0

        # check contact matrix only contains valid entries 
        assert all([i in [-1,0,1] for i in contact_matrix.flatten()]), 'Contact matrix must contain only 0, 1, or -1 in entries'
        # assert the matrix is square and symmetric 
        shape = contact_matrix.shape 
        assert len(shape) == 2 
        assert shape[0] == shape[1]
        for i in range(shape[0]):
            for j in range(shape[1]):
                assert contact_matrix[i,j] == contact_matrix[j,i]
        self.nchain=shape[0]

         
    def _get_idx(self,i,L):
        """
        Returns the zero-indexed indices of the residues in chain i
        """
        assert L%self.nchain == 0
        Lchain = L//self.nchain
        return i*Lchain + torch.arange(Lchain)


    def compute(self, xyz):
        """
        Iterate through the contact matrix, compute contact potentials between chains that need it,
        and negate contacts for any 
        """
        L = xyz.shape[0]

        all_contacts = 0
        start = 0
        for i in range(self.nchain):
            for j in range(self.nchain):
                # only compute for upper triangle, disregard zeros in contact matrix 
                if (i <= j) and (self.contact_matrix[i,j] != 0):

                    # get the indices for these two chains 
                    idx_i = self._get_idx(i,L)
                    idx_j = self._get_idx(j,L)

                    Ca_i = xyz[idx_i,1]  # slice out crds for this chain 
                    Ca_j = xyz[idx_j,1]  # slice out crds for that chain 
                    dgram           = torch.cdist(Ca_i[None,...].contiguous(), Ca_j[None,...].contiguous(), p=2) # [1,Lb,Lb]

                    divide_by_r_0   = (dgram - self.d_0) / self.r_0
                    numerator       = torch.pow(divide_by_r_0,6)
                    denominator     = torch.pow(divide_by_r_0,12)
                    ncontacts       = (1 - numerator) / (1 - denominator)

                    # weight, don't double count intra 
                    scalar = (i==j)*self.weight_intra/2 + (i!=j)*self.weight_inter

                    #                 contacts              attr/repuls          relative weights 
                    all_contacts += ncontacts.sum() * self.contact_matrix[i,j] * scalar 

        return all_contacts 

          
class hetero_olig(Potential):
    """
    Applies PV's num contacts potential within/between chains in hetero-oligomers

    Author: MH&NZR
    """

    def __init__(self,
                 chains,
                 interactions,
                 chain_lengths,
                 weight_intra=1,
                 weight_inter=1,
                 r_0=8, d_0=2):
        """
        Parameters:
            chains (list, required): List of chain labels, e.g., ['A', 'B', 'C', 'D']
            interactions (list, required): List of interaction specifications, e.g., ['ACD', 'ACD-B']
            weight_intra (int/float, optional): Scaling/weighting factor for intra-chain interactions
            weight_inter (int/float, optional): Scaling/weighting factor for inter-chain interactions
            r_0 (int/float, optional): Distance scaling parameter
            d_0 (int/float, optional): Distance offset parameter
        """
        chains = chains.replace('[', '').replace(']', '').split(';')
        interactions = interactions.replace('[', '').replace(']', '').split(';')
        self.chain_lengths = chain_lengths
        print('chain_lengths', self.chain_lengths)
        
        # Generate contact matrix based on chains and interactions
        self.contact_matrix = self.create_contact_matrix(chains, interactions)
        self.weight_intra = weight_intra
        self.weight_inter = weight_inter
        self.r_0 = r_0
        self.d_0 = d_0
        print('contact_matrix', self.contact_matrix)
        # Check that the contact matrix only contains valid entries
        assert all([i in [-1, 0, 1] for i in self.contact_matrix.flatten()]), 'Contact matrix must contain only 0, 1, or -1 in entries'
        # Assert the matrix is square
        shape = self.contact_matrix.shape
        assert len(shape) == 2
        assert shape[0] == shape[1]

        self.nchain = shape[0]

    def create_contact_matrix(self, chains, interactions):
        """
        Creates a contact matrix based on specified chain interactions.

        Parameters:
            chains (list): List of chain labels (e.g., ['A', 'B', 'C', 'D']).
            interactions (list): List of interaction specifications (e.g., ['ACD', 'ACD-B']).

        Returns:
            np.array: Contact matrix representing the interactions.
        """
        n_chains = len(chains)
        contact_matrix = np.zeros((n_chains, n_chains))

        # Create a mapping from chain labels to indices
        chain_to_index = {chain: i for i, chain in enumerate(chains)}

        for interaction in interactions:
            if '-' in interaction:
                # Handle interactions of the form "ACD-B"
                attract_part, repel_part = interaction.split('-')
                for i in range(len(attract_part) - 1):
                    chain1 = attract_part[i]
                    chain2 = attract_part[i + 1]
                    idx1 = chain_to_index[chain1]
                    idx2 = chain_to_index[chain2]
                    contact_matrix[idx1, idx2] = 1
                    contact_matrix[idx2, idx1] = 1

                for chain in repel_part:
                    for attract_chain in attract_part:
                        idx_attract = chain_to_index[attract_chain]
                        idx_repel = chain_to_index[chain]
                        contact_matrix[idx_attract, idx_repel] = -1
                        contact_matrix[idx_repel, idx_attract] = -1

            else:
                # Handle interactions of the form "ACD"
                for i in range(len(interaction) - 1):
                    chain1 = interaction[i]
                    chain2 = interaction[i + 1]
                    idx1 = chain_to_index[chain1]
                    idx2 = chain_to_index[chain2]
                    contact_matrix[idx1, idx2] = 1
                    contact_matrix[idx2, idx1] = 1

        return contact_matrix

    def _get_idx(self, i, L):
        """
        Returns the zero-indexed indices of the residues in chain i.
        """
        Lchain = self.chain_lengths[i]
        offset = sum(self.chain_lengths[:i])
        return offset + torch.arange(Lchain)

    def compute(self, xyz):
        """
        Iterate through the contact matrix, compute contact potentials between chains that need it.
        """
        L = xyz.shape[0]
        all_contacts = 0

        for i in range(self.nchain):
            for j in range(self.nchain):
                # Disregard zeros in contact matrix
                if self.contact_matrix[i, j] != 0:
                    # Get the indices for these two chains
                    idx_i = self._get_idx(i, L)
                    idx_j = self._get_idx(j, L)

                    Ca_i = xyz[idx_i, 1]  # slice out coordinates for this chain
                    Ca_j = xyz[idx_j, 1]  # slice out coordinates for that chain
                    dgram = torch.cdist(Ca_i[None, ...].contiguous(), Ca_j[None, ...].contiguous(), p=2)  # [1, Lb, Lb]

                    divide_by_r_0 = (dgram - self.d_0) / self.r_0
                    numerator = torch.pow(divide_by_r_0, 6)
                    denominator = torch.pow(divide_by_r_0, 12)
                    ncontacts = (1 - numerator) / (1 - denominator)

                    # Weight intra-chain interactions separately
                    scalar = self.weight_intra if i == j else self.weight_inter

                    # Calculate contacts, apply scalar, and sum up
                    all_contacts += ncontacts.sum() * self.contact_matrix[i, j] * scalar

        return all_contacts


def get_damped_lj(r_min, r_lin,p1=6,p2=12):
    
    y_at_r_lin = lj(r_lin, r_min, p1, p2)
    ydot_at_r_lin = lj_grad(r_lin, r_min,p1,p2)
    
    def inner(dgram):
        return (dgram < r_lin) * (ydot_at_r_lin * (dgram - r_lin) + y_at_r_lin) + (dgram >= r_lin) * lj(dgram, r_min, p1, p2)
    return inner

def lj(dgram, r_min,p1=6, p2=12):
    return 4 * ((r_min / (2**(1/p1) * dgram))**p2 - (r_min / (2**(1/p1) * dgram))**p1)

def lj_grad(dgram, r_min,p1=6,p2=12):
    return -p2 * r_min**p1*(r_min**p1-dgram**p1) / (dgram**(p2+1))

def mask_expand(mask, n=1):
    mask_out = mask.clone()
    assert mask.ndim == 1
    for i in torch.where(mask)[0]:
        for j in range(i-n, i+n+1):
            if j >= 0 and j < len(mask):
                mask_out[j] = True
    return mask_out

def contact_energy(dgram, d_0, r_0):
    divide_by_r_0 = (dgram - d_0) / r_0
    numerator = torch.pow(divide_by_r_0,6)
    denominator = torch.pow(divide_by_r_0,12)
    
    ncontacts = (1 - numerator) / ((1 - denominator)).float()
    return - ncontacts

def poly_repulse(dgram, r, slope, p=1):
    a = slope / (p * r**(p-1))

    return (dgram < r) * a * torch.abs(r - dgram)**p * slope

#def only_top_n(dgram


class substrate_contacts(Potential):
    '''
    Implicitly models a ligand with an attractive-repulsive potential.
    '''

    def __init__(self, weight=1, r_0=8, d_0=2, s=1, eps=1e-6, rep_r_0=5, rep_s=2, rep_r_min=1):

        self.r_0       = r_0
        self.weight    = weight
        self.d_0       = d_0
        self.eps       = eps
        
        # motif frame coordinates
        # NOTE: these probably need to be set after sample_init() call, because the motif sequence position in design must be known
        self.motif_frame = None # [4,3] xyz coordinates from 4 atoms of input motif
        self.motif_mapping = None # list of tuples giving positions of above atoms in design [(resi, atom_idx)]
        self.motif_substrate_atoms = None # xyz coordinates of substrate from input motif
        r_min = 2
        self.energies = []
        self.energies.append(lambda dgram: s * contact_energy(torch.min(dgram, dim=-1)[0], d_0, r_0))
        if rep_r_min:
            self.energies.append(lambda dgram: poly_repulse(torch.min(dgram, dim=-1)[0], rep_r_0, rep_s, p=1.5))
        else:
            self.energies.append(lambda dgram: poly_repulse(dgram, rep_r_0, rep_s, p=1.5))


    def compute(self, xyz):
        
        # First, get random set of atoms
        # This operates on self.xyz_motif, which is assigned to this class in the model runner (for horrible plumbing reasons)
        self._grab_motif_residues(self.xyz_motif)
        
        # for checking affine transformation is corect
        first_distance = torch.sqrt(torch.sqrt(torch.sum(torch.square(self.motif_substrate_atoms[0] - self.motif_frame[0]), dim=-1))) 

        # grab the coordinates of the corresponding atoms in the new frame using mapping
        res = torch.tensor([k[0] for k in self.motif_mapping])
        atoms = torch.tensor([k[1] for k in self.motif_mapping])
        new_frame = xyz[self.diffusion_mask][res,atoms,:]
        # calculate affine transformation matrix and translation vector b/w new frame and motif frame
        A, t = self._recover_affine(self.motif_frame, new_frame)
        # apply affine transformation to substrate atoms
        substrate_atoms = torch.mm(A, self.motif_substrate_atoms.transpose(0,1)).transpose(0,1) + t
        second_distance = torch.sqrt(torch.sqrt(torch.sum(torch.square(new_frame[0] - substrate_atoms[0]), dim=-1)))
        assert abs(first_distance - second_distance) < 0.01, "Alignment seems to be bad" 
        diffusion_mask = mask_expand(self.diffusion_mask, 1)
        Ca = xyz[~diffusion_mask, 1]

        #cdist needs a batch dimension - NRB
        dgram = torch.cdist(Ca[None,...].contiguous(), substrate_atoms.float()[None], p=2)[0] # [Lb,Lb]

        all_energies = []
        for i, energy_fn in enumerate(self.energies):
            energy = energy_fn(dgram)
            all_energies.append(energy.sum())
        return - self.weight * sum(all_energies)

        #Potential value is the average of both radii of gyration (is avg. the best way to do this?)
        return self.weight * ncontacts.sum()

    def _recover_affine(self,frame1, frame2):
        """
        Uses Simplex Affine Matrix (SAM) formula to recover affine transform between two sets of 4 xyz coordinates
        See: https://www.researchgate.net/publication/332410209_Beginner%27s_guide_to_mapping_simplexes_affinely

        Args: 
        frame1 - 4 coordinates from starting frame [4,3]
        frame2 - 4 coordinates from ending frame [4,3]
        
        Outputs:
        A - affine transformation matrix from frame1->frame2
        t - affine translation vector from frame1->frame2
        """

        l = len(frame1)
        # construct SAM denominator matrix
        B = torch.vstack([frame1.T, torch.ones(l)])
        D = 1.0 / torch.linalg.det(B) # SAM denominator

        M = torch.zeros((3,4), dtype=torch.float64)
        for i, R in enumerate(frame2.T):
            for j in range(l):
                num = torch.vstack([R, B])
                # make SAM numerator matrix
                num = torch.cat((num[:j+1],num[j+2:])) # make numerator matrix
                # calculate SAM entry
                M[i][j] = (-1)**j * D * torch.linalg.det(num)

        A, t = torch.hsplit(M, [l-1])
        t = t.transpose(0,1)
        return A, t

    def _grab_motif_residues(self, xyz) -> None:
        """
        Grabs 4 atoms in the motif.
        Currently random subset of Ca atoms if the motif is >= 4 residues, or else 4 random atoms from a single residue
        """
        idx = torch.arange(self.diffusion_mask.shape[0])
        idx = idx[self.diffusion_mask].float()
        if torch.sum(self.diffusion_mask) >= 4:
            rand_idx = torch.multinomial(idx, 4).long()
            # get Ca atoms
            self.motif_frame = xyz[rand_idx, 1]
            self.motif_mapping = [(i,1) for i in rand_idx]
        else:
            rand_idx = torch.multinomial(idx, 1).long()
            self.motif_frame = xyz[rand_idx[0],:4]
            self.motif_mapping = [(rand_idx, i) for i in range(4)]


class binder_RMSD(Potential):
    """ Potential to reduce RMSD of binder towards specific shape. """

    def __init__(self, weight=1, ref_pdb='/nas/longleaf/home/bkuhlman/snap/fold_conditioning/backbone_305.pdb', basis=0, squared=True, roi= list(range(0,180))):
        self.weight = weight
        if basis == 0:
            self.basis = 'ca'
        elif basis == 1:
            self.basis = 'bb'
        else:
            raise ValueError(f'Unrecognized basis mode: {basis}. Use either 0 for "ca" or 1 for "bb"')
        self.squared = squared
        self.roi = roi

        # Parse the reference pdb and get appropriate coordinates.
        from rfdiffusion.inference.utils import parse_pdb
        ref_pdb = parse_pdb(ref_pdb)
        self.ref_xyz = torch.from_numpy(self.get_basis_xyz(ref_pdb['xyz'], self.basis, roi))
        self.ref_xyz = self.ref_xyz - self.ref_xyz.mean(0)        

    def compute(self, xyz):
        # Grab only the binder residues.
        binder_xyz = xyz

        # Get the atoms used for RMSD calculation.
        binder_xyz = self.get_basis_xyz(binder_xyz, self.basis, self.roi) # [num_atom, 3]

        # Compute alignment.
        R = self._calc_kabsch_R(binder_xyz, self.ref_xyz)

        # Compute MSD after centering binder_xyz and rotating ref_xyz.
        binder_xyz = binder_xyz - binder_xyz.mean(0)
        ref_xyz_aligned = self.ref_xyz @ R
        
        msd = torch.mean(torch.square(binder_xyz - ref_xyz_aligned).sum(-1))
        print('RMSD:', torch.sqrt(msd))
        
        # Return either MSD or RMSD
        if self.squared:
            return -1 * self.weight * msd
        else:
            return -1 * self.weight * torch.sqrt(msd)

    def get_basis_xyz(self, xyz, basis, roi):
        """ Gets the coordinates to be used in the RMSD calculation.
        
        Args:
            xyz (torch.tensor, size: [L,27,3] or [L,14,3]): Complete set of coordinates.
            basis (str, ["ca", "bb"]): Either 'ca' or 'bb', indicating to use only CA atom or all backbone atoms.
            roi (list(int)): The residues of interest as a list of the index (0-based) of the residues.

        Returns:
            basis_xyz (torch.tensor, size: [num_atom, 3]): Atomic coordinates of the residues of interest.
        """

        # Get the roi coordinates.
        roi_xyz = xyz[roi]

        # Get the appropriate atomic coordinates.
        if basis == 'ca':
            roi_xyz = roi_xyz[:, 1]
        elif basis == 'bb':
            roi_xyz = roi_xyz[:, :3]
        else:
            raise ValueError(f"Unrecognized basis {basis}. Use either 'ca' or 'bb'.")
        
        # Reshape and return
        return roi_xyz.reshape(-1, 3)

    @torch.no_grad()
    def _calc_kabsch_R(self, xyz1, xyz2):
        """
        Calculates optimal rotation matrix for alignment, following the Kabsch algorithm.
        """
        # Center coordinates on origin.
        xyz1 = xyz1 - xyz1.mean(0)
        xyz2 = xyz2 - xyz2.mean(0)

        # Computate the covariance matrix
        C = xyz2.T @ xyz1

        # Compute optimal rotation matrix using SVD
        try:
            U, S, Vh = torch.linalg.svd(C)
        except:
            print('Hit SVD exception.')
            U, S, Vh = torch.linalg.svd(C + 1e-2 * C.mean() * torch.rand_like(C))

        # Get the sign to ensure right handedness
        d = torch.ones([3,3])
        d[:,-1] = torch.sign(torch.linalg.det(U) * torch.linalg.det(Vh))

        # Rotation matrix R
        R = (d * U) @ Vh

        return R


class res_pair_constraints(Potential):
    """ Potential to force residue pair distances to match a reference pdb """

    def __init__(self, weight=10, ref_pdb='/nas/longleaf/home/bkuhlman/snap/fold_conditioning/backbone_305.pdb', basis=0, squared=True, roi= list(range(0,180))):
        #self.binderlen = binderlen
        self.weight = weight
        if basis == 0:
            self.basis = 'ca'
        elif basis == 1:
            self.basis = 'bb'
        else:
            raise ValueError(f'Unrecognized basis mode: {basis}. Use either 0 for "ca" or 1 for "bb"')
        self.squared = squared
        self.roi = roi

        # Parse the reference pdb and get appropriate coordinates.
        from rfdiffusion.inference.utils import parse_pdb
        ref_pdb = parse_pdb(ref_pdb)
        self.ref_xyz = torch.from_numpy(self.get_basis_xyz(ref_pdb['xyz'], self.basis, roi))
        self.dgram_ref = torch.cdist(self.ref_xyz[None, ...].contiguous(),self.ref_xyz[None, ...].contiguous(), p=2)
        
    def compute(self, xyz):
        # Grab only the binder residues.
        #binder_xyz = xyz[:self.binderlen] # [L_binder, 27, 3]
        binder_xyz = xyz
        # Get the atoms used for the distance calculations.
        binder_xyz = self.get_basis_xyz(binder_xyz, self.basis, self.roi) # [num_atom, 3]
        dgram_model = torch.cdist(binder_xyz[None, ...].contiguous(),binder_xyz[None, ...].contiguous(), p=2)
        #dgram_model = torch.cdist(binder_xyz,binder_xyz)
        dgram_diff = torch.abs(self.dgram_ref - dgram_model)
        dgram_diff = dgram_diff**1.5

        
        #print('dgram_model distance',dgram_model[0,0,1])
        #print('dgram_model',dgram_model)
        #print('dgram_model_shape',dgram_model.shape)
        #print('dgram_model other side distance',dgram_model[0,1,0])
        #print('dgram_model diagonal',dgram_model[0,0,0])
        #print('dgram_ref distance',self.dgram_ref[0,0,1])
        #print('dgram_diff distance',dgram_diff[0,0,1])

        mask = torch.ones_like(dgram_diff).triu(diagonal=1).bool()
        score = dgram_diff[mask].mean()
        
        #score = (dgram_diff ** 1.5).sum()
        #score = dgram_diff.sum()
        print('distance score:',score)

        return -1 * self.weight * score

    def get_basis_xyz(self, xyz, basis, roi):
        """ Gets the coordinates to be used in the distance calculations.
        
        Args:
            xyz (torch.tensor, size: [L,27,3] or [L,14,3]): Complete set of coordinates.
            basis (str, ["ca", "bb"]): Either 'ca' or 'bb', indicating to use only CA atom or all backbone atoms.
            roi (list(int)): The residues of interest as a list of the index (0-based) of the residues.

        Returns:
            basis_xyz (torch.tensor, size: [num_atom, 3]): Atomic coordinates of the residues of interest.
        """

        # Get the roi coordinates.
        roi_xyz = xyz[roi]

        # Get the appropriate atomic coordinates.
        if basis == 'ca':
            roi_xyz = roi_xyz[:, 1]
        elif basis == 'bb':
            roi_xyz = roi_xyz[:, :2]
        else:
            raise ValueError(f"Unrecognized basis {basis}. Use either 'ca' or 'bb'.")
        
        # Reshape and return
        return roi_xyz.reshape(-1, 3)


# Dictionary of types of potentials indexed by name of potential. Used by PotentialManager.
# If you implement a new potential you must add it to this dictionary for it to be used by
# the PotentialManager
implemented_potentials = { 'monomer_ROG':          monomer_ROG,
                           'binder_ROG':           binder_ROG,
                           'dimer_ROG':            dimer_ROG,
                           'binder_ncontacts':     binder_ncontacts,
                           'interface_ncontacts':  interface_ncontacts,
                           'monomer_contacts':     monomer_contacts,
                           'olig_contacts':        olig_contacts,
                           'substrate_contacts':    substrate_contacts,
                           'binder_RMSD':           binder_RMSD,
                           'loop_contacts':         loop_contacts,
                           'res_pair_constraints':  res_pair_constraints,
                           'hetero_olig':           hetero_olig}

require_binderlen      = { 'binder_ROG',
                           'binder_distance_ReLU',
                           'binder_any_ReLU',
                           'dimer_ROG',
                           'binder_ncontacts',
                           'interface_ncontacts'}

