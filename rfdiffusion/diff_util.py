import torch 
import numpy as np 
import random 

from rfdiffusion.chemical import INIT_CRDS
from icecream import ic 


def th_min_angle(start, end, radians=False):
    """
    Finds the angle you would add to <start> in order to get to <end>
    on the shortest path.
    """
    a,b,c = (np.pi, 2*np.pi, 3*np.pi) if radians else (180, 360, 540)
    shortest_angle = ((((end - start) % b) + c) % b) - a
    return shortest_angle


def th_interpolate_angles(start, end, T, n_diffuse,mindiff=None, radians=True):
    """
    
    """
    # find the minimum angle to add to get from start to end
    angle_diffs = th_min_angle(start, end, radians=radians)
    if mindiff is not None:
        assert torch.sum(mindiff.flatten()-angle_diffs) == 0.
    if n_diffuse is None:
        # default is to diffuse for max steps 
        n_diffuse = torch.full((len(angle_diffs)), T)


    interps = []
    for i,diff in enumerate(angle_diffs):
        N = int(n_diffuse[i])
        actual_interp = torch.linspace(start[i], start[i]+diff, N)
        whole_interp = torch.full((T,), float(start[i]+diff))
        temp=torch.clone(whole_interp)
        whole_interp[:N] = actual_interp
        
        interps.append(whole_interp)

    return torch.stack(interps, dim=0)


def th_interpolate_angle_single(start, end, step, T, mindiff=None, radians=True):
    """
    
    """
    # find the minimum angle to add to get from start to end
    angle_diffs = th_min_angle(start, end, radians=radians)
    if mindiff is not None:
        assert torch.sum(mindiff.flatten()-angle_diffs) == 0.

    # linearly interpolate between x = [0, T-1], y = [start, start + diff]
    x_range = T-1
    interps = step / x_range * angle_diffs + start

    return interps


def get_aa_schedule(T, L, nsteps=100):
    """
    Returns the steps t when each amino acid should be decoded, 
    as well as how many steps that amino acids chi angles will be diffused 
    
    Parameters:
        T (int, required): Total number of steps we are decoding the sequence over 
        
        L (int, required): Length of protein sequence 
        
        nsteps (int, optional): Number of steps over the course of which to decode the amino acids 

    Returns: three items 
        decode_times (list): List of times t when the positions in <decode_order> should be decoded 

        decode_order (list): List of lists, each element containing which positions are going to be decoded at 
                             the corresponding time in <decode_times> 

        idx2diffusion_steps (np.array): Array mapping the index of the residue to how many diffusion steps it will require 

    """
    # nsteps can't be more than T or more than length of protein
    if (nsteps > T) or (nsteps > L):
        nsteps = min([T,L])

    
    decode_order = [[a] for a in range(L)]
    random.shuffle(decode_order)
    
    while len(decode_order) > nsteps:
        # pop an element and then add those positions randomly to some other step
        tmp_seqpos = decode_order.pop()
        decode_order[random.randint(0,len(decode_order)-1)] += tmp_seqpos
        random.shuffle(decode_order)
    
    decode_times = np.arange(nsteps)+1
    
    # now given decode times, calculate number of diffusion steps each position gets
    aa_masks = np.full((200,L), False)
    
    idx2diffusion_steps = np.full((L,),float(np.nan))
    for i,t in enumerate(decode_times):
        decode_pos = decode_order[i]    # positions to be decoded at this step 
        
        for j,pos in enumerate(decode_pos):
            # calculate number of diffusion steps this residue gets 
            idx2diffusion_steps[pos] = int(t)
            aa_masks[t,pos] = True
    
    aa_masks = np.cumsum(aa_masks, axis=0)
    
    return decode_times, decode_order, idx2diffusion_steps, ~(aa_masks.astype(bool))

####################
### for SecStruc ###
####################

def ss_to_tensor(ss_dict):
    """
    Function to convert ss files to indexed tensors
    0 = Helix
    1 = Strand
    2 = Loop
    3 = Mask/unknown
    4 = idx for pdb
    """
    ss_conv = {'H':0,'E':1,'L':2}
    ss_int = np.array([int(ss_conv[i]) for i in ss_dict['ss']])
    return ss_int

def mask_ss(ss, min_mask = 0, max_mask = 0.75):
    """
    Function to take ss array, find the junctions, and randomly mask these until a random proportion (up to 75%) is masked
    Input: numpy array of ss (H=0,E=1,L=2,mask=3)
    output: tensor with some proportion of junctions masked
    """
    mask_prop = random.uniform(min_mask, max_mask)
    transitions = np.where(ss[:-1] - ss[1:] != 0)[0] #gets last index of each block of ss
    counter = 0
    #TODO think about masking whole ss elements
    while len(ss[ss == 3])/len(ss) < mask_prop and counter < 100: #very hacky - do better
        try:
            width = random.randint(1,9)
            start = random.choice(transitions)
            offset = random.randint(-8,1)
            ss[start+offset:start+offset+width] = 3
            counter += 1
        except:
            counter += 1
    ss = torch.tensor(ss)
    mask = torch.where(ss == 3, True, False)
    ss = torch.nn.functional.one_hot(ss, num_classes=4)
    return ss, mask

def construct_block_adj_matrix( sstruct, xyz, nan_mask, cutoff=6, include_loops=False ):
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

    # Remove nans at this stage, as ss doesn't consider nans
    xyz_nonan = xyz[nan_mask]
    L = xyz_nonan.shape[0]
    assert L == sstruct.shape[0]
    # three anchor atoms
    N  = xyz_nonan[:,0]
    Ca = xyz_nonan[:,1]
    C  = xyz_nonan[:,2]

    # recreate Cb given N,Ca,C
    Cb = generate_Cbeta(N,Ca,C)

    dist = get_pair_dist(Cb,Cb) # [L,L]
    dist[torch.isnan(dist)] = 999.9
    assert torch.sum(torch.isnan(dist)) == 0
    dist += 999.9*torch.eye(L,device=xyz.device)

    # Now we have dist matrix and sstruct specification, turn this into a block adjacency matrix

    # First: Construct a list of segments and the index at which they begin and end
    in_segment = True
    segments = []

    begin = -1
    end = -1
    # need to expand ss out to size L


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

    # Second: Using segments and dgram, determine adjacent blocks
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
