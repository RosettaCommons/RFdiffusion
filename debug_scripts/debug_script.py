import os
import numpy as np
import torch

from rfdiffusion.inference import utils as iu
from rfdiffusion.inference.model_runners import Sampler
from omegaconf import OmegaConf
from hydra import compose, initialize


WORKDIR = '/work/gllab/noyes.m/RFdiffusion'

pdb_path = f"{WORKDIR}/examples/input_pdbs/5an7.pdb"

# My current conf script
with initialize(config_path="config/inference"):
    conf = compose(config_name="base")

conf.inference.input_pdb = pdb_path
conf.inference.ckpt_override_path = f"{WORKDIR}/models/ActiveSite_ckpt.pt"
conf.contigmap.contigs=['10-100/A1083-1083/10-100/A1051-1051/10-100/A1180-1180/10-100']
conf.potentials.guide_scale=1 
conf.potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1"]

conf.potentials.substrate = 'LLK'

# Declare sampler
sampler = iu.sampler_selector(conf)
x_init, seq_init = sampler.sample_init()

x_init.shape    # torch.Size([166, 14, 3])
seq_init.shape  # torch.Size([166, 22])

# Just process target
target_feats = iu.process_target(pdb_path, parse_hetatom=True, center=False)

# Deeper into process target
target_struct = iu.parse_pdb(pdb_path, parse_hetatom=True)

# Deeper into parse_pdb

# How to get coords for each CA
target_struct['xyz'][0][1]  # array([-14.383,   6.009, -20.015], dtype=float32)
target_struct['xyz'][1][1]  # array([-13.563,   4.09 , -16.86 ], dtype=float32)
target_struct['xyz'][2][1]  # array([-15.97 ,   1.969, -14.875], dtype=float32)

ca_array = target_struct['xyz'][:, 1, :]
center = np.mean(ca_array, axis=0)
cmax = np.max(ca_array, axis=0)
cmin = np.min(ca_array, axis=0)
cdiff = cmax - cmin  # array([36.774998, 44.698997, 37.512997], dtype=float32)

# Walkthrough of the monomer_ROG potential derived from 5an7 structure
xyz_np = target_struct['xyz']
xyz = torch.from_numpy(xyz_np)
xyz.requires_grad = True

min_dist = 15
weight = 1

Ca = xyz[:,1] # [L,3]
centroid = torch.mean(Ca, dim=0, keepdim=True) # [1,3]
dgram = torch.cdist(Ca[None,...].contiguous(), centroid[None,...].contiguous(), p=2) # [1,L,1,3]
dgram = torch.maximum(min_dist * torch.ones_like(dgram.squeeze(0)), dgram.squeeze(0)) # [L,1,3]
rad_of_gyration = torch.sqrt( torch.sum(torch.square(dgram)) / Ca.shape[0] ) # [1]
weighted_rog = -1 * weight * rad_of_gyration

# Compute all potentials
potential_list = [weighted_rog]
potential_stack = torch.stack(potential_list, dim=0)
potential = torch.sum(potential_stack, dim=0)

# get_potential_gradients
potential.backward()
Ca_grads = xyz.grad[:, 1, :]

# get_next_pose

# For tonight, come up with hypothesis for active site prediction
#  - From centroid, create potential w/radius
#  - Figure out how to design a cleft
