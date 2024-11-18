import os
import numpy as np
import torch

from rfdiffusion.contigs import ContigMap
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

sampler = iu.sampler_selector(conf)
x_init, seq_init = sampler.sample_init()

x_init.shape    # torch.Size([177, 14, 3])
seq_init.shape  # torch.Size([177, 22])

# This processes the PDB, need to figure out what xyz_het is
# xyz_het is the xyz coordinates of the heteroatoms
# TODO: Test with 7yke. It should get everything, maybe check the MG atom
# NOTE: This will likely need to be updated to take in a dual substrate
target_feats = iu.process_target(pdb_path, parse_hetatom=True, center=False)
target_struct = iu.parse_pdb(pdb_path, parse_hetatom=True)

# Processing through line
het_names = np.array([i['name'].strip() for i in target_feats['info_het']])
xyz_het = target_feats['xyz_het'][het_names == 'LLK']

# Diffusion mask
target_feats = iu.process_target(pdb_path, parse_hetatom=True, center=False)
# contig_map = ContigMap(target_feats, conf.contigmap)


# cleft w substrate
# ../scripts/run_inference.py inference.output_prefix=example_outputs/design_enzyme inference.input_pdb=input_pdbs/5an7.pdb 'contigmap.contigs=[10-100/A1083-1083/10-100/A1051-1051/10-100/A1180-1180/10-100]' potentials.guide_scale=1 'potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1"]' potentials.substrate=LLK potentials.cleft_w_centered_substrate_potential=True inference.ckpt_override_path=../models/ActiveSite_ckpt.pt