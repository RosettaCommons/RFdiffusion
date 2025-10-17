#!/bin/bash
# Here, we're asking diffusion to make a TIM barrel, by providing course-grained specification of the fold
# We specify the output path
# We tell RFdiffusion that we want to do scaffoldguided design, and that we are not making a binder to a target (just a monomer)
# We provide a path to a directory of TIM barrel scaffolds, generated with the helper script
# We generate 10 designs, with a reduced noise scale during inference of 0.5
# We sample additional length to increase diversity of the outputs. Specifically, we mask the loops and insert 0-5 residues (randomly sampled per-loop) into each loop
# We add 0-5 residues (randomly sampled) to the N and the C-terminus
# This will allow the generation of diverse TIM barrels with slightly different length helices and strands

../scripts/run_inference.py inference.output_prefix=example_outputs/design_tim_barrel scaffoldguided.scaffoldguided=True scaffoldguided.target_pdb=False scaffoldguided.scaffold_dir=tim_barrel_scaffold/ inference.num_designs=10 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 scaffoldguided.sampled_insertion=0-5 scaffoldguided.sampled_N=0-5 scaffoldguided.sampled_C=0-5
