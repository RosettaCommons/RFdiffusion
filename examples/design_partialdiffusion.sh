#!/bin/bash
# Here, we're partially noising and denoising from a starting structure, to generate diversity around this fold
# We first specify the output path, and the path to the 79aa input pdb, 2KL8
# We then specify the contig corresponding to what we want to noise.
# In this case, we want to noise the whole structure, so we specify [79-79], i.e. noise all 79 residues
# We could instead just noise the latter 59 residues, by specifying [A1-19/59]
# But, in either case, the contig length must sum to the length of the input pdb file
# We generate 10 designs, and noise and denoise 10 steps (20% of the full trajectory)

../scripts/run_inference.py inference.output_prefix=example_outputs/design_partialdiffusion inference.input_pdb=input_pdbs/2KL8.pdb 'contigmap.contigs=[79-79]' inference.num_designs=10 diffuser.partial_T=10
