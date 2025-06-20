#!/bin/bash 

prefix=./outputs/uncond_cycpep
# Note that the indices in this pdb file have been 
# shifted by +2 in chain A relative to pdbID 7zkr.
pdb='./input_pdbs/7zkr_GABARAP.pdb'

num_designs=10
script="../scripts/run_inference.py"
$script --config-name base \
inference.output_prefix=$prefix \
inference.num_designs=$num_designs \
'contigmap.contigs=[12-18]' \
inference.input_pdb=$pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a'
