#!/bin/bash 

# Note that in the example below the indices in the 
# input_pdbs/7zkr_GABARAP.pdb file have been shifted
# by +2 in chain A relative to pdbID 7zkr.

../scripts/run_inference.py \
--config-name base \
inference.output_prefix=example_outputs/uncond_cycpep \
inference.num_designs=2 \
'contigmap.contigs=[12-18]' \
inference.input_pdb=input_pdbs/7zkr_GABARAP.pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a'
