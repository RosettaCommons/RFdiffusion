#!/bin/bash

prefix=./outputs/diffused_binder_cyclic2

# Note that the indices in this pdb file have been 
# shifted by +2 in chain A relative to pdbID 7zkr.
pdb='./input_pdbs/7zkr_GABARAP.pdb'

num_designs=10
script="../scripts/run_inference.py"
$script --config-name base \
inference.output_prefix=$prefix \
inference.num_designs=$num_designs \
'contigmap.contigs=[12-18 A3-117/0]' \
inference.input_pdb=$pdb \
inference.cyclic=True \
diffuser.T=50 \
inference.cyc_chains='a' \
ppi.hotspot_res=[\'A51\',\'A52\',\'A50\',\'A48\',\'A62\',\'A65\'] \
