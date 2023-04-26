#!/bin/bash
# This example is a copy of design_partialdiffusion_withseq.sh
# The only difference is we specify multiple disjoint sequences to hold constant.
# In this case, we provide the residues of the beginning and end of the peptide, instead of the whole sequence.
# We can provide comma-separated ranges to specify this: provide_seq=[172-177,200-205]
# Note the ranges do not necessarily need to lie on the same chain as in this example.
# However, positions are 0-indexed over the whole sequence--not per-chain-- so care must be taken when providing ranges to provide_seq.

../scripts/run_inference.py inference.output_prefix=example_outputs/design_partialdiffusion_peptidewithmultiplesequence inference.input_pdb=input_pdbs/peptide_complex_ideal_helix.pdb 'contigmap.contigs=["172-172/0 34-34"]' diffuser.partial_T=10 inference.num_designs=10 'contigmap.provide_seq=[172-177,200-205]'
