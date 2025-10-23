#!/bin/bash
# Here, we're generating a C4-symmetric oligomer scaffolding a symmetric nickel-binding site
# This is described in the manuscript
# First, we specify we want a C4 symmetric oligomer, and 15 designs, to a specific output path
# We then apply potentials to promote compact structures
# This external potential promotes contacts both within (with a relative weight of 1) and between chains (relative weight 0.06)
# We specify that we want to apply these potentials to all chains, with a guide scale of 2.0 (a sensible starting point)
# We decay this potential with quadratic form, so that it is applied more strongly initially 
# We then describe the protein we want to make, with a contig input
# The input pdb contains four copies of the motif, arranged in a C4 symmetric arrangement around the symmetry (Z) axis
# These four copies are residues A2-4, A7-9, A12-14, and A17-19
# The contig therefore defines the following protein:
#   - 50 residues, then the first motif (A2-4), followed by 50 residues
#   - a chainbreak (this is not strictly necessary in the symmetry case, but we leave it here for clarity
#   - Three more "chains" of this. Note that the order of the motifs doesn't matter, but their symmetric arrangement around the symmetry axis does
# Based on some preliminary tests, we found that a different model checkpoint performs better on this task. We therefore override the default with these weights

ckpt='../models/Base_epoch8_ckpt.pt'

python ../scripts/run_inference.py inference.symmetry="C4" inference.num_designs=15 inference.output_prefix=example_outputs/design_nickel 'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.06"]' potentials.olig_intra_all=True potentials.olig_inter_all=True potentials.guide_scale=2 potentials.guide_decay="quadratic" inference.input_pdb=input_pdbs/nickel_symmetric_motif.pdb 'contigmap.contigs=[50/A2-4/50/0 50/A7-9/50/0 50/A12-14/50/0 50/A17-19/50/0]' inference.ckpt_override_path=$ckpt
