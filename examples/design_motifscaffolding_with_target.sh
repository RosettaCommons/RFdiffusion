#!/bin/bash
# Here, we're running one of the motif-scaffolding benchmark examples, in the presence of a target
# Specifically, we're scaffolding the Mdm2-interacting peptide from p53
# We specify the output path and input pdb (the p53-Mdm2 complex)
# We specify the protein we want to build, with the contig input:
#   - the Mdm2 target protein (residues A25-109), with a chain break.
#   - 0-70 residues (randomly sampled)
#   - residues 17-29 (inclusive) on the B chain of the input (the p53 helix)
#   - 0-70 residues (randomly sampled)
# We also constrain the total length of the diffused chain to be within 70 and 120 residues
# We generate 10 designs
# As in the paper (at least for some of the designs we tested), we use the complex-finetuned model

python ../scripts/run_inference.py inference.output_prefix=example_outputs/design_motifscaffolding_with_target inference.input_pdb=input_pdbs/1YCR.pdb 'contigmap.contigs=[A25-109/0 0-70/B17-29/0-70]' contigmap.length=70-120 inference.num_designs=10 inference.ckpt_override_path=../models/Complex_base_ckpt.pt
