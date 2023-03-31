#!/bin/bash
# In this example, we scaffold an active site from a retroaldolase, and implicitly model the substrate with an external potential
# We specify the output path and the input pdb (a native retroaldolase)
# We provide a contig describing the protein we want to build:
#   - 10-100 residues (randomly sampled), then
#   - residue 1083 on the A chain of the pdb file, then
#   - 10-100 residues (randomly sampled), then
#   - residue 1051 on the A chain of the pdb file, then
#   - 10-100 residues (randomly sampled), then
#   - residue 1180 of the A chain of the pdb file, then
#   - 10-100 residues (randomly sampled)
# We then set the potential guide scale to 1, and use the "substrate_contact" potential, with parameters described in the manuscript
# We specify the identity of the substrate molecule (from which to apply the potential)
# We then specify the we want to use the model weights fine tuned to scaffold small motifs (i.e. three single residues)

../scripts/run_inference.py inference.output_prefix=example_outputs/design_enzyme inference.input_pdb=input_pdbs/5an7.pdb 'contigmap.contigs=[10-100/A1083-1083/10-100/A1051-1051/10-100/A1180-1180/10-100]' potentials.guide_scale=1 'potentials.guiding_potentials=["type:substrate_contacts,s:1,r_0:8,rep_r_0:5.0,rep_s:2,rep_r_min:1"]' potentials.substrate=LLK inference.ckpt_override_path=../models/ActiveSite_ckpt.pt
