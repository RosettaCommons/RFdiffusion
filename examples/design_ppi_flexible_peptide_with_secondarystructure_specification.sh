#!/bin/bash
# Here, we're designing binders to a tau peptide. We want to specify that this peptide should adopt a helical conformation.
# We first provide the output path and input pdb of the target protein (tau_peptide.pdb)
# We then describe the protein we want with the contig input:
#   - A 70-100 residue binder to be diffused (the exact length is sampled each iteration of diffusion)
#   - a chainbreak (as we don't want the binder fused to the target!)
#   - residues 165-178 of the B chain of the target protein
# We make 10 designs
# We mask (diffuse) the structure of the peptide using the inpaint_str flag. This has the effect of having RFdiffusion simultaneously design a binder and predict the structure of the peptide within the complex.
# We then specify that we want this to adopt a helical conformation
# If you wanted to specify that it should adopt a strand conformation, you would specify `contigmap.inpaint_str_strand`

../scripts/run_inference.py inference.output_prefix=example_outputs/design_ppi_flexible_peptide_with_secondarystructure inference.input_pdb=input_pdbs/tau_peptide.pdb 'contigmap.contigs=[70-100/0 B165-178]' inference.num_designs=10 'contigmap.inpaint_str=[B165-178]' scaffoldguided.scaffoldguided=True 'contigmap.inpaint_str_helix=[B165-178]'
