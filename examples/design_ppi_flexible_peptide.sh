#!/bin/bash
# Here, we're designing binders to the glp1 helical peptide, without specifying the topology of the binder a priori, and without specifying the structure of the peptide (we know peptides can be flexible).
# We first provide the output path and input pdb of the target protein (5uul)
# We then describe the protein we want with the contig input:
#   - residues 10-35 of the B chain of the target protein
#   - a chainbreak (as we don't want the binder fused to the target!)
#   - A 70-100 residue binder to be diffused (the exact length is sampled each iteration of diffusion)
# We tell diffusion to target two specific residues on the target, specifically residues 28 and 29 of the B chain
# We make 10 designs
# We mask (diffuse) the structure of the peptide using the inpaint_str flag. This has the effect of having RFdiffusion simultaneously design a binder and predict the structure of the peptide within the complex.

../scripts/run_inference.py inference.output_prefix=example_outputs/design_ppi_flexible_peptide inference.input_pdb=input_pdbs/3IOL.pdb 'contigmap.contigs=[B10-35/0 70-100]' 'ppi.hotspot_res=[B28,B29]' inference.num_designs=10 'contigmap.inpaint_str=[B10-35]'
