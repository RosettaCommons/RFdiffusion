#!/bin/bash
# Here, we're making a binder, and we're specifying the course-grained topology of the binder, from a set of processed inputs
# We specify the path to the target protein, in this case the insulin receptor, along with an output path
# We tell RFdiffusion that we want to do "scaffoldguided" diffusion (i.e. we want to specify the fold of the protein)
# We tell RFdiffusion where on the (cropped) input protein we want to bind, in this case to residues 59, 83 and 91 on the A chain
# We tell RFdiffusion that we're wanting to make a binder to a target, and provide the secondary structure and block adjacency input for these. This may not be necessary
# We then provide a path to a directory of different scaffolds (we've provided some for you to use, from Cao et al., 2022)
# We generate 10 designs, and reduce the noise added during inference to 0 (which improves the quality of designs)

../scripts/run_inference.py scaffoldguided.target_path=input_pdbs/insulin_target.pdb inference.output_prefix=example_outputs/design_ppi_scaffolded scaffoldguided.scaffoldguided=True 'ppi.hotspot_res=[A59,A83,A91]' scaffoldguided.target_pdb=True scaffoldguided.target_ss=target_folds/insulin_target_ss.pt scaffoldguided.target_adj=target_folds/insulin_target_adj.pt scaffoldguided.scaffold_dir=./ppi_scaffolds/ inference.num_designs=10 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
