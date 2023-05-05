#!/bin/bash
# Here, we're making some unconditional designs
# We are also applying a Radius Of Gyration potential which encourages the design to be more compact
# We specify the path for the outputs
# We tell RFdiffusion that designs should be 100-200 residues in length (randomly sampled each design)
# We generate 10 such designs
# We use the monomer_ROG potential, with guide scale 2 and quadratic decay
# Note that this potential is probably not necessary in this kind of case, but is provided as an example

../scripts/run_inference.py inference.output_prefix=example_outputs/design_monomer_ROG_unconditional 'contigmap.contigs=[100-200]' inference.num_designs=10 'potentials.guiding_potentials=["type:monomer_ROG,weight:1,min_dist:5"]' potentials.guide_scale=2 potentials.guide_decay="quadratic"
