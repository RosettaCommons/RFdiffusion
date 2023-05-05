#!/bin/bash
# Here, we're making some unconditional designs, and also using the 'monomer_contacts' auxillary potential function
# We specify the path for the outputs
# We tell RFdiffusion that designs should be 100-200 residues in length (randomly sampled each design)
# We generate 10 such designs

../scripts/run_inference.py inference.output_prefix=example_outputs/design_unconditional_w_contact_potential 'contigmap.contigs=[100-200]' inference.num_designs=10 'potentials.guiding_potentials=["type:monomer_contacts,weight:0.05"]'
