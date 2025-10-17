#!/bin/bash
# Here, we're making some unconditional designs
# We specify the path for the outputs
# We tell RFdiffusion that designs should be 100-200 residues in length (randomly sampled each design)
# We generate 10 such designs

../scripts/run_inference.py inference.output_prefix=example_outputs/design_unconditional 'contigmap.contigs=[100-200]' inference.num_designs=10
