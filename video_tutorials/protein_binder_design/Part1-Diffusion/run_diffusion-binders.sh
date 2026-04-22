#!/bin/bash

script='/home/dlopezma/RFdiffusion/scripts/run_inference.py'


python $script \
inference.output_prefix=outputs-test/binder_ \
inference.input_pdb=channel-toxin.pdb 'contigmap.contigs=[A1-110/0 B111-220/0 C221-330/0 D331-440/0 30-50/H462-463/20-50]' \
'ppi.hotspot_res=[A68,A70,B178,B180,C288,C290,D398,D400]' \
'potentials.guiding_potentials=["type:binder_ROG,min_dist:5,weight:10"]' potentials.guide_scale=10 potentials.guide_decay="quadratic" \
denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 \
inference.num_designs=200