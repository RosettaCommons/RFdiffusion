#!/bin/bash
# In this example, we generate a C6 symmetric oligomer
# We use the symmetry config, and specify we want a C6 symmetric oligomer, with 10 designs generated
# We specify the output prefix, and also the potential we want to apply
# This external potential promotes contacts both within (with a relative weight of 1) and between chains (relative weight 0.1)
# We specify that we want to apply these potentials to all chains, with a guide scale of 2.0 (a sensible starting point)
# We decay this potential with quadratic form, so that it is applied more strongly initially
# We specify a total length of 480aa, so each chain is 80 residues long

python ../scripts/run_inference.py --config-name=symmetry inference.symmetry="C6" inference.num_designs=10 inference.output_prefix="example_outputs/C6_oligo" 'potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.1"]' potentials.olig_intra_all=True potentials.olig_inter_all=True potentials.guide_scale=2.0 potentials.guide_decay="quadratic" 'contigmap.contigs=[480-480]'
