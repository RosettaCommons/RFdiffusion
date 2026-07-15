"""Named constants for RFdiffusion.

Centralizes magic numbers used across the codebase. Each constant
documents its meaning and where it originates.
"""

# ===== Amino Acid Encoding =====
NUM_AA_CLASSES = 22        # 20 standard amino acids + UNK + MASK
AA_MASK_TOKEN = 21         # Index for masked/unknown residue
AA_GLYCINE = 7             # Index for glycine in the alphabet

# ===== Atom Counts =====
N_BACKBONE_ATOMS = 4       # N, CA, C, O
N_HEAVY = 14               # Heavy atoms per residue (backbone + sidechain)
N_ALLATOM = 27             # All atoms per residue including hydrogens

# ===== Virtual Cbeta Reconstruction =====
# Coefficients for computing virtual Cbeta from backbone N, CA, C atoms.
# Derived from the cross product (CA-N) x (C-CA) basis vectors.
# Used in generate_Cbeta() across util.py, Embeddings.py, coords6d.py.
CBETA_A = -0.58273431      # coefficient for cross product vector
CBETA_B = 0.56802827       # coefficient for (CA - N) vector
CBETA_C = -0.54067466      # coefficient for (C - CA) vector

# ===== Distance Sentinel =====
# Large distance value used to indicate "no contact" or to mask
# self-interactions in distance matrices and top-k graphs.
NO_CONTACT_DIST = 999.9

# ===== Chain/Contig =====
# Index jump inserted between chains in the residue index array.
# Used by ContigMap and positional encodings to detect chain boundaries.
CHAIN_BREAK_INDEX_JUMP = 200

# Chain break detection threshold: gaps in residue index larger than
# this value are treated as chain breaks.
CHAIN_BREAK_DETECTION_THRESH = 35

# Maximum number of attempts when randomly sampling a valid contig
# length from a specified range.
CONTIG_MAX_SAMPLE_ATTEMPTS = 100_000

# ===== SE(3) Prediction Scaling =====
# Divisors applied to raw SE(3) transformer outputs to bring
# translations and rotations into physical scale.
SE3_TRANSLATION_SCALE = 10.0
SE3_ROTATION_SCALE = 100.0

# ===== Diffusion Schedule =====
# Reference number of timesteps for beta schedule scaling.
# When T != 200, betas are rescaled: beta *= 200/T.
BETA_SCHEDULE_REF_T = 200

# Minimum number of diffusion steps required for the schedule
# approximation to remain valid.
MIN_DIFFUSION_STEPS = 15

# ===== IGSO3 (Rotation Diffusion) =====
# Number of discrete sigma values for the IGSO3 distribution.
IGSO3_NUM_SIGMA = 500

# Truncation level L for the power series expansion of the
# IGSO3 probability density.
IGSO3_TRUNCATION_LEVEL = 2000

# ===== Peptide Geometry =====
# Ideal C-N peptide bond length in Angstroms.
PEPTIDE_BOND_LENGTH = 1.33
