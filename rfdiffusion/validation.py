"""Input validation for RFdiffusion inference.

Catches common configuration and input errors early, before model loading
and GPU allocation, so users get clear error messages instead of cryptic
tensor shape mismatches deep in the forward pass.
"""

import os
import re
import logging

logger = logging.getLogger(__name__)


class ValidationError(ValueError):
    """Raised when input validation fails with a user-friendly message."""
    pass


def validate_pdb_path(pdb_path: str) -> None:
    """Validate that a PDB file exists and contains parseable ATOM records.

    Args:
        pdb_path: Path to input PDB file.

    Raises:
        ValidationError: If file doesn't exist or has no ATOM records.
    """
    if not os.path.isfile(pdb_path):
        raise ValidationError(
            f"Input PDB file not found: {pdb_path}"
        )

    has_atoms = False
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                has_atoms = True
                try:
                    float(line[30:38])
                    float(line[38:46])
                    float(line[46:54])
                except ValueError:
                    raise ValidationError(
                        f"Invalid coordinates in PDB line: {line.rstrip()}"
                    )
                break

    if not has_atoms:
        raise ValidationError(
            f"PDB file contains no ATOM/HETATM records: {pdb_path}"
        )


def validate_contig_string(contigs: list) -> None:
    """Validate contig string syntax before parsing.

    Args:
        contigs: List of contig specification strings.

    Raises:
        ValidationError: If contig syntax is invalid.
    """
    if not contigs or not isinstance(contigs, (list, tuple)):
        raise ValidationError(
            "contigs must be a non-empty list of strings. "
            "Example: ['10-20/A5-50/0 30-40']"
        )

    contig_str = contigs[0]
    if not isinstance(contig_str, str) or not contig_str.strip():
        raise ValidationError(
            f"Contig string must be a non-empty string, got: {contig_str!r}"
        )

    for segment in contig_str.strip().split():
        for part in segment.split("/"):
            part = part.strip()
            if not part:
                continue
            # Chain break marker
            if part == "0":
                continue
            # Numeric range: "10-20" or "10"
            if part[0].isdigit():
                if "-" in part:
                    pieces = part.split("-")
                    if len(pieces) != 2:
                        raise ValidationError(
                            f"Invalid contig range format: '{part}'. "
                            f"Expected 'N-M' (e.g., '10-20')."
                        )
                    try:
                        lo, hi = int(pieces[0]), int(pieces[1])
                    except ValueError:
                        raise ValidationError(
                            f"Non-integer values in contig range: '{part}'"
                        )
                    if lo < 0 or hi < 0:
                        raise ValidationError(
                            f"Negative value in contig range: '{part}'"
                        )
                    if lo > hi:
                        raise ValidationError(
                            f"Invalid contig range: '{part}' (start > end)"
                        )
            # Chain-residue range: "A5-50" or "A5"
            elif part[0].isalpha():
                if not re.match(r"^[A-Za-z]\d+(-\d+)?$", part):
                    logger.warning(f"Unusual contig segment: '{part}'")


def validate_checkpoint_path(ckpt_path: str) -> None:
    """Validate that a model checkpoint file exists.

    Args:
        ckpt_path: Path to model checkpoint.

    Raises:
        ValidationError: If checkpoint file doesn't exist.
    """
    if not os.path.isfile(ckpt_path):
        raise ValidationError(
            f"Model checkpoint not found: {ckpt_path}. "
            f"Please download models following the README instructions."
        )


def validate_hotspot_res(hotspot_res: list) -> None:
    """Validate hotspot residue format (e.g., ['A50', 'B123']).

    Args:
        hotspot_res: List of hotspot residue strings.

    Raises:
        ValidationError: If format is invalid.
    """
    if hotspot_res is None:
        return

    for res in hotspot_res:
        if not isinstance(res, str) or len(res) < 2:
            raise ValidationError(
                f"Invalid hotspot residue format: {res!r}. "
                f"Expected format like 'A50' (chain letter + residue number)."
            )
        if not res[0].isalpha():
            raise ValidationError(
                f"Hotspot residue must start with a chain letter: {res!r}"
            )
        try:
            int(res[1:])
        except ValueError:
            raise ValidationError(
                f"Hotspot residue number must be an integer: {res!r}"
            )


def validate_diffuser_config(diffuser_conf) -> None:
    """Validate diffuser configuration parameters.

    Args:
        diffuser_conf: Diffuser configuration object.

    Raises:
        ValidationError: If parameters are out of valid range.
    """
    T = getattr(diffuser_conf, "T", None)
    partial_T = getattr(diffuser_conf, "partial_T", None)

    if T is not None and T < 1:
        raise ValidationError(
            f"diffuser.T must be >= 1, got {T}"
        )
    if partial_T is not None:
        if partial_T < 1:
            raise ValidationError(
                f"diffuser.partial_T must be >= 1, got {partial_T}"
            )
        if T is not None and partial_T > T:
            raise ValidationError(
                f"diffuser.partial_T ({partial_T}) cannot exceed "
                f"diffuser.T ({T})"
            )
