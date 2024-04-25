from typing import Any, Literal


class ContigSettingsConfig(dict[str, Any]):
    ref_idx: Any
    hal_idx: Any
    idx_rf: Any


class PreprocessConfig(dict[str, Any]):
    sidechain_input: bool
    motif_sidechain_input: bool
    d_t1d: int
    d_t2d: int


class LoggingConfig(dict[str, Any]):
    inputs: Any
    level: str


class InferenceConfig(dict[str, Any]):
    input_pdb: Any
    num_designs: int
    design_startnum: int
    ckpt_override_path: Any
    symmetry: Any
    recenter: bool
    radius: float
    model_only_neighbors: bool
    output_prefix: str
    write_trajectory: bool
    scaffold_guided: bool
    model_runner: Literal["SelfConditioning", "ScaffoldedSampler", "default"]
    cautious: bool
    align_motif: bool
    symmetric_self_cond: bool
    final_step: int
    deterministic: bool
    seed: int
    trb_save_ckpt_path: Any
    schedule_directory_path: None | str
    model_directory_path: None | str
    device_name: str


class ContigMapConfig(dict[str, Any]):
    contigs: Any
    inpaint_seq: Any
    inpaint_str: Any
    provide_seq: Any
    length: str | None


class ScaffoldguidedConfig(dict[str, Any]):
    scaffoldguided: bool
    target_pdb: bool
    target_path: Any
    scaffold_list: Any
    scaffold_dir: Any
    sampled_insertion: Any
    sampled_N: Any
    sampled_C: Any
    ss_mask: Any
    systematic: Any
    target_ss: Any
    target_adj: Any
    mask_loops: Any
    contig_crop: Any


class DiffuserConfig(dict[str, Any]):
    T: int
    partial_T: int


class RFDiffusionConfig(dict[str, Any]):
    inference: InferenceConfig
    contigmap: ContigMapConfig
    model: Any
    diffuser: DiffuserConfig
    denoiser: Any
    ppi: Any
    potentials: Any
    contig_settings: ContigSettingsConfig
    preprocess: PreprocessConfig
    logging: LoggingConfig
    scaffoldguided: ScaffoldguidedConfig
