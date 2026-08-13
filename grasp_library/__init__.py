"""GRASP library designer: masked codon optimization with Pareto trade-offs."""

from .ligation_fidelity import (
    LigationFidelityCalculator,
    fidelity_calculator_for_level,
    score_grasp_cloning_stages,
)
from .assembly_interfaces import (
    build_order_fragment,
    custom_interface_preset,
    deposited_grasp_interface_preset,
    resolve_assembly_interfaces,
    validate_assembly_interfaces,
)
from .arelf import (
    ARELF_MOTIF,
    ARELF_OFFSETS,
    build_arelf_candidates,
    materialize_arelf_parts,
    parse_cut_token,
)
from .objectives import ObjectiveScores, evaluate_design
from .pareto import (
    ParetoPoint,
    dominates,
    knee_point,
    optimize_pareto_overhangs,
    pareto_front,
)
from .sample_codon_tables import (
    CUSTOM_FILE,
    FETCH_FROM_KAZUSA,
    SAMPLE_CODON_TABLES,
    UPLOAD_OWN_TABLE,
    builtin_sample_names,
    codon_table_dataframe,
    parse_frequency_block,
    sample_names,
    write_sample_codon_table,
)
from .kazusa import (
    KAZUSA_HOME,
    fetch_kazusa_codon_table,
    kazusa_table_url,
    normalize_kazusa_species_id,
)
from .codon_upload import (
    parse_codon_table_text,
    prompt_colab_codon_upload,
    write_uploaded_codon_table,
)
from .import_grasp import (
    compile_target_gap,
    import_grasp_profile,
    pick_parts_for_target,
)
from .restriction_sites import (
    COMMON_RESTRICTION_SITES,
    DEFAULT_SITE_BLACKLIST,
    apply_site_blacklist_to_config,
    parse_site_blacklist,
    resolve_site_blacklist,
    restriction_site_names,
    restriction_site_options,
)
from .synthesis_vendors import (
    ASSEMBLY_ENZYMES,
    GRASP_LIGATION_BY_LEVEL,
    GRASP_STAGE_MATCHED_LIGATION,
    LEVEL0_LIGATION,
    LEVEL1_LIGATION,
    LEVEL_MINUS1_LIGATION,
    LIGATION_TABLES,
    SYNTHESIS_VENDORS,
    apply_enzyme_to_config,
    apply_ligation_table_to_config,
    apply_vendor_to_config,
    enzyme_names,
    grasp_ligation_by_level,
    ligation_protocol_for_level,
    ligation_table_names,
    redesign_ligation_table_names,
    resolve_dashboard_ligation_choice,
    twist_length_advice,
    vendor_names,
)
from .control_panel import GraspControlPanel, build_default_config, wire_reload_button
from .workflows import (
    run_overhang_redesign,
    run_library_optimize,
    rescore_pareto_front_after_anneal,
    load_and_validate_parts,
    ensure_grasp_imported,
    export_optimized_library,
    compile_and_assemble_target,
    write_oligo_fasta,
    plot_library_pareto_after_anneal,
    select_from_rescored_front,
    run_library_redesign_and_anneal,
)
from .codon_tables import load_codon_usage, apply_organism_codon_table, validate_parts_for_organism
from .codon_validation import (
    analyze_cut_site_aa_risks,
    cds_matches_organism,
    validate_codon_table_aas,
    verify_cds_for_organism,
)
from .plotting import plot_pareto_front
from .optimizer import (
    optimize_coding_sequence,
    optimize_library,
    synthesis_qc,
    simulate_assembled_cds,
)
from .oneshot import run_oneshot_design, sanitize_rna_name
from .binder import (
    rna_to_binder_aa,
    describe_binder,
    describe_part_id,
    annotate_module_roles,
    normalize_target_rna,
)
from .gga_split import plan_gga_from_optimized_cds, suggest_fragment_count
from .paths import (
    bundled_profile_genbank,
    materialize_project,
    project_paths,
)
from .colab import write_notebook
from .genbank_export import write_annotated_genbank

try:
    from importlib.metadata import PackageNotFoundError, version

    __version__ = version("grasp-library-designer")
except PackageNotFoundError:  # pragma: no cover - editable/source tree
    __version__ = "0.1.12"

__all__ = [
    "LigationFidelityCalculator",
    "fidelity_calculator_for_level",
    "score_grasp_cloning_stages",
    "build_order_fragment",
    "custom_interface_preset",
    "deposited_grasp_interface_preset",
    "resolve_assembly_interfaces",
    "validate_assembly_interfaces",
    "ARELF_MOTIF",
    "ARELF_OFFSETS",
    "build_arelf_candidates",
    "materialize_arelf_parts",
    "parse_cut_token",
    "ObjectiveScores",
    "evaluate_design",
    "ParetoPoint",
    "dominates",
    "pareto_front",
    "optimize_pareto_overhangs",
    "knee_point",
    "SAMPLE_CODON_TABLES",
    "CUSTOM_FILE",
    "FETCH_FROM_KAZUSA",
    "UPLOAD_OWN_TABLE",
    "builtin_sample_names",
    "codon_table_dataframe",
    "parse_frequency_block",
    "sample_names",
    "write_sample_codon_table",
    "KAZUSA_HOME",
    "fetch_kazusa_codon_table",
    "kazusa_table_url",
    "normalize_kazusa_species_id",
    "parse_codon_table_text",
    "prompt_colab_codon_upload",
    "write_uploaded_codon_table",
    "import_grasp_profile",
    "compile_target_gap",
    "pick_parts_for_target",
    "SYNTHESIS_VENDORS",
    "COMMON_RESTRICTION_SITES",
    "DEFAULT_SITE_BLACKLIST",
    "apply_site_blacklist_to_config",
    "parse_site_blacklist",
    "resolve_site_blacklist",
    "restriction_site_names",
    "restriction_site_options",
    "ASSEMBLY_ENZYMES",
    "LIGATION_TABLES",
    "GRASP_LIGATION_BY_LEVEL",
    "GRASP_STAGE_MATCHED_LIGATION",
    "LEVEL_MINUS1_LIGATION",
    "LEVEL0_LIGATION",
    "LEVEL1_LIGATION",
    "vendor_names",
    "enzyme_names",
    "ligation_table_names",
    "redesign_ligation_table_names",
    "grasp_ligation_by_level",
    "ligation_protocol_for_level",
    "resolve_dashboard_ligation_choice",
    "apply_vendor_to_config",
    "apply_enzyme_to_config",
    "apply_ligation_table_to_config",
    "twist_length_advice",
    "GraspControlPanel",
    "build_default_config",
    "wire_reload_button",
    "run_overhang_redesign",
    "run_library_optimize",
    "rescore_pareto_front_after_anneal",
    "load_and_validate_parts",
    "ensure_grasp_imported",
    "export_optimized_library",
    "compile_and_assemble_target",
    "write_oligo_fasta",
    "plot_library_pareto_after_anneal",
    "select_from_rescored_front",
    "run_library_redesign_and_anneal",
    "load_codon_usage",
    "apply_organism_codon_table",
    "validate_parts_for_organism",
    "analyze_cut_site_aa_risks",
    "validate_codon_table_aas",
    "verify_cds_for_organism",
    "cds_matches_organism",
    "plot_pareto_front",
    "optimize_coding_sequence",
    "optimize_library",
    "synthesis_qc",
    "simulate_assembled_cds",
    "run_oneshot_design",
    "sanitize_rna_name",
    "rna_to_binder_aa",
    "describe_binder",
    "describe_part_id",
    "annotate_module_roles",
    "normalize_target_rna",
    "plan_gga_from_optimized_cds",
    "suggest_fragment_count",
    "bundled_profile_genbank",
    "materialize_project",
    "project_paths",
    "write_notebook",
    "write_annotated_genbank",
    "__version__",
]
