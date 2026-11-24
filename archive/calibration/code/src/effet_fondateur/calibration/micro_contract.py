"""Topologie technique fermée du micro-panel transversal R10."""

from __future__ import annotations

from effet_fondateur.calibration.contracts import STAGE_IDS


CORE_STAGES = tuple(stage for stage in STAGE_IDS if stage not in {"18", "19"})

STAGE_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "00": (),
    "01": ("00",),
    "02": ("01",),
    "03": ("02",),
    "04": ("03",),
    "05": ("03",),
    "06": ("05",),
    "07": ("06",),
    "08": ("07",),
    "09": ("04", "07", "08"),
    "10": ("09",),
    "10A": ("10",),
    "11": ("10",),
    "12": ("11",),
    "12A": ("10A",),
    "13": ("12",),
    "14": ("13",),
    "15": ("11",),
    "16": ("10", "11"),
    "16A": ("08", "12", "16"),
    "16B": ("13", "16A"),
    "16C": ("12", "12A", "13", "16A"),
    "17": (
        "03", "05", "07", "08", "10", "10A", "11", "12", "12A",
        "13", "14", "15", "16", "16A", "16B", "16C",
    ),
    "18": tuple(stage for stage in STAGE_IDS if stage not in {"18", "19"}),
    "19": tuple(stage for stage in STAGE_IDS if stage != "19"),
}

STAGE_PRIMARY_OUTPUT: dict[str, str] = {
    "00": "initialization_summary.json",
    "01": "source_qc.tsv",
    "02": "sample_registry_review.tsv",
    "03": "conversion_report.json",
    "04": "target_variant_report.json",
    "05": "qc_preliminary_report.json",
    "06": "kinship_panel_report.json",
    "07": "kinship_report.json",
    "08": "population_structure_report.json",
    "09": "cohort_freeze_report.json",
    "10": "qc_final_report.json",
    "10A": "autosomal_phasing_panel_summary.json",
    "11": "target_region_report.json",
    "12": "phase_target_region_summary.json",
    "12A": "autosomal_phasing_manifest.json",
    "13": "founder_analysis_summary.json",
    "14": "variant_age_summary.json",
    "15": "local_ld_analysis_summary.json",
    "16": "roh_analysis_summary.json",
    "16A": "reference_ancestry_summary.json",
    "16B": "founder_haplotype_enrichment_summary.json",
    "16C": "explicit_ibd_summary.json",
    "17": "sensitivity_analysis_summary.json",
    "18": "visualization_completeness.json",
    "19": "interpretation_facts.json",
}

STAGE_VIEWS: dict[str, tuple[str, ...]] = {
    "00": (),
    "01": ("V-SOURCE",),
    "02": ("V-SOURCE", "V-PEDIGREE"),
    "03": ("V-SOURCE", "V-DENSITY", "V-MEASUREMENT"),
    "04": ("V-TARGET",),
    "05": ("V-DENSITY", "V-MEASUREMENT", "V-AUTOSOMAL"),
    "06": ("V-DENSITY", "V-AUTOSOMAL"),
    "07": ("V-PEDIGREE", "V-AUTOSOMAL"),
    "08": ("V-DENSITY", "V-AUTOSOMAL"),
    "09": ("V-TARGET", "V-PEDIGREE"),
    "10": ("V-TARGET", "V-DENSITY", "V-MEASUREMENT"),
    "10A": ("V-AUTOSOMAL",),
    "11": ("V-REGIONAL",),
    "12": ("V-REGIONAL", "V-PHASING-REFERENCE", "V-PEDIGREE"),
    "12A": ("V-AUTOSOMAL", "V-PHASING-REFERENCE"),
    "13": ("V-REGIONAL", "V-FOUNDER-NULL"),
    "14": ("V-REGIONAL",),
    "15": ("V-REGIONAL", "V-DENSITY"),
    "16": ("V-REGIONAL", "V-ROH-IBD"),
    "16A": ("V-REGIONAL", "V-PHASING-REFERENCE", "V-AUTOSOMAL"),
    "16B": ("V-REGIONAL", "V-PHASING-REFERENCE", "V-FOUNDER-NULL"),
    "16C": ("V-AUTOSOMAL", "V-REGIONAL", "V-ROH-IBD"),
    "17": ("V-SENSITIVITY",),
    "18": (),
    "19": (),
}


def validate_micro_topology() -> None:
    """Refuse une étape absente, un cycle ou une dépendance vers l'avenir."""
    if tuple(STAGE_DEPENDENCIES) != STAGE_IDS:
        raise ValueError("micro_stage_set_mismatch")
    if set(STAGE_PRIMARY_OUTPUT) != set(STAGE_IDS) or set(STAGE_VIEWS) != set(STAGE_IDS):
        raise ValueError("micro_stage_contract_incomplete")
    position = {stage: index for index, stage in enumerate(STAGE_IDS)}
    for stage, dependencies in STAGE_DEPENDENCIES.items():
        if any(dependency not in position for dependency in dependencies):
            raise ValueError(f"unknown_micro_dependency:{stage}")
        if any(position[dependency] >= position[stage] for dependency in dependencies):
            raise ValueError(f"non_acyclic_micro_dependency:{stage}")
