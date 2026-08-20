from effet_fondateur.orchestrator.pipeline import BUILD_REPORT_STAGE, DEFAULT_STAGE_DEFINITIONS


def test_founder_enrichment_is_between_16a_and_17() -> None:
    stages = [(stage.stage_id, stage.stage_name) for stage in DEFAULT_STAGE_DEFINITIONS]
    names = [name for _, name in stages]
    enrichment_index = names.index("evaluate_founder_haplotype_enrichment")
    assert stages[enrichment_index] == ("16B", "evaluate_founder_haplotype_enrichment")
    assert names[enrichment_index - 1] == "analyze_reference_ancestry"
    assert names[enrichment_index + 1] == "call_explicit_ibd"
    assert stages[enrichment_index + 1] == ("16C", "call_explicit_ibd")
    assert names[enrichment_index + 2] == "run_sensitivity_analyses"


def test_report_direct_fact_producers_are_explicit_dependencies() -> None:
    assert {
        "infer_founder_haplotype",
        "estimate_variant_age",
        "analyze_roh",
        "analyze_reference_ancestry",
        "call_explicit_ibd",
        "run_sensitivity_analyses",
    } <= set(BUILD_REPORT_STAGE.dependencies)
