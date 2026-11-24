"""Noyau contractuel de la calibration unifiée R00–R14.

Les anciens bancs exploratoires restent dans le paquet pour audit historique,
mais ne sont pas exportés et ne sont pas des producteurs autorisés de
``panel_master_generic_v1``.
"""

from effet_fondateur.calibration.contracts import (
    SCHEMA_FILES,
    STAGE_IDS,
    TRUTH_CONSUMERS,
    VIEW_CONSUMERS,
    CalibrationContractError,
)
from effet_fondateur.calibration.evaluation import (
    CalibrationEvaluationError,
    EvaluationStatus,
    StageDecision,
    decide_stage,
)
from effet_fondateur.calibration.readiness import (
    ReadinessAssessment,
    ReadinessError,
    assess_r09_repository,
    validate_readiness_receipt,
)
from effet_fondateur.calibration.release_guard import (
    ReleaseGuardError,
    partition_pseudonym,
    validate_blinded_inventory,
    validate_partition_separation,
    validate_single_factor_pair,
)
from effet_fondateur.calibration.truth_router import (
    TruthGrant,
    TruthRoutingError,
    authorize_truth_access,
)


__all__ = [
    "SCHEMA_FILES",
    "STAGE_IDS",
    "TRUTH_CONSUMERS",
    "VIEW_CONSUMERS",
    "CalibrationContractError",
    "CalibrationEvaluationError",
    "EvaluationStatus",
    "ReadinessAssessment",
    "ReadinessError",
    "ReleaseGuardError",
    "StageDecision",
    "TruthGrant",
    "TruthRoutingError",
    "assess_r09_repository",
    "authorize_truth_access",
    "decide_stage",
    "partition_pseudonym",
    "validate_blinded_inventory",
    "validate_partition_separation",
    "validate_readiness_receipt",
    "validate_single_factor_pair",
]
