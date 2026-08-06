"""Inférence conservatrice d'un haplotype fondateur local."""

from effet_fondateur.founder.local_ibs import (
    FounderAnalysisError,
    FounderAnalysisResult,
    infer_target_centered_ibs,
)

__all__ = [
    "FounderAnalysisError",
    "FounderAnalysisResult",
    "infer_target_centered_ibs",
]
