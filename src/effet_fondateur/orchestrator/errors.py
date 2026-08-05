"""Erreurs contrôlées de l'orchestrateur V2."""


class PipelineError(RuntimeError):
    """Erreur attendue empêchant la poursuite sûre d'un run."""


class IntegrityError(PipelineError):
    """Une sortie publiée ne correspond plus à son empreinte auditée."""


class StageExecutionError(PipelineError):
    """Un sous-processus d'étape a échoué ou publié des sorties invalides."""

    def __init__(self, message: str, return_code: int = 5) -> None:
        super().__init__(message)
        self.return_code = return_code
