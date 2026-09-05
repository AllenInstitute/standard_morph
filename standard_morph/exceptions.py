"""Exception types for the QC framework."""


class QCError(Exception):
    """Base class for all QC framework errors."""


class IncompatibleMetricContextError(QCError):
    """Raised when a requested metric is incompatible with the run context.

    Per the architecture spec, incompatibility is terminal: the run halts
    immediately rather than skipping the offending metric.
    """


class MissingPolicyValueError(QCError):
    """Raised by ``policy[metric_name, key]`` when a required threshold is absent."""


class MissingPolicyValuesError(QCError):
    """Raised during pre-flight when one or more metrics lack required policy keys."""
