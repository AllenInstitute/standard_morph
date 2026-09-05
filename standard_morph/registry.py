"""Metric registration, metadata, and lookup.

Metrics self-register on import (see each module in ``qc/metrics``). The engine
resolves metric names to registered instances via the module-level ``REGISTRY``.
Metrics are stateless with respect to a run -- ``evaluate`` receives all inputs --
so a single shared instance per metric is sufficient.
"""


class MetricRegistry:
    def __init__(self):
        self._metrics = {}

    def register(self, metric):
        """Register a metric instance keyed by its ``name``."""
        name = metric.name
        if name is None:
            raise ValueError(f"Metric {metric!r} has no name and cannot be registered.")
        if name in self._metrics:
            raise ValueError(f"A metric named '{name}' is already registered.")
        self._metrics[name] = metric
        return metric

    def get(self, name):
        """Return the registered metric instance for ``name``."""
        if name not in self._metrics:
            raise KeyError(
                f"Unknown metric '{name}'. Registered metrics: {self.names()}"
            )
        return self._metrics[name]

    def names(self):
        """Return all registered metric names, sorted."""
        return sorted(self._metrics)

    def all(self):
        """Return all registered metric instances."""
        return list(self._metrics.values())

    def __contains__(self, name):
        return name in self._metrics


#: Process-wide registry populated as metric modules are imported.
REGISTRY = MetricRegistry()


def register(metric):
    """Register a metric instance with the global registry (convenience)."""
    return REGISTRY.register(metric)
