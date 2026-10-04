import numpy as np
import copy

from line_solver import MetricType


class Event:
    """Represents a queueing network event (arrival, departure, etc.)."""

    def __init__(self, event_type, node, jobclass=None):
        self.event = event_type
        self.node = node
        self.jobclass = jobclass


class SampledMetric:
    """Observed data for a metric.

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """

    def __init__(self, metric_type, ts, data, node, jobclass=None):
        self.t = np.asarray(ts, dtype=float).flatten()
        self.data = np.asarray(data, dtype=float).flatten()
        self.type = metric_type
        self.node = node
        self.jobclass = jobclass
        self.cond = None
        self.format = 'timeseries'

    def set_conditional(self, event):
        self.cond = event

    def set_trace(self):
        self.format = 'trace'

    def is_aggregate(self):
        return self.jobclass is None

    def is_conditional(self):
        return self.cond is not None

    def is_trace(self):
        return self.format == 'trace'

    def copy(self):
        return copy.deepcopy(self)

    # MATLAB-compatible aliases
    setConditional = set_conditional
    setTrace = set_trace
    isAggregate = is_aggregate
    isConditional = is_conditional
    isTrace = is_trace
