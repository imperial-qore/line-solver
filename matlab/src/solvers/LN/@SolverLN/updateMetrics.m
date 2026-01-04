function updateMetrics(self, it)
switch self.options.method
    case 'moment3'
        % this method propagates through the layers 3 moments of the
        % response time distribution computed from the CDF obtained by the
        % solvers of the individual layers. In the present implementation,
        % calls are still assumed to be exponentially distributed.
        updateMetricsMomentBased(self,it)
    otherwise
        % default method for 'default', 'mva', 'nc', etc.
        updateMetricsDefault(self,it)
end
end