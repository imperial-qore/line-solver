function updateMetrics(self, it)
switch self.options.method
    case 'moment3'
        % see _kb/06-solver-catalog.md (LN section) for rationale
        updateMetricsMomentBased(self,it)
    otherwise
        % default method for 'default', 'mva', 'nc', etc.
        updateMetricsDefault(self,it)
end
end