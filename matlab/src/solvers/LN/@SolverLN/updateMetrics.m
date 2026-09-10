function updateMetrics(self, it)
switch self.lnmethod
    case 'moment3'
        % see _kb/06-solver-catalog.md (LN section) for rationale
        updateMetricsMomentBased(self,it)
    case {'srvn.ph','flat.ph'}
        % see _kb/06-solver-catalog.md (LN section) for rationale
        updateMetricsPH(self,it)
    otherwise
        % 'srvn.cs': the routing encoding of the activity graph
        updateMetricsDefault(self,it)
end
end
