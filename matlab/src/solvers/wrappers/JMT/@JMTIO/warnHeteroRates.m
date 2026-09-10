function warnHeteroRates(self, ind)
% WARNHETERORATES(IND)
% Warns that per-server-type service rates cannot reach the JMT engine.
%
% JMT keys the ServiceStrategy array of a station by refClass, so its loader
% (jmt.engine.simEngine.SimLoader) keeps one strategy per class however many
% (type, class) entries are written, and every pool of a station ends up
% serving at the class rate. Pool sizes, class compatibilities and the
% assignment policy do cross; setHeteroService rates do not.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
np = sn.nodeparam{ind};
if ~isstruct(np) || ~isfield(np, 'heterorates') || isempty(np.heterorates)
    return
end
rates = np.heterorates;
if ~any(rates(:) > 0)
    return
end
% One rate row, or identical rows, is the homogeneous case and crosses intact
distinct = rates(rates > 0);
if size(rates, 1) == 1 || all(abs(distinct - distinct(1)) < 1e-12)
    return
end
line_warning(mfilename, 'JMT keys service strategies by job class, so the per-server-type service rates of station %s cannot be exported; every pool will serve at the class service rate. Use the LDES or CTMC solver for per-type rates.', sn.nodenames{ind});
end
