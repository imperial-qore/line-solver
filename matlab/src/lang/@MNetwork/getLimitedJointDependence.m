function eta = getLimitedJointDependence(self)
% eta = GETLIMITEDJOINTDEPENDENCE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
eta = cell(M,1);
hasAny = false;
for i=1:M
    if ~isempty(self.stations{i}.ljdScaling)
        eta{i} = self.stations{i}.ljdScaling; % function handle
        hasAny = true;
    end
end
if ~hasAny
    % No joint-dependent station: return {} so that ISEMPTY(sn.jdscaling)
    % remains the test for "the model has no joint dependence at all".
    eta = {};
    return
end
% Non joint-dependent stations are deliberately left EMPTY rather than filled
% with @(nvec) 1. The entry is the joint-dependent service-rate scaling
% eta_i(n) (non-product-form: eta may read the joint vector arbitrarily).
% Consumers treat an empty entry as "no joint dependence": PFQN_JDFUN skips it.
end
