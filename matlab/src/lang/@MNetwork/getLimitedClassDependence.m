function gamma = getLimitedClassDependence(self)
% gamma = GETLIMITEDCLASSDEPENDENCE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
gamma = cell(M,1);
hasAny = false;
for i=1:M
    if ~isempty(self.stations{i}.lcdScaling)
        gamma{i} = self.stations{i}.lcdScaling; % function handle
        hasAny = true;
    end
end
if ~hasAny
    % No class-dependent station: return {} so that ISEMPTY(sn.cdscaling)
    % remains the test for "the model has no class dependence at all".
    gamma = {};
    return
end
% Non class-dependent stations are deliberately left EMPTY rather than
% filled with @(nvec) 1. The entry is the class-dependent service rate
% beta_{i,r}(n) (Sauer 1983, eq. (40)), for which a constant 1 would assert
% "every class completes at rate 1" -- not load independence (an LI station is
% beta_{i,r}(n) = mu_i * n_r/|n|, whose n_r/|n| factor is what regenerates the
% multinomial). Consumers treat an empty entry as "no class dependence":
% PFQN_CDFUN skips it, and PFQN_CONV keeps the station on the load-independent
% recurrence.
end
