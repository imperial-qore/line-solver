function f = cdclassfactor(beta, nirmat, rows, class)
% F = CDCLASSFACTOR(BETA, NIRMAT, ROWS, CLASS)
%
% Per-row class-dependence factor for event rates. BETA is a class-dependence
% handle mapping a 1xR per-class population vector n to the 1xR vector of
% dimensionless rate scalings beta_r(n) (see fes_beta_handle). NIRMAT is the
% per-state per-class population matrix (nstates x R, one row per enabled
% state). ROWS selects the enabled rows (logical mask or index vector over the
% rows of NIRMAT) and CLASS is the class completing service.
%
% For each selected row j the handle is evaluated on that row's population
% vector and the component of the completing class is returned:
%   F(j) = v(min(CLASS, numel(v))),  v = BETA(NIRMAT(ROWS(j),:))
% so a neutral scaling @(n) 1 (scalar) yields F(j) = 1 for every class.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if islogical(rows)
    rows = find(rows);
end
rows = rows(:);
f = ones(numel(rows),1);
for j = 1:numel(rows)
    v = beta(nirmat(rows(j),:));
    f(j) = v(min(class, numel(v)));
end
end
