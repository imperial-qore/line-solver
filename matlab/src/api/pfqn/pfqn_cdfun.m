%{
%{
 % @file pfqn_cdfun.m
 % @brief AMVA-QD class-dependence function for queue-dependent scaling.
%}
%}

%{
%{
 % @brief AMVA-QD class-dependence function for queue-dependent scaling.
 % @fn pfqn_cdfun(nvec, cdscaling, classIdx)
 % @param nvec Population state vector.
 % @param cdscaling Cell array of class-dependent scaling functions.
 % @param classIdx Optional class index selecting beta_{i,r} (default: 1).
 % @return r Scaling factor vector for each station.
%}
%}
function r = pfqn_cdfun(nvec, cdscaling, classIdx)
% R = PFQN_CDFUN(NVEC, CDSCALING, CLASSIDX)
%
% AMVA-QD class-dependence function. Returns, for every station i, the
% reciprocal of the class-dependent scaling
%   beta_{i,r}(n_i1, ..., n_iR)
% evaluated at the per-class population vector NVEC(i,:), for class r=CLASSIDX.
%
% CDSCALING{i} is a function handle of the per-class population vector at
% station i. It may return either
%   - a scalar, i.e. a chain-independent scaling beta_i(n) shared by every
%     class (the common case, and the historical contract), or
%   - a vector of length R, i.e. the per-class scalings
%     [beta_{i,1}(n), ..., beta_{i,R}(n)], of which element CLASSIDX is taken.
% The per-class form expresses Sauer's chain-dependent service rates
% mu_{r,i}(n) (Sauer 1983, "Computational Algorithms for State-Dependent
% Queueing Networks", eq. (40)), so a single class-dependence mechanism covers
% both the chain-independent and the chain-specific cases.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin < 3 || isempty(classIdx)
    classIdx = 1;
end
M = size(nvec,1);
r = ones(M,1);
if ~isempty(cdscaling)
    for i = 1:M
        if isempty(cdscaling{i})
            continue
        end
        v = cdscaling{i}(nvec(i,:));
        if numel(v) > 1
            % per-class beta_{i,r}: select the requested class
            v = v(classIdx);
        end
        r(i) = 1 / v;
    end
end
end
