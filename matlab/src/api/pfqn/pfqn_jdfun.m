%{
%{
 % @file pfqn_jdfun.m
 % @brief AMVA joint-dependence function for non-product-form scaling.
%}
%}

%{
%{
 % @brief AMVA joint-dependence function for non-product-form scaling.
 % @fn pfqn_jdfun(nvec, jdscaling, classIdx)
 % @param nvec Population state vector.
 % @param jdscaling Cell array of joint-dependent scaling functions.
 % @param classIdx Optional class index selecting eta_{i,r} (default: 1).
 % @return r Scaling factor vector for each station.
%}
%}
function r = pfqn_jdfun(nvec, jdscaling, classIdx)
% R = PFQN_JDFUN(NVEC, JDSCALING, CLASSIDX)
%
% AMVA joint-dependence function. Returns, for every station i, the
% reciprocal of the joint-dependent scaling
%   eta_i(n_i1, ..., n_iR)
% evaluated at the per-class population vector NVEC(i,:), for class r=CLASSIDX.
%
% JDSCALING{i} is a function handle of the joint per-class population vector at
% station i. It may return either
%   - a scalar, i.e. a scaling eta_i(n) shared by every class (broadcast, as
%     in the flagship min(ni(1),c)), or
%   - a vector of length R, i.e. per-class scalings [eta_{i,1}(n), ...,
%     eta_{i,R}(n)], of which element CLASSIDX is taken (Sauer chain-dependent
%     rate mu_{r,i}(n)).
% Unlike PFQN_CDFUN (product-form beta_{i,r} depending on the own-class
% marginal n_{i,r}), eta may read the joint vector arbitrarily and is
% therefore NON-product-form: the AMVA result is an approximation with no
% exactness/uniqueness guarantee. The numerical evaluation matches PFQN_CDFUN;
% the distinction is semantic (product-form vs joint) and is carried by the
% separate sn.jdscaling field.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin < 3 || isempty(classIdx)
    classIdx = 1;
end
M = size(nvec,1);
r = ones(M,1);
if ~isempty(jdscaling)
    for i = 1:M
        if isempty(jdscaling{i})
            continue
        end
        v = jdscaling{i}(nvec(i,:));
        if numel(v) > 1
            % per-class eta_{i,r}: select the requested class
            v = v(classIdx);
        end
        r(i) = 1 / v;
    end
end
end
