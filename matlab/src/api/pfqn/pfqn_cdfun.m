%{
%{
 % @file pfqn_cdfun.m
 % @brief AMVA-QD class-dependence function for queue-dependent scaling.
%}
%}

%{
%{
 % @brief AMVA-QD class-dependence function for queue-dependent scaling.
 % @fn pfqn_cdfun(nvec, cdscaling)
 % @param nvec Population state vector.
 % @param cdscaling Cell array of class-dependent scaling functions.
 % @return r Scaling factor vector for each station.
%}
%}
function r = pfqn_cdfun(nvec,cdscaling)
% R = PFQN_CDFUN(NVEC,PHI)

% AMVA-QD class-dependence function

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
M = size(nvec,1);
r = ones(M,1);
if ~isempty(cdscaling)
    for i = 1:M
        r(i) = 1 / cdscaling{i}(nvec(i,:));
    end
end
end
