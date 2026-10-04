%{
%{
 % @file pfqn_lldfun.m
 % @brief AMVA-QD load and queue-dependent scaling function.
%}
%}

%{
%{
 % @brief AMVA-QD load and queue-dependent scaling function.
 % @fn pfqn_lldfun(n, lldscaling, nservers)
 % @param n Queue population vector.
 % @param lldscaling Load-dependent scaling matrix.
 % @param nservers Number of servers per station.
 % @return r Scaling factor vector.
%}
%}
function r = pfqn_lldfun(n,lldscaling, nservers)
% R = PFQN_LLDFUN(N,MU,C)

% AMVA-QD queue-dependence function

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = length(n);
r = ones(M,1);
smax = size(lldscaling,2);
alpha = 20; % softmin parameter
for i = 1:M
    %% handle servers
    if nargin>=3
        if isinf(nservers(i)) % delay server
            r(i) = 1; %1 / n(i); % handled in the main code differently so not needed
        else
            r(i) = r(i) / softmin(n(i),nservers(i),alpha);
            if isnan(r(i)) % if numerical problems in soft-min
                r(i) = 1 / min(n(i),nservers(i));
            end
        end
    end
    %% handle generic lld
    % A CONSTANT ROW IS A NO-OP ONLY WHEN IT IS ONE. The test here used to be
    % range(lldscaling(i,:))>0, which skipped alpha(n) = c for EVERY c: a station
    % declaring a uniform rate multiplier was then solved at its UNSCALED rate,
    % while the exact recursions (solver_mvald, pfqn_mvaldmx) and SolverCTMC
    % applied it, so AMVA contradicted them on the same model -- silently, and by
    % the whole factor c. Only a row of ones divides by 1 and may be skipped.
    % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
    if ~isempty(lldscaling) && any(lldscaling(i,1:smax) ~= 1)
        if smax == 1
            % interp1 needs two sample points, and a one-column lattice is the
            % constant alpha(n) = lldscaling(i,1): there is nothing to interpolate.
            r(i) = r(i) / lldscaling(i,1);
        else
            r(i) = r(i) / interp1(1:smax, lldscaling(i,1:smax), min(max(n(i),1),smax), 'linear');
        end
    end
end
end