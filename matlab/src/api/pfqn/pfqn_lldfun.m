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
    if ~isempty(lldscaling) && range(lldscaling(i,:))>0
        r(i) = r(i) / interp1(1:smax, lldscaling(i,1:smax), n(i), 'spline');
    end
end
end