%{
%{
 % @file pfqn_mvams_ilock.m
 % @brief MVA entry point for models carrying the interlocked-flow correction.
%}
%}

%{
%{
 % @brief MVA entry point for models carrying the interlocked-flow correction.
 % @fn pfqn_mvams_ilock(lambda, L, N, Z, mi, S, IL)
 % @param lambda Arrival rate vector.
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mi Queue replication factors (default: ones).
 % @param S Number of servers per station (default: ones).
 % @param IL Interlock matrix (R x R), see PFQN_MVA_ILOCK.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Residence times (M x R).
 % @return lG Always NaN.
%}
%}
function [XN,QN,UN,CN,lG]=pfqn_mvams_ilock(lambda,L,N,Z,mi,S,IL)
% [XN,QN,UN,CN,LOGG]=PFQN_MVAMS_ILOCK(LAMBDA,L,N,Z,MI,S,IL)
%
% The interlock is defined only for closed single-server models, so this is the
% one shape accepted here; anything else is refused rather than served without
% the correction. Models with no interlock go to PFQN_MVAMS.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,R]=size(L); % get number of queues (M) and classes (R)
hasOpenClasses = false;
for r = 1:R
    if isinf(N(r))
        hasOpenClasses = true;
    end
end
if nargin<7 || isempty(IL)
    line_error(mfilename,'an interlock matrix is required; use pfqn_mvams for the standard arrival theorem');
end
if nargin<6 %~exist('S','var')
    S = ones(M,1);
end
if nargin<5 %~exist('mi','var')
    mi = ones(M,1);
end
if nargin<4 %~exist('Z','var')
    Z = zeros(1,R);
end

hasMultiServer = false;
for ist = 1:M
    if isfinite(S(ist)) && S(ist) > 1
        hasMultiServer = true;
        break;
    end
end
if hasOpenClasses || hasMultiServer
    line_error(mfilename,'the interlock correction is available in exact MVA for closed single-server models only; use an AMVA method for this model.');
end
if any(lambda ~= 0)
    line_error(mfilename,'the interlock correction is available in exact MVA for closed single-server models only; use an AMVA method for this model.');
end

[XN,QN,UN,CN,lG] = pfqn_mva_ilock(L,N,Z,mi,IL);
end
