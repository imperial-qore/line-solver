function [QN,UN,RN,TN,CN,XN,lG,runtime,iter,method] = solver_nc_sdr_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,LG,RUNTIME,ITER,METHOD] = SOLVER_NC_SDR_ANALYZER(SN, OPTIONS)
%
% Exact product-form analysis of a closed multiclass network with the
% state-dependent routing of Krzesinski (1987), "Multiclass Queueing Networks
% with State-Dependent Routing", Performance Evaluation 7:125-143. The joint
% distribution is eq. (16); the coefficients xi are those of Section 3.2,
% obtained from the state-independent part of the routing matrix.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
iter = 1;
% 'sdr' evaluates the product form (16) exactly by state enumeration, which is
% general in the branch topology; 'sdr.mva' runs the paper's Section 4 MVA and
% convolution, which costs O(J T M (V_1...V_J)^2) instead of the state-space
% size but requires single-center branches. Both are exact.
method = 'sdr';
if isfield(options,'method') && strcmpi(options.method,'sdr.mva')
    method = 'sdr.mva';
end

M = sn.nstations;
R = sn.nclasses;
sdr = sn.sdr;

% The premises of the product form (closed, no class switching, every stateful
% node a station, one FCFS rate per chain) are decided by NC_SDR_REFUSAL, which
% SolverNC's support gate asks too, so a pair the report offers is a pair this
% analyzer accepts.
sdrReason = nc_sdr_refusal(sn);
if ~isempty(sdrReason)
    line_error(mfilename, sdrReason);
end

% Mean service times and the load-dependent rate scalings alpha_i(n)
S = zeros(M,R);
for i = 1:M
    for r = 1:R
        if sn.rates(i,r) > 0 && isfinite(sn.rates(i,r))
            S(i,r) = 1/sn.rates(i,r);
        end
    end
end
Ntot = sum(sn.njobs);
alpha = ones(M, max(1,Ntot));
for i = 1:M
    if sn.sched(i) == SchedStrategy.INF
        alpha(i,:) = 1:size(alpha,2);
    elseif sn.nservers(i) > 1 && isfinite(sn.nservers(i))
        alpha(i,:) = min(1:size(alpha,2), sn.nservers(i));
    end
end
if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling)
    for i = 1:M
        for k = 1:min(size(alpha,2), size(sn.lldscaling,2))
            if sn.lldscaling(i,k) > 0
                alpha(i,k) = alpha(i,k) * sn.lldscaling(i,k);
            end
        end
    end
end

% State-independent routing probabilities, station indexed, one page per chain
P = zeros(M,M,R);
for r = 1:R
    for i = 1:M
        isf = sn.stationToStateful(i);
        for j = 1:M
            jsf = sn.stationToStateful(j);
            P(i,j,r) = sn.rt((isf-1)*R+r, (jsf-1)*R+r);
        end
    end
end
xi = pfqn_sdrvisits(sdr, P);

if strcmp(method,'sdr.mva')
    [QN, TN, UN, RN, lG] = pfqn_sdrmva(S, xi, sn.njobs(:)', sdr, alpha);
else
    [QN, TN, UN, RN, ~, lG] = pfqn_sdr(S, xi, sn.njobs(:)', sdr, alpha);
end

XN = zeros(1,R);
CN = zeros(1,R);
for r = 1:R
    ref = sn.refstat(r);
    XN(r) = TN(ref,r);
    if XN(r) > 0
        CN(r) = sn.njobs(r) / XN(r);
    end
end

runtime = toc(Tstart);
end
