function [rtst, Vst] = sn_rt_stations(sn)
% [RTST, VST] = SN_RT_STATIONS(SN)
%
% Station-to-station routing probabilities and per-station visits.
%
% sn.rt and sn.visits are indexed by STATEFUL node, so a solver that writes
% traffic equations over stations and indexes them by station index silently
% reads the wrong rows as soon as the model owns a stateful node that is not a
% station (Router, Cache, stateful class switch). RTST is obtained from sn.rt
% by absorbing those nodes,
%
%   Pst = P_AA + P_AB * (I - P_BB)^-1 * P_BA,
%
% with A the station rows in station order and B the remaining stateful rows,
% which is exact because a non-station stateful node holds no jobs: it passes
% every arrival on instantaneously. When every stateful node is a station,
% RTST is sn.rt unchanged.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = sn.nclasses;
S = sn.nstateful;
M = sn.nstations;

isSt = false(1,S);
isSt(sn.stationToStateful) = true;

A = zeros(1, M*K);
for ist=1:M
    isf = sn.stationToStateful(ist);
    A((ist-1)*K + (1:K)) = (isf-1)*K + (1:K);
end
B = [];
for isf=1:S
    if ~isSt(isf)
        B = [B, (isf-1)*K + (1:K)]; %#ok<AGROW>
    end
end

P = full(sn.rt);
if isempty(B)
    rtst = P(A,A);
else
    rtst = P(A,A) + P(A,B) * ((eye(numel(B)) - P(B,B)) \ P(B,A));
end

if nargout > 1
    Vall = cellsum(sn.visits);
    Vst = Vall(sn.stationToStateful, :);
end
end
