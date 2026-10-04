function [ok, QN, UN, RN, TN, CN, XN] = solver_mam_mapmap1_exact(sn)
% SOLVER_MAM_MAPMAP1_EXACT Exact MAP/MAP/1 solution for a single-class,
% single-server, open Source -> FCFS Queue -> Sink model.
%
% The decomposition methods (dec.source/dec.mmap) approximate this queue: they
% either fit the arrival to a simpler process or treat the service as a renewal
% phase-type, discarding autocorrelation. When the arrival or the service is a
% genuinely correlated (non-renewal) MAP, this routine returns the exact mean
% queue length via the matrix-geometric MAP/MAP/1 solver (Q_CT_MAP_MAP_1), which
% carries both phase processes across arrivals and departures.
%
% Returns ok=false (and empty metrics) when the model is not an exactly-solvable
% single MAP/MAP/1 queue, so the caller falls back to the decomposition methods.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
ok = false; QN = []; UN = []; RN = []; TN = []; CN = []; XN = [];

M = sn.nstations;
K = sn.nclasses;
if K ~= 1
    return;
end
if ~isinf(sn.njobs)
    return;   % open model only
end

sourceIdx = find(sn.sched == SchedStrategy.EXT);
queueIdx  = find(sn.sched == SchedStrategy.FCFS);
if numel(sourceIdx) ~= 1 || numel(queueIdx) ~= 1
    return;
end
if sn.nservers(queueIdx) ~= 1
    return;
end
% Require source-adjacent queue (only two stations); see _kb/06-solver-catalog.md for rationale
if any(~ismember(1:M, [sourceIdx queueIdx]))
    return;
end

arrProc = sn.proc{sourceIdx}{1};
svcProc = sn.proc{queueIdx}{1};
if numel(arrProc) < 2 || numel(svcProc) < 2
    return;
end
Da0 = arrProc{1}; Da1 = arrProc{2};
Ds0 = svcProc{1}; Ds1 = svcProc{2};

% Genuine MAPs only (RAP/ME left to fallback); see _kb/06-solver-catalog.md for rationale
if ~is_markovian_map(Da0, Da1) || ~is_markovian_map(Ds0, Ds1)
    return;
end

% Fire only for non-renewal (correlated) processes; see _kb/06-solver-catalog.md for rationale
if is_renewal_map(Da0, Da1) && is_renewal_map(Ds0, Ds1)
    return;
end

lambda = map_lambda({Da0, Da1});
mu = map_lambda({Ds0, Ds1});
if ~(lambda < mu)
    return;   % unstable or degenerate: leave to the fallback path
end

ql = Q_CT_MAP_MAP_1(Da0, Da1, Ds0, Ds1, 'MaxNumComp', 100000);
ql = ql(:);
EN = sum((0:numel(ql)-1)' .* ql);
rho = lambda / mu;

QN = zeros(M, K); UN = zeros(M, K); RN = zeros(M, K); TN = zeros(M, K);
TN(sourceIdx, 1) = lambda;
TN(queueIdx, 1) = lambda;
QN(queueIdx, 1) = EN;
UN(queueIdx, 1) = rho;
RN(queueIdx, 1) = EN / lambda;
CN = sum(RN, 1);
XN = lambda;
ok = true;
end

function tf = is_markovian_map(D0, D1)
% Single definition in mam_is_markovian_map, shared with solver_mam_ag, which
% needs the same predicate to keep a RAP/ME out of the RCAT phase construction.
tf = mam_is_markovian_map(D0, D1);
end

function tf = is_renewal_map(D0, D1)
% Single definition in mam_is_renewal_map, shared with solver_mam_basic, which
% needs the same predicate to keep a correlated arrival off the marginal-only
% PH/M/c closed form.
tf = mam_is_renewal_map(D0, D1);
end
