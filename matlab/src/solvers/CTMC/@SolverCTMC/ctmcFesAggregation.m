function [QN,UN,RN,TN,CN,XN] = ctmcFesAggregation(self, sn, options)
% [QN,UN,RN,TN,CN,XN] = CTMCFESAGGREGATION(SN, OPTIONS)
%
% Solve the model with a station subset replaced by a FLOW-EQUIVALENT SERVER,
% then recover the collapsed stations' own metrics by conditioning. Selected by
% options.config.fes_stations, the 1-based indices of the stations to collapse.
%
% WHY THIS EXISTS. ModelAdapter.aggregateFES has existed in all four codebases
% with no solver consumer at all: it was exercised by examples and tests only,
% so nothing in the solver stack depended on it. Flow-equivalent aggregation is
% the standard route to HIERARCHICAL DECOMPOSITION -- a subnetwork is solved in
% isolation and enters the outer chain as a single load-dependent station, which
% is what makes an otherwise intractable state space tractable. This is that
% consumer.
%
% HOW THE COLLAPSED STATIONS ARE REPORTED. The reduced model answers for the
% surviving stations directly. For a collapsed station the answer is the
% Chandy-Herzog-Woo conditional sum
%
%   E[Q_i] = sum_n P(N_fes = n) * Q_i(n),
%
% with P(N_fes = n) read off the reduced chain's stationary law and Q_i(n) the
% isolated subnetwork's queue length at population n (FES_COMPUTE_METRICS).
% Throughput needs no conditioning: flow through a station is fixed by the
% routing, T_i = X * V_i, and an exact reduction leaves X unchanged.
%
% WHAT IS EXACT AND WHAT IS NOT. The decomposition is EXACT when the collapsed
% subnetwork is product-form, which is the condition aggregateFES already
% imposes. It is an approximation when it is not, and the state-space saving is
% the reason to accept that.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

subsetIdx = options.config.fes_stations(:)';
M = sn.nstations;
K = sn.nclasses;

if numel(subsetIdx) < 2
    line_error(mfilename, ['options.config.fes_stations must name at least two stations: ' ...
        'collapsing one station into a flow-equivalent server saves nothing.']);
end
if any(subsetIdx < 1) || any(subsetIdx > M) || numel(unique(subsetIdx)) ~= numel(subsetIdx)
    line_error(mfilename, sprintf(['options.config.fes_stations must be distinct station ' ...
        'indices in 1..%d.'], M));
end
if numel(subsetIdx) >= M
    line_error(mfilename, ['options.config.fes_stations names every station: there is no ' ...
        'complement left to solve.']);
end

stationSubset = cell(1, numel(subsetIdx));
for a = 1:numel(subsetIdx)
    stationSubset{a} = self.model.stations{subsetIdx(a)};
end
[fesModel, ~, info] = ModelAdapter.aggregateFES(self.model, stationSubset);

% The reduced solve. solver_ctmc_analyzer returns the generator and both state
% spaces alongside the metrics, so the stationary law below costs no second
% enumeration.
subopts = options;
subopts.config.fes_stations = [];
snRed = fesModel.getStruct();
[Qr,Ur,~,Tr,~,Xr,InfGen,SS,SSq,~,~,~,snRedCopy] = solver_ctmc_analyzer(snRed, subopts);
piRed = ctmc_stationary(InfGen, SS, snRedCopy, subopts);
piRed(piRed < GlobalConstants.Zero) = 0;

% P(N_fes = n): the aggregate state space carries K columns per stateful node,
% so the FES's block is the one at its stateful index.
fesIst = snRed.nodeToStation(info.fesNodeIdx);
fesIsf = snRed.nodeToStateful(info.fesNodeIdx);
cols = (fesIsf-1)*K + (1:K);
cutoffs = info.cutoffs;
Pn = zeros(1, prod(cutoffs+1));
for s = 1:size(SSq,1)
    nvec = SSq(s, cols);
    idx = ljd_linearize(nvec, cutoffs);
    Pn(idx) = Pn(idx) + piRed(s);
end

% The isolated subnetwork's metrics at every population it can hold.
[Qtab, Utab] = fes_compute_metrics(info.isolatedDemands, info.isolatedServers, ...
    info.isolatedIsDelay, cutoffs, struct('verbose', false));

QN = zeros(M,K); UN = zeros(M,K); RN = zeros(M,K); TN = zeros(M,K);

% Surviving stations: the reduced model lists them in the order aggregateFES
% walked the complement, before the FES it appends.
comp = info.complementIndices(:)';
for a = 1:numel(comp)
    QN(comp(a),:) = Qr(a,:);
    UN(comp(a),:) = Ur(a,:);
    TN(comp(a),:) = Tr(a,:);
end

% Collapsed stations: the conditional sum. The subnetwork's stations are listed
% in the order of info.subsetIndices, which is the order fes_build_isolated
% used, so row a of the isolated tables is station subsetIndices(a).
sub = info.subsetIndices(:)';
Qsub = zeros(numel(sub), K);
Usub = zeros(numel(sub), K);
for idx = 1:numel(Pn)
    if Pn(idx) <= 0
        continue
    end
    Qsub = Qsub + Pn(idx) * Qtab{idx};
    Usub = Usub + Pn(idx) * Utab{idx};
end
for a = 1:numel(sub)
    QN(sub(a),:) = Qsub(a,:);
    UN(sub(a),:) = Usub(a,:);
end

% Throughput needs no conditioning: flow through a station is fixed by the
% routing, so it is the chain throughput times the ORIGINAL visit ratio. The
% reduction leaves the chain throughput unchanged, and the FES carries exactly
% the flow that enters the subnetwork.
for c = 1:sn.nchains
    inchain = sn.inchain{c};
    for k = inchain(:)'
        for i = sub
            if sn.visits{c}(sn.stationToStateful(i), k) > 0
                TN(i,k) = Tr(fesIst, k) * ...
                    (sn.visits{c}(sn.stationToStateful(i), k) / ...
                     max(sn.visits{c}(sn.stationToStateful(sn.refstat(k)), k), GlobalConstants.Zero)) / ...
                    max(snRed.visits{c}(snRed.stationToStateful(fesIst), k) / ...
                        max(snRed.visits{c}(snRed.stationToStateful(snRed.refstat(k)), k), GlobalConstants.Zero), ...
                        GlobalConstants.Zero);
            end
        end
    end
end

RN = QN ./ TN;
RN(~isfinite(RN)) = 0;
XN = Xr;
CN = sum(RN, 1);
end
