function [QN,UN,RN,TN,CN,XN,lG,runtime,method,pf] = solver_nc_spn_analyzer(model, sn, options)
% [QN,UN,RN,TN,CN,XN,LG,RUNTIME,METHOD,PF] = SOLVER_NC_SPN_ANALYZER(MODEL, SN, OPTIONS)
% Stationary analysis of a PRODUCT-FORM stochastic Petri net by MDD-rec: the
% normalising constant is obtained from one memoised walk of the decision
% diagram holding the reachable set, and every reported measure is a masked
% walk of the same diagram.
%
% This is the 'rec' method of SolverNC, and the first analytical route LINE
% offers for a Petri net -- CTMC solves the explicit generator, SSA and LDES
% simulate, FLD fluidises. Three functions do the work and each is the subject
% of its own reference:
%
%   SPN_PF      decides the product form and derives the per-place factors g_l
%               (Coleman-Henderson-Taylor complex balance)
%   MDD_REC     G = sum_S prod_l g_l(s_l) in O(sum_l nodes_l * |S_l|) rather
%               than O(|S|)                              (Balsamo-Marin-Stojic)
%   SPN_METRICS mean tokens, place and mode utilisation, and throughputs, all
%               from masked walks of the same diagram    (Balsamo-Marin-Stojic)
%
% WHAT THIS REACHES THAT THE EXPLICIT GENERATOR DOES NOT. The diagram stores the
% reachable set, never the generator, so the cost is set by the number of
% diagram nodes and not by |S|. It also does not need the marking to be a
% conserved job population: a mode may consume two tokens and produce one, or
% consume one and produce two, which is the fork-join and batch case that the
% MDD-rec paper exists to serve.
%
% -- Output
% QN,UN,RN,TN : per (Place station, class): mean tokens, utilisation, Little
%               response time and token throughput
%
% UN FOLLOWS LINE, NOT THE PAPER. A Place is an INF station, and LINE reports
% U = Q at an infinite server, which is what SolverCTMC returns for the same
% net. The paper's place utilisation u(P_j) = 1 - P(m_j = 0) is a different
% quantity and is reported separately, as PF.metrics.placeUtil.
% CN,XN       : per class, system response time and throughput at the reference
% LG          : log of the normalising constant
% PF          : the SPN_PF certificate (kind, complexes, deficiency, y, ...)
%
% See also: spn_pf, mdd_rec, spn_metrics, spn_mdd, SolverNC.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
method = 'rec';

pfopt = struct('verbose', options.verbose > 1);
if isfield(options, 'config') && isstruct(options.config) && isfield(options.config, 'spn_bound')
    pfopt.bound = options.config.spn_bound;
end
if isfield(options, 'tol') && ~isempty(options.tol) && options.tol > 0
    pfopt.tol = max(options.tol, 1e-12);
end
pf = spn_pf(model, pfopt);
met = spn_metrics(pf.mdds, pf.g, pf.info);

M = sn.nstations;
R = sn.nclasses;
QN = zeros(M, R); UN = zeros(M, R); RN = zeros(M, R); TN = zeros(M, R);

places = pf.info.places;
for pp = 1:numel(places)
    ist = sn.nodeToStation(places(pp));
    if ist < 1, continue; end
    for k = 1:R
        l = (pp - 1) * R + k;
        QN(ist, k) = met.tokens(l);
        % INF station: LINE charges one server per resident token, so U = Q.
        % The paper's 1 - P(m = 0) is MET.placeUtil, carried on the certificate.
        UN(ist, k) = met.tokens(l);
        TN(ist, k) = met.placeTput(l);
        if TN(ist, k) > 0
            RN(ist, k) = QN(ist, k) / TN(ist, k);       % Little's law at the place
        end
    end
end

% System throughput at the reference station of each class, per unit visit, and
% the response time that Little's law then fixes. A Petri net whose class
% population is not conserved has no meaningful N/X, so CN is left at zero
% there rather than reporting a ratio against a moving population.
XN = zeros(1, R); CN = zeros(1, R);
for k = 1:R
    ref = sn.refstat(k);
    if ref >= 1 && ref <= M
        XN(k) = TN(ref, k);
    end
    Nk = sum(QN(:, k));
    if XN(k) > 0 && Nk > 0
        CN(k) = Nk / XN(k);
    end
end

pf.metrics = met;
lG = log(met.G);
runtime = toc(Tstart);
end
