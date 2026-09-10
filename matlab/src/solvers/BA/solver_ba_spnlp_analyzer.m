function [Q,U,R,T,C,X,lG,runtime] = solver_ba_spnlp_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_BA_SPNLP_ANALYZER(SN, OPTIONS)
%
% Linear-programming bounds on the mean marking and the throughputs of a
% stochastic timed Petri net. The polytope and the LP are SPN_LPBND; this
% analyzer maps the LINE model onto them and reads one side of the bracket back
% per place and class.
%
% METHOD NAMES. Four, in two families:
%
%   spnlp.upper     Markovian LP, upper side          exponential firing
%   spnlp.lower     Markovian LP, lower side          exponential firing
%   spnlp.op.upper  operational LP, upper side        any phase-type firing
%   spnlp.op.lower  operational LP, lower side        any phase-type firing
%
% The operational variant drops the second-moment, covariance and Little's-law
% families and the whole E[X_p e_t] block with them, which is what removes the
% exponential requirement. It is much looser, and is the reference's own
% "without Markovian assumption" column.
%
% BOUND CONVENTION. Q(i,r) is the reported side of the bracket on the mean
% number of class-r tokens in place i. T(i,r) is the same side of the bracket
% on the token throughput of that place, and R follows by Little's law from the
% two. U(i,r) = Q(i,r) DELIBERATELY: a Place is an INF station and LINE reports
% U = Q at an infinite server, which is what SolverCTMC and
% SOLVER_NC_SPN_ANALYZER both do on the same net. The reference's place
% utilization 1 - P(m = 0) is a different quantity and is not this column.
%
% A TRANSITION GETS NO ROW. It is a StatefulNode and not a Station, so it has
% no station index; the mode throughputs and enabling probabilities the LP also
% brackets stay inside SPN_LPBND's return value, the same way SPN_METRICS keeps
% modeTput and modeUtil off the table.
%
% THE TWO SIDES ARE NOT ONE SOLVE. Each reported cell is its own linear form,
% minimised and maximised over the same polytope, so a '.upper' run and a
% '.lower' run cost the same and neither is derivable from the other.
%
% HOW TIGHT. The reference's own Table 2 measures it on a four-server
% production line: the upper side lands 2% to 11% above simulation and the
% lower side 30% to 40% below it, both comfortably inside the operational
% bounds it also reports. Expect a usable upper bound and a weak lower one.
%
% Reference: Z. Liu (1998). Performance analysis of stochastic timed Petri nets
% using linear programming approach. IEEE Transactions on Software Engineering
% 24(11), 1014-1030.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
M = sn.nstations;
K = sn.nclasses;
Q = zeros(M,K); U = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
C = zeros(1,K); X = zeros(1,K); lG = NaN;

method = options.method;
switch method
    case 'spnlp.upper',    markovian = true;  side = 2;
    case 'spnlp.lower',    markovian = true;  side = 1;
    case 'spnlp.op.upper', markovian = false; side = 2;
    case 'spnlp.op.lower', markovian = false; side = 1;
    otherwise
        line_error(mfilename, ['Unknown SPN bound method ''%s''. Valid: spnlp.upper, ' ...
            'spnlp.lower, spnlp.op.upper, spnlp.op.lower.'], method);
end

% ----- model gates -----
% BA_SPNLP_REFUSAL is the predicate SolverBA.supportsModelMethod reports, so
% the run raises the sentence the report gave: a Petri net at all, no queueing
% place, and every firing mode timed with one finite constant rate (phase-type
% admitted by the operational pair alone). SPN_LPBND keeps its own copies of
% the mode tests as the library's guard.
reason = ba_spnlp_refusal(sn, method);
if ~isempty(reason)
    line_error(mfilename, '%s', reason);
end

lpopt = struct();
lpopt.markovian = markovian;
lpopt.verbose = options.verbose > VerboseLevel.STD;
if isfield(options,'config') && isfield(options.config,'spnlp_assumelive') && ...
        ~isempty(options.config.spnlp_assumelive)
    % The reference's liveness rows hold only on a live net, which SPN_LPBND
    % cannot certify, so they are opt-in. See its header.
    lpopt.assumelive = options.config.spnlp_assumelive;
end
if isfield(options,'config') && isfield(options.config,'spnlp_init') && ...
        ~isempty(options.config.spnlp_init)
    lpopt.init = options.config.spnlp_init;
end

bnd = spn_lpbnd(sn, lpopt);

% ----- read the reported side back per place and class -----
Rn = bnd.nclasses;
for pp = 1:numel(bnd.places)
    ist = sn.nodeToStation(bnd.places(pp));
    if ist < 1, continue; end
    for k = 1:Rn
        l = (pp - 1) * Rn + k;
        Q(ist,k) = bnd.tokens(side, l);
        U(ist,k) = Q(ist,k);
        T(ist,k) = bnd.placeTput(side, l);
        if T(ist,k) > 0
            R(ist,k) = Q(ist,k) / T(ist,k);
        end
    end
end

for k = 1:K
    ref = sn.refstat(k);
    if ref >= 1 && ref <= M
        X(k) = T(ref,k);
    end
    Nk = sum(Q(:,k));
    if X(k) > 0 && Nk > 0
        C(k) = Nk / X(k);
    end
end

line_debug(options, sprintf(['SPNLP bound: %s, %d place levels, %d modes, %d variables, ' ...
    '%d rows'], method, bnd.nplacelevels, numel(bnd.modes), bnd.nvars, bnd.nrows));
runtime = toc(T0);
end
