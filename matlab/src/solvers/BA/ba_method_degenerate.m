function reason = ba_method_degenerate(sn, method)
% REASON = BA_METHOD_DEGENERATE(SN, METHOD)
%
% Whether METHOD APPLIES to SN but its bound carries no information there, and
% why. Empty when the bound is informative, and empty for every method that has
% no such regime.
%
% THIS IS A DIFFERENT QUESTION FROM BA_METHOD_REFUSAL, which is why it is a
% different function. That one answers "is this model outside the method's
% domain", and its answer is what the analyzer raises. This one answers "inside
% the domain, does the formula still say anything", and its answer is NOT
% raised: a degenerate bound is a VALID bound, just a vacuous one, so an
% analyzer asked for it by name is entitled to publish it. What must not happen
% is OFFERING it: findSolver and listValidMethods exist to name the pairs a
% caller can act on, and a table of zeros over a network with jobs circulating
% in it is not something anyone can act on.
%
% THE ONE METHOD WITH SUCH A REGIME IS 'ldbcmp.lower'. The Anselmi-Cremonesi
% bound is built from the population SURPLUS a = N - Qhat, where Qhat is the
% occupancy the non-bottleneck stations and the think time would hold in the
% open network fed at the bottleneck's saturation rate. pfqn_ldbcmp returns NaN
% below the regime (a < 0), which solver_ba_analyzer already refuses by name;
% AT the boundary a = 0 it returns Xlo = 0, which is formally the trivial bound
% X >= 0 and propagates into a table whose queue lengths, utilizations and
% throughputs are all zero. Every entry of that table is a true lower bound and
% none of them is usable, and a caller cannot tell it from a real answer of
% zero. So the bound is computed here and the name withheld when it degenerates.
%
% The cost is one pfqn_ldbcmp evaluation, a closed form plus a scalar fixed
% point, and only for the single method that has the regime -- every other name
% returns immediately.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if ~strcmp(ba_resolve_method(method), 'ldbcmp.lower')
    return
end
% The applicability rules come first and are not restated: a model this method
% is outside the domain of has no bound to be degenerate about.
if ~isempty(ba_method_refusal(sn, method))
    return
end

% TWO INDEX SPACES MEET HERE, and conflating them is what made this predicate
% misread a Petri net. sn.visits is STATEFUL-indexed -- sn_refresh_visits builds
% it as zeros(sn.nstateful,K) -- while sn.sched and sn.rates are STATION-indexed.
% Every station is stateful, but not every stateful node is a station: a
% Transition is stateful without being one, and a Place is both. On a fork-join
% SPN that is 4 stations against 7 stateful nodes. SN.STATIONTOSTATEFUL is the
% converter, and it is what sn_refresh_visits itself uses.
Vsf = sn.visits{1}(:);
V = Vsf(sn.stationToStateful(1:sn.nstations));
isinf_ = (sn.sched == SchedStrategy.INF);
Zt = sum(V(isinf_) ./ sn.rates(isinf_));
D  = V(~isinf_) ./ sn.rates(~isinf_);
N  = sn.nclosedjobs;
% NO QUEUEING STATION, NO BOUND. Every station is a delay (or a Place, on a
% Petri net, which is an INF station too), so there is no bottleneck to build
% Qhat on and pfqn_ldbcmp has no demand vector to read. A PREDICATE MUST NOT
% ERROR -- this one is asked once per name by listValidMethods, before the Petri
% sieve has had a chance to drop anything -- so the case is answered here.
if isempty(D)
    reason = ['Method ''ldbcmp.lower'' has no queueing station to bound here: every ' ...
        'station is an infinite server, so the bottleneck the open-network occupancy ' ...
        'is built on does not exist.'];
    return
end
[Xlo, ~, Qhat] = pfqn_ldbcmp(D, N, Zt, zeros(numel(D),1));

if ~isfinite(Qhat)
    reason = ['Method ''ldbcmp.lower'' does not apply here: the open-network occupancy ' ...
        'Qhat the bound is built from does not exist, because a non-bottleneck station ' ...
        'saturates at the bottleneck''s arrival rate.'];
    return
end
if ~isfinite(Xlo) || Xlo <= 0
    reason = sprintf(['Method ''ldbcmp.lower'' needs a population strictly above the ' ...
        'open-network occupancy the bound is built from (Qhat=%.4f, N=%d): with no ' ...
        'surplus it degenerates to the trivial bound X >= 0 and reports a table of ' ...
        'zeros.'], Qhat, N);
end
end
