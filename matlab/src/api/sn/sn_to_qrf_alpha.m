%{ @file sn_to_qrf_alpha.m
 %  @brief Per-station load-dependent rate scaling alpha(i,n) for the QRF bounds
 %
 %  @author LINE Development Team
%}

%{
 % @brief Returns the QRF load-dependent scaling alpha, and why it may not exist
 %
 % @details
 % The load-dependent QRF arms carry a scaling alpha(i,n) that multiplies EVERY
 % rate out of station i while it holds n jobs, completions mu and background
 % phase changes v alike -- see the q construction in QRF_NOBLO_MMI_LD, which
 % writes r(i,j)*mu*alpha off the diagonal and (v + r(i,i)*mu)*alpha on it.
 % That is exactly the rate law of
 %
 %   an infinite server           alpha(i,n) = n
 %   a c-server station           alpha(i,n) = min(n, c_i)
 %   limited load dependence      alpha(i,n) = sn.lldscaling(i,n)
 %
 % so the three COMPOSE BY MULTIPLICATION and not one of them is an
 % approximation: the relaxed chain is the model's own, and the QRF answer
 % keeps whatever status it had on a single-server model.
 %
 % WHERE IT STOPS BEING THE MODEL'S OWN IS PHASE-TYPE SERVICE AT A STATION THAT
 % SERVES SEVERAL JOBS AT ONCE. The QRF local state carries ONE phase per
 % station, a faithful description of one job in service and of nothing else:
 % min(n,c) jobs served in parallel each advance through a phase of their own,
 % and no scaling of a single-phase process reproduces that joint motion.
 % A multiserver or delay station must therefore be exponential. Scaling a PH
 % server by min(n,c) would answer a DIFFERENT chain, so the relaxation would
 % no longer contain the model's stationary distribution and the number would
 % not bound anything -- which is the one failure this file exists to prevent.
 %
 % Limited load dependence at a SINGLE server is exempt and admits PH freely:
 % one job is in service whatever the rate, so alpha rescales that job's whole
 % phase process and the local state still describes it exactly.
 %
 % The scaling is tabulated at n = 1..N. sn.lldscaling may be narrower than N,
 % and is read clamped at its last column, the convention PFQN_LLDFUN uses.
 %
 % THE UTILIZATION NORMALIZER IS THE DECLARED PEAK, NOT max(alpha). LINE reports
 % U = T*S/peak at every station whose rate scales with the population, one
 % convention shared by multiserver, lld and class dependence (see
 % CD_PEAK_SCALING, "the same convention as solver_ncld does for lldscaling,
 % U/max(lldscaling)"). PEAK is therefore nservers(i) times the largest lld
 % scaling the model can REACH, and not max(alpha(i,:)): at c = 3 with N = 2 the
 % reachable alpha peaks at 2 while the station still has three servers, and
 % normalizing by 2 would report a utilization the model never attains. Inf at a
 % delay, where LINE reports U = QN instead.
 %
 % @par Syntax:
 % @code
 % [alpha, msg, ld, peak] = sn_to_qrf_alpha(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>alpha<td>(nstations x N) rate scaling at population n = 1..N
 % <tr><td>msg<td>Empty on success, otherwise why alpha is not defined for this model
 % <tr><td>ld<td>true when alpha is not identically 1, i.e. the model needs a load-dependent arm
 % <tr><td>peak<td>(nstations x 1) utilization normalizer nservers(i)*max reachable lld scaling; Inf at a delay
 % </table>
%}
function [alpha, msg, ld, peak] = sn_to_qrf_alpha(sn)
msg = '';
ld = false;
M = sn.nstations;
N = sum(sn.njobs);
peak = ones(M, 1);

if ~isfinite(N) || N < 1
    alpha = ones(M, 1);
    msg = 'the QRF bounds need a closed model with a finite population.';
    return
end

alpha = ones(M, N);
lld = sn.lldscaling;
smax = size(lld, 2);

for i = 1:M
    c = sn.nservers(i);
    isDelay = (sn.sched(i) == SchedStrategy.INF) || isinf(c);
    servesMany = isDelay || c > 1;
    ki = qrf_alpha_phases(sn, i);
    if servesMany && ki > 1
        % ld stays TRUE through the refusal: the model IS load dependent, and
        % the caller has to tell "no arm serves this" from "the arm you asked
        % for does not". Clearing it here would report the latter for both.
        ld = true;
        msg = sprintf(['station %d serves %s jobs at once with %d-phase service, and the QRF ' ...
            'local state carries one phase per station, which describes one job in service ' ...
            'and no more. Give that station exponential service, or use a single-server ' ...
            'model.'], i, qrf_alpha_many(isDelay, c), ki);
        return
    end
    lldpeak = 1;
    for n = 1:N
        if isDelay
            alpha(i, n) = n;
        elseif c > 1
            alpha(i, n) = min(n, c);
        end
        if smax > 0
            lldpeak = max(lldpeak, lld(i, min(n, smax)));
            alpha(i, n) = alpha(i, n) * lld(i, min(n, smax));
        end
    end
    if isDelay
        peak(i) = Inf;
    else
        peak(i) = c * lldpeak;
    end
end

ld = any(alpha(:) ~= 1);
end

function ki = qrf_alpha_phases(sn, i)
% Phases of station i's service process, read from sn.proc as the QRF adapter
% reads them: that {D0,D1} pair is what sizes the local state, so testing it
% rather than sn.phases keeps the refusal and the formulation on one quantity.
ki = 1;
if ~isempty(sn.proc) && ~isempty(sn.proc{i}) && ~isempty(sn.proc{i}{1})
    ki = size(sn.proc{i}{1}{1}, 1);
end
end

function s = qrf_alpha_many(isDelay, c)
if isDelay
    s = 'unboundedly many';
else
    s = sprintf('up to %d', c);
end
end
