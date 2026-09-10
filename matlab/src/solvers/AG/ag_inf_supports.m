function [bool, reason] = ag_inf_supports(sn)
% [BOOL, REASON] = AG_INF_SUPPORTS(SN)
%
% @brief Can the RCAT decomposition serve the infinite-server stations of
%        this model?
%
% EVERY RCAT COMPONENT IS A SINGLE-SERVER QBD: build_rcat (solver_ag) gives
% each (station, class) pair a level-INDEPENDENT completion block
% kron(I, D1^s) and never reads sn.sched or sn.nservers, so a Delay, or any
% station scheduled INF, is decomposed as an M/PH/1 queue. That is exact
% only when the component never holds more than one job, i.e. when the class
% is closed with a population of one: then the level is 0 or 1 and a
% single-server chain IS the infinite-server chain. With two jobs or more, or
% with an open class, the component serves at mu where the station serves at
% n*mu, and the answer is that of a different model.
%
% A self-looping class is exempt: it never leaves its reference station, and
% solver_ag reports it from sn.njobs directly (QN = N, UN = N at an INF
% station), which is exact.
%
% ONE PREDICATE, TWO CALLERS: SolverAG.supportsModelMethod answers with it and
% solver_ag raises with it, so model.help cannot offer what the run refuses.
% The registry cannot state the rule because 'Delay' and 'SchedStrategy_INF'
% ARE declared; what decides is the population behind them.
%
% @param sn NetworkStruct of the model
% @return bool true when every INF station may be decomposed
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = true;
reason = '';
isslc = false(1, sn.nclasses);
if isfield(sn, 'isslc') && ~isempty(sn.isslc)
    isslc = logical(sn.isslc(:)');
end
for ist = 1:sn.nstations
    if sn.sched(ist) ~= SchedStrategy.INF
        continue
    end
    for r = 1:sn.nclasses
        if isslc(r) || isnan(sn.rates(ist, r)) || sn.rates(ist, r) <= 0
            continue
        end
        if isinf(sn.njobs(r)) || sn.njobs(r) > 1
            bool = false;
            reason = sprintf(['SolverAG decomposes every station into a single-server component ' ...
                '(build_rcat never reads sn.nservers), so the infinite-server station %d is served ' ...
                'exactly only by a closed class of population one; class %d has population %g there. ' ...
                'Use SolverMVA, SolverNC or SolverCTMC.'], ist, r, sn.njobs(r));
            return
        end
    end
end
end
