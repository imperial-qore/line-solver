function [bool, reason] = solver_ctmc_state_supports(sn, solverName)
% [BOOL, REASON] = SOLVER_CTMC_STATE_SUPPORTS(SN, SOLVERNAME)
%
% @brief Can the State machinery (lang/+State) serve the impatience laws
%        this model declares?
%
% The rules of State.afterEventStation that a feature name cannot state,
% asked as a predicate rather than raised. The explicit-generator methods of
% SolverCTMC and BOTH engines of SolverSSA drive the same afterEvent code
% (the serial engine directly, the NRM through its own propensities, which
% read sn.retrialMu as a memoryless rate), so each solver used to carry a
% copy of these tests at the top of its analyzer (solver_ctmc.m,
% solver_ssa.m) and nothing above them asked: model.help offered
% ctmc.default and every ssa.* row on a model with a phase-type patience,
% and each then raised. ONE BODY, FOUR CALLERS: the two analyzers, which
% raise, and the two supportsModelMethod gates, which answer.
%
% What it refuses, and why the registry cannot: Reneging, Balking and Retrial
% are DECLARED by both solvers, and the rules are about the LAW behind them,
% i.e. reneging with a non-exponential patience, balking outside
% QUEUE_LENGTH, a non-exponential retrial delay, a finite maxAttempts, a
% retrial station serving more than one class.
%
% The multi-server DPS/GPS station afterEventStation also refuses is NOT a
% rule here: Queue.setNumServers refuses it at construction, so no built
% model carries one and the analyzer's own raise is the only guard it needs.
%
% @param sn NetworkStruct of the model
% @param solverName 'SolverCTMC' or 'SolverSSA', named in the refusal
% @return bool true when the model may run
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = false;
if nargin < 2 || isempty(solverName)
    solverName = 'SolverCTMC';
end
% see _kb/06-solver-catalog.md (CTMC section, support gates) for rationale
if isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass)
    badRenege = (sn.impatienceClass==ImpatienceType.RENEGING) & (sn.impatienceType~=ProcessType.EXP);
    if any(badRenege(:))
        reason = sprintf('%s supports only exponential (memoryless) patience for reneging. Use SolverLDES or SolverJMT for phase-type patience.', solverName);
        return
    end
end
if isfield(sn,'balkingStrategy') && ~isempty(sn.balkingStrategy)
    badBalk = (sn.balkingStrategy~=0) & (sn.balkingStrategy~=BalkingStrategy.QUEUE_LENGTH);
    if any(badBalk(:))
        reason = sprintf('%s supports only QUEUE_LENGTH balking. Use SolverLDES or SolverJMT for wait-time-based balking.', solverName);
        return
    end
end
if isfield(sn,'retrialProc') && ~isempty(sn.retrialProc)
    hasRetrial = ~cellfun(@isempty, sn.retrialProc);
    if any(hasRetrial(:))
        if any(hasRetrial(:) & (sn.retrialType(:)~=ProcessType.EXP))
            reason = sprintf('%s supports only exponential (memoryless) retrial delay. Use SolverLDES or SolverMAM for phase-type retrials.', solverName);
            return
        end
        if any(hasRetrial(:) & (sn.retrialMaxAttempts(:)>=0))
            reason = sprintf('%s supports only unlimited retrials (maxAttempts=-1). Use SolverLDES for finite max-attempts.', solverName);
            return
        end
        % A retrial orbit is enumerated per single populated class; reject a
        % retrial station that serves more than one class.
        for ii = find(any(hasRetrial,2))'
            served = 0;
            for rr = 1:sn.nclasses
                if ~isempty(sn.proc{ii}{rr}) && ~any(any(isnan(sn.proc{ii}{rr}{1})))
                    served = served + 1;
                end
            end
            if served > 1
                reason = sprintf('%s supports retrial only for single-class stations. Use SolverLDES for multi-class retrial.', solverName);
                return
            end
        end
    end
end
bool = true;
reason = '';
end
