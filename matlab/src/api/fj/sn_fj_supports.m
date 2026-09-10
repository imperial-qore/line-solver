function [bool, reason] = sn_fj_supports(sn)
% [BOOL, REASON] = SN_FJ_SUPPORTS(SN)
%
% @brief Can the exact fork-join construction be asked for this model?
%
% The fork-join model class SN_FJ_VALIDATE admits, asked as a predicate rather
% than raised. SolverCTMC.supportsModelMethod and SolverSSA.supportsModelMethod
% call it so that a CALLER (model.help, findSolver, SolverAUTO) sees the verdict
% before paying for a run, and BOTH analyzers reach the SAME rules through
% MODELADAPTER.FJTAG -> SN_FJ_VALIDATE. The sentence the validator raises names
% 'the native CTMC/SSA fork-join implementation', which is why this predicate
% lives beside it rather than inside either solver.
%
% IT WRAPS THE VALIDATOR RATHER THAN RESTATING IT, and that is the point: the
% rules are eight and they move (pairing, join strategy, tasks-per-link, branch
% probability, open classes through a fork), so a second copy here would be a
% second thing to keep in step. There is exactly one body of rules and two ways
% in -- one that raises, for the run, and this one, which answers.
%
% WHAT IT REFUSES AND WHY THE ANALYZER IS RIGHT TO. The fork-join PAIRING is a
% declaration carried by the Join (`joinOf` here, `_fork` in native python, the
% third constructor argument in all four codebases), not a derivation from the
% routing: a nested model such as examples/basic/forkJoin/fj_basic_nesting has
% two forks and two joins whose pairing the routing alone does not determine.
% So a Join built without naming its fork leaves SN.FJ empty, and 'Fork nodes
% without a matched Join' is the honest answer to a model that declares none --
% not a topology test that failed to see one.
%
% @param sn NetworkStruct of the model
% @return bool true when the fork-join construction may run
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = true;
reason = '';
if ~any(sn.nodetype == NodeType.Fork) && ~any(sn.nodetype == NodeType.Join)
    return
end
try
    sn_fj_validate(sn);
catch ME
    bool = false;
    % LINE_ERROR stamps '[caller.m @ line N] ' on the front, which is a
    % diagnostic for a thrown error and noise in a report column; the sentence
    % after it is the rule.
    reason = regexprep(ME.message, '^\[[^\]]*\]\s*', '');
end
end
