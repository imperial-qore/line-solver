function [ok, reason] = qns_multiserver_refusal(sn, method)
% [OK, REASON] = QNS_MULTISERVER_REFUSAL(SN, METHOD)
% Whether qnsolver's own -m switch offers this multiserver approximation.
%
% ASKED ON THE MULTISERVER APPROXIMATION, i.e. the name after 'qns.' that
% SolverLQNS.qnsMultiserver extracts, or options.config.multiserver at run time.
% runAnalyzerNetwork sets the config value from the method name before
% SOLVER_QNS runs, so SolverLQNS(model,'qns.suri') reaches qnsolver asking for
% 'suri', which qnsolver has no flag for. SOLVER_QNS used to fall off its inner
% switch there, leaving CMD unassigned, so the failure surfaced as an
% undefined-variable error about a temporary rather than as a diagnosis.
%
% 'qnsolver -m' accepts conway, reiser, rolia and zhou. 'suri' and 'schmidt' are
% LQNS approximations, reachable only on the non-product-form closed branch of
% SolverLQNS, where QN2LQN hands the model to lqns, and qnsolver has no flag for either.
%
% THE RULE IS INSIDE THE MULTISERVER BRANCH, and that is not a detail. Without a
% multiserver station the reference emits no -m at all and answers under the
% caller's method name, so refusing 'suri' there would refuse a model this
% solver does solve. With one, SOLVER_QNS used to fall off its inner switch and
% leave CMD unassigned, so the failure surfaced as an undefined-variable error
% about a temporary rather than as a diagnosis naming the method. The C++ port
% has diagnosed this since it was written (solver_qns.h, is_qnsolver_multiserver);
% this is the same rule, stated where MATLAB, the JAR and python can all reach it.

ok = true;
reason = '';
if nargin < 2 || isempty(method)
    return
end
if ~any(sn.nservers > 1 & sn.nservers < Inf)
    % No multiserver station, so no -m flag is emitted and every method name is
    % served by the plain invocation.
    return
end
ms = lower(char(method));
if any(strcmp(ms, {'default','conway','reiser','rolia','zhou'}))
    return
end
ok = false;
reason = sprintf(['SolverLQNS: the multiserver approximation ''%s'' is one LQNS offers and ' ...
    'qnsolver does not: ''qnsolver -m'' accepts conway, reiser, rolia and zhou only; ' ...
    'qns.suri and qns.schmidt are available only on a closed non-product-form Network, ' ...
    'which SolverLQNS solves with lqns.'], ms);
end
