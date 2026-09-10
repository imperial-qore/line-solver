function [ok, reason] = qns_multiserver_refusal(sn, method)
% [OK, REASON] = QNS_MULTISERVER_REFUSAL(SN, METHOD)
% Whether qnsolver's own -m switch offers this multiserver approximation.
%
% ASKED ON THE CONFIG VALUE, NOT ON THE METHOD NAME, and the distinction is the
% whole reason this is not a gate. The outer switch in SOLVER_QNS sets
% options.config.multiserver only for conway/reiser/rolia/zhou, so
% SolverQNS(model,'suri') leaves it at 'default', takes the {'default','conway'}
% arm and answers with Conway -- reporting 'conway' as the actual method. That
% is an honest substitution and the row is genuinely runnable, so
% @SolverQNS/supportsModelMethod deliberately does NOT withdraw it.
%
% What this refuses is the direct options.config.multiserver='suri', which no
% arm matches: CMD was never assigned and the failure surfaced as an
% undefined-variable error about a temporary rather than as a diagnosis. A
% condition on the OPTIONS is not something a feature set or a model gate can
% state, which is why it lives here.
%
% 'qnsolver -m' accepts conway, reiser, rolia and zhou. 'suri' and 'schmidt' are
% LQNS approximations, reachable only on the non-product-form closed SolverLQNS
% branch, and qnsolver has no flag for either.
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
reason = sprintf(['SolverQNS: the multiserver approximation ''%s'' is one LQNS offers and ' ...
    'qnsolver does not: ''qnsolver -m'' accepts conway, reiser, rolia and zhou only; ' ...
    'suri and schmidt are available only on the non-product-form closed SolverLQNS branch.'], ms);
end
