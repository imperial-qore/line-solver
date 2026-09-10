function requested = lqn_ln_method(method)
% REQUESTED = LQN_LN_METHOD(METHOD)
%
% Normalise a SolverLN method name onto one the solver dispatches on. A method
% name here carries TWO decisions: the LAYERING, which fixes what a submodel is,
% and the ENCODING, which fixes how an activity graph is written into it.
%
%   'srvn.cs'  srvn layering, activity graph as ROUTING: one class per task,
%              entry, activity and call, plus the Fork/Join/Router nodes and the
%              ClassSwitch that an off-diagonal P{r,s} mints. Serves every
%              layered feature this solver implements.
%   'srvn.ph'  srvn layering, activity graph as a composed PHASE-TYPE server
%              law, so a layer is a two-station cycle with one class per caller
%              task. Far smaller and faster, but it cannot represent second
%              phases, forwarding, cache tasks, quorum joins, routed call groups
%              or per-entry admission constraints.
%   'srvn'     the ALIAS, and the default: 'srvn.ph' where it can serve the
%              model, 'srvn.cs' otherwise. The choice is made once, at layer
%              build time, and reported in SolverLN.lnmethod.
%   'flat.cs'  the squashed layering, in which every server is a station of ONE
%              submodel, again with the routing encoding.
%   'flat.ph'  the squashed layering with the PHASE-TYPE encoding: the same
%              single submodel, but a caller visits each server once per
%              invocation under the composed law rather than once per activity
%              or call under a class chain. Squashing does not conflict with the
%              composition, because a task station's service law is ALREADY the
%              inflated entry law that the outer fixed point updates, exactly as
%              lqns does under --squashed-layering. It cannot represent what
%              'srvn.ph' cannot, and additionally loses the routed call groups,
%              whose dispatch only the routing encoding carries.
%   'flat'     the ALIAS for 'flat.cs'. It resolves unconditionally rather than
%              probing 'flat.ph', because a model is squashed in order to express
%              what only the routing encoding carries.
%   'moment3'  the three-moment response-time distribution pass, which builds
%              the routing layers and adds a distribution propagation.
%
% 'default' is the srvn alias, so a model solved without naming a method takes
% the better of the two srvn encodings rather than always the routing one. An
% unrecognised method name takes 'srvn.cs', as every name other than 'moment3' did
% before the alias existed.
%
% A method that names a layering SETS it: 'flat'/'flat.cs' squash the model even
% when options.config.layering is unset. Conversely options.config.layering is
% still honoured on its own, and SolverLN.lnmethod then reports what was built,
% so config.layering='flat' with the default method reports 'flat.cs'.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(method) || ~ischar(method)
    requested = 'srvn';
    return
end
switch lower(method)
    case {'srvn.ph','ph'}
        requested = 'srvn.ph';
    case {'srvn.cs','srvncs','cs'}
        requested = 'srvn.cs';
    case {'srvn','default','auto',''}
        requested = 'srvn';
    case {'flat.cs','flatcs','flat','squashed'}
        requested = 'flat.cs';
    case {'flat.ph','flatph','squashed.ph'}
        requested = 'flat.ph';
    case 'moment3'
        requested = 'moment3';
    otherwise
        % An unrecognised method name takes the routing encoding, which is what every
        % method name other than 'moment3' resolved to before the alias existed.
        requested = 'srvn.cs';
end
end
