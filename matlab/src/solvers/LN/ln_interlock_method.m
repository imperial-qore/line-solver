function [ilm, why] = ln_interlock_method(lqn, config, lnmethod)
% [ILM, WHY] = LN_INTERLOCK_METHOD(LQN, CONFIG, LNMETHOD)
% Interlock tracking method a run will actually use, asked as a query.
%
% ONE QUERY, TWO CALLERS, for the reason LN_METHOD_REFUSAL exists: BUILDLAYERS
% asks it on the run path and stores the answer in SolverLN.interlockMethod, and
% SolverLN.interlockMethodFor asks it before a solve, so a caller reading the
% option back can never be told something the builder did not do.
%
% IT DOWNGRADES, IT DOES NOT REFUSE. An encoding or a model feature refpath
% cannot carry is not a reason to decline the model: 'ilrate' analyses it exactly
% as it did before the option existed. That is why this is not a clause of
% LN_METHOD_REFUSAL, whose OK=false means the METHOD cannot encode the model at
% all -- 'srvn.ph' is fully supported on a model that merely cannot have refpath.
%
% CONFIG is options.config, read unvalidated: a user option struct reaches
% setOptions without passing through SolverOptions, so every key is isfield-
% guarded here rather than assumed present. LNMETHOD is the resolved encoding
% ('srvn.cs', 'srvn.ph', ...), which BUILDLAYERS has already settled by the time
% it asks.
%
% ILM is one of 'ilrate' (the rate-based discount of Franks' Eq. 4.7, the
% default and what LINE did before this option existed), 'refpath' (the merged
% reference-path chain in the client Delay) or 'none'. WHY is empty exactly when
% the request was honoured, and names the downgrade otherwise.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

why = '';
ilm = 'ilrate';
if nargin < 2 || isempty(config) || ~isstruct(config)
    config = struct();
end
if nargin < 3
    lnmethod = '';
end
if isfield(config,'interlock_method') && ~isempty(config.interlock_method)
    ilm = lower(char(config.interlock_method));
end
if ~any(strcmp(ilm, {'ilrate','refpath','none'}))
    line_error(mfilename, sprintf(['Unknown config.interlock_method ''%s'', use ' ...
        '''ilrate'', ''refpath'' or ''none''.'], ilm));
end

% The master switch wins over the tracking method: interlocking=false means no
% correction of any kind, which is what 'none' is.
if ~isfield(config,'interlocking') || ~config.interlocking
    ilm = 'none';
    return
end
if ~strcmp(ilm, 'refpath')
    return
end

% refpath writes class switches into the routing matrix of the client Delay,
% which is what the '.cs' encodings are. A composed phase-type entry law has no
% routing for them to be written into. 'flat.cs' is a .cs encoding but is
% refused too: under it every server is in one submodel and every non-REF task
% is a caller with a chain of its own, so the re-expansion would rewire a large
% fraction of the model at once and the hop stages would re-charge work that
% other chains in the same submodel already carry as call classes.
if ~strcmp(lnmethod, 'srvn.cs')
    ilm = 'ilrate';
    why = sprintf(['config.interlock_method=''refpath'' is carried by the ''srvn.cs'' ' ...
        'routing encoding only; ''%s'' falls back to ''ilrate''.'], lnmethod);
    return
end

% A phase-2 entry replies to its caller BEFORE phase 2 runs, so the residence the
% caller sees is servt_ph1 + prOvertake*servt_ph2, while the summed activity
% residence a hop stage charges counts phase 2 in full. Until the stage mean is
% taught the overtaking form, refpath must not be applied to such a model.
if isfield(lqn,'actphase') && any(full(lqn.actphase(:)) > 1)
    ilm = 'ilrate';
    why = ['config.interlock_method=''refpath'' does not yet carry second-phase ' ...
        'activities, whose caller-visible residence is servt_ph1 + prOvertake*servt_ph2 ' ...
        'and not the summed activity residence a hop stage charges; falling back to ''ilrate''.'];
    return
end

% A routed call group states the order in which one caller visits several
% callees. The re-expansion rewrites exactly those visits, so the two would be
% writing the same client routing.
if isfield(lqn,'callgroups') && ~isempty(lqn.callgroups)
    ilm = 'ilrate';
    why = ['config.interlock_method=''refpath'' does not carry routed call groups, ' ...
        'whose dispatch owns the client routing the re-expansion rewrites; ' ...
        'falling back to ''ilrate''.'];
    return
end
end
