function policy = nc_multiserver_policy(options)
% POLICY = NC_MULTISERVER_POLICY(OPTIONS)
%
% Resolves options.config.multiserver into the multiserver handling SolverNC
% implements. SolverNC represents a finite multiserver station in one of two
% ways: Seidmann's approximation (demand L/c plus a delay L(c-1)/c, applied in
% SOLVER_NC) or the exact load-dependent lattice mu(n)=min(n,c), which routes
% the model to SOLVER_NCLD. Returns one of:
%
%   'default'   the historical dispatch: Seidmann on method 'default', the
%               load-dependent lattice on method 'exact'/'is'/'panald', and
%               the lattice on the 2-station Delay+multiserver topology
%   'seidmann'  Seidmann everywhere, including on 'exact'
%   'lld'       the load-dependent lattice everywhere it is admissible,
%               including on method 'default'
%
% options.config.multiserver is a field of the GENERAL SolverOptions, shared
% with SolverMVA, which implements approximations SolverNC has no counterpart
% for ('softmin', 'conway', 'krzesinski', 'suri', 'erlang'). Those are not
% errors here -- one options struct is commonly reused across solvers -- but
% they are not silently honoured either: they warn and fall back to 'default'.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

policy = 'default';
if ~isfield(options,'config') || ~isstruct(options.config) ...
        || ~isfield(options.config,'multiserver') || isempty(options.config.multiserver)
    return
end

requested = lower(char(options.config.multiserver));
switch requested
    case {'default',''}
        policy = 'default';
    case 'seidmann'
        policy = 'seidmann';
    case {'lld','exact','loaddep','load-dependent'}
        policy = 'lld';
    otherwise
        % see _kb/06-solver-catalog.md (NC section, multiserver handling)
        line_warning(mfilename, ...
            ['SolverNC does not implement config.multiserver=''%s'' (it is a SolverMVA ' ...
             'approximation); using ''default''. SolverNC accepts ''default'', ''seidmann'' ' ...
             'and ''lld''.\n'], requested);
        policy = 'default';
end
end
