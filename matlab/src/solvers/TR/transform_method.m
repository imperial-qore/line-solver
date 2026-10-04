function [token, strategy] = transform_method(method)
% [TOKEN, STRATEGY] = TRANSFORM_METHOD(METHOD)
%
% Normalise an options.config.transform method name onto one TRANSFORMSOLVE
% dispatches on, and return the strategy function that implements it. A
% transformation rewrites the model into one or more subproblems, solves those
% with a real solver, and maps the metrics back onto the original classes and
% stations.
%
%   'none'     no transformation; the caller solves its own model directly.
%   'chains'   collapse every chain onto a single class
%              (ModelAdapter.aggregateChains) and map the chain metrics back
%              through alpha (sn_deaggregate_chain_results). Single pass, and
%              EXACT on a product-form model.
%   'lc'       load concealment (Birman-Kogan Algorithm 2): chain aggregation,
%              then a Gauss-Seidel sweep in which chain l is solved on its own
%              against the residual capacity the other chains leave it. The
%              first ITERATED strategy. Note that SolverNC's 'lc' METHOD method name
%              is a different thing: it selects the pfqn_bklc KERNEL, which
%              stays the fast path and is untouched by this.
%
% An unrecognised token is an error rather than a silent 'none', because a
% mistyped transform that quietly solved the untransformed model would report a
% plausible number for the wrong problem.
%
% METHOD NAME is the canonical name, recorded on the result so the reported method
% names the transformation. STRATEGY is the phase-dispatched function handle,
% empty for 'none'.
%
% NOT every transformation in the tree is reachable here. The fork-join tag
% augmentation (solver_tr_fjtag_analyzer) is selected STRUCTURALLY, from the
% presence of a Fork or Join node, and never by a token; the fork-join fixed
% point keeps its own options.config.fork_join. See _kb/05-solvers-overview.md.
%
% See also TRANSFORMSOLVE, LQN_LN_METHOD.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(method)
    token = 'none';
    strategy = [];
    return
end
if ~ischar(method) && ~isstring(method)
    line_error(mfilename, 'options.config.transform must be a character token.');
end
switch lower(char(method))
    case {'none','off',''}
        token = 'none';
        strategy = [];
    case {'chains','chain','chainaggr','chain_aggregation'}
        token = 'chains';
        strategy = @solver_tr_chains_analyzer;
    case {'lc','loadconceal','load_concealment','thinning'}
        token = 'lc';
        strategy = @solver_tr_lc_analyzer;
    otherwise
        line_error(mfilename, sprintf(['options.config.transform=''%s'' is not a known model ' ...
            'transformation. Use ''none'', ''chains'' or ''lc''.'], char(method)));
end
end
