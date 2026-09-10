function out = transformSolve(self, options)
% OUT = TRANSFORMSOLVE(OPTIONS)
%
% Solver-agnostic driver of a model TRANSFORMATION, and the sibling of
% FJFIXEDPOINT. Where fjFixedPoint drives one fixed transformation, this drives
% whichever one options.config.transform names: the strategy rewrites the model
% into subproblems, each subproblem is solved by a REAL solver, and the strategy
% maps the metrics back onto the original classes and stations.
%
%   expand -> [ for e = 1:nsub: solve(e); couple(e) ] -> converged -> lift
%
% SINGLE PASS BY DEFAULT. The strategy's expand phase returns a context; unless
% it sets ctx.iterated the loop runs one sweep and lifts, so the common case
% never pays for a convergence test it does not need.
%
% THE INNER SOLVER IS THE OUTER SOLVER. Each subproblem is solved by
% feval(class(self), submodel, innerOptions), the reflection already used by
% MAPENVAPPROX, so a transformation written once serves every NetworkSolver
% rather than the one it was first written for.
%
% COUPLING IS GAUSS-SEIDEL BY CONSTRUCTION. couple(e) runs immediately after
% solve(e) inside the sweep, so subproblem e+1 sees the updated state of
% 1..e. A strategy that wants Jacobi must accumulate in couple and apply in
% converged; nothing here offers a parallel sweep, deliberately.
%
% RECURSION is cut by options.config.transform_depth: the inner options carry
% transform='none' and depth+1, so a transformed submodel cannot re-enter the
% driver. A KERNEL selection inside the inner solver is not a transform and is
% correctly not cut: SolverNC picking an estimator for the aggregated model is
% the ordinary method choice, not a second transformation.
%
% OUT carries QN, UN, RN, TN, CN, XN, lG, runtime, iter and method, the
% fjFixedPoint contract minus actualmethod (the inner solvers report their own).
%
% See also TRANSFORM_METHOD, FJFIXEDPOINT, MAPENVAPPROX.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
requested = '';
if isfield(options.config,'transform') && ~isempty(options.config.transform)
    requested = options.config.transform;
end
[methodName, strategy] = transform_method(requested);
if strcmp(methodName,'none')
    line_error(mfilename, 'transformSolve called with options.config.transform=''none''.');
end

depth = 0;
if isfield(options.config,'transform_depth') && ~isempty(options.config.transform_depth)
    depth = options.config.transform_depth;
end
if depth > 0
    line_error(mfilename, sprintf(['a model transformation (''%s'') cannot be nested inside another ' ...
        'one; options.config.transform_depth is %d.'], methodName, depth));
end

sn = getStruct(self);
[submodels, ctx] = strategy(self, 'expand', sn, options);
if ~isfield(ctx,'iterated') || isempty(ctx.iterated)
    ctx.iterated = false;
end

% The inner solve must not re-enter this driver, and must not reprint the
% banner or rerun the feature gate on every sweep of an iterated strategy.
innerOptions = options;
innerOptions.config.transform = 'none';
innerOptions.config.transform_depth = depth + 1;

nsub = numel(submodels);
res = cell(1, nsub);
iter = 0;
for it = 1:max(1, options.iter_max)
    iter = it;
    for e = 1:nsub
        % A FRESH handle per sweep, deliberately. getAvg caches on the solver,
        % so an iterated strategy that re-solved a mutated model through the
        % same handle would read the previous sweep and converge to the wrong
        % point; constructing here makes that unrepresentable rather than
        % relying on a reset() call being remembered.
        inner = feval(class(self), submodels{e}, innerOptions);
        res{e} = collectInner(inner);
        if ctx.iterated
            [submodels, ctx] = strategy(self, 'couple', ctx, submodels, res, e);
        end
    end
    if ~ctx.iterated
        break
    end
    [done, ctx] = strategy(self, 'converged', ctx, res, it);
    if done
        break
    end
end

[QN,UN,RN,TN,CN,XN] = strategy(self, 'lift', ctx, res);

out = struct('QN',QN, 'UN',UN, 'RN',RN, 'TN',TN, 'CN',CN, 'XN',XN, ...
    'lG',res{1}.lG, 'runtime',toc(T0), 'iter',iter, 'method',methodName);
end

function r = collectInner(inner)
% The inner solve answers on FOUR channels, not one. getAvg returns
% (Q,U,R,T,A,W) whose sixth output is the RESIDENCE time, so the system
% throughput and the normalizing constant have to be asked for separately.
[QN,UN,RN,TN,AN,WN] = inner.getAvg();
r = struct('QN',QN, 'UN',UN, 'RN',RN, 'TN',TN, 'AN',AN, 'WN',WN);
r.XN = inner.getAvgSysTput();
r.lG = NaN;
if isprop(inner,'result') && isstruct(inner.result) && isfield(inner.result,'Prob') ...
        && isfield(inner.result.Prob,'logNormConstAggr')
    r.lG = inner.result.Prob.logNormConstAggr;
end
r.method = '';
if isprop(inner,'result') && isstruct(inner.result) && isfield(inner.result,'Avg') ...
        && isfield(inner.result.Avg,'method')
    r.method = inner.result.Avg.method;
end
end
