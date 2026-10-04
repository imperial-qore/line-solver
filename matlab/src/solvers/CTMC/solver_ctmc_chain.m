function [pi, infGen, stateSpace, runtime] = solver_ctmc_chain(chain, options)
% [PI, INFGEN, STATESPACE, RUNTIME] = SOLVER_CTMC_CHAIN(CHAIN, OPTIONS)
%
% Steady-state analysis of a user-supplied Markov chain, i.e. a MarkovProcess
% (CTMC, generator Q) or a MarkovChain (DTMC, transition matrix P). INFGEN is
% Q for a CTMC and the uniformized generator P-I for a DTMC, which carries the
% same stationary vector. STATESPACE is the state space attached to the chain,
% or the state indices when the chain carries none.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = SolverCTMC.defaultOptions();
end

T0 = tic;
if isa(chain,'MarkovChain')
    P = chain.getTransMat();
    n = size(P,1);
    if issym(P)
        pi = dtmc_solve(P);
    else
        pi = dtmc_solve(P);
        if ~solver_ctmc_chain_isvalid(pi, n)
            pi = dtmc_solve_reducible(P);
        end
    end
    if issparse(P)
        infGen = P - speye(n);
    else
        infGen = P - eye(n);
    end
elseif isa(chain,'MarkovProcess')
    infGen = chain.getGenerator();
    n = size(infGen,1);
    if issym(infGen)
        pi = ctmc_solve(infGen);
    else
        pi = ctmc_solve(infGen, options);
        if ~solver_ctmc_chain_isvalid(pi, n)
            pi = ctmc_solve_reducible(infGen);
        end
    end
else
    line_error(mfilename,'solver_ctmc_chain requires a MarkovProcess or a MarkovChain.');
end

pi = reshape(pi, 1, n);
stateSpace = chain.stateSpace;
if isempty(stateSpace)
    stateSpace = (1:n)';
end
runtime = toc(T0);
end

function bool = solver_ctmc_chain_isvalid(pi, n)
% BOOL = SOLVER_CTMC_CHAIN_ISVALID(PI, N)
% Reject a solution that the primary solver could not produce on a reducible
% chain, so that the reducible fallback is used instead.
bool = numel(pi) == n && all(isfinite(pi(:))) && ...
    all(pi(:) >= -GlobalConstants.FineTol) && ...
    abs(sum(pi(:)) - 1) <= sqrt(GlobalConstants.FineTol);
end
