%{ @file dtmc_solve.m
 %  @brief Equilibrium distribution of the discrete-time Markov chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Equilibrium distribution of the discrete-time Markov chain
 %
 % @details
 % Calculates the equilibrium distribution of a discrete-time Markov chain given its transition matrix.
 %
 % @par Syntax:
 % @code
 % PROB = dtmc_solve(P)
 % PROB = dtmc_solve(P, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Stochastic transition matrix of the discrete-time Markov chain
 % <tr><td>options<td>(Optional) Solver options
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>PROB<td>Equilibrium distribution vector
 % </table>
 %
 % @par Examples:
 % @code
 % PROB = dtmc_solve([0.5,0.5;0.2,0.8]);
 % @endcode
%}
function PROB=dtmc_solve(P,options)

% dtmc_solve is a pure function of P, and the layered fixed point asks it the
% same question over and over: a layer's routing does not change between
% SolverLN iterations, only its rates do, and visits do not depend on rates, so
% sn_refresh_visits re-solves one identical chain per layer per iteration.
% Native Python measured 2410 calls carrying ONE distinct matrix on lqn_ofbiz.
% The twins are Dtmc_solve.java, python/line_solver/api/mc/dtmc.py and
% cpp/include/line/api/mc/dtmc_solve.h; keep the four in step.
%
% The key is the matrix's exact BYTES plus its shape and sparsity, and a hit
% re-compares them rather than trusting a digest, so no collision can return
% another matrix's answer. Bytes also keep -0 apart from +0 and each NaN payload
% apart, which makes the cache conservative rather than clever: those are misses,
% never wrong hits. Only the no-options path is cached, which is the one
% sn_refresh_visits takes; a call carrying options recomputes.
persistent cacheKeys cacheVals
if isempty(cacheKeys)
    cacheKeys = {};
    cacheVals = {};
end
CACHE_MAX = 32;

if nargin<2
    key = [typecast(full(double(P(:))).','uint8'), ...
           typecast(uint32([size(P,1), size(P,2), uint32(issparse(P))]),'uint8')];
    for k = 1:numel(cacheKeys)
        if isequal(cacheKeys{k}, key)
            PROB = cacheVals{k};
            return
        end
    end
    if issparse(P)
        PROB=ctmc_solve(P-speye(size(P)));
    else
        PROB=ctmc_solve(P-eye(size(P)));
    end
    cacheKeys{end+1} = key;
    cacheVals{end+1} = PROB;
    if numel(cacheKeys) > CACHE_MAX
        cacheKeys(1) = [];
        cacheVals(1) = [];
    end
else
    if issparse(P)
        PROB=ctmc_solve(P-speye(size(P)),options);
    else
        PROB=ctmc_solve(P-eye(size(P)),options);
    end
end
end