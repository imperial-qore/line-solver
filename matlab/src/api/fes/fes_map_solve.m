%{ @file fes_map_solve.m
 %  @brief Solves the closed model left by a MAP flow-equivalent server
 %
 %  @author LINE Development Team
%}

%{
 % @brief Solves the reduced model made of a delay and a load-dependent MAP
 % flow-equivalent server
 %
 % @details
 % Closes the aggregation of Section 5.2.1 of Casale, Mi, Cherkasova and
 % Smirni, IEEE Trans. Soft. Eng. 37(5), 2011. Once a subnetwork has been
 % replaced by the load-dependent MAP of fes_map_aggregate, the model left
 % is a delay holding the think times and one station, which is a finite
 % level-dependent quasi birth-death process: level k is the number of jobs
 % held by the flow-equivalent server and N-k jobs are thinking. The chain
 % is the same block bidiagonal pair used to measure the inter-departure
 % times, now read as a generator rather than as a MAP, so the delay is a
 % station whose process is scaled by the number of jobs it holds and the
 % marked transitions are the arrivals into the flow-equivalent server.
 %
 % The think time may itself be a MAP, which is how Section 5.3.1 models a
 % bounded flash crowd: burstiness in the stream of requests is carried by
 % (Z0,Z1) and the rates scale with the population at the delay.
 %
 % @par Syntax:
 % @code
 % [XN,RN,QN,pk] = fes_map_solve(FES, thinkMAP, N)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>FES<td>Flow-equivalent server, one MAP per level, from fes_map_aggregate
 % <tr><td>thinkMAP<td>Think time process {Z0,Z1}; use map_exponential(Z) for
 %                     an exponential think time and a MAP for a flash crowd
 % <tr><td>N<td>Number of jobs in the closed model
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>XN<td>System throughput
 % <tr><td>RN<td>Mean response time of the aggregated subnetwork, N/XN-E[Z]
 % <tr><td>QN<td>Mean number of jobs held by the flow-equivalent server
 % <tr><td>pk<td>Distribution of the jobs held by the flow-equivalent server,
 %               pk(k+1) = P(k jobs), k = 0..N
 % </table>
 %
 % @see fes_map_aggregate, fes_map_deaggregate, fes_map_interdeparture
%}
function [XN,RN,QN,pk] = fes_map_solve(FES, thinkMAP, N)

if N < 1
    line_error(mfilename,'The population N must be at least 1.');
end

FESlev = fes_map_levels(FES, N);
[T0,T1] = fes_map_interdeparture(thinkMAP, FESlev, N, [Inf 1]);
Q = T0 + T1;
dim = size(Q,1);

phi = ctmc_solve(Q);
phi = reshape(phi, 1, dim);
XN = phi*T1*ones(dim,1);

blk = dim/(N+1);
pk = zeros(1,N+1);
for k = 0:N
    pk(k+1) = sum(phi((k*blk+1):((k+1)*blk)));
end

QN = (0:N)*pk';
RN = N/XN - map_moment(thinkMAP,1);
end
