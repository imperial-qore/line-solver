function percentileRT_K = mainFJ( arrival, service, pers, K, Cs, T_mode )

% If you use the scripts available here in your work, please cite our paper
% entitled 'Beyond the Mean in Fork-Join Queues: Efficient Approximation
% for Response-Time Tails' (IFIP Performance 2015).
%
% This function implements the approximation method proposed in the paper.
% It returns the approximated values of the RT percentiles for a K-node FJ queue
% according to Section 6 of the paper.
% The RT percentiles for 1-node queue is exact, while the results for a
% 2-node FJ queue are based on the approximation proposed in Section 4,
%
% % arrival is a structure that must contain arrival.lambda, arrival.lambda0
% and arrival.lambda1, where lambda is the mean arrival rate, lambda0 is
% the intensity-matrix for state changes in the arrival process not
% accompanied by arrivals and lambda1 is the intensity-matrix for state
% changes accompanied by arrivals.
%
% % service is a structure that must contain service.mu, service.ST,
% service.tau_st, where mu is the mean service rate, and tau_st and ST
% corresponds to the PH representation (tau_st,ST) for the service time
% of each of the homogeneous K servers.
%
% % pers are the targeted percentiles, e.g., 90-th, 95-th and 99-th.
%
% % K is the number of nodes in the FJ queue, e.g., K=16 represents a
% 16-node FJ queue.
%
% % Cs: limits of the value C as introduced in Section 4.
% The larger the C, the more accurate the results.
%
% % This function returns a cell, in which each element is a structure
% with the following properties:
% K: number of nodes, i.e., K-node FJ queue
% percentiles: the targeted percentiles, e.g., 95-th percentile
% RTp: response time percentiles according to percentiles
%
% Script example.m shows an example of how to use this function.
%
%
% % Two methods are provided to solve the T matrix: (default: 'NARE')
%   'Sylves' : solves a Sylvester matrix equation at each step
%                          using a Hessenberg algorithm
%   'NARE  ' : solves the Sylvester matrix equation by defining 
%            a non-symmetric algebraic Riccati equation (NARE) 
%
%


% check utilization
load = arrival.lambda/service.mu;
if (load >= 1)
    error('System not stable: mean arrival rate %d > mean service rate %d',arrival.lambda,service.mu);
end

for c = 1 : length(Cs)
    
    C = Cs(c);
    
    % the RT percentiles for 1-node queue
    percentileRT_1  = returnRT1(arrival, service, pers);
    
    % the RT percentiles for 2-node FJ queues
    percentileRT_2 = returnRT2(arrival, service, pers, C, T_mode);
    
end

% predict the response time percentiles based on the results of the 1-node
% and 2-node FJ queue with the same setting
percentileRT_K = cell(1,length(K));
for k = 1 : length(K)
    percentileRT_K{k}.K = K(k);
    percentileRT_K{k}.percentiles = 100*pers;
    percentileRT_K{k}.RTp = zeros(1,length(pers));
    for p = 1 : length(pers)
        percentileRT_K{k}.RTp(p) = percentileRT_1(p,2)+(percentileRT_2(p,2)-percentileRT_1(p,2))*log(K(k))/log(2);
    end
end

end


