% This script provides an example use of the returnPercentiles.m function,
% which implements the approximation method proposed in paper
% 'Beyond the Mean in Fork-Join Queues: Efficient Approximation for
% Response-Time Tails', which is accepted in IFIP Performance 2015.
% 
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
% % pers are the targeted percentiles, e.g., 0.90 = 90-th, 0.95 = 95-th.
%
% % K is the number of nodes in the FJ queue, e.g., K=16 represents a
% 16-node FJ queue.
%
% % Cs: limits of the value C as introduced in Section 4.
% The larger the C, the more accurate the results.
%
%
% % Two methods are provided to solve the T matrix: (default: 'NARE')
%   'Sylves' : solves a Sylvester matrix equation at each step
%                          using a Hessenberg algorithm
%   'NARE  ' : solves the Sylvester matrix equation by defining 
%            a non-symmetric algebraic Riccati equation (NARE) 
%
%

clear
clc

%% Arrival
% Four possibilities of arrival process are provided
% 1=Exp; 2=HE2 (2-phase Hyper-exponential);  
% 3=ER2 (2-phase Erlang); 4 = MAP2 (2-phase MAP)
ArrChoice = 2;

arrival.lambda = 0.5; % mean arrival rate
CX2 = 10; % squared coefficient of variation of the arrival process
decay = 0.5; % decay rate of the auto-correlation function, only applies to MAP arrivals
[arrival.lambda0, arrival.lambda1] = get_distribution(arrival.lambda, ArrChoice, CX2, decay);
arrival.ma = size(arrival.lambda0,2);
arrival.Ia = eye(arrival.ma);

%% Service
% Three possibilities of arrival process are provided
service.SerChoice = 1; % 1: Exp;  2: HE2;  3: ER2;

service.mu = 1; % mean service rate
service_CX2 = 10; % squared coefficient of variation of the services
[service.ST, service.St, service.tau_st] = get_serviceDis(service.mu, service.SerChoice, service_CX2);

%% Targeted percentiles
% 0.95 = 95-th percentile of the response times
pers = [ 0.90, 0.95, 0.99 ];

%% System settings
K = [10, 512, 1024]; % the targeted K-node queues
Cs = 100; % the limit of the difference
T_Mode = 'NARE'; % the method to compute T

%% Response time
percentileRT = mainFJ(arrival, service, pers, K, Cs, T_Mode); 
for k = 1:length(K)
    disp(['K: ', int2str(K(k))]);
    for j = 1:length(percentileRT{k}.percentiles)
        disp(['RT(', num2str(percentileRT{k}.percentiles(j)), '): ', num2str(percentileRT{k}.RTp(j))]);
    end
    disp(' ');
end
    
