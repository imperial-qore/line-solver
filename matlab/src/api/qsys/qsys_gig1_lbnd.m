function [W,rhohat]=qsys_gig1_lbnd(lambda,mu,ca,cs)
% [W,RHOHAT]=QSYS_GIG1_LBND(LAMBDA,MU,CA,CS)
%
% Computes fundamental theoretical lower bounds for G/G/1 queues.
% These are the minimum possible values that performance measures
% cannot fall below for any realization of the arrival and service
% processes.
%
% Inputs:
%   LAMBDA - Arrival rate
%   MU     - Service rate
%   CA     - Coefficient of variation of inter-arrival time
%   CS     - Coefficient of variation of service time
%
% Returns:
%   W      - Lower bound on average time in system (= 1/mu)
%   RHOHAT - Modified utilization (so that M/M/1 formulas still hold)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

W = 1/mu;  % At least the mean service time
rhohat = W*lambda/(1+W*lambda);

end
