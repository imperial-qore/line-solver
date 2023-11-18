function [ode_h,q_indices] = solver_fluid_odes(sn, N, Mu, phi, PH, P, nservers, sched, schedparam, options)
% [ODE_H,Q_INDICES] = SOLVER_FLUID_ODES(sn, N, MU, PHI, PH, P, NSERVERS, SCHED, SCHEDPARAM)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

M = length(nservers);    % number of stations
K = length(Mu{1});   % number of classes
w = ones(M,K);

enabled = false(M,K); % indicates whether a class is served at a station
% for i = 1:M
%     for c = 1:K
%         %enabled(i,c) = sum( P(:,(i-1)*K+c) ) > 0;
%         %changed to consider rows instead of columns, for response time
%         %analysis (first articial class never returns to delay node)
%         enabled(i,c) = sum( P((i-1)*K+c,:) ) > 0;
%     end
% end

q_indices = zeros(M,K);
Kic = zeros(M,K);
cumsum = 1;
for i = 1 : M
    for c = 1:K       
        if isnan(Mu{i}{c})
            numphases = 0;
            enabled(i,c) = false;
            q_indices(i,c) = cumsum;
        elseif isempty(Mu{i}{c})
            enabled(i,c) = false;
            numphases = 0;
            q_indices(i,c) = cumsum;
        else
            numphases = length(Mu{i}{c});            
            q_indices(i,c) = cumsum;
            enabled(i,c) = true;
        end
        Kic(i,c) = numphases;
        cumsum = cumsum + numphases;
    end
end

% to speed up convert sched strings in numerical values
schedid = zeros(1,M);
for i = 1 : M
    schedid(i) = SchedStrategy.toId(sched{i});
    switch schedid(i) % source
        case SchedStrategy.ID_DPS
            w(i,:) = schedparam(i,:);
    end
end

%[rateBase, eventIdx] = ode_hybrid_rate_base(Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, sched_id, all_jumps);

%% define ODE system to be returned
switch options.method
    case 'softmin'
        alpha = 20; % softmin parameter
        ode_sm_h = @(t,x) ode_softmin(x, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, schedid, alpha);
        ode_h = ode_sm_h;
    case 'statedep'
        ode_sd_h = @(t,x) ode_statedep(x, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, schedid);
        ode_h = ode_sd_h;
        %ode_sd_h = @(t,x) ode_hybrid(x, Phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched_id, all_jumps, rateBase, eventIdx);
    otherwise
        % determine all the jumps, and saves them for later use
        all_jumps = ode_jumps_new(M, K, enabled, q_indices, P, Kic);
        % determines a vector with the fixed part of the rates,
        % and defines the indexes that correspond to the events that occur
        [rateBase, eventIdx] = ode_rate_base(sn, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, schedid, all_jumps);
        ode_si_h = @(t,x) all_jumps * ode_rates_closing(x, M, K, enabled, q_indices, Kic, nservers, w, schedid, rateBase, eventIdx);
        ode_h = ode_si_h;
end
end