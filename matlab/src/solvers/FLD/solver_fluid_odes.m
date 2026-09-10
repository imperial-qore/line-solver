function [ode_h,q_indices,rt_breaks,absorb] = solver_fluid_odes(sn, N, Mu, phi, PH, P, nservers, sched, schedparam, options)
% [ODE_H,Q_INDICES,RT_BREAKS,ABSORB] = SOLVER_FLUID_ODES(sn, N, MU, PHI, PH, P, NSERVERS, SCHED, SCHEDPARAM)
%
% RT_BREAKS are the instants at which the time-varying rate multiplier jumps,
% i.e. where the returned drift is DISCONTINUOUS in t. Empty unless an NHPP
% schedule is configured. The caller must integrate up to each of them and
% restart there rather than step across; see SOLVER_FLUID_RATEMULT.

% Copyright (c) 2012-2026, Imperial College London
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
% Only the closing family builds a time-varying multiplier, so every other
% method leaves this empty and integrates one uninterrupted window.
rt_breaks = [];
% Projector onto the coordinates that survive the immediate elimination, empty
% when nothing was eliminated. The caller must apply it to the initial point:
% mass parked on an eliminated coordinate has no event left to move it.
absorb = [];
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
for i = 1 : M
    switch sched(i) % source
        case {SchedStrategy.DPS, SchedStrategy.GPS}
            w(i,:) = schedparam(i,:);
    end
end

%% define ODE system to be returned
switch options.method
    case 'softmin'
        % options.config.alpha overrides the smoothing parameter, scalar or
        % one entry per station. It cannot be calibrated against the closure
        % variance: this Boltzmann operator lies ABOVE min() for n>c whereas
        % E[min(X,c)] lies below it by concavity, so no alpha reproduces the
        % Gaussian closure. Use options.method='minnormal' for that.
        alpha = 20;
        if isfield(options,'config') && isfield(options.config,'alpha') && ~isempty(options.config.alpha) ...
                && isnumeric(options.config.alpha)
            alpha = options.config.alpha;
        end
        ode_sm_h = @(t,x) ode_softmin(x, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched, alpha);
        ode_h = ode_sm_h;
    case 'pnorm'
        % p-norm smoothing method based on Ruuskanen et al., PEVA 151 (2021)
        if isfield(options, 'pstar') && ~isempty(options.pstar)
            pstar = options.pstar;
        else
            pstar = 20; % default p-norm parameter (similar behavior to softmin at alpha=20)
        end
        ode_pn_h = @(t,x) ode_pnorm(x, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched, pstar);
        ode_h = ode_pn_h;
    case 'statedep'
        ode_sd_h = @(t,x) ode_statedep(x, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, nservers, w, sched);
        ode_h = ode_sd_h;
    otherwise
        % Gaussian moment closure: SOLVER_FLUID_MOMENTS parks the converged
        % station population variances here, so the same closing ODE serves
        % both the first-order closure (absent or zero) and the second-order
        % one (see FLUID_MIN_CLOSURE)
        moment_sigma2 = [];
        if isfield(options,'config') && isfield(options.config,'moment_sigma2')
            moment_sigma2 = options.config.moment_sigma2;
        end

        % the DPS capacity share is a ratio of populations, so its closure
        % needs the covariance BETWEEN station coordinates, parked here by
        % SOLVER_FLUID_MOMENTS alongside the station variances
        moment_cov = {};
        if isfield(options,'config') && isfield(options.config,'moment_cov')
            moment_cov = options.config.moment_cov;
        end

        % limited load dependence: the station rate multiplier alpha(n_i)
        % multiplies whatever share the scheduling policy already applies, so
        % it composes with the closure rather than replacing it. Only the
        % closing family reads it; the other methods reject LD at the featset
        % gate (see SolverFLD.getMethodFeatureSet).
        lldscaling = [];
        if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling)
            lldscaling = sn.lldscaling;
        end

        % determine all the jumps, and saves them for later use
        all_jumps = ode_jumps_new(M, K, enabled, q_indices, P, Kic);
        % determines a vector with the fixed part of the rates,
        % and defines the indexes that correspond to the events that occur
        [rateBase, eventIdx] = ode_rate_base(sn, phi, Mu, PH, M, K, enabled, q_indices, P, Kic, sched, all_jumps);

        % Stochastic-complement the instantaneous coordinates out of the event
        % set, so no integrator has to step through an InfRate mode.
        if fluid_hide_immediate(sn, options)
            [all_jumps, rateBase, eventIdx, ~, ~, absorb] = ...
                ode_eliminate_immediate(all_jumps, rateBase, eventIdx, sn, options);
        end

        % see _kb/06-solver-catalog.md for rationale
        numEvents = numel(rateBase);
        [rt_tgrid, rt_Mmat, rt_breaks] = solver_fluid_ratemult(numEvents, M, K, enabled, ...
            q_indices, Kic, Mu, eventIdx, options);
        if isempty(rt_Mmat)
            ode_si_h = @(t,x) all_jumps * ode_rates_closing(x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx, moment_sigma2, lldscaling, moment_cov);
        else
            ode_si_h = @(t,x) all_jumps * ( fluid_interpcols(rt_tgrid, rt_Mmat, t) .* ode_rates_closing(x, M, K, enabled, q_indices, Kic, nservers, w, sched, rateBase, eventIdx, moment_sigma2, lldscaling, moment_cov) );
        end
        ode_h = ode_si_h;
end
end