function [QN,xvec_it,QNt,UNt,xvec_t,t,iters,runtime] = solver_fluid(sn, options)
% [QN,XVEC_IT,QNT,UNT,XVEC_T,T,ITERS,RUNTIME] = SOLVER_FLUID(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
N = sn.nclosedjobs;    %population
Mu = sn.mu;
Phi = sn.phi;
PH = sn.proc;
sched = sn.sched;
rt = sn.rt;
S = sn.nservers;
NK = sn.njobs';  %initial population

match = zeros(M,K); % indicates whether a class is served at a station
for ist = 1:M
    for k = 1:K
        if isnan(Mu{ist}{k})
            Mu{ist}{k} = [];
            Phi{ist}{k}=[];
        end
        match(ist,k) = sum(rt(:, (ist-1)*K+k)) > 0;
    end
    %Set number of servers in delay station = population
    if isinf(S(ist))
        S(ist) = N;
    end
end

%% initialization
max_time = Inf;
Tstart = tic;

%phases = sn.phases;
phases = zeros(M,K);
for ist = 1:M
    for k = 1:K
        phases(ist,k) = length(Mu{ist}{k});
    end
end

slowrate = zeros(M,K);
for ist = 1:M
    for k = 1:K
        if ~isempty(Mu{ist}{k})
            slowrate(ist,k) = min(Mu{ist}{k}(:)); %service completion (exit) rates in each phase
        else
            slowrate(ist,k) = Inf;
        end
    end
end

%NK(isinf(NK))=1e6;
iters = 0;
xvec_it = {};
y0 = [];
assigned = zeros(1,K); %number of jobs of each class already assigned
for ist = 1:M
    for k = 1:K
        if match(ist,k) > 0 % indicates whether a class is served at a station
            if isinf(NK(k))
                if sched(ist)==SchedStrategy.EXT
                    toAssign = 1; % open job pool
                else
                    toAssign = 0; % set to zero open jobs everywhere
                end
            else
                toAssign = floor(NK(k)/sum(match(:,k)));
                if sum( match(ist+1:end, k) ) == 0 % if this is the last station for this job class
                    toAssign = NK(k) - assigned(k);
                end
            end
            y0 = [y0, toAssign, zeros(1,phases(ist,k)-1)]; % this is just for PS
            assigned(k) = assigned(k) + toAssign;
        else
            y0 = [y0, zeros(1,phases(ist,k))];
        end
    end
end

if isempty(options.init_sol)
    xvec_it{iters +1} = y0(:)'; % average state embedded at stage change transitions out of e
    ydefault = y0(:)'; % not used in this case
else
    xvec_it{iters +1} = options.init_sol(:)';
    ydefault = y0(:)'; % default solution if init_sol fails
end

%% solve ode
[xvec_it, xvec_t, t, iters] = solver_fluid_iteration(sn, N, Mu, Phi, PH, rt, S, xvec_it, ydefault, slowrate, Tstart, max_time, options);

runtime = toc(Tstart);
% if options.verbose >= 2
%     if iters==1
%         line_printf('Fluid analysis iteration completed in %0.6f sec [%d iteration]\n',runtime,iters);
%     else
%         line_printf('Fluid analysis iteration completed in %0.6f sec [%d iterations]\n',runtime,iters);
%     end
% end

% this part assumes PS, DPS, GPS scheduling
QN = zeros(M,K);
QNt = cell(M,K);
Qt = cell(M,1);
UNt = cell(M,K);
for ist=1:M
    Qt{ist} = 0;
    for k = 1:K
        shift = sum(sum(phases(1:ist-1,:))) + sum( phases(ist,1:k-1) ) + 1;
        QN(ist,k) = sum(xvec_it{end}(shift:shift+phases(ist,k)-1));
        QNt{ist,k} = sum(xvec_t(:,shift:shift+phases(ist,k)-1),2);
        Qt{ist} = Qt{ist} + QNt{ist,k};
        % computes queue length in each node and stage, summing the total
        % number in service and waiting in that station
        % results are weighted with the stat prob of the stage
    end
end

for ist=1:M
    if sn.nservers(ist) > 0 % not INF
        for k = 1:K
            UNt{ist,k} = min(QNt{ist,k} / S(ist), QNt{ist,k} ./ Qt{ist}); % if not an infinite server then this is a number between 0 and 1
            UNt{ist,k}(isnan(UNt{ist,k})) = 0; % fix cases where qlen is 0
        end
    else % infinite server
        for k = 1:K
            UNt{ist,k} = QNt{ist,k};
        end
    end
end
return
end
