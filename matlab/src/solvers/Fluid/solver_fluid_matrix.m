function [QN,UN,RN,TN,xvec_it,QNt,UNt,TNt,xvec_t,t,iters,runtime] = solver_fluid_matrix(sn, options)

% [QN,UN,RN,TN,CN,RUNTIME] = SOLVER_FLUID_MATRIX(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
pie = sn.pie;
PH = sn.proc;
P = sn.rt;
NK = sn.njobs';  %initial population
S = sn.nservers;
infServers = isinf(S);
S(infServers) = sum(NK);
nphases = sn.phases;
%refstat = sn.refstat; % reference station
weights = ones(M,K);

% ODE building as per Ruuskanen et al., PEVA 151 (2021).
Psi = [];
A = [];
B = [];
for i=1:M
    for r=1:K
        if nphases(i,r)==0
            Psi = blkdiag(Psi,0);
            B = blkdiag(B,0);
            A = blkdiag(A,NaN);        
        else
            Psi = blkdiag(Psi,PH{i}{r}{1});
            B = blkdiag(B,sum(PH{i}{r}{2},2));
            A = blkdiag(A,pie{i}{r}');        
        end
    end
end
W = Psi + B*P*A';

% remove disabled transitions
keep = find(~isnan(sum(W,1)));
W = W(keep,:);
W = W(:,keep);

Qa = []; % state mapping to queues (called Q(a) in Ruuskanen et al.)
SQC = zeros(M*K,0); % to compute per-class queue length at the end
SUC = zeros(M*K,0); % to compute per-class utilizations at the end
STC = zeros(M*K,0); % to compute per-class throughput at the end
x0 = options.init_sol(:);
%x0 = []; % initial state
state = 0;
for i=1:M
    for r=1:K
        for k=1:nphases(i,r)
            state = state + 1;
            Qa(1,state)=i; %#ok<*AGROW>
            SQC((i-1)*K+r,state) = 1;
            SUC((i-1)*K+r,state) = 1/S(i);
            STC((i-1)*K+r,state) = sum(sn.proc{i}{r}{2}(k,:));
        end
        % code to initialize all jobs at ref station
        %if i == refstat(r)
        %    x0 = [x0; NK(r)*pie{i}{r}']; % initial probability of PH
        %else
        %    x0 = [x0; zeros(nphases(i,r),1)];
        %end
    end
end

SQ = zeros(0,length(x0)); % to compute total queue length in ODEs
for i=1:M
    for r=1:K
        for k=1:nphases(i,r)
            SQ(end+1,find(Qa==i)) = weights(i,r); %#ok<FNDSB>
        end
    end
end

%x0

tol = options.tol;
timespan = options.timespan;
itermax = options.iter_max;
odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(x0));
nonZeroRates = abs(W(abs(W)>0)); nonZeroRates=nonZeroRates(:);
trange = [timespan(1),min(timespan(2),abs(10*itermax/min(nonZeroRates)))];

T0 = tic;
iters = 1;
if options.stiff
    [t, xvec_t] = ode_solve_stiff(@(t,x) W'*(x./(GlobalConstants.FineTol+SQ*x).*min(S(Qa),GlobalConstants.FineTol+SQ*x)), trange, x0, odeopt, options);
else
    [t, xvec_t] = ode_solve(@(t,x) W'*(x./(GlobalConstants.FineTol+SQ*x).*min(S(Qa),GlobalConstants.FineTol+SQ*x)), trange, x0, odeopt, options);
end
runtime = toc(T0);

Tmax = size(xvec_t,1);
QNtmp = cell(1,Tmax);
UNtmp = cell(1,Tmax);
RNtmp = cell(1,Tmax);
TNtmp = cell(1,Tmax);
Sa = S(Qa);
S = repmat(S,1,K)'; S=S(:);
for j=1:Tmax
    x = xvec_t(j,:)';
    QNtmp{j} = zeros(K,M);
    TNtmp{j} = zeros(K,M);
    UNtmp{j} = zeros(K,M);
    RNtmp{j} = zeros(K,M);

    QNtmp{j}(:) = SQC*x;
    TNtmp{j}(:) = STC*(x./(GlobalConstants.FineTol+SQ*x).*min(Sa,GlobalConstants.FineTol+SQ*x));
    UNtmp{j}(:) = SUC*(x./(GlobalConstants.FineTol+SQ*x).*min(Sa,GlobalConstants.FineTol+SQ*x));
    % Little's law is invalid in transient so this vector is not returned
    % except the last element as an approximation of the actual RN
    RNtmp{j}(:) = QNtmp{j}(:)./TNtmp{j}(:);

    QNtmp{j} = QNtmp{j}';
    UNtmp{j} = UNtmp{j}';
    RNtmp{j} = RNtmp{j}';
    TNtmp{j} = TNtmp{j}';
end
% steady state metrics
for j=1:Tmax
    QNtmp{j} = QNtmp{j}(:);
    UNtmp{j} = UNtmp{j}(:);
    RNtmp{j} = RNtmp{j}(:);
    TNtmp{j} = TNtmp{j}(:);
end

% compute cell array with time-varying metrics for stations and classes
QNtmp = cell2mat(QNtmp)';
UNtmp = cell2mat(UNtmp)';
RNtmp = cell2mat(RNtmp)';
TNtmp = cell2mat(TNtmp)';
QNt = cell(M,K);
UNt = cell(M,K);
RNt = cell(M,K);
TNt = cell(M,K);
for i=1:M
    for r=1:K
        QNt{i,r} = QNtmp(:,(r-1)*M+i);
        UNt{i,r} = UNtmp(:,(r-1)*M+i);
        RNt{i,r} = RNtmp(:,(r-1)*M+i);
        TNt{i,r} = TNtmp(:,(r-1)*M+i);
    end
end
QN = reshape(QNtmp(end,:),M,K);
UN = reshape(UNtmp(end,:),M,K);
RN = reshape(RNtmp(end,:),M,K);
TN = reshape(TNtmp(end,:),M,K);
xvec_it = {xvec_t(end,:)};
end