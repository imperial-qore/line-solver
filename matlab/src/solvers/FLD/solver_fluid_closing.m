function [QN, UN, RN, TN, xvec_it, QNt, UNt, TNt, xvec_t, t, iters, runtime] = solver_fluid_closing(sn, options)
% [QN, UN, RN, TN, XVEC, QNt, UNt, TNt, XVEC_T, T, ITERS, RUNTIME] = SOLVER_FLUID_CLOSING(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
chains = sn.chains;

% inner iteration of fluid analysis
[Q,xvec_it,QNt,UNt,xvec_t,t,iters,runtime] = solver_fluid(sn, options);

%% assumes the existence of a delay node through which all classes pass
delayNodes = zeros(1,sn.nstations);
delayrefstat = zeros(1,sn.nstations); % delay nodes that are also reference nodes
for i = 1:sn.nstations
    if sn.sched(i) == SchedStrategy.INF
        delayNodes(i) = 1;
    end
end
for k = 1:sn.nclasses
    if sn.refstat(k) > 0 % artificial classes do not need to have a reference node
        delayrefstat(sn.refstat(k)) = 1;
    end
end

%% simpler version: no chain analysis
M = sn.nstations;
K = sn.nclasses;
refstat = sn.refstat;
Lambda = sn.mu;
phi = sn.phi;
phases = sn.phases;

%Qlength for all stations, classes,
QN = xvec_it{end};

%% Tput for all classes in each station
Rlittle = zeros(sn.nstations,sn.nclasses); % throughput of every class at each station
TN = zeros(sn.nstations,sn.nclasses); % throughput of every class at each station
TNt = cell(size(TN));
for i=1:size(TNt,1)
    for j=1:size(TNt,2)
        TNt{i,j} = 0;
    end
end
Xservice = cell(M,K); %throughput per class, station and phase
for i = 1:sn.nstations
    if delayNodes(i) == 1
        for k = 1:sn.nclasses
            idx = sum(sum(phases(1:i-1,:))) + sum( phases(i,1:k-1) );
            Xservice{i,k} = zeros(phases(i,k),1);
            for f = 1:phases(i,k)
                TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f);
                TNt{i,k} = TNt{i,k} + QNt{i,k}*Lambda{i}{k}(f)*phi{i}{k}(f);
                Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f);
            end
        end
    else
        xi = sum(Q(i,:)); %number of jobs in the station
        xi_t = QNt{i,1};
        for r=2:size(QNt,2)
            xi_t = xi_t + QNt{i,r};
        end

        if xi>0 || sn.sched(i)==SchedStrategy.EXT % step in to compute source tput in all cases
            switch sn.sched(i)
                case SchedStrategy.FCFS
                    wni = GlobalConstants.FineTol;
                    w = zeros(1,sn.nclasses);
                    for k = 1:sn.nclasses
                        idx = sum(sum(phases(1:i-1,:))) + sum( phases(i,1:k-1) );
                        wi(k) = map_mean(sn.proc{i}{k});
                        wni = wni + wi(k)*sum(QN((idx+1):(idx+phases(i,k))));
                    end
                    wni_t = 0*QNt{i,1};
                    for r=1:sn.nclasses
                        wni_t = wni_t + wi(r)* QNt{i,r};
                    end

            end

            for k = 1:sn.nclasses
                idx = sum(sum(phases(1:i-1,:))) + sum( phases(i,1:k-1) );
                Xservice{i,k} = zeros(phases(i,k),1);
                for f = 1:phases(i,k)
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            if f==1
                                TN(i,k) = TN(i,k) + (1-sum(QN(idx+(2:phases(i,k)))))*Lambda{i}{k}(f)*phi{i}{k}(f);
                                TNt{i,k} = TNt{i,k} + (1-sum(xvec_t(:,idx+(2:phases(i,k))),2))*Lambda{i}{k}(f)*phi{i}{k}(f);
                                Xservice{i,k}(f) = (1-sum(QN(idx+(2:phases(i,k)))))*Lambda{i}{k}(f);
                            else
                                TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f);
                                TNt{i,k} = TNt{i,k} + xvec_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f);
                                Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f);
                            end
                        case {SchedStrategy.INF, SchedStrategy.PS}
                            TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)/xi*min(xi,sn.nservers(i));
                            TNt{i,k} = TNt{i,k} + xvec_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)./xi_t.*min(xi_t,sn.nservers(i));
                            Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f)/xi*min(xi,sn.nservers(i));
                        case SchedStrategy.DPS
                            w = sn.schedparam(i,:);
                            wxi = w*Q(i,:)'; %number of jobs in the station
                            wxi_t = w(1)*QNt{i,1};
                            for r=2:size(QNt,2)
                                wxi_t = wxi_t + w(r)*QNt{i,r};
                            end
                            TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)*w(k)/wxi*min(xi,sn.nservers(i));
                            TNt{i,k} = TNt{i,k} + xvec_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)*w(k)./wxi_t.*min(xi_t,sn.nservers(i));
                            Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f)*w(k)/wxi*min(xi,sn.nservers(i));
                        case {SchedStrategy.FCFS, SchedStrategy.SIRO}
                            switch options.method
                                case {'default','closing','pnorm','softmin'}
                                    TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)/xi*min(xi,sn.nservers(i));
                                    TNt{i,k} = TNt{i,k} + xvec_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)./xi_t.*min(xi_t,sn.nservers(i));
                                    Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f)/xi*min(xi,sn.nservers(i));
                                case 'statedep'
                                    TN(i,k) = TN(i,k) + QN(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)*wi(k)/wni*min(xi,sn.nservers(i));
                                    TNt{i,k} = TNt{i,k} + xvec_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)./xi_t.*min(xi_t,sn.nservers(i));
                                    Xservice{i,k}(f) = QN(idx+f)*Lambda{i}{k}(f)*wi(k)/wni*min(xi,sn.nservers(i));

                                    %Tfull(i,k) = Tfull(i,k) + Qfull(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)*min(xi,sn.nservers(i))*wi(k);
                                    %Tfull_t{i,k} = Tfull_t{i,k} + ymean_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f).*min(xi_t,sn.nservers(i)) * wi(k);
                                    %Xservice{i,k}(f) = Qfull(idx+f)*Lambda{i}{k}(f)*min(xi,sn.nservers(i))*wi(k);

                                    %                                    Tfull(i,k) = Tfull(i,k) + Qfull(idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f)*min(xi,sn.nservers(i))*wi(k)/wni;
                                    %                                    Tfull_t{i,k} = Tfull_t{i,k} + ymean_t(:,idx+f)*Lambda{i}{k}(f)*phi{i}{k}(f).*min(xi_t,sn.nservers(i)) * wi(k)/wni;
                                    %                                    Xservice{i,k}(f) = Qfull(idx+f)*Lambda{i}{k}(f)*min(xi,sn.nservers(i))*wi(k)/wni;
                            end
                        otherwise
                            line_error(mfilename,'Unsupported scheduling policy');
                    end
                end
            end
        end
    end
end
for i=1:size(TNt,1)
    for j=1:size(TNt,2)
        if isscalar(TNt{i,j})
            TNt{i,j} = TNt{i,j} * ones(size(xvec_t(:,1),1),1);
        end
    end
end
%% response times
origK = size(chains,1);
% This is approximate, Little's law does not hold in transient
R = zeros(sn.nstations, sn.nclasses);
for i = 1:sn.nstations
    R(i, TN(i,:)>0) = Q(i,TN(i,:)>0) ./ TN(i,TN(i,:)>0);
end
%R = Rlittle;
%Tfull = Tlittle;
newR = zeros(sn.nstations,origK);
X = zeros(1,origK);
newQ = zeros(sn.nstations,origK);
% determine eventual probability of visiting each station in each class (expected number of visits)
% weight response time of each original class with the expected number of
% visits to each station in each associated artificial class
idxNR = reshape(1:sn.nstations*sn.nclasses, sn.nclasses, sn.nstations); %indices of non-delay (non-reference) nodes
idxNR = reshape(idxNR(:,delayrefstat==0),1,[]);
rtTrans = sn.rt(idxNR,idxNR); %transient transition matrix for non-reference nodes
eventualVisit = pinv( eye(size(rtTrans)) - rtTrans ); % pinv handles better self-looping customers setting the relevant entries to zero instead than infinity
idxOrigClasses = zeros(origK,1);

for k = 1:origK
    idxOrigClasses(k) = find(chains(k,:),1);
    refNode = refstat(idxOrigClasses(k));

    eventualVisitProb = reshape( sn.rt((refNode-1)*sn.nclasses+k,idxNR)*eventualVisit, sn.nclasses , sn.nstations-sum(delayrefstat>0) )'; %#ok<MINV> %probability of eventual visit
    eventualVisitProb = eventualVisitProb(:,chains(k,:)==1);

    newR(refNode,k) = sum( R(refNode,chains(k,:)==1),2 );
    newR(delayrefstat==0,k) = sum( R(delayrefstat==0,chains(k,:)==1).*eventualVisitProb,  2 );
    newR(refNode,chains(k,:)==1) = R(refNode,chains(k,:)==1);
    newR(delayrefstat==0,chains(k,:)==1) = R(delayrefstat==0,chains(k,:)==1).*eventualVisitProb;

    %X(k,find(chains(k,:)))= sum(Tfull(refNode,find(chains(k,:))));
    X(1,k) = sum( TN(refNode,chains(k,:)==1) );
    newQ(:,k) = sum( Q(:,chains(k,:)==1),2 );
end
QN = Q;
RN = R;

%% Utilization

UN = zeros(M,K);
for i =1:M
    for k = 1:K
        idx = Xservice{i,k}>0;
        UN(i,k) = sum(Xservice{i,k}(idx)./ Lambda{i}{k}(idx));
        switch sn.sched(i)
            case SchedStrategy.FCFS
                switch options.method
                    case 'statedep'
                        TN(i,k) = sum(Xservice{i,k}(idx));
                end
        end
    end
end

UN(delayNodes==0,:) = UN(delayNodes==0,:)./repmat(sn.nservers(delayNodes==0),1,K);

end
