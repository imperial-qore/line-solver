function [Q,U,R,T,C,X,lG,totiter] = solver_qna(sn, options)
% [Q,U,R,T,C,X,lG,totiter] = SOLVER_QNA(QN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
%
% Implementation as per Section 7.2.3 of N. Gautaum, Analysis of Queues, CRC Press, 2012.
% Minor corrections applied in discussion with the author.

config = options.config;
config.space_max = 1;

K = sn.nclasses;
rt = sn.rt;
S = 1./sn.rates;
scv = sn.scv; scv(isnan(scv))=0;

%% immediate feedback elimination
% this is adapted for class-switching, an alternative implementation would
% rescale by sum(rt((i-1)*K+r,(i-1)*K+1:K)) rather than rt((i-1)*K+r,(i-1)*K+r)
% for i=1:size(rt,1)
%     for r=1:K
%         for j=1:size(rt,2)
%             for s=1:K
%                 if i~=j
%                     rt((i-1)*K+r, (j-1)*K+s) = rt((i-1)*K+r, (j-1)*K+s) / (1-rt((i-1)*K+r,(i-1)*K+r));
%                 end
%             end
%         end
%         S(i,r) = S(i,r) / (1-rt((i-1)*K+r,(i-1)*K+r));
%         scv(i,r) = rt((i-1)*K+r,(i-1)*K+r) + (1-rt((i-1)*K+r,(i-1)*K+r))*scv(i,r);
%         rt((i-1)*K+r,(i-1)*K+r) = 0;
%     end
% end

%% generate local state spaces
I = sn.nnodes;
M = sn.nstations;
C = sn.nchains;
V = cellsum(sn.visits);
Q = zeros(M,K);
Q_1 = Q+Inf;

U = zeros(M,K);
R = zeros(M,K);
T = zeros(M,K);
X = zeros(1,K);

lambda = zeros(1,C);

it = 0;

if any(isfinite(sn.njobs))
    %    line_error(mfilename,'QNA does not support closed classes.');
end

%isMixed = isOpen & isClosed;
%if isMixed
% treat open as having higher priority than closed
%sn.classprio(~isfinite(sn.njobs)) = 1 + max(sn.classprio(isfinite(sn.njobs)));
%end

%% compute departure process at source
a1 = zeros(M,K);
a2 = zeros(M,K);
d2 = zeros(M,1);
f2 = zeros(M*K,M*K); % scv of each flow pair (i,r) -> (j,s)
for ist=1:M
    for jst=1:M
        if sn.nodetype(sn.stationToNode(jst)) ~= NodeType.Source
            for r=1:K
                for s=1:K
                    if rt((ist-1)*K+r, (jst-1)*K+s)>0
                        f2((ist-1)*K+r, (jst-1)*K+s) = 1; % C^2ij,r
                    end
                end
            end
        end
    end
end
lambdas_inchain = cell(1,C);
scvs_inchain = cell(1,C);
d2c = [];

for c=1:C
    inchain = sn.inchain{c};
    refStatIdx = sn.refstat(inchain(1));
    if isfinite(sn.njobs(c))
        Q(refStatIdx,:)=State.toMarginal(sn,sn.stationToNode(refStatIdx));
    end
    lambdas_inchain{c} = sn.rates(refStatIdx,inchain);
    scvs_inchain{c} = scv(refStatIdx,inchain);
    lambda(c) = sum(lambdas_inchain{c}(isfinite(lambdas_inchain{c})));
    d2c(c) = qna_superpos(lambdas_inchain{c},scvs_inchain{c});
    if isinf(sum(sn.njobs(inchain))) % if open chain
        T(refStatIdx,inchain') = lambdas_inchain{c};
    end
end

for c=1:C
    inchain = sn.inchain{c};
    refStatIdx = sn.refstat(inchain(1));
    d2(refStatIdx)=d2c*lambda'/sum(lambda);
end

%% main iteration

while max(max(abs(Q-Q_1))) > options.iter_tol && it <= options.iter_max %#ok<max>
    it = it + 1;
    for c=1:C
        inchain = sn.inchain{c};
        Q(:,c) = sn.njobs(c) .* Q(:,c) / sum(Q(:,c));
    end
    Q_1 = Q;

    for k=1:K
        if sn.isslc(k)
            Q(:,k) = 0;
            Q(sn.refstat(k),k) = sn.njobs(k);
        end
    end

    % update throughputs at all stations
    if it==1
        for c=1:C
            inchain = sn.inchain{c};
            for m=1:M
                T(m,inchain) = V(m,inchain) .* lambda(c);
            end
        end
    end

    % superposition
    for ist=1:M
        a1(ist,:) = 0;
        a2(ist,:) = 0;
        lambda_i = sum(T(ist,:));
        for jst=1:M
            for r=1:K
                for s=1:K
                    a1(ist,r) = a1(ist,r) + T(jst,s)*rt((jst-1)*K+s, (ist-1)*K+r);
                    a2(ist,r) = a2(ist,r) + (1/lambda_i) * f2((jst-1)*K+s, (ist-1)*K+r)*T(jst,s)*rt((jst-1)*K+s, (ist-1)*K+r);
                end
            end
        end
    end

    % update flow trhough queueing station
    for ind=1:I
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            switch sn.nodetype(ind)
                case NodeType.Join
                    % no-op
                    %                     for c=1:C
                    %                         inchain = sn.inchain{c};
                    %                         for k=inchain
                    %                             fanin = nnz(sn.rtnodes(:, (ind-1)*K+k));
                    %                             TN(ist,k) = lambda(c)*V(ist,k)/fanin;
                    %                             UN(ist,k) = 0;
                    %                             QN(ist,k) = 0;
                    %                             RN(ist,k) = 0;
                    %                         end
                    %                     end
                otherwise
                    switch sn.sched(ist)
                        case SchedStrategy.INF
                            for r=1:K
                                for s=1:K
                                    d2(ist,s) = a2(ist,s);
                                end
                            end
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    T(ist,k) = a1(ist,k);
                                    Q(ist,k) = T(ist,k).*S(ist,k)*V(ist,k);
                                    U(ist,k) = Q(ist,k);
                                    R(ist,k) = Q(ist,k)/T(ist,k);
                                end
                            end
                        case SchedStrategy.PS
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    T(ist,k) = lambda(c)*V(ist,k);
                                    U(ist,k) = S(ist,k)*T(ist,k);
                                end
                                Nc = sum(sn.njobs(inchain)); % closed population
                                Uden = min([1-options.tol,sum(U(ist,:))]);
                                for k=inchain
                                    Q(ist,k) = (U(ist,k)-U(ist,k)^(Nc+1))/(1-Uden); % geometric bound type approximation
                                    %Q(ist,k) = UN(ist,k)/(1-Uden);
                                    R(ist,k) = Q(ist,k)/T(ist,k);
                                end
                            end
                        case {SchedStrategy.FCFS}
                            mu_ist = sn.rates(ist,1:K);
                            mu_ist(isnan(mu_ist))=0;
                            rho_ist_class = a1(ist,1:K)./(GlobalConstants.FineTol+sn.rates(ist,1:K));
                            rho_ist_class(isnan(rho_ist_class))=0;
                            lambda_ist = sum(a1(ist,:));
                            mi = sn.nservers(ist);
                            rho_ist = sum(rho_ist_class) / mi;
                            if rho_ist < 1-options.tol
                                for k=1:K
                                    if rho_ist > 0.7
                                        alpha_mi = (rho_ist^mi+rho_ist) / 2;
                                    else
                                        alpha_mi = rho_ist^((mi+1)/2);
                                    end
                                    mubar(ist) = lambda_ist ./ rho_ist;
                                    c2(ist) = -1;
                                    for r=1:K
                                        if mu_ist(r)>0
                                            c2(ist) = c2(ist) + a1(ist,r)/lambda_ist * (mubar(ist)/mi/mu_ist(r))^2 * (scv(ist,r)+1 );
                                        end
                                    end
                                    Wiq(ist) = (alpha_mi / mubar(ist)) * 1/ (1-rho_ist) * (sum(a2(ist,:))+c2(ist))/2;
                                    Q(ist,k) = a1(ist,k) / mu_ist(k) + a1(ist,k)*Wiq(ist);
                                end
                                d2(ist) = 1 + rho_ist^2*(c2(ist)-1)/sqrt(mi) + (1 - rho_ist^2) *(sum(a2(ist,:))-1);
                            else
                                for k=1:K
                                    Q(ist,k) = sn.njobs(k);
                                end
                                d2(ist) = 1;
                            end
                            for k=1:K
                                T(ist,k) = a1(ist,k);
                                U(ist,k) = T(ist,k) * S(ist,k) /sn.nservers(ist);
                                R(ist,k) = Q(ist,k) ./ T(ist,k);
                            end
                    end
            end
        else % not a station
            switch sn.nodetype(ind)
                case NodeType.Fork
                    line_error(mfilename,'Fork nodes not supported yet by QNA solver.');
            end
        end
    end


    % splitting - update flow scvs
    for ist=1:M
        for jst=1:M
            if sn.nodetype(sn.stationToNode(jst)) ~= NodeType.Source
                for r=1:K
                    for s=1:K
                        if rt((ist-1)*K+r, (jst-1)*K+s)>0
                            f2((ist-1)*K+r, (jst-1)*K+s) = 1 + rt((ist-1)*K+r, (jst-1)*K+s) * (d2(ist)-1); % C^2ij,r
                        end
                    end
                end
            end
        end
    end

end

for k=1:K
    if sn.isslc(k)
        Q(:,k) = 0;
        ist = sn.refstat(k);
        Q(ist,k) = sn.njobs(k);
        T(ist,k) = sn.njobs(k)*sn.rates(ist,k);
        R(ist,k) = Q(ist,k) ./ T(ist,k);
        U(ist,k) = S(ist,k)*T(ist,k);
    end
end


for c=1:C
    inchain = sn.inchain{c};
    if isfinite(sn.njobs(c))
        Q(:,c) = sn.njobs(c) .* Q(:,c) / sum(Q(:,c));
    end
end
for ist=1:sn.nstations
    switch sn.sched(ist)
        case SchedStrategy.INF
            U(ist,:) = Q(ist,:);
    end
end

C = sum(R,1);
Q = abs(Q);
Q(isnan(Q))=0;
U(isnan(U))=0;
R(isnan(R))=0;
C(isnan(C))=0;
X(isnan(X))=0;
lG = 0;
totiter = it;
end

function [d2]=qna_superpos(lambda,a2)
a2 = a2(isfinite(lambda));
lambda = lambda(isfinite(lambda));
d2 = a2(:)'*lambda(:) / sum(lambda);
end