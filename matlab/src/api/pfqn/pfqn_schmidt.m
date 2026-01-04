%{
%{
 % @file pfqn_schmidt.m
 % @brief Schmidt's exact MVA for networks with general scheduling disciplines.
%}
%}

%{
%{
 % @brief Schmidt's exact MVA for networks with general scheduling disciplines.
 % @fn pfqn_schmidt(D, N, S, sched, v)
 % @param D Service demand matrix.
 % @param N Population vector.
 % @param S Number of servers per station (matrix or vector).
 % @param sched Scheduling discipline per station.
 % @param v Visit ratio matrix (optional, defaults to ones).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Cycle times.
 % @return T Results table.
%}
%}
function [XN,QN,UN,CN,T] = pfqn_schmidt(D,N,S,sched,v)
% [XN,QN,UN,CN] = PFQN_SCHMIDT(D,N,S,SCHED,V)

% utilization in general ld case does not work
[M,R] = size(D);
closedClasses = 1:R;
XN = zeros(1,R);
UN = zeros(M,R);
CN = zeros(M,R);
QN = zeros(M,R);

% Default visit ratios to ones if not provided
if nargin < 5 || isempty(v)
    v = ones(M,R);
end

% Compute pure service times from demands: S_pure = D / v
% D is demands (= visit_ratio * service_time), so S_pure = D / v
% Bcn calculations need pure service times, not demands
S_pure = D ./ max(v, 1e-12);  % Protect against division by zero

C = length(closedClasses); % number of closed classes
Dc = D(:,closedClasses);
Nc = N(closedClasses);
prods = zeros(1,C); % needed for fast hashing
for r=1:C
    prods(r)=prod(Nc(1:r-1)+1);
end
% Start at nc=(0,...,0)
kvec = pprod(Nc);
% Initialize L and Pc
L = cell(1, M); % mean queue-length
Pc = cell(1, M); % state probabilities (empty cells for stations that don't need them)
for ist=1:M
    % Always initialize L for ALL stations (queue lengths needed everywhere)
    L{ist} = zeros(R, prod(1+Nc)); % mean queue-length
    switch sched(ist)
        case SchedStrategy.INF
            % No Pc needed for infinite server
        case SchedStrategy.PS
            if ~all(S(ist,:) == 1)
                % Multi-server PS needs state probabilities
                Pc{ist} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
            end
        case SchedStrategy.FCFS
            if all(D(ist,:)==D(ist,1)) % class-independent
                if ~all(S(ist,:) == 1) % multi server
                    Pc{ist} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
                end
            else % class-dependent - needs vector probabilities
                Pc{ist} = zeros(prod(1+Nc), prod(1+Nc)); % Pr(jvec|N)
            end
    end
end
% Per-station throughput like Java: x(ist,c,hkvec) = v(ist,c) * kvec(c) / denom
x = zeros(M,C,prod(1+Nc));
w = zeros(M,C,prod(1+Nc));
for ist=1:M
    if ~isempty(Pc{ist})
        Pc{ist}(1 + 0, hashpop(kvec,Nc,C,prods)) = 1.0; %Pj(0|0) = 1
    end
end
u = zeros(M,C);
% Population recursion
while all(kvec>=0) && all(kvec <= Nc)
    nc = sum(kvec);
    kprods = zeros(1,C); % needed for fast hashing
    for r=1:C
        kprods(r)=prod(kvec(1:r-1)+1);
    end
    for ist=1:M
        for c=1:C
            hkvec = hashpop(kvec,Nc,C,prods);
            hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
            if size(S(ist,:)) == 1
                ns = S(ist);
            else
                ns = S(ist,c);
            end
            if kvec(c) > 0
                switch sched(ist)
                    case SchedStrategy.INF
                        w(ist,c,hkvec) = D(ist,c);
                    case SchedStrategy.PS
                        if ns == 1
                            % Sum queue lengths of ALL classes (arrival theorem)
                            totalQueueLength = sum(L{ist}(:, hkvec_c));
                            w(ist,c,hkvec) = Dc(ist,c) * (1 + totalQueueLength);
                        else
                            % Sum queue lengths of ALL classes (arrival theorem)
                            totalQueueLength = sum(L{ist}(:, hkvec_c));
                            w(ist,c,hkvec) = (Dc(ist,c) / ns) * (1 + totalQueueLength);
                            for j=1:ns-1
                                % Pc{ist}(j,...) = Pr(j-1 jobs) due to MATLAB 1-indexing
                                w(ist,c,hkvec) = w(ist,c,hkvec) + (ns-1-(j-1))*Pc{ist}(j, hkvec_c) * (Dc(ist,c) / ns);
                            end
                        end
                    case SchedStrategy.FCFS
                        if all(D(ist,:)==D(ist,1)) % product-form case
                            if ns == 1
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist,c,hkvec) = Dc(ist,c) * (1 + totalQueueLength);
                            else
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist,c,hkvec) = (Dc(ist,c) / ns) * (1 + totalQueueLength);
                                for j=1:ns-1
                                    % Pc{ist}(j,...) = Pr(j-1 jobs) due to MATLAB 1-indexing
                                    w(ist,c,hkvec) = w(ist,c,hkvec) + (ns-1-(j-1))*Pc{ist}(j, hkvec_c) * (Dc(ist,c) / ns);
                                end
                            end
                        else
                            if ns == 1
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist,c,hkvec) = Dc(ist,c) * (1 + totalQueueLength);
                            else
                                % Multi-server FCFS class-dependent case
                                nvec = pprod(kvec);
                                while nvec >= 0
                                    if nvec(c) > 0
                                        % Use Nc and prods for hash (matches Java)
                                        hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                                        % Bcn calculation for multi-server using pure service times
                                        if sum(nvec) <= ns
                                            Bcn = S_pure(ist,c);
                                        else
                                            % Weighted average service time based on queue composition
                                            % Use epsilon protection for division (matches Java)
                                            Bcn = S_pure(ist,c) + max(0,sum(nvec)-ns)/max(ns*(sum(nvec)-1), 1e-12) * (nvec*S_pure(ist,:)' - S_pure(ist,c));
                                        end
                                        w(ist,c,hkvec) = w(ist,c,hkvec) + Bcn * Pc{ist}(hnvec_c, hkvec_c);
                                    end
                                    nvec = pprod(nvec, kvec);
                                end
                            end
                        end
                end
            end
        end
    end
    % Compute throughputs (per-station like Java)
    for c=1:C
        % denom = sum over i of (v(i,c) * w(i,c,hkvec))
        denom = 0;
        for ist=1:M
            denom = denom + v(ist,c) * w(ist,c,hkvec);
        end
        % Per-station throughput: x(ist,c) = v(ist,c) * N(c) / denom
        for ist=1:M
            if denom > 0
                x(ist,c,hkvec) = v(ist,c) * kvec(c) / denom;
            else
                x(ist,c,hkvec) = 0;
            end
        end
    end
    % Update queue lengths: L = x * w
    for ist=1:M
        for c=1:C
            L{ist}(c,hkvec) = x(ist,c,hkvec) * w(ist,c,hkvec);
        end
        if size(S(ist,:)) == 1
            ns = S(ist);
        else
            ns = S(ist,c);
        end
        switch sched(ist)
            case SchedStrategy.PS
                if ns > 1
                    for n=1:min(S(ist),sum(kvec))
                        for c=1:C
                            if kvec(c) > 0
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                Pc{ist}(1 + n, hkvec) = Pc{ist}(1 + n, hkvec) + Dc(ist,c) * (1/n) * x(ist,c,hkvec) * Pc{ist}(1+(n-1), hkvec_c);
                            end
                        end
                        Pc{ist}(1 + 0, hkvec) = max(eps,1-sum(Pc{ist}(1 + (1:min(S(ist),sum(kvec))), hkvec)));
                    end
                end
            case SchedStrategy.FCFS
                if all(D(ist,:)==D(ist,1))
                    if ns > 1
                        for n=1:(min(ns,sum(kvec))-1)
                            for c=1:C
                                if kvec(c) > 0
                                    hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                    Pc{ist}(1 + n, hkvec) = Pc{ist}(1 + n, hkvec) + Dc(ist,c) * (1/n) * x(ist,c,hkvec) * Pc{ist}(1+(n-1), hkvec_c);
                                end
                            end
                            Pc{ist}(1 + 0, hkvec) = max(eps,1-sum(Pc{ist}(1 + (1:min(ns,sum(kvec))), hkvec)));
                        end
                    end
                else
                    nvec = pprod(kvec);
                    nvec = pprod(nvec, kvec); % Skip zero vector like Java
                    sumOfAllProbs = 0;
                    while nvec >= 0
                        % Use Nc and prods for hash (matches Java)
                        hnvec = hashpop(nvec,Nc,C,prods);
                        prob = 0;
                        for c=1:C
                            if nvec(c)>0
                                % Use Nc and prods for hash (matches Java)
                                hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                % Bcn calculation using pure service times (matches Java getBcn)
                                Bcn = S_pure(ist,c);
                                if sum(nvec) > 1
                                    sumVal = nvec*S_pure(ist,:)';
                                    Bcn = Bcn + max(0,sum(nvec)-ns)/max(ns*(sum(nvec)-1), 1e-12) * (sumVal - S_pure(ist,c));
                                end
                                % Use 1/sum(nvec) instead of 1/nvec(c) to match Java
                                capacity_inv = 1/sum(nvec);
                                classProb = Bcn * capacity_inv * x(ist,c,hkvec) * Pc{ist}(hnvec_c, hkvec_c);
                                prob = prob + classProb;
                            end
                        end
                        Pc{ist}(hnvec, hkvec) = prob;
                        sumOfAllProbs = sumOfAllProbs + prob;
                        nvec = pprod(nvec, kvec);
                    end
                    Pc{ist}(1 + 0, hkvec) = max(1e-12, 1 - sumOfAllProbs);
                end
        end
    end
    kvec = pprod(kvec, Nc);
end

% Throughput - compute system throughput as N(c) / sum(w) like Java
for c=1:C
    totalResponseTime = 0;
    for ist=1:M
        totalResponseTime = totalResponseTime + w(ist,c,hkvec);
    end
    if totalResponseTime > 0
        XN(c) = Nc(c) / totalResponseTime;
    else
        XN(c) = 0;
    end
end
if M>1
    XN = repmat(XN,M,1);
end
% Utilization
UN(1:M,closedClasses) = u(1:M,1:C); % this will return 0 
% Response time
CN(1:M,closedClasses) = w(1:M,1:C,hkvec);
for ist=1:M
    QN(ist,closedClasses) = L{ist}(closedClasses,hkvec);
end

T = table(XN,CN,QN,UN); % for display purposes

end