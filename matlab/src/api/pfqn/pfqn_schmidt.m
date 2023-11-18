function [XN,QN,UN,CN,T] = pfqn_schmidt(D,N,S,sched)
% [XN,QN,UN,CN] = PFQN_SCHMIDT(D,N,S,SCHED)

% utilization in general ld case does not work
[M,R] = size(D);
closedClasses = 1:R;
XN = zeros(1,R);
UN = zeros(M,R);
CN = zeros(M,R);
QN = zeros(M,R);

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
L = {}; % mean queue-length
Pc = {}; % state probabilities
for i=1:M
    switch sched(i)
        case SchedStrategy.ID_INF
            L{i} = zeros(R, prod(1+Nc)); % mean queue-length
        case SchedStrategy.ID_PS
            if all(S(i,:) == 1)
                L{i} = zeros(R, prod(1+Nc)); % mean queue-length
            else
                Pc{i} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
            end
        case SchedStrategy.ID_FCFS
            if all(D(i,:)==D(i,1)) % class-independent
                if all(S(i,:) == 1) % single server
                    L{i} = zeros(R, prod(1+Nc)); % mean queue-length
                else % multi server
                    Pc{i} = zeros(1 + sum(Nc), prod(1+Nc)); % Pr(j|N)
                end
            else % all general product nodes
                L{i} = zeros(R, prod(1+Nc));
                Pc{i} = zeros(prod(1+Nc), prod(1+Nc)); % Pr(jvec|N)
            end
    end
end
x = zeros(C,prod(1+Nc));
w = zeros(M,C,prod(1+Nc));
for i=1:M
    Pc{i}(1 + 0, hashpop(kvec,Nc,C,prods)) = 1.0; %Pj(0|0) = 1
end
u = zeros(M,C);
% Population recursion
while all(kvec>=0) && all(kvec <= Nc)
    nc = sum(kvec);
    kprods = zeros(1,C); % needed for fast hashing
    for r=1:C
        kprods(r)=prod(kvec(1:r-1)+1);
    end
    for i=1:M
        for c=1:C
            hkvec = hashpop(kvec,Nc,C,prods);
            hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
            if size(S(i,:)) == 1
                ns = S(i);
            else
                ns = S(i,c);
            end
            if kvec(c) > 0
                switch sched(i)
                    case SchedStrategy.ID_INF
                        w(i,c,hkvec) = D(i,c);
                    case SchedStrategy.ID_PS
                        if ns == 1
                            w(i,c,hkvec) = Dc(i,c) * (1 + L{i}(c, hkvec_c));
                        else
                            w(i,c,hkvec) = (Dc(i,c) / ns) * (1 + L{i}(c,hkvec_c));
                            for j=1:ns-1
                                w(i,c,hkvec) = w(i,c,hkvec) + (ns-1-(j-1))*Pc{i}(j, hkvec_c) * (Dc(i,c) / ns);
                            end
                        end
                    case SchedStrategy.ID_FCFS
                        if all(D(i,:)==D(i,1)) % product-form case
                            if ns == 1
                                w(i,c,hkvec) = Dc(i,c) * (1 + L{i}(c, hkvec_c));
                            else
                                w(i,c,hkvec) = (Dc(i,c) / ns) * (1 + L{i}(c,hkvec_c));
                                for j=1:ns-1
                                    w(i,c,hkvec) = w(i,c,hkvec) + (ns-1-(j-1))*Pc{i}(j, hkvec_c) * (Dc(i,c) / ns);
                                end
                            end
                        else 
                            if ns == 1
                                w(i,c,hkvec) = Dc(i,c) * (1 + L{i}(c, hkvec_c));
                            else
                                nvec = pprod(kvec);
                                while nvec >= 0
                                    if nvec(c) > 0
                                        hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods);
                                        if ns == 1
                                            Bcn = norm(nvec) * ns;
                                        else 
                                            if norm(nvec) <= ns
                                                Bcn = D(i,c);
                                            else
                                                Bcn = D(i,c) + max(0,norm(nvec)-ns)/(ns*(norm(nvec)-1)) * (nvec*D(i,:)' - D(i,c));
                                            end
                                        end
                                        w(i,c,hkvec) = w(i,c,hkvec) + Bcn * Pc{i}(hnvec_c, hkvec_c);
                                    end
                                    nvec = pprod(nvec, kvec);
                                end
                            end
                        end
                end
            end
        end
    end
    % Compute tput
    for c=1:C
        x(c,hkvec) = kvec(c) / sum(w(1:M,c,hkvec));
        x(isnan(x))=0; % avoid nan for base case
    end
    for i=1:M
        for c=1:C
            L{i}(c,hkvec) = x(c,hkvec) * w(i,c,hkvec);
        end
        if size(S(i,:)) == 1
            ns = S(i);
        else
            ns = S(i,c);
        end
        switch sched(i)
            case SchedStrategy.ID_PS
                if ns > 1
                    for n=1:min(S(i),sum(kvec))
                        for c=1:C
                            if kvec(c) > 0
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                Pc{i}(1 + n, hkvec) = Pc{i}(1 + n, hkvec) + Dc(i,c) * (1/n) * x(c,hkvec) * Pc{i}(1+(n-1), hkvec_c);
                            end 
                        end
                        Pc{i}(1 + 0, hkvec) = max(eps,1-sum(Pc{i}(1 + (1:min(S(i),sum(kvec))), hkvec)));
                    end
                end
            case SchedStrategy.ID_FCFS
                if all(D(i,:)==D(i,1))
                    if ns > 1
                        for n=1:(min(ns,sum(kvec))-1)
                            for c=1:C
                                if kvec(c) > 0
                                    hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                    Pc{i}(1 + n, hkvec) = Pc{i}(1 + n, hkvec) + Dc(i,c) * (1/n) * x(c,hkvec) * Pc{i}(1+(n-1), hkvec_c);
                                end 
                            end
                            Pc{i}(1 + 0, hkvec) = max(eps,1-sum(Pc{i}(1 + (1:min(ns,sum(kvec))), hkvec)));
                        end
                    end
                else
                    nvec = pprod(kvec);
                    while nvec >= 0
                        hnvec = hashpop(nvec,kvec,C,kprods);
                        for c=1:C
                            if nvec(c)>0
                                hnvec_c = hashpop(oner(nvec,c),kvec,C,kprods);
                                hkvec_c = hashpop(oner(kvec,c),Nc,C,prods);
                                if ns == 1
                                    Bcn = sum(nvec) * ns;
                                else 
                                    if sum(nvec) <= ns
                                        Bcn = D(i,c);
                                    else
                                        Bcn = D(i,c) + max(0,sum(nvec)-ns)/(ns*(sum(nvec)-1)) * (nvec*D(i,:)' - D(i,c));
                                    end
                                end
                                Pc{i}(hnvec, hkvec) = (1/nvec(c))*x(c,hkvec)*Bcn*Pc{i}(hnvec_c, hkvec_c);
                            end
                        end
                        Pc{i}(1 + 0, hkvec) = max(eps,1-sum(Pc{i}(1 + (1:min(ns,sum(kvec))), hkvec)));
                        nvec = pprod(nvec, kvec);
                    end
                end
        end
    end
    kvec = pprod(kvec, Nc);
end

% Throughput
XN(closedClasses) = x(1:C,hkvec);
if M>1
    XN = repmat(XN,M,1);
end
% Utilization
UN(1:M,closedClasses) = u(1:M,1:C); % this will return 0 
% Response time
CN(1:M,closedClasses) = w(1:M,1:C,hkvec);
for i=1:M
    QN(i,closedClasses) = L{i}(closedClasses,hkvec);
end

T = table(XN,CN,QN,UN); % for display purposes

end