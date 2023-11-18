function [XN,QN,UN,CN,lGN,Pc] = pfqn_mvaldmx(lambda,D,N,Z,mu,S)
% [XN,QN,UN,CN,lGN,Pc] = PFQN_MVALDMX(LAMBDA,D,N,Z,MU,S)

if nargin<5
    mu=ones(size(D,1),sum(N(isfinite(N))));
    S=ones(size(D,1),1);
end
if size(mu,2) < sum(N(isfinite(N)))
    line_error(mfilename,'MVALDMX requires to specify the load-dependent rates with one job more than the maximum closed population.');
end
if any(N(find(lambda))>0 & isfinite(N(find(lambda))))
    line_error(mfilename,'Arrival rate cannot be specified on closed classes.');
end
[M,R] = size(D);
openClasses = find(isinf(N));
closedClasses = setdiff(1:length(N), openClasses);

XN = zeros(1,R);
UN = zeros(M,R);
CN = zeros(M,R);
QN = zeros(M,R);
lGN = 0;
mu(:,end+1) = mu(:,end); % we need up to sum(N)+1, but there is limited load dep
[EC,E,Eprime] = pfqn_mvaldmx_ec(lambda,D,mu);
C = length(closedClasses); % number of closed classes
Dc = D(:,closedClasses);
Nc = N(closedClasses);
Zc = Z(closedClasses);
prods = zeros(1,C); % needed for fast hashing
for r=1:C
    prods(r)=prod(Nc(1:r-1)+1);
end
% Start at nc=(0,...,0)
nvec = pprod(Nc);
% Initialize Pc
Pc = zeros(M,1+sum(Nc),prod(1+Nc));
x = zeros(C,prod(1+Nc));
w = zeros(M,C,prod(1+Nc));
for i=1:M
    Pc(i, 1 + 0, hashpop(nvec,Nc,C,prods)) = 1.0;
end
u = zeros(M,C);
% Population recursion
while nvec>=0
    hnvec = hashpop(nvec,Nc,C,prods);
    nc = sum(nvec);
    for i=1:M
        for c=1:C
            if nvec(c)>0
                hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                % Compute mean residence times
                for n=1:nc
                    w(i,c,hnvec) = w(i,c,hnvec) + Dc(i,c) * n * EC(i,n) * Pc(i, 1+(n-1), hnvec_c);
                end
            end
        end
    end
    % Compute tput
    for c=1:C
        x(c,hnvec) = nvec(c) / (Zc(c)+sum(w(1:M,c,hnvec)));
    end
    for i=1:M
        for n=1:nc
            for c=1:C
                if nvec(c)>0
                    hnvec_c = hashpop(oner(nvec,c),Nc,C,prods);
                    Pc(i, 1 + n, hnvec) = Pc(i, 1 + n, hnvec) + Dc(i,c) * EC(i,n) * x(c,hnvec) * Pc(i, 1+(n-1), hnvec_c);
                end
            end
        end
        Pc(i, 1 + 0, hnvec) = max(eps,1-sum(Pc(i, 1 + (1:nc), hnvec)));
    end
    
    % now compute the normalizing constant
    last_nnz = find(nvec>0, 1, 'last' );
    if sum(nvec(1:last_nnz-1)) == sum(Nc(1:last_nnz-1)) && sum(nvec((last_nnz+1):C))==0
        logX = log(XN(last_nnz));
        if ~isempty(logX)
            lGN = lGN - logX;
        end
    end
    
    nvec = pprod(nvec, Nc);
end

% compute performance indexes at Nc for closed classes
hnvec = hashpop(Nc,Nc,C,prods);
for c=1:C
    hnvec_c = hashpop(oner(Nc,c),Nc,C,prods);
    for i=1:M
        u(i,c) = 0;
        for n=1:sum(Nc) % closed class utilization
            u(i,c) = u(i,c) + Dc(i,c) * x(c,hnvec) * Eprime(i,1+n-1) / E(i,1+n-1) * Pc(i, 1+n-1, hnvec_c);
        end
    end
end

% Throughput
XN(closedClasses) = x(1:C,hnvec);
% Utilization
UN(1:M,closedClasses) = u(1:M,1:C);
% Response time
CN(1:M,closedClasses) = w(1:M,1:C,hnvec);
% Queue-length
QN(1:M,closedClasses) = repmat(XN(closedClasses),M,1) .* CN(1:M,closedClasses);

% Compute performance indexes at Nc for open classes
for r=openClasses
    % Throughput
    XN(r) = lambda(r);
    for i=1:M
        % Queue-length
        QN(i,r) = 0;   
        for n=0:sum(Nc)             
            QN(i,r) = QN(i,r) + lambda(r) * D(i,r) * (n+1) * EC(i,n+1) * Pc(i, 1+n, hnvec);
        end
        % Response time
        CN(i,r) = QN(i,r) / lambda(r);
        % Utilization - the formula from Bruell-Balbo-Ashfari does not
        % match simulation, this appears to be simly lambda_r*D_{ir}
        UN(i,r) = 0;
        for n=0:sum(Nc)
            UN(i,r) = UN(i,r) + lambda(r) * Eprime(i,1+n+1) / E(i,1+n+1) * Pc(i, 1+n, hnvec);
        end
    end
end

Pc = Pc(1:M,:,hnvec);
end
