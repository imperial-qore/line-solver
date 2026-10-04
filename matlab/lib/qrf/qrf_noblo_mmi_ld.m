function [UN,QN,BN,p2opt]=qrf_noblo_mmi_ld(MAPs,N,rt,alpha)
%%%  PARAMETERS  %%%
%  f; % finite capacity queue
%  M, integer, > 0; % number of queues
%  MR, integer, > 0; % number of independent blocking configurations
%  MM {m in 1:MR, i in 1:M} >=0; % blocking order
%  MM1 {m in 1:MR, i in 1:M}; % blocking order
%  ZZ {m in 1:MR} >=0; % nonzeros in independent blocking configurations
%  ZM, integer, >=0; % max of ZZ
%  BB {m in 1:MR, i in 1:M} >=0; % blocking state
%  K {i in 1:M}, integer, > 0; % number of phases for each queue
%  F {i in 1:M}, integer, > 0; % capacity
%  N, integer, >0;  % population
%  mu {i in 1:M, k in 1:K(i), h in 1:K(i)} >=0; % completion transition rates
%  v {i in 1:M, k in 1:K(i), h in 1:K(i)} >=0; % background transition rates
%  r {i in 1:M, j in 1:M} >=0; % routing probabilities

%%%  VARIABLES  %%%
%var p2 {j = 1:M, nj = 1+(0:N), k = 1:K(j), i = 1:M, ni = 1+(0:N), h = 1:K(i), m = 1:MR} >= 0;
%var e {i = 1:M, k = 1:K(i)} >=0;

M = length(MAPs);
if nargin < 4 || isempty(alpha)
    alpha = ones(M, N); % load independent, as the three ports already default it
end
MR = 1;
MM = 0;
F = repmat(N,M,1);
MM1 = ones(1,3);
ZZ = zeros(1,M);
ZM = 0;
BB = zeros(MR,M);

for i=1:M
    K(i)=size(MAPs{i}{1},1);
end

r = zeros(M);
for i=1:M
    for j=1:M
        r(i,j)=rt(i,j);
    end
end

for i=1:M
    for h=1:size(MAPs{i}{1},1)
        for k=1:size(MAPs{i}{1},1)
            mu(i,h,k)=MAPs{i}{2}(h,k);
        end
    end
end

for i=1:M
    for h=1:size(MAPs{i}{1},1)
        for k=1:size(MAPs{i}{1},1)
            if h==k
                v(i,h,k)=0; 
            else
                v(i,h,k)=MAPs{i}{1}(h,k); % (from,to), as mu is: v(i,k,h) reverses an Erlang
            end
        end
    end
end

% q {i in 1:M, j in 1:M, k in 1:K(i), h in 1:K(i), n in 0:N}, with the
% load-dependent scaling alpha(i,n) folded in, and q(...,n=0)=0. The
% population index n is that of station i, the station emitting the
% transition; MATLAB stores it shifted as 1+n.
q = zeros(M,M,max(K),max(K),N+1);
for i = 1:M
    for j = 1:M
        for k = 1:K(i)
            for h = 1:K(i)
                for ni = 1:N
                    if j ~= i
                        q(i,j,k,h,1+ni) = rt(i,j)*mu(i,k,h)*alpha(i,ni);
                    else
                        q(i,j,k,h,1+ni) = v(i,k,h)*alpha(i,ni)+rt(i,i)*mu(i,k,h)*alpha(i,ni);
                    end
                end
            end
        end
    end
end

n = M*(N+1)*max(K)*M*(N+1)*max(K)*MR + M*max(K);

options = optimset('fmincon');
options.Display = 'off';
%options.LargeScale = 'off';
options.MaxIter =  100;
options.MaxFunEvals = 1e3;
options.MaxSQPIter = 50;
options.TolCon = 1e-8;
options.Algorithm = 'sqp';
%options.OutputFcn =  @outfun;


% Start on the polytope: from an infeasible start fmincon exhausts its
% iteration budget restoring feasibility and returns a point outside the
% polytope. See qrf_noblo_start.
x = qrf_noblo_start(@(z) sub_qrfcon(z,q,M,MR,MM,MM1,ZZ,ZM,BB,F,N), n);

[xopt, fopt] = fmincon(@(x) mmi(x),x,[],[],[],[],x*0,x*0+1,@(x) sub_qrfcon(x,q,M,MR,MM,MM1,ZZ,ZM,BB,F,N),options);
[p2opt,~] = sub_qrfvar(xopt);

for ti=1:M
    UN(ti) = 0;
    QN(ti) = 0;
    BN(ti) = 0;
    for m=1:MR
        for ni=1+(1:F(ti))
            for ki=1:K(ti)
                UN(ti) = UN(ti) + p2opt(ti,ni,ki,ti,ni,ki, m);
                QN(ti) = QN(ti) + (ni-1)*p2opt(ti,ni,ki,ti,ni,ki, m); % rescaled back ni
                % BN is the alpha-weighted marginal mean: E[min(n,c)] at a
                % c-server station, E[n] at a delay, P(n>=1) where alpha is 1.
                % It is what the departure rate is proportional to, so the
                % station throughput is BN/stime exactly -- see sn_to_qrf_alpha.
                BN(ti) = BN(ti) + alpha(ti,ni-1)*p2opt(ti,ni,ki,ti,ni,ki, m);
            end
        end
    end
end

% add LB>=0 to both p2 and e
% LINEAR PROGRAMMING: UTILIZATION UPPER BOUND AT QUEUE 1
%minimize U1min: sum {m in 1..MR} sum {k in 1..K[1]} sum {n1 in 1..F[1]} p2[1,n1,k,1,n1,k,m]; 


    function fobj = mmi(x)
    % MMI
    %minimize MI: sum {m in 1..MR} sum {i in 1..M, j in 1..M, ki in 1..K[i], kj in  1..K[j]: i<>j} sum {nj in 0..F[j]} sum {ni in 0..F[j]} p2[i,ni,ki,j,nj,kj,m]*(log(1e-6+p2[i,ni,ki,j,nj,kj,m])-log(1e-6+p2[i,ni,ki,i,ni,ki,m])-log(1e-6+p2[j,nj,kj,j,nj,kj,m]));
        [p2,~] = sub_qrfvar(x);
        fobj = 0;
        for m = 1:MR
            for i = 1:M
                for ki = 1:K(i)
                    for j = 1:M
                        if i~=j
                            for kj = 1:K(j)
                                for ni = 1+(0:F(i))
                                    for nj = 1+(0:F(j))
                                         fobj = fobj + p2(i,ni,ki,j,nj,kj,m)*(log(1e-6+p2(i,ni,ki,j,nj,kj,m))-log(1e-6+p2(i,ni,ki,i,ni,ki,m))-log(1e-6+p2(j,nj,kj,j,nj,kj,m)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    function fobj = mem(x)
        % MEM
        %maximize H: -sum {m in 1..MR} sum {i in 1..M} sum {k in 1..K[i]} sum {ni in 1..F[i]} p2[i,ni,k,i,ni,k,m]*log(1e-6+p2[i,ni,k,i,ni,k,m]);
        %
        % Returned NEGATED, i.e. as -H, because fmincon MINIMIZES and the AMPL
        % objective above is a MAXIMIZE. Returning +H picks the minimum-entropy
        % point of the polytope, the opposite face, under a method documented as
        % maximum-entropy; that is what every port did until 2026-08-29.
        [p2,~] = sub_qrfvar(x);
        fobj = 0;
        for m = 1:MR
            for i = 1:M
                for k = 1:K(i)
                    for ni = 1+(1:F(i))
                        fobj = fobj + p2(i,ni,k,i,ni,k,m)*log(1e-6 + p2(i,ni,k,i,ni,k,m));
                    end
                end
            end
        end        
    end

    function fobj = umin(x, i)
        [p2,~] = sub_qrfvar(x);
        fobj = sum(p2(i,1+(1:F(i)),1:K(i),i,1+(1:F(i)),1:K(i), 1:MR),'all');
    end

    function [p2,e] = sub_qrfvar(x)
        ctr = 1;
        p2 = zeros(M,N+1,max(K),M,N+1,max(K),MR);
        for j = 1:M
            for nj = 1+(0:N)
                for k = 1:K(j)
                    for i = 1:M
                        for ni = 1+(0:N)
                            for h = 1:K(i)
                                for m = 1:MR
                                    p2(j,nj,k,i,ni,h,m) = x(ctr);
                                    ctr = ctr + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        e = zeros(M,max(K));
        for i=1:M
            for k=1:K(i)
                e(i,k) = x(ctr);
                ctr = ctr + 1;
            end
        end
    end

    function [c,ceq] = sub_qrfcon(x,q,M,MR,MM,MM1,ZZ,ZM,BB,F,N)
        c=zeros(0,1);
        ceq=zeros(0,1);
        
        %%%  VARIABLES  %%%
        [p2,e] = sub_qrfvar(x);
        
        %% DEFINITIONS
        % subject to ONE {j in 1..M}: sum {nj in 0..N, k in 1..K[j], m in 1..MR} p2[j,nj,k,j,nj,k,m]=1;
        for j = 1:M
            % LHS            
            ceq(end+1) = 0;
            for nj = 1+(0:N), for k = 1:K(j), for m = 1:MR
                ceq(end) =  ceq(end) + p2(j,nj,k,j,nj,k,m);
            end, end, end
            % RHS
            ceq(end) = ceq(end) -1;
        end
        
        % subject to ZERO1 {j in 1..M, k in 1..K[j], nj in 0..N, i in 1..M, h in 1..K[i], ni in 0..N, m in 1..MR: i==j and nj==ni and h<>k}: p2[j,nj,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+(0:N), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 1:MR 
            if i==j && (nj-1)==(ni-1) && h~=k % rescaled back nj and ni
                ceq(end+1) = p2(j,nj,k,i,ni,h,m); %=0
            end
        end, end, end, end, end, end, end
        
        % subject to ZERO2 {j in 1..M, k in 1..K[j], nj in 0..N, i in 1..M, h in 1..K[i], ni in 0..N, m in 1..MR: i==j and nj<>ni}: p2[j,nj,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+(0:N), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 1:MR 
            if i==j && (nj-1)~=(ni-1) % rescaled back nj and ni
                ceq(end+1) = p2(j,nj,k,i,ni,h,m); %=0
            end
        end, end, end, end, end, end, end
        
        % subject to ZERO3 {j in 1..M, k in 1..K[j], nj in 0..N, i in 1..M, h in 1..K[i], ni in 0..N, m in 1..MR: i<>j and nj+ni>N}: p2[j,nj,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+(0:N), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 1:MR
        	if i~=j && (nj-1)+(ni-1)>N  % rescaled back nj and ni
                ceq(end+1) = p2(j,nj,k,i,ni,h,m);
            end
        end, end, end, end, end, end, end        
        
        % subject to ZERO5 {j in 1..M, k in 1..K[j], i in 1..M, h in 1..K[i], ni in 0..F[i], m in 2..MR: BB[m,j]==1}: p2[j,0,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for i = 1:M, for h = 1:K(i), for ni = 1+(0:F(i)), for m = 2:MR
            if BB(m,j)==1
                ceq(end+1) = p2(j,1+0,k,i,ni,h,m);
            end
        end, end, end, end, end, end
        
        % subject to ZERO6 {j in 1..M, k in 1..K[j], nj in F[j]+1..N, i in 1..M, h in 1..K[i], ni in 0..N, m in 1..MR}: p2[j,nj,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+((F(j)+1):N), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 1:MR
            ceq(end+1) = p2(j,nj,k,i,ni,h,m);
        end, end, end, end, end, end, end
        
        % subject to ZERO7 {j in 1..M, k in 1..K[j], nj in 1..F[j], i in 1..M, h in 1..K[i], ni in 0..N, m in 2..MR: BB[m,j]==1 and i<>j and i<>f and ni+nj+F[f]>N}: p2[j,nj,k,i,ni,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+(1:F(j)), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 2:MR 
            if BB(m,j)==1 && i~=j && i~=f && (ni-1)+(nj-1)+F(f)>N % rescaled back ni and nj
                ceq(end+1) = p2(j,nj,k,i,ni,h,m);
            end
        end, end, end, end, end, end, end                                                           
        
        % subject to SIMMETRY {j in 1..M, nj in 0..N, k in 1..K[j], i in 1..M, ni in 0..N, h in 1..K[i], m in 1..MR}: p2[i,ni,h,j,nj,k,m] = p2[j,nj,k,i,ni,h,m];
        for j = 1:M, for nj = 1+(0:N), for k =1:K(j), for i = 1:M, for ni = 1+(0:N), for h = 1:K(i), for m = 1:MR
            ceq(end+1) = p2(i,ni,h,j,nj,k,m) - p2(j,nj,k,i,ni,h,m);
        end, end, end, end, end, end, end                                
        
        % subject to MARGINALS {j in 1..M, k in 1..K[j], nj in 0..N, i in 1..M, m in 1..MR: i<>j}: p2[j,nj,k,j,nj,k,m]= sum {ni in 0..N-nj} sum {h in 1..K[i]} p2[j,nj,k,i,ni,h,m];
        for j = 1:M, for k =1:K(j), for nj = 1+(0:N), for i = 1:M, for m = 1:MR
            if i~=j
                % LHS
                ceq(end+1) = p2(j,nj,k,j,nj,k,m);
                % RHS
                for ni = 1+(0:N) % full range per AMPL MARGINALS; ZERO3 zeroes nj+ni>N                    
                    for h = 1:K(i) 
                        ceq(end) = ceq(end) - p2(j,nj,k,i,ni,h,m);
                    end
                end
            end
        end, end, end, end, end
        
        %subject to UEFF {j in 1..M, i in 1..M, ki in 1..K[i]}: e[i,ki] = sum {nj in 0..N, kj in 1..K[j], m in 1..MR, ni in 1..N: BB[m,i]==0} p2[j,nj,kj,i,ni,ki,m];
        for j = 1:M, for i = 1:M, for ki = 1:K(i)
            % LHS
            ceq(end+1) = e(i,ki);
            % RHS
            for nj = 1+(0:N), for kj = 1:K(j), for m = 1:MR, for ni = 1+(1:N)
                if BB(m,i)==0
                    ceq(end) = ceq(end) - p2(j,nj,kj,i,ni,ki,m);
                end
            end, end, end, end
        end, end, end       
        
        %subject to THM1 {i in 1..M, k in 1..K[i]}: sum {j in 1..M, h in 1..K[i]} q[i,j,k,h]*e[i,k] =sum {j in 1..M, h in 1..K[i]} q[i,j,h,k]*e[i,h];
        % Load-dependent form (cf. THM2 of qrboundsrsrd_skel.mod): q carries the
        % population index of station i, so the aggregated e[i,k] must be resolved
        % per population. UEFF plus MARGINALS give e[i,k] = sum_{ni>=1} of the
        % station-i marginal p2[i,ni,k,i,ni,k,m], so that marginal is the
        % population-resolved summand. For alpha(i,:) constant this collapses
        % exactly onto the e-based constraint above.
        for i = 1:M, for k =1:K(i)
                ceq(end+1) = 0;
                for ni = 1+(1:F(i)), for m = 1:MR
                    % LHS
                    for j = 1:M, for h = 1:K(i)
                       ceq(end) = ceq(end) + q(i,j,k,h,ni)*p2(i,ni,k,i,ni,k,m);
                    end, end
                    % RHS
                    for j = 1:M, for h = 1:K(i)
                       ceq(end) = ceq(end) - q(i,j,h,k,ni)*p2(i,ni,h,i,ni,h,m);
                    end, end
                end, end
        end, end
        
        %subject to THM2 {j in 1..M, k in 1..K[j], nj in 0..F[j], m in 1..MR}: sum {i in 1..M, ni in 1..F[i], ki in 1..K[i]} ni*p2[j,nj,k,i,ni,ki,m]= N*p2[j,nj,k,j,nj,k,m];
        for j = 1:M, for k =1:K(j), for nj = 1+(0:F(j)), for m = 1:MR
            ceq(end+1) = 0;
            % LHS
            for i = 1:M, for ni = 1+(1:F(i)), for ki = 1:K(i)
                ceq(end) = ceq(end) + (ni-1)*p2(j,nj,k,i,ni,ki,m); % recaled back ni
            end, end, end
            % RHS
            ceq(end) = ceq(end) - N*p2(j,nj,k,j,nj,k,m);
        end, end, end, end
        
        %subject to COR1 : sum {m in 1..MR, i in 1..M, j in 1..M, nj in 1..F[j], ni in 1..F[i], ki in 1..K[i], kj in 1..K[j]} ni*nj*p2[j,nj,kj,i,ni,ki,m]= N^2;        
        ceq(end+1) =  0;
        % LHS
        for m = 1:MR, for i = 1:M, for j = 1:M, for nj = 1+(1:F(j)), for ni = 1+(1:F(i)), for ki = 1:K(i), for kj = 1:K(j) 
            ceq(end) = ceq(end) + (ni-1)*(nj-1)*p2(j,nj,kj,i,ni,ki,m); % rescaled back ni and nj
        end, end, end, end, end, end, end
        % RHS
        ceq(end) = ceq(end) - N^2;
        
        %subject to THM30 {i in 1..M, u in 1..K[i]: i<>f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,j]==0} 
        % q[j,i,k,h]*p2[j,nj,k, i,0,u, m] +  sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], m in 1..MR: j<>i and j==f and MM[m,1]<>i} q[j,i,k,h]*p2[j,nj,k, i,0,u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], h in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,i]==0} q[i,j,k,u]*p2[j,nj,h, i,0+1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], k in 1..K[i], h in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==0 and nj<F[j]} q[i,j,k,u]*p2[j,nj,h, i,1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], y in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==1 and nj==F[j] and MM[m,1]==i} sum {p in 1..K[f], w in 1..M: w<>f and w<>i} q[f,w,y,p]*p2[f,nj,y, i,1,u, m] ;
        
        f = -1; % never equal to FC station
        for i = 1:M, for u = 1:K(i)
                if i~=f
                    ceq(end+1)=0;
                    % LHS
                    for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for m = 1:MR 
                        if j~=i && j~=f && BB(m,j)==0 
                            ceq(end)= ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,1+0,u, m);
                        end
                    end, end, end, end, end     
    
                    for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for m = 1:MR 
                        if j~=i && j==f && MM(m,1)~=i 
                            ceq(end)= ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,1+0,u, m);
                        end
                    end, end, end, end, end
                    % RHS
                    for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(j), for m = 1:MR 
                        if j~=i && j~=f && BB(m,i)==0 
                            ceq(end)= ceq(end) - q(i,j,k,u,1+1)*p2(j,nj,h, i,1+1,k, m); 
                        end
                    end, end, end, end, end
    
                    for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(j), for m = 1:MR 
                        if j~=i && j==f && BB(m,i)==0 && (nj-1)<F(j)  % rescaled back nj 
                            ceq(end)= ceq(end) - q(i,j,k,u,1+1)*p2(j,nj,h, i,1+1,k, m); 
                        end
                    end, end, end, end, end
    
                    for j = 1:M, for nj = 1+(0:F(j)), for y = 1:K(j), for m = 1:MR 
                        if j~=i && j==f && BB(m,i)==1 && (nj-1)==F(j) && MM(m,1)==i  % rescaled back nj 
                            for p = 1:K(f), for w = 1:M
                                if w~=f && w~=i 
                                    ceq(end)= ceq(end) - q(f,w,y,p,nj)*p2(f,nj,y, i,1+1,u, m);
                                end
                            end, end
                        end
                    end, end, end, end
                    end % if
        end, end
                
        % subject to THM3 {i in 1..M, ni in 0..(F[i]-1): i<>f}: 
        % sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j<>f and BB[m,j]==0} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % +  sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j==f and MM[m,1]<>i} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,i]==0} sum {h in 1..K[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==0 and nj<F[j]} sum {h in 1..K[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==1 and nj==F[j] and MM[m,1]==i} sum {p in 1..K[f], w in 1..M: w<>f and w<>i} q[f,w,u,p]*p2[f,nj,u, i,ni+1,k, m] ;
        for i = 1:M, for ni = 1+(0:(F(i)-1))
                if i~=f
                ceq(end+1)=0;
                % LHS                
                for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR 
                    if j~=i && j~=f && BB(m,j)==0 
                        ceq(end)= ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,ni,u, m); 
                    end                   
                end, end, end, end, end, end
                for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR 
                    if j~=i && j==f && MM(m,1)~=i 
                        ceq(end)= ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,ni,u, m);
                    end
                end, end, end, end, end, end
                % RHS
                for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
                    if j~=i && j~=f && BB(m,i)==0 
                        for h = 1:K(i) 
                            ceq(end)= ceq(end) - q(i,j,k,h,ni+1)*p2(j,nj,u, i,ni+1,k, m); 
                        end
                    end
                end, end, end, end, end
                for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
                    if j~=i && j==f && BB(m,i)==0 && (nj-1)<F(j)   % rescaled back nj
                        for h = 1:K(i) 
                            ceq(end)= ceq(end) - q(i,j,k,h,ni+1)*p2(j,nj,u, i,ni+1,k, m);
                        end
                    end
                end, end, end, end, end
                for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
                    if j~=i && j==f && BB(m,i)==1 && (nj-1)==F(j) && MM(m,1)==i  % rescaled back nj 
                        for p = 1:K(f), for w = 1:M 
                            if w~=f && w~=i 
                                ceq(end)= ceq(end) - q(f,w,u,p,nj)*p2(f,nj,u, i,ni+1,k, m);
                            end
                        end, end
                    end
                end, end, end, end, end
                end
        end, end
                
        % subject to THM3f {i in 1..M, ni in 0..(F[i]-1): i==f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j<>f and BB[m,j]==0 and ni < F[i]} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % +  sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j==f and ni==F[i] and MM[m,1]==j} sum {w in 1..M: w<>f} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], h in 1..K[i], u in 1..K[j]: j<>i and ni < F[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, 1];
        for i = 1:M, for ni = 1+(0:(F(i)-1))
                if i==f
                    ceq(end+1) = 0;
                    % LHS
                    for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR
                        if j~=i && j~=f && BB(m,j)==0 && (ni-1) < F(i)  % rescaled back ni
                            ceq(end) = ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,ni,u, m);
                        end
                    end, end, end, end, end, end
                    for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR
                        if j~=i && j==f && (ni-1)==F(i) && MM(m,1)==j  % rescaled back ni
                            for w = 1:M
                                if w~=f
                                    ceq(end) = ceq(end) + q(j,i,k,h,nj)*p2(j,nj,k, i,ni,u, m);
                                end
                            end
                        end
                    end, end, end, end, end, end
                    % RHS
                    for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(i), for u = 1:K(j)
                        if j~=i && (ni-1) < F(i)  % rescaled back ni
                            ceq(end) = ceq(end) - q(i,j,k,h,ni+1)*p2(j,nj,u, i,ni+1,k, 1);
                        end
                    end, end, end, end, end
                end % if
        end, end
    
        % subject to THM4 {j in 1..M, k in 1..K[j], i in 1..M, m in 1..MR}: sum{t in 1..M} sum {h in 1..K[t]}  sum {nj in 0..N} sum {nt in 0..N} nt*p2[j,nj,k,t,nt,h,m] 
        % >= N*sum {h in 1..K[i]}  sum {nj in 0..N} sum {ni in 1..N} (p2[j,nj,k,i,ni,h,m]);
        for j = 1:M, for k = 1:K(j), for i = 1:M, for m = 1:MR
            c(end+1) = 0; % <= inequality
            % LHS with sign swapped since >= in GLPK
            for t = 1:M, for h = 1:K(t), for nj = 1+(0:N), for nt = 1+(0:N)
                c(end) = c(end) - (nt-1)*p2(j,nj,k,t,nt,h,m);  % rescaled back nt
            end, end, end, end
            % RHS with sign swapped since >= in GLPK
            for h = 1:K(i), for nj = 1+(0:N), for ni = 1+(1:N)
                c(end) = c(end) + N*(p2(j,nj,k,i,ni,h,m));        
            end, end, end 
        end, end, end, end
	end
end
        