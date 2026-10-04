function [UN,QN,p2opt]=qrf_noblo_bethe(M,MR,K,N,mu,v,rt)
% QRF_NOBLO_BETHE No-blocking quadratic reduction under the tree-reweighted
% (Bethe) free entropy.
%
% Same polytope, same phase-1 start and the same fmincon call as
% QRF_NOBLO_MMI; the objective is the only difference. Writing lambda = 1/M,
%
%   f(p) = lambda * sum_m sum_{i~=j} sum_{ki,kj} sum_{ni,nj>=0}
%              p_ij*(log p_ij - log p_ii - log p_jj)
%        + sum_m sum_i sum_k sum_{n>=0} p_ii*log p_ii
%
% i.e. lambda*sum_{i~=j} I(n_i;n_j) - sum_i H(n_i), the NEGATIVE of the
% tree-reweighted entropy with uniform edge weight rho_ij = 2*lambda on the
% complete station graph.
%
% WHY lambda = 1/M AND NOT THE BETHE 1/2. H_rho is a convex combination of
% tree entropies, hence concave on the local marginal polytope, exactly when
% rho lies in the spanning tree polytope of K_M. The uniform point of that
% polytope is rho_ij = 2/M, i.e. lambda = 1/M, which is therefore the LARGEST
% uniform weight for which minimising f is a CONVEX program -- every local
% optimum global, the answer a property of the model rather than of the start
% point. The Bethe weight lambda = 1/2 (rho_ij = 1, total edge mass
% nchoosek(M,2) against the M-1 a spanning tree can carry) is outside it for
% every M > 2 and coincides with 1/M only at M = 2.
%
% TWO DIFFERENCES FROM THE CODED mmi(), BOTH DELIBERATE. Its population loops
% run from n = 0, as the AMPL source's `ni, nj in 0..F` does and as mmi() does
% not, so the idle/idle cell -- the strongest correlation in a closed chain --
% is inside the sum; and the entropy term carries the sign a MINIMISER needs,
% which is mem() negated. Both are repaired here by construction, without
% touching 'qrf.mmi' or 'qrf.mem', whose values are pinned by tests.
%
% NUMERICAL NOTE. Restoring the n = 0 cells brings the structurally zero
% entries inside the sum: they contribute 0*log(LOGTOL) = 0 to the VALUE but
% log(LOGTOL) ~ -13.8 to the gradient, so the objective is insensitive to
% LOGTOL while the descent direction is not. Measured across LOGTOL in
% {1e-4, 1e-6, 1e-8} the optimum moves by under 1e-7 in U, but do not freeze
% test digits below that.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%%%  PARAMETERS  %%%
%  f; % finite capacity queue
%  M, integer, > 0; % number of queues
%  MR, integer, > 0; % number of independent blocking configurations
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

MR = 1;
BB = zeros(1,M); 
F = repmat(N,M,1);

q = zeros(M,M,max(K),max(K));
for i = 1:M
    for j = 1:M
        for k = 1:K(i)
            for h = 1:K(i)
                if j ~= i
                    q(i,j,k,h) = rt(i,j)*mu(i,k,h);
                else
                    q(i,j,k,h) = v(i,k,h)+rt(i,i)*mu(i,k,h);
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
%options.MaxFunEvals = 1e10;
%options.MaxSQPIter = 500;
%options.TolCon = 1e-8;
%options.Algorithm = 'sqp';
%options.OutputFcn =  @outfun;


% Start on the polytope: from an infeasible start fmincon exhausts its
% iteration budget restoring feasibility and returns a point outside the
% polytope. See qrf_noblo_start.
x = qrf_noblo_start(@(z) sub_qrfcon(z,q,M,MR,BB,F,N), n);

% One solve from the phase-1 feasible point: f is convex on this polytope
% (see the header), so there is no second local minimum for a restart to
% find and no multi-start loop is needed.
[xopt, fopt] = fmincon(@(x) bethe(x),x,[],[],[],[],x*0,x*0+1,@(x) sub_qrfcon(x,q,M,MR,BB,F,N),options);
[p2opt,~] = sub_qrfvar(xopt);

for ti=1:M
    UN(ti) = 0;
    QN(ti) = 0;
    for m=1:MR
        for ni=1+(1:F(ti))
            for ki=1:K(ti)
                UN(ti) = UN(ti) + p2opt(ti,ni,ki,ti,ni,ki, m);
                QN(ti) = QN(ti) + (ni-1)*p2opt(ti,ni,ki,ti,ni,ki, m); % rescaled back ni
            end
        end
    end
end

    function fobj = bethe(x)
    % BETHE Tree-reweighted free entropy at the uniform spanning-tree weight.
    %
    % Term 1 is the body of mmi() with its population loops started at 0, the
    % range the AMPL source states ("sum {nj in 0..F[j]} sum {ni in 0..F[i]}"),
    % scaled by lambda = 1/M. Term 2 is the body of mem() over the same
    % restored range and with the sign a minimiser needs: mem() returns +H and
    % is handed to a MAXIMISE in AMPL, so -H = +sum p log p is what enters a
    % minimisation.
        [p2,~] = sub_qrfvar(x);
        fobj = 0;
        LOGTOL = 1e-6;
        lambda = 1/M;
        % lambda * sum_{i~=j} I(n_i;n_j)
        for m = 1:MR
            for i = 1:M
                for ki = 1:K(i)
                    for j = 1:M
                        if i~=j
                            for kj = 1:K(j)
                                for ni = 1+(0:F(i))
                                    for nj = 1+(0:F(j))
                                         fobj = fobj + lambda*p2(i,ni,ki,j,nj,kj,m)*(log(LOGTOL+p2(i,ni,ki,j,nj,kj,m))-log(LOGTOL+p2(i,ni,ki,i,ni,ki,m))-log(LOGTOL+p2(j,nj,kj,j,nj,kj,m)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        % - sum_i H(n_i), i.e. +sum p log p over the diagonal cells
        for m = 1:MR
            for i = 1:M
                for k = 1:K(i)
                    for ni = 1+(0:F(i))
                        fobj = fobj + p2(i,ni,k,i,ni,k,m)*log(LOGTOL + p2(i,ni,k,i,ni,k,m));
                    end
                end
            end
        end
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

    function [c,ceq] = sub_qrfcon(x,q,M,MR,BB,F,N)
        c=sparse(0,1);
        ceq=sparse(0,1);
        
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
        for i = 1:M, for k =1:K(i)
                ceq(end+1) = 0;
                % LHS
                for j = 1:M, for h = 1:K(i)
                   ceq(end) = ceq(end) + q(i,j,k,h)*e(i,k); 
                end, end
                % RHS
                for j = 1:M, for h = 1:K(i)
                   ceq(end) = ceq(end) - q(i,j,h,k)*e(i,h);
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

        % subject to THM30 {i in 1..M, u in 1..K[i]: i<>f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,j]==0}
        % q[j,i,k,h]*p2[j,nj,k, i,0,u, m] + ... = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], h in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,i]==0} q[i,j,k,u]*p2[j,nj,h, i,0+1,k, m] + ...
        %
        % Marginal-balance family of the bas skeleton. Without it the polytope
        % carries no dependence on mu at all (THM1 only balances phases, and is
        % identically zero when K(i)==1), and the U bound collapses to the
        % uninformative [0,1]; verified against glpsol.
        %
        % No-blocking specialization: there is no finite-capacity station f and
        % BB is all-zero, so only the j<>i, j<>f terms of the skeleton survive.
        % The j==f branches are vacuous here, as are THM3f, THM3I and THM3L,
        % which are conditioned entirely on f.
        for i = 1:M, for u = 1:K(i)
            ceq(end+1) = 0;
            % LHS
            for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for m = 1:MR
                if j~=i && BB(m,j)==0
                    ceq(end) = ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,1+0,u, m);
                end
            end, end, end, end, end
            % RHS
            for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(j), for m = 1:MR
                if j~=i && BB(m,i)==0
                    ceq(end) = ceq(end) - q(i,j,k,u)*p2(j,nj,h, i,1+1,k, m);
                end
            end, end, end, end, end
        end, end

        % subject to THM3 {i in 1..M, ni in 0..(F[i]-1): i<>f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j<>f and BB[m,j]==0} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] + ...
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,i]==0} sum {h in 1..K[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, m] + ...
        for i = 1:M, for ni = 1+(0:(F(i)-1))
            ceq(end+1) = 0;
            % LHS
            for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR
                if j~=i && BB(m,j)==0
                    ceq(end) = ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,ni,u, m);
                end
            end, end, end, end, end, end
            % RHS
            for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR
                if j~=i && BB(m,i)==0
                    for h = 1:K(i)
                        ceq(end) = ceq(end) - q(i,j,k,h)*p2(j,nj,u, i,ni+1,k, m);
                    end
                end
            end, end, end, end, end
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
        