function [UN,QN,p2opt]=qrf_bas_bethe(f,M,MR,MM,MM1,ZZ,ZM,BB,K,F,N,mu,v,rt)
% QRF_BAS_BETHE BAS-blocking quadratic reduction under the tree-reweighted
% (Bethe) free entropy.
%
% Same polytope, same phase-1 start and the same fmincon call as
% QRF_BAS_MEM; the objective is the only difference. Writing lambda = 1/M,
%
%   f(p) = lambda * sum_m sum_{i~=j} sum_{ki,kj} sum_{ni,nj>=0}
%              p_ij*(log p_ij - log p_ii - log p_jj)
%        + sum_m sum_i sum_k sum_{n>=0} p_ii*log p_ii
%
% i.e. lambda*sum_{i~=j} I(n_i;n_j) - sum_i H(n_i), the NEGATIVE of the
% tree-reweighted entropy with uniform edge weight rho_ij = 2*lambda on the
% complete station graph. It is QRF_NOBLO_BETHE's objective evaluated over the
% BAS decision vector, so the blocking configurations m and the per-station
% capacities F enter only through the ranges.
%
% WHY lambda = 1/M AND NOT THE BETHE 1/2. H_rho is a convex combination of
% tree entropies, hence concave on the local marginal polytope, exactly when
% rho lies in the spanning tree polytope of K_M. The uniform point of that
% polytope is rho_ij = 2/M, i.e. lambda = 1/M, the LARGEST uniform weight for
% which minimising f is convex there. The BAS polytope adds the blocking
% families on top of the marginal ones, so it is a convex SUBSET and the
% objective's convexity carries; what is not inherited is a proof that the
% per-configuration marginals stay consistent under ZERO5, so start-point
% independence is MEASURED (python/tests/test_qrf_bas_bethe.py), not assumed.
%
% ENTROPY RANGE. The population sums of BOTH blocks run from n = 0, as they do
% in QRF_NOBLO_BETHE and as QRF_BAS_MEM's entropy does not: the idle cell is
% the strongest correlation in a closed chain, and an entropy taken over a
% different range than the mutual information it is combined with is not a free
% entropy of anything.
%
% THE POLYTOPE IS THE FULL FAMILY SET, THE SAME ONE THE PORTS ASSEMBLE. All
% three MATLAB BAS routines implement the five THM30-group families once, in
% the phase-aggregated form near the top of sub_qrfcon. QRF_BAS_MMI used to
% carry a SECOND, stale AMPL-literal transcription of each on top of that, so
% it emitted every one of them twice under different index conventions and its
% phase 1 could not reach feasibility (maximum equality residual 2.469e-01 on
% cqn_bas_blocking); QRF_BAS_MEM had commented those duplicates out, which is
% why it worked. The duplicates were deleted from QRF_BAS_MMI on 2026-09-03,
% so all three objectives now share one polytope.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
ZM = max(ZZ(:)); % max(ZZ) by definition; a larger one empties the polytope via THM3I
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

MR=size(MM,1);


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

% Start ON the polytope where that is affordable. With THM30/THM3/THM3f in the
% constraint set fmincon cannot recover from the all-zero point: it stops at a
% maximum residual of 1.2e-01, which the feasibility check after the solve then
% refuses. qrf_noblo_start recovers the affine form by calling the constraint
% callback once per variable, so it costs O(n) callback evaluations and is only
% usable on small instances; above the threshold the zero start is kept and the
% same check decides whether the result is a bound.
nVarsNlp = M*(N+1)*max(K)*M*(N+1)*max(K)*MR + M*max(K);
if nVarsNlp <= 5000
    x = qrf_noblo_start(@(z) sub_qrfcon(z,q,f,M,MR,MM,MM1,ZZ,ZM,BB,F,N), nVarsNlp);
else
    x = zeros(nVarsNlp, 1);
end

options = optimset('fmincon');
options.Display = 'off';
options.LargeScale = 'off';
options.MaxIter =  100;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 500;
%options.TolCon = 1e-8;
options.Algorithm = 'sqp';
%options.OutputFcn =  @outfun;

[xopt, fopt] = fmincon(@(x) bethe(x),x,[],[],[],[],x*0,x*0+1,@(x) sub_qrfcon(x,q,f,M,MR,MM,MM1,ZZ,ZM,BB,F,N),options);

% fmincon returns whatever point it reached, converged or not. On a 3-station
% N=4 no-blocking instance it stops at U = 1 everywhere with a total queue
% length of 7 for a population of 4, i.e. THM2 violated: a point of that shape
% is not a bound, it is a failed solve, and reporting it as one is the defect
% this check exists to prevent.
[cChk, ceqChk] = sub_qrfcon(xopt,q,f,M,MR,MM,MM1,ZZ,ZM,BB,F,N);
resNlp = 0;
if ~isempty(ceqChk), resNlp = max(resNlp, full(max(abs(ceqChk(:))))); end
if ~isempty(cChk), resNlp = max(resNlp, full(max(cChk(:)))); end
if resNlp > 1e-6
    line_error(mfilename, ...
        ['fmincon did not reach a feasible point (max constraint residual ' ...
         '%.3e). The qrf_bas_bethe bound is not defined for this instance.'], resNlp);
end
[p2opt,eopt] = sub_qrfvar(xopt);

% UTILIZATION, not occupancy. UN used to sum the diagonal p2 over ALL blocking
% configurations, i.e. P(n_i >= 1) with the BLOCKED ones included. A blocked BAS
% server holds a job it has already finished and does no work, so that is
% occupancy: on the M=2, N=3, F=[2 3] cyclic model it reported U2 = 1 where the
% exact utilization is 7/15, which the LP twin qrf_bas already returns. The e
% variables carry the right quantity -- UEFF pins e(i,ki) to the mass with
% n_i >= 1 in the configurations where i is NOT blocked. NO 1/M here, unlike
% qrf_bas.m: this formulation emits UEFF as one row per (j,i,ki), leaving e
% unscaled, where the LP sums j inside a single row and leaves e scaled by M.
for ti=1:M
    UN(ti) = sum(eopt(ti,1:K(ti)));
    QN(ti) = 0;
    for m=1:MR
        for ni=1+(1:F(ti))
            for ki=1:K(ti)
                QN(ti) = QN(ti) + (ni-1)*p2opt(ti,ni,ki,ti,ni,ki, m); % rescaled back ni: the loop index is population + 1
            end
        end
    end
end

% add LB>=0 to both p2 and e
% LINEAR PROGRAMMING: UTILIZATION UPPER BOUND AT QUEUE 1
%minimize U1min: sum {m in 1..MR} sum {k in 1..K[1]} sum {n1 in 1..F[1]} p2[1,n1,k,1,n1,k,m]; 


    function fobj = bethe(x)
    % BETHE
    % lambda*sum_{i~=j} I(n_i;n_j) - sum_i H(n_i) at lambda = 1/M, over the
    % BAS decision vector. The MI block is qrf_bas_mmi's mmi() body scaled by
    % lambda; the entropy block is mem()'s with its population sum from n = 0.
        [p2,~] = sub_qrfvar(x);
        lam = 1/M;
        fobj = 0;
        for m = 1:MR
            for i = 1:M
                for ki = 1:K(i)
                    for j = 1:M
                        if i~=j
                            for kj = 1:K(j)
                                for ni = 1+(0:F(i))
                                    for nj = 1+(0:F(j))
                                        fobj = fobj + lam*p2(i,ni,ki,j,nj,kj,m)*(log(1e-6+p2(i,ni,ki,j,nj,kj,m))-log(1e-6+p2(i,ni,ki,i,ni,ki,m))-log(1e-6+p2(j,nj,kj,j,nj,kj,m)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        for m = 1:MR
            for i = 1:M
                for k = 1:K(i)
                    for ni = 1+(0:F(i))
                        fobj = fobj + p2(i,ni,k,i,ni,k,m)*log(1e-6 + p2(i,ni,k,i,ni,k,m));
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

    function [c,ceq] = sub_qrfcon(x,q,f,M,MR,MM,MM1,ZZ,ZM,BB,F,N)
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
                
        % subject to ZERO4 {j in 1..M, k in 1..K[j], nj in 0..N, m in 2..MR: f<>j}: sum {nf in 0..F[f]-1, h in 1..K[f]} p2[j,nj,k,f,nf,h,m]=0;
        for j = 1:M, for k =1:K(j), for nj = 1+(0:N), for m = 2:MR
            if f~=j
                ceq(end+1) = 0;
                for nf =1+ (0:(F(f)-1)), for h = 1:K(f) 
                        ceq(end) = ceq(end) + p2(j,nj,k,f,nf,h,m);                    
                end, end
            end
        end, end, end, end
        
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
        
        % subject to ZERO8 {nf in 1..F[f]-1, k in 1..K[f], i in 1..M, h in 1..K[i], ni in 0..N, m in 2..MR}: p2[f,nf,k,i,ni,h,m]=0;
        for nf = 1+(1:(F(f)-1)), for k =1:K(f), for i = 1:M, for h = 1:K(i), for ni = 1+(0:N), for m = 2:MR
            ceq(end+1) = p2(f,nf,k,i,ni,h,m);
        end, end, end, end, end, end                                
        
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
                for ni = 1+(0:(N-(nj-1))) % nj is the 1-based index, so the population is nj-1                    
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
        
        %subject to THM1old {i in 1..M, k in 1..K[i]}: sum {j in 1..M, h in 1..K[i]: h<>k and j==i} q[i,j,k,h]*e[i,k] =sum {j in 1..M, h in 1..K[i]:h<>k and j==i} q[i,j,h,k]*e[i,h];
%         for i = 1:M, for k =1:K(i)
%                 ceq(end+1) = 0;
%                 % LHS
%                 for j = 1:M, for h = 1:K(i)
%                         if h~=k && j==i
%                             ceq(end) = ceq(end) + q(i,j,k,h)*e(i,k);
%                         end
%                 end, end
%                 % RHS
%                 for j = 1:M, for h = 1:K(i)
%                     if h~=k && j==i 
%                         ceq(end) = ceq(end) - q(i,j,h,k)*e(i,h);
%                     end
%                 end, end
%         end, end
        
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
        
        %% THM30, THM3, THM3f: the marginal-balance families of
        %% qrboundsbas_skel.mod. WITHOUT THEM the only other family mentioning q
        %% is THM1, which vanishes when every station has one phase, so nothing
        %% ties the utilization to the service rates and the relaxation is the
        %% [0,1] box. Ported from the LP implementation in qrf_bas.m, which is
        %% validated against the AMPL model.

        %subject to THM30 {i in 1..M, u in 1..K[i]: i<>f}: arrivals into i at
        %population 0 balance departures from i at population 1.
        for i = 1:M
            if i == f, continue; end
            for ui = 1:K(i)
                ceq(end+1) = 0;
                for j = 1:M
                    if j == i || j == f, continue; end
                    for nj = 1:F(j), for kj = 1:K(j), for hj = 1:K(j), for m = 1:MR
                        if BB(m,j) == 0
                            ceq(end) = ceq(end) + q(j,i,kj,hj)*p2(j,nj+1,kj,i,1,ui,m);
                        end
                    end, end, end, end
                end
                for nj = 1:F(f), for kj = 1:K(f), for hj = 1:K(f), for m = 1:MR
                    if MM(m,1) ~= i
                        ceq(end) = ceq(end) + q(f,i,kj,hj)*p2(f,nj+1,kj,i,1,ui,m);
                    end
                end, end, end, end
                for j = 1:M
                    if j == i || j == f, continue; end
                    for nj = 0:F(j), for ki = 1:K(i), for hj = 1:K(j), for m = 1:MR
                        if BB(m,i) == 0
                            ceq(end) = ceq(end) - q(i,j,ki,ui)*p2(j,nj+1,hj,i,2,ki,m);
                        end
                    end, end, end, end
                end
                for nj = 0:(F(f)-1), for ki = 1:K(i), for hj = 1:K(f), for m = 1:MR
                    if BB(m,i) == 0
                        ceq(end) = ceq(end) - q(i,f,ki,ui)*p2(f,nj+1,hj,i,2,ki,m);
                    end
                end, end, end, end
                for m = 1:MR
                    if BB(m,i) == 1 && MM(m,1) == i
                        for kf = 1:K(f), for pf = 1:K(f), for w = 1:M
                            if w ~= f && w ~= i
                                ceq(end) = ceq(end) - q(f,w,kf,pf)*p2(f,F(f)+1,kf,i,2,ui,m);
                            end
                        end, end, end
                    end
                end
            end
        end

        %subject to THM3 {i in 1..M, ni in 0..(F[i]-1): i<>f}: the same balance
        %at population ni, phase-aggregated.
        for i = 1:M
            if i == f, continue; end
            for ni = 1:(F(i)-1)
                ceq(end+1) = 0;
                for j = 1:M
                    if j == i || j == f, continue; end
                    for nj = 1:F(j), for kj = 1:K(j), for hj = 1:K(j), for ui = 1:K(i), for m = 1:MR
                        if BB(m,j) == 0
                            ceq(end) = ceq(end) + q(j,i,kj,hj)*p2(j,nj+1,kj,i,ni+1,ui,m);
                        end
                    end, end, end, end, end
                end
                for nj = 1:F(f), for kj = 1:K(f), for hj = 1:K(f), for ui = 1:K(i), for m = 1:MR
                    if MM(m,1) ~= i
                        ceq(end) = ceq(end) + q(f,i,kj,hj)*p2(f,nj+1,kj,i,ni+1,ui,m);
                    end
                end, end, end, end, end
                for j = 1:M
                    if j == i || j == f, continue; end
                    for nj = 0:F(j), for ki = 1:K(i), for hi = 1:K(i), for uj = 1:K(j), for m = 1:MR
                        if BB(m,i) == 0
                            ceq(end) = ceq(end) - q(i,j,ki,hi)*p2(j,nj+1,uj,i,ni+2,ki,m);
                        end
                    end, end, end, end, end
                end
                for nj = 0:(F(f)-1), for ki = 1:K(i), for hi = 1:K(i), for uj = 1:K(f), for m = 1:MR
                    if BB(m,i) == 0
                        ceq(end) = ceq(end) - q(i,f,ki,hi)*p2(f,nj+1,uj,i,ni+2,ki,m);
                    end
                end, end, end, end, end
                for m = 1:MR
                    if BB(m,i) == 1 && MM(m,1) == i
                        for ki = 1:K(i), for kf = 1:K(f), for pf = 1:K(f), for w = 1:M
                            if w ~= f && w ~= i
                                ceq(end) = ceq(end) - q(f,w,kf,pf)*p2(f,F(f)+1,kf,i,ni+2,ki,m);
                            end
                        end, end, end, end
                    end
                end
            end
        end

        %subject to THM3f {i in 1..M, ni in 0..(F[i]-1): i==f}: the same balance
        %at the finite-capacity station itself.
        for ni = 0:(F(f)-1)
            ceq(end+1) = 0;
            for j = 1:M
                if j == f, continue; end
                for nj = 1:F(j), for kj = 1:K(j), for hj = 1:K(j), for uf = 1:K(f), for m = 1:MR
                    if BB(m,j) == 0
                        ceq(end) = ceq(end) + q(j,f,kj,hj)*p2(j,nj+1,kj,f,ni+1,uf,m);
                    end
                end, end, end, end, end
            end
            for j = 1:M
                if j == f, continue; end
                for nj = 0:F(j), for kf = 1:K(f), for hf = 1:K(f), for uj = 1:K(j)
                    ceq(end) = ceq(end) - q(f,j,kf,hf)*p2(j,nj+1,uj,f,ni+2,kf,1);
                end, end, end, end
            end
        end

        %subject to THM3I {i in 1..M, z in 0..(ZM-1): i==f}: balance across
        %blocking depth at the finite station, transcribed from qrf_bas.m.
        for z = 0:(ZM-1)
            ceq(end+1) = 0;
            for j = 1:M
                if j == f, continue; end
                for nj = 1:F(j), for kj = 1:K(j), for hj = 1:K(j), for uf = 1:K(f), for m = 1:MR
                    if BB(m,j) == 0 && ZZ(m) == z
                        ceq(end) = ceq(end) + q(j,f,kj,hj)*p2(j,nj+1,kj,f,F(f)+1,uf,m);
                    end
                end, end, end, end, end
            end
            for j = 1:M
                if j == f, continue; end
                for nj = 0:F(j), for kf = 1:K(f), for hf = 1:K(f), for uj = 1:K(j), for m = 1:MR
                    if ZZ(m) == z+1
                        ceq(end) = ceq(end) - q(f,j,kf,hf)*p2(j,nj+1,uj,f,F(f)+1,kf,m);
                    end
                end, end, end, end, end
            end
        end

        %subject to THM3L {m in 1..MR: ZZ[m]==ZM-1}: the deepest blocking
        %configuration closes onto the one MM1 names.
        for m = 1:MR
            if ZZ(m) ~= ZM - 1, continue; end
            ceq(end+1) = 0;
            for j = 1:M
                if j == f || BB(m,j) ~= 0 || MM1(m,j) <= 0, continue; end
                for nj = 1:F(j), for kj = 1:K(j), for hj = 1:K(j), for uf = 1:K(f)
                    ceq(end) = ceq(end) + q(j,f,kj,hj)*p2(j,nj+1,kj,f,F(f)+1,uf,m);
                end, end, end, end
            end
            for j = 1:M
                if j == f || BB(m,j) ~= 0 || MM1(m,j) <= 0, continue; end
                mp = MM1(m,j);
                for kf = 1:K(f), for uf = 1:K(f), for w = 1:M
                    if w ~= f
                        ceq(end) = ceq(end) - q(f,w,kf,uf)*p2(f,F(f)+1,kf,f,F(f)+1,kf,mp);
                    end
                end, end, end
            end
        end

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
%         for i = 1:M, for u = 1:K(i)
%                 if i~=f
%                     ceq(end+1)=0;
%                     % LHS
%                     for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for m = 1:MR 
%                         if j~=i && j~=f && BB(m,j)==0 
%                             ceq(end)= ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,1+0,u, m);
%                         end
%                     end, end, end, end, end     
%     
%                     for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for m = 1:MR 
%                         if j~=i && j==f && MM(m,1)~=i 
%                             ceq(end)= ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,1+0,u, m);
%                         end
%                     end, end, end, end, end
%                     % RHS
%                     for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(j), for m = 1:MR 
%                         if j~=i && j~=f && BB(m,i)==0 
%                             ceq(end)= ceq(end) - q(i,j,k,u)*p2(j,nj,h, i,1+0,k, m); 
%                         end
%                     end, end, end, end, end
%     
%                     for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(j), for m = 1:MR 
%                         if j~=i && j==f && BB(m,i)==0 && (nj-1)<F(j)  % rescaled back nj 
%                             ceq(end)= ceq(end) - q(i,j,k,u)*p2(j,nj,h, i,1+1,k, m); 
%                         end
%                     end, end, end, end, end
%     
%                     for j = 1:M, for nj = 1+(0:F(j)), for y = 1:K(j), for m = 1:MR 
%                         if j~=i && j==f && BB(m,i)==1 && (nj-1)==F(j) && MM(m,1)==i  % rescaled back nj 
%                             for p = 1:K(f), for w = 1:M
%                                 if w~=f && w~=i 
%                                     ceq(end)= ceq(end) - q(f,w,y,p)*p2(f,nj,y, i,1+1,u, m);
%                                 end
%                             end, end
%                         end
%                     end, end, end, end
%                     end % if
%         end, end
                
        % subject to THM3 {i in 1..M, ni in 0..(F[i]-1): i<>f}: 
        % sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j<>f and BB[m,j]==0} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % +  sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j==f and MM[m,1]<>i} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j<>f and BB[m,i]==0} sum {h in 1..K[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==0 and nj<F[j]} sum {h in 1..K[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, m] 
        % + sum {j in 1..M, nj in 0..F[j], k in 1..K[i], u in 1..K[j], m in 1..MR: j<>i and j==f and BB[m,i]==1 and nj==F[j] and MM[m,1]==i} sum {p in 1..K[f], w in 1..M: w<>f and w<>i} q[f,w,u,p]*p2[f,nj,u, i,ni+1,k, m] ;
%         for i = 1:M, for ni = 1+(0:(F(i)-1))
%                 if i~=f
%                 ceq(end+1)=0;
%                 % LHS                
%                 for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR 
%                     if j~=i && j~=f && BB(m,j)==0 
%                         ceq(end)= ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,ni,u, m); 
%                     end                   
%                 end, end, end, end, end, end
%                 for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR 
%                     if j~=i && j==f && MM(m,1)~=i 
%                         ceq(end)= ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,ni,u, m);
%                     end
%                 end, end, end, end, end, end
%                 % RHS
%                 for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
%                     if j~=i && j~=f && BB(m,i)==0 
%                         for h = 1:K(i) 
%                             ceq(end)= ceq(end) - q(i,j,k,h)*p2(j,nj,u, i,ni+1,k, m); 
%                         end
%                     end
%                 end, end, end, end, end
%                 for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
%                     if j~=i && j==f && BB(m,i)==0 && (nj-1)<F(j)   % rescaled back nj
%                         for h = 1:K(i) 
%                             ceq(end)= ceq(end) - q(i,j,k,h)*p2(j,nj,u, i,ni+1,k, m);
%                         end
%                     end
%                 end, end, end, end, end
%                 for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for u = 1:K(j), for m = 1:MR 
%                     if j~=i && j==f && BB(m,i)==1 && (nj-1)==F(j) && MM(m,1)==i  % rescaled back nj 
%                         for p = 1:K(f), for w = 1:M 
%                             if w~=f && w~=i 
%                                 ceq(end)= ceq(end) - q(f,w,u,p)*p2(f,nj,u, i,ni+1,k, m);
%                             end
%                         end, end
%                     end
%                 end, end, end, end, end
%                 end
%         end, end
                
        % subject to THM3f {i in 1..M, ni in 0..(F[i]-1): i==f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j<>f and BB[m,j]==0 and ni < F[i]} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % +  sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[i], m in 1..MR: j<>i and j==f and ni==F[i] and MM[m,1]==j} sum {w in 1..M: w<>f} q[j,i,k,h]*p2[j,nj,k, i,ni,u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[i], h in 1..K[i], u in 1..K[j]: j<>i and ni < F[i]} q[i,j,k,h]*p2[j,nj,u, i,ni+1,k, 1];
%         for i = 1:M, for ni = 1+(0:(F(i)-1))
%                 if i==f
%                     ceq(end+1) = 0;
%                     % LHS
%                     for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR
%                         if j~=i && j~=f && BB(m,j)==0 && (ni-1) < F(i)  % rescaled back ni
%                             ceq(end) = ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,ni,u, m);
%                         end
%                     end, end, end, end, end, end
%                     for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(i), for m = 1:MR
%                         if j~=i && j==f && (ni-1)==F(i) && MM(m,1)==j  % rescaled back ni
%                             for w = 1:M
%                                 if w~=f
%                                     ceq(end) = ceq(end) + q(j,i,k,h)*p2(j,nj,k, i,ni,u, m);
%                                 end
%                             end
%                         end
%                     end, end, end, end, end, end
%                     % RHS
%                     for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(i), for h = 1:K(i), for u = 1:K(j)
%                         if j~=i && (ni-1) < F(i)  % rescaled back ni
%                             ceq(end) = ceq(end) - q(i,j,k,h)*p2(j,nj,u, i,ni+1,k, 1);
%                         end
%                     end, end, end, end, end
%                 end % if
%         end, end
    
        % subject to THM3I {i in 1..M, z in 0..(ZM-1): i==f}: sum {j in 1..M, nj in 1..F[j], k in 1..K[j], h in 1..K[j], u in 1..K[f], m in 1..MR: j<>f and BB[m,j]==0 and ZZ[m]==z} q[j,f,k,h]*p2[j,nj,k, f,F[f],u, m] 
        % = sum {j in 1..M, nj in 0..F[j], k in 1..K[f], h in 1..K[f], u in 1..K[j], m in 1..MR: j<>f and ZZ[m]=z+1 } q[f,j,k,h]*p2[j,nj,u, f,F[f],k, m];
%         for i = 1:M, for z = 1+(0:(ZM-1))
%                 if i==f
%                 ceq(end+1)=0;
%                 % LHS                
%                     for j = 1:M, for nj = 1+(1:F(j)), for k = 1:K(j), for h = 1:K(j), for u = 1:K(f), for m = 1:MR 
%                         if j~=f && BB(m,j)==0 && ZZ(m)==z 
%                             ceq(end) = ceq(end) + q(j,f,k,h)*p2(j,nj,k, f,F(f),u, m);
%                         end
%                     end, end, end, end, end, end
%                 % RHS                            
%                     for j = 1:M, for nj = 1+(0:F(j)), for k = 1:K(f), for h = 1:K(f), for u = 1:K(j), for m = 1:MR 
%                         if j~=f && ZZ(m)==z+1  
%                             ceq(end) = ceq(end) - q(f,j,k,h)*p2(j,nj,u, f,F(f),k, m);
%                         end
%                     end, end, end, end, end, end
%                 end % if
%         end, end
                
        % subject to THM3L {m in 1..MR: ZZ[m]==ZM-1}: sum {j in 1..M, k in 1..K[j], u in 1..K[f], nj in 1..F[j], h in 1..K[j]: j<>f and BB[m,j]==0 and MM1[m,j]>0} q[j,f,k,h]*p2[j,nj,k, f,F[f],u, m] 
        % = sum {j in 1..M: j<>f and BB[m,j]==0 and MM1[m,j]>0}  sum {w in 1..M, k in 1..K[f], u in 1..K[f]:w<>f} q[f,w,k,u]*p2[f,F[f],k, f,F[f],k, MM1[m,j]];
%         for m = 1:MR 
%             if ZZ(m)==ZM-1
%                 ceq(end+1)=0;
%             % LHS
%             	for j = 1:M, for k = 1:K(j), for u = 1:K(f), for nj = 1+(1:F(j)), for h = 1:K(j) 
%                     if j~=f && BB(m,j)==0 && MM1(m,j)>0 
%                         ceq(end) = ceq(end) + q(j,f,k,h)*p2(j,nj,k, f,F(f),u, m);                             
%                     end
%                 end, end, end, end, end
%             % RHS                        
%                 for j = 1:M
%                     if j~=f && BB(m,j)==0 && MM1(m,j)>0  
%                         for w = 1:M, for k = 1:K(f), for u = 1:K(f)
%                             if w~=f 
%                                 ceq(end) = ceq(end) - q(f,w,k,u)*p2(f,F(f),k, f,F(f),k, MM1(m,j));
%                             end
%                         end, end, end
%                     end 
%                 end
%             end
%         end
        
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
        