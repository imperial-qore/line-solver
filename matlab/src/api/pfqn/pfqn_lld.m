%{
%{
 % @file pfqn_lld.m
 % @brief Normalizing constant of a multiclass LIMITED load-dependent closed model.
%}
%}

%{
%{
 % @brief Normalizing constant of a multiclass LIMITED load-dependent closed model.
 %
 % Same recursion, same arithmetic and the same result as pfqn_gld, but with the
 % rate shift saturated at the LIMITED LOAD-DEPENDENCE threshold, which turns the
 % recursion's state space finite and lets it be memoised. This is what
 % pfqn_lldsingle does to pfqn_gldsingle, one level up: there the rate offset is
 % an index into a table, here it is the shift that pfqn_mushift applies.
 %
 % pfqn_gld peels station M and advances its rate lattice one job at a time,
 %
 %   g(m,n,j) = g(m-1,n,0) + sum_r L(m,r)/alpha_m(j+1) * g(m,n-e_r,j+1)
 %
 % with j the number of shifts row m has taken, so that pfqn_mushift's leading
 % element is alpha_m(j+1). Once j >= s_m-1, where s_m is the population past
 % which alpha_m stays constant, every remaining entry of the row is
 % alpha_m(s_m) and a further shift LEAVES THE ROW UNCHANGED over the columns
 % the recursion can still read. Saturating j at s_m-1 therefore returns the
 % same value and makes the state (m, n, j) repeat, at which point one memo
 % answers what pfqn_gld recomputes down an exponential tree.
 %
 % COST. The state space is M * prod_r(N_r+1) * max_k s_k, against the
 % unmemoised binary recursion of pfqn_gld, which revisits the same states
 % exponentially often. Without the saturation a memo would still be bounded,
 % but by M * prod_r(N_r+1) * (|N|+1): the threshold is what replaces the
 % population by the server count, exactly as in pfqn_lldsingle.
 %
 % There is no gain on a station whose rates never settle, an infinite server
 % alpha(n)=n being the usual case: it gets s_k = |N| and the memo keeps its
 % population axis. The saving is over the OTHER stations.
 %
 % Every shortcut of pfqn_gld is kept and evaluated on the same materialised
 % arguments, so the two agree to the last bit rather than merely to a
 % tolerance: the node reached at (m,n,j) sees the demand block L(1:m,:) and a
 % rate block whose read range 1..sum(n) is entrywise the one pfqn_gld would
 % have built there.
 %
 % @fn pfqn_lld(L, N, mu, options)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param mu Load-dependent rate matrix (Mx sum(N)), alpha_i(j) = mu(i,j); default all ones.
 % @param options Solver options.
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
 % @return s Detected per-station thresholds (Mx1), alpha_i(n)=alpha_i(s_i) for n>=s_i.
%}
%}
function [G,lG,s]=pfqn_lld(L,N,mu,options)
% [G,LG]=PFQN_LLD(L,N,MU,OPTIONS)

[M,R]=size(L);

% As in pfqn_gld, the mu default has to precede every branch that reads mu.
if nargin<3 || isempty(mu)
    mu=ones(M,sum(N));
end
if nargin<4
    options = SolverNC.defaultOptions;
end

Ntot = sum(N);
ncols = size(mu,2);

% ---- s_k: the smallest offset past which the rate row is constant ----
% Equality is tested first so that an infinite rate, which the recursion admits
% and zeroes through a division by Inf, ties with itself instead of producing
% Inf-Inf. The tolerance is then confined to FINITE pairs: at tail=Inf the bound
% eps*max(abs(tail),1) is itself Inf and abs(prev-Inf)<=Inf would tie every
% finite rate to it, collapsing the row on a false tie.
s = ones(M,1);
if ncols >= 1
    for m=1:M
        s(m) = ncols;
        tail = mu(m,ncols);
        for n=ncols:-1:2
            if mu(m,n-1)==tail || (isfinite(tail) && isfinite(mu(m,n-1)) ...
                    && abs(mu(m,n-1)-tail) <= eps*max(abs(tail),1))
                s(m) = n-1;
            else
                break
            end
        end
    end
end
s = max(s,1);

memo = containers.Map('KeyType','char','ValueType','double');
G = node(M, N, 0);
lG = log(G);

% ------------------------------------------------------------------
    function g = node(m, n, j)
        % The value pfqn_gld would return at demands L(1:m,:), population n and
        % row m shifted j times. j is saturated by the caller.
        key = char([m+1, j+1, n+1]);
        if isKey(memo,key)
            g = memo(key);
            return
        end
        [Lm, mum] = materialize(m, n, j);
        g = gld_node(Lm, n, mum, m, j);
        memo(key) = g;
    end

% ------------------------------------------------------------------
    function [Lm, mum] = materialize(m, n, j)
        % Rows 1..m-1 are never shifted, row m is shifted j times, and every row
        % is truncated by one column per job already placed, exactly as
        % pfqn_mushift leaves them.
        cols = ncols - (Ntot - sum(n));
        Lm = L(1:m,:);
        mum = zeros(m, max(cols,0));
        if cols > 0
            if m > 1
                mum(1:(m-1),:) = mu(1:(m-1),1:cols);
            end
            mum(m,:) = mu(m,(j+1):(j+cols));
        end
    end

% ------------------------------------------------------------------
    function g = gld_node(Lm, n, mum, m, j)
        % pfqn_gld's own cascade, evaluated on the materialised block.
        if m==1
            if any(n>0 & abs(Lm(1,:))==0)
                g = 0;
                return
            end
            g = exp(factln(sum(n)) - sum(factln(n)) + n(Lm(1,:)>0)*log(Lm(1,Lm(1,:)>0))' ...
                - sum(log(mum(1,1:sum(n)))));
            return
        end
        if R==1
            g = exp(pfqn_lldsingle(Lm,n,mum));
            return
        end
        if isempty(Lm)
            g = 0;
            return
        end
        isLoadDep = false;
        isInfServer = false(1,m);
        for ist=1:m
            row = mum(ist,1:sum(n));
            % ELEMENTWISE &, as pfqn_gld has it: a node with no jobs left leaves
            % row empty, min([])==1 is [] rather than a logical scalar, and &&
            % rejects that where & yields [] and the branch falls through to the
            % all([]) test, which is the arm pfqn_gld takes there
            if min(row)==1 & max(row)==1 %#ok<AND2>
                isInfServer(ist) = false;
            elseif all(row == 1:sum(n))
                isInfServer(ist) = true;
            else
                isLoadDep = true;
            end
        end
        if ~isLoadDep
            Lli = Lm(~isInfServer,:);
            if isempty(Lli)
                Lli = 0*n;
            end
            Zli = Lm(isInfServer,:);
            if isempty(Zli)
                Zli = 0*n;
            end
            opts = options;
            opts.method = 'exact';
            g = exp(pfqn_nc(zeros(1,R), Lli, n, sum(Zli,1), opts));
            return
        end
        if m==0
            g = 0;
            return
        end
        if sum(n==zeros(1,R))==R
            g = 1;
            return
        end
        % The recursion itself, memoised through node(). The shift saturates at
        % s(m)-1, past which the row the child would see is the one it sees now.
        g = node(m-1, n, 0);
        for r=1:R
            if n(r)>0
                g = g + (Lm(m,r)/mum(m,1)) * node(m, oner(n,r), min(j+1, s(m)-1));
            end
        end
    end
end
