%{
%{
 % @file pfqn_amvasjn.m
 % @brief Approximate MVA of closed networks with shortest-job-next stations.
%}
%}

function [XN,QN,UN,CN,WX,it] = pfqn_amvasjn(L,N,Z,scv,sjnset,V,options)
%{
%{
 % @brief Fixed-point (Bard-Schweitzer) counterpart of pfqn_mvasjn for closed
 %        networks with non-preemptive shortest-job-next (SJN/SJF) stations.
 %
 %        pfqn_mvasjn carries the conditional waiting time profile W(x,n) over
 %        the whole population lattice, which costs prod(N+1) steps and rules
 %        the method out for large populations. The closure used here rests on
 %        the observation that
 %
 %          lam_k(n) W_k(x,n) f_k(x) dx
 %
 %        is the mean number of queued class-k customers whose service
 %        requirement lies in (x, x+dx), that is, the queue length resolved by
 %        job size. Schweitzer's assumption is applied to that density rather
 %        than to its integral: removing one customer of class r scales the
 %        class-r size-resolved queue length by (N_r-1)/N_r and leaves the
 %        other classes unchanged. Integrating over x recovers the usual
 %        Schweitzer rule for the aggregate queue lengths, so the closure is
 %        the exact analogue of the one applied at the ordinary stations, and
 %        the two are used together consistently.
 %
 %        The unknowns are therefore the profiles W_r(x,N) on the quadrature
 %        grid together with the queue lengths, and they are found by
 %        successive substitution. Cost per iteration is O(M R ns), against
 %        the prod(N+1) M R ns of the exact recursion, and the population may
 %        be arbitrarily large. What is given up is the population dependence
 %        of the *shape* of W(x): the closure lets its level scale but keeps
 %        its shape fixed, whereas the true profile stiffens with the load
 %        because the denominator 1 - sum_k lam_k theta_k(x) sharpens. The
 %        error therefore concentrates at high utilization, where the SJN
 %        approximation is already at its weakest, and pfqn_mvasjn should be
 %        preferred whenever the lattice is affordable.
 %
 %        The response time equation, the two-moment branching-Erlang fit of
 %        the size distribution, the quadrature grid, the analytic tail and
 %        the pooled and priority multiclass readings are all shared with
 %        pfqn_mvasjn; see that function and the reference for their
 %        derivation.
 %
 %        Reference: K. Kant, "MVA approximations for SJN scheduling",
 %        Performance Evaluation 15(1):41-61, 1992. The bidirectional use of
 %        the priority equations, of which this is the limiting form, is
 %        discussed in section 3.1 of that paper.
 % @fn pfqn_amvasjn(L, N, Z, scv, sjnset, V, options)
 % @param L Service demand matrix (M x R) of the queueing stations.
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R). Default: zeros.
 % @param scv Squared coefficient of variation of the service times (M x R). Default: ones.
 % @param sjnset Indices of the stations scheduling by SJN. Default: none.
 % @param V Visit ratios (M x R), so that the per-visit service time is L./V. Default: ones.
 % @param options Struct with fields ns (grid subdivisions, default 32),
 %        Lfactor (grid extent in mean service times, default 8), prio
 %        (1 x R priority levels, default [] for the pooled reading), tol
 %        (default 1e-8) and iter_max (default 1000).
 % @return XN System throughput (1 x R).
 % @return QN Mean queue length (M x R).
 % @return UN Utilization (M x R).
 % @return CN Residence time (M x R).
 % @return WX Struct array with the converged conditional waiting times:
 %        WX.station, WX.x (grid), WX.W (ns+1 x R) and WX.tail (R x 3 tail
 %        parameters a, b, c).
 % @return it Number of iterations performed.
%}
%}
% [XN,QN,UN,CN,WX,IT] = PFQN_AMVASJN(L,N,Z,SCV,SJNSET,V,OPTIONS)

if nargin < 3, Z = []; end
if nargin < 4, scv = []; end
if nargin < 5, sjnset = []; end
if nargin < 6, V = []; end
if nargin < 7, options = struct(); end
[M,R,N,Z,scv,sjnset,V,S,options] = sjn_args(mfilename,L,N,Z,scv,sjnset,V,options);
prio = options.prio;
useprio = ~isempty(prio);

ns = options.ns;
ngrid = ns + 1;
nsjn = length(sjnset);
G = cell(1,nsjn);
for q = 1:nsjn
    G{q} = sjn_setup(S(sjnset(q),:), scv(sjnset(q),:), ns, options.Lfactor);
end

% start from the product-form Schweitzer solution: a light-load guess would put the
% deflated utilization above one and the SJN denominator has no solution there
[XN,QN,UN,CN] = pfqn_bs(L,N,Z,options.tol,options.iter_max);
XN = reshape(XN,1,R);
W = cell(1,nsjn);
P = cell(1,nsjn);
Iinf = cell(1,nsjn);
T = cell(1,nsjn);
for q = 1:nsjn
    W{q} = zeros(ngrid,R);
    P{q} = zeros(ngrid,R);
    Iinf{q} = zeros(1,R);
    T{q} = zeros(R,3);
end
it = 0;
capped = false;
converged = false;
while ~converged && it < options.iter_max
    it = it + 1;
    Xit = XN; Qit = QN; Uit = UN; Cit = CN;
    Wit = W; Pit = P; Iit = Iinf; Tit = T;
    for r = 1:R
        if N(r) == 0
            continue
        end
        beta = ones(1,R);
        beta(r) = (N(r) - 1) / N(r);
        Cr = zeros(M,1);
        for m = 1:M
            q = find(sjnset == m, 1);
            if isempty(q)
                Cr(m) = L(m,r) * (1 + sum(beta .* QN(m,:)));
                continue
            end
            st = struct('lam', XN .* V(m,:), 'U', UN(m,:), 'Q', QN(m,:), ...
                'W', W{q}, 'phi', P{q}, 'phiinf', Iinf{q});
            [Cr(m), Wprof, phiprof, phiinf, tailpar] = sjn_station(mfilename, m, r, G{q}, ...
                S(m,:), scv(m,:), V(m,:), st, beta, useprio, prio);
            Wit{q}(:,r) = Wprof;
            Pit{q}(:,r) = phiprof;
            Iit{q}(r) = phiinf;
            Tit{q}(r,:) = tailpar;
        end
        Cit(:,r) = Cr;
    end
    [Cit, Xit, kappa, bound] = sjn_cap(mfilename, Cit, L, N, Z, sjnset, options.umax);
    if bound
        capped = true;
        for q = 1:nsjn
            Wit{q} = kappa(q) * Wit{q};
            Pit{q} = kappa(q) * Pit{q};
            Iit{q} = kappa(q) * Iit{q};
            Tit{q}(:,1:2) = kappa(q) * Tit{q}(:,1:2);
        end
    end
    Qit = repmat(Xit,M,1) .* Cit;
    Uit = repmat(Xit,M,1) .* L;
    delta = max(abs(Qit(:) - QN(:)));
    for q = 1:nsjn
        delta = max(delta, max(abs(Wit{q}(:) - W{q}(:))));
    end
    XN = Xit; QN = Qit; UN = Uit; CN = Cit;
    W = Wit; P = Pit; Iinf = Iit; T = Tit;
    converged = delta < options.tol;
end
if ~converged
    line_warning(mfilename,sprintf('the SJN fixed point did not converge in %d iterations, residual %g',options.iter_max,delta));
end
if capped
    line_warning(mfilename,sprintf(['the utilization cap of %g was binding at an SJN station: the station is in the\n' ...
        'starvation regime, where long jobs are held back and the arrival theorem is badly violated.\n' ...
        'The results are stable but their accuracy is not warranted, use SolverCTMC or SolverLDES there.'],options.umax));
end

WX = struct('station',{},'x',{},'W',{},'tail',{});
for q = 1:nsjn
    WX(q).station = sjnset(q);
    WX(q).x = G{q}.x;
    WX(q).W = W{q};
    WX(q).tail = T{q};
end
end
