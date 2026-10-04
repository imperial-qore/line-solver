%{
%{
 % @file pfqn_mvasjn.m
 % @brief Mean value analysis of closed networks with shortest-job-next stations.
%}
%}

function [XN,QN,UN,CN,WX] = pfqn_mvasjn(L,N,Z,scv,sjnset,V,options)
%{
%{
 % @brief Approximate MVA for closed queueing networks in which a subset of
 %        the single-server stations schedules non-preemptively by shortest
 %        job next (SJN/SJF), the job size being known on arrival.
 %
 %        The SJN station is modelled by the conditional waiting time W(x,n)
 %        of a tagged customer whose service requirement is x, obtained from
 %        the arrival theorem as the sum of the residual life of the job in
 %        service, the work of the queued jobs that will be served before the
 %        tagged one, and the work of the jobs that overtake it while it waits
 %        (Kant 1992, eqs. 1-7):
 %
 %          W(x,n) = [ (1+CV^2) s U(n-1)/2 + X(n-1) phi(x,n-1) ]
 %                   / [ 1 - X(n-1) theta(x) ],
 %          theta(x) = int_0^x  t f(t) dt,
 %          phi(x,n) = int_0^x  W(t,n) t f(t) dt,
 %          R(n)     = s + int_0^inf W(x,n) f(x) dx.
 %
 %        The recursion is explicit: W(.,n) needs only phi(.,n-1), so it is
 %        carried alongside the population recursion of exact MVA. This is the
 %        unidirectional scheme of the reference and steps over the whole
 %        population lattice; pfqn_amvasjn is the fixed-point counterpart that
 %        trades the lattice for a Schweitzer closure on the same profile.
 %
 %        Two multiclass readings of "shortest job" are supported and selected
 %        by options.prio:
 %
 %        - pooled (options.prio empty, the default): a job of any class with
 %          service requirement below x precedes the tagged job, so the sums
 %          over classes run over all of them. This is the direct multiclass
 %          reading of the derivation above and collapses to eq. (6) of the
 %          reference for a single class.
 %        - priority (options.prio a vector of distinct levels, 1 = highest):
 %          classes are non-preemptively prioritised at the SJN station and
 %          SJN applies only within a class, which is method A of Kant 1992,
 %          eq. (21). Method B of that paper, which evaluates the denominator
 %          at the non-integral population n - Q(n), is not implemented here;
 %          pfqn_nintmva provides the fractional-population recursion it needs.
 %
 %        The service time density at an SJN station is not an input: only its
 %        mean and squared coefficient of variation are, and the density is
 %        reconstructed by the two-moment branching-Erlang fit the reference
 %        prescribes (see sjn_fit). This makes theta(x) and the tail integrals
 %        closed form.
 %
 %        The x-integrals are evaluated on a uniform grid of options.ns
 %        subdivisions of [0, Lx] with Lx = options.Lfactor times the largest
 %        mean service time at the station, by composite Simpson, exactly as
 %        the reference does; W(.,n) is needed at the next population step so
 %        quadrature rules that sample at arbitrary abscissae cannot be used.
 %        Beyond Lx the profile is closed by the analytic tail
 %        W(x,n) = a_n - b_n exp(-c_n (x - Lx)) of eqs. (11)-(14).
 %
 %        SJN starves long jobs as the station saturates, and the arrival
 %        theorem then fails badly: when the fraction of work due to jobs no
 %        longer than x reaches one the recursion has no solution and the
 %        function throws 'LINE:SjnStarvation' rather than returning a value,
 %        which is the behaviour reported in the reference. Capping the
 %        utilization instead is not an option here: the cap rescales the
 %        conditional waiting time profile that the next population step reads
 %        back, so the correction compounds along the lattice and the
 %        recursion oscillates. pfqn_amvasjn solves the same profile by a
 %        fixed point and does cap, the iteration being self-consistent.
 %
 %        Reference: K. Kant, "MVA approximations for SJN scheduling",
 %        Performance Evaluation 15(1):41-61, 1992.
 % @fn pfqn_mvasjn(L, N, Z, scv, sjnset, V, options)
 % @param L Service demand matrix (M x R) of the queueing stations.
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R). Default: zeros.
 % @param scv Squared coefficient of variation of the service times (M x R). Default: ones.
 % @param sjnset Indices of the stations scheduling by SJN. Default: none.
 % @param V Visit ratios (M x R), so that the per-visit service time is L./V. Default: ones.
 % @param options Struct with fields ns (grid subdivisions, default 32),
 %        Lfactor (grid extent in mean service times, default 8) and prio
 %        (1 x R priority levels, default [] for the pooled reading).
 % @return XN System throughput (1 x R).
 % @return QN Mean queue length (M x R).
 % @return UN Utilization (M x R).
 % @return CN Residence time (M x R).
 % @return WX Struct array with the conditional waiting times at population N:
 %        WX.station, WX.x (grid), WX.W (ns+1 x R) and WX.tail (R x 3 tail
 %        parameters a, b, c).
%}
%}
% [XN,QN,UN,CN,WX] = PFQN_MVASJN(L,N,Z,SCV,SJNSET,V,OPTIONS)

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

stride = cumprod([1, N(1:end-1)+1]);
npop = prod(N+1);
Xp = zeros(npop,R);
Qp = zeros(M,R,npop);
Up = zeros(M,R,npop);
Cp = zeros(M,R,npop);
Wp = cell(1,nsjn);   % conditional waiting profiles
Pp = cell(1,nsjn);   % phi on the grid
Ip = cell(1,nsjn);   % phi at infinity
Tp = cell(1,nsjn);   % tail parameters a, b, c
for q = 1:nsjn
    Wp{q} = zeros(ngrid,R,npop);
    Pp{q} = zeros(ngrid,R,npop);
    Ip{q} = zeros(R,npop);
    Tp{q} = zeros(R,3,npop);
end

for idx = 2:npop
    n = local_decode(idx, stride, N);
    Call = zeros(M,R);
    for r = 1:R
        if n(r) == 0
            continue
        end
        nprev = n;
        nprev(r) = nprev(r) - 1;
        iprev = 1 + sum(nprev .* stride);
        for m = 1:M
            q = find(sjnset == m, 1);
            if isempty(q)
                Call(m,r) = L(m,r) * (1 + sum(Qp(m,:,iprev)));
                continue
            end
            % the population step already supplies the neighbouring profile, so no deflation is needed
            beta = ones(1,R);
            st = struct('lam', Xp(iprev,:) .* V(m,:), 'U', reshape(Up(m,:,iprev),1,R), ...
                'Q', reshape(Qp(m,:,iprev),1,R), 'W', Wp{q}(:,:,iprev), 'phi', Pp{q}(:,:,iprev), ...
                'phiinf', Ip{q}(:,iprev)');
            [Call(m,r), Wprof, phiprof, phiinf, tailpar] = sjn_station(mfilename, m, r, G{q}, ...
                S(m,:), scv(m,:), V(m,:), st, beta, useprio, prio);
            Wp{q}(:,r,idx) = Wprof;
            Pp{q}(:,r,idx) = phiprof;
            Ip{q}(r,idx) = phiinf;
            Tp{q}(r,:,idx) = tailpar;
        end
    end
    [Call, Xn, ~, bound] = sjn_cap(mfilename, Call, L, n, Z, sjnset, options.umax);
    if bound
        % the cap has invalidated the profile the next population step reads back
        throw(MException('LINE:SjnStarvation', ...
            ['[%s.m] the utilization cap of %g was binding at an SJN station at population %s: ' ...
            'the station is in the starvation regime, where the conditional waiting time equation ' ...
            'has no solution and the population lattice no valid continuation. Use the Schweitzer ' ...
            'fixed point (pfqn_amvasjn, method ''amva''), SolverCTMC or SolverLDES.'], ...
            mfilename, options.umax, mat2str(n)));
    end
    Xp(idx,:) = Xn;
    Cp(:,:,idx) = Call;
    Qp(:,:,idx) = repmat(Xn,M,1) .* Call;
    Up(:,:,idx) = repmat(Xn,M,1) .* L;
end
XN = Xp(npop,:);
QN = Qp(:,:,npop);
UN = Up(:,:,npop);
CN = Cp(:,:,npop);

WX = struct('station',{},'x',{},'W',{},'tail',{});
for q = 1:nsjn
    WX(q).station = sjnset(q);
    WX(q).x = G{q}.x;
    WX(q).W = Wp{q}(:,:,npop);
    WX(q).tail = reshape(Tp{q}(:,:,npop),R,3);
end
end

function n = local_decode(idx, stride, N)
% Mixed-radix decoding of a linear lattice index into a population vector.
R = length(N);
n = zeros(1,R);
for r = 1:R
    n(r) = mod(floor((idx-1)/stride(r)), N(r)+1);
end
end
