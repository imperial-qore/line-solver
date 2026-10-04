%{
%{
 % @file pfqn_nre.m
 % @brief Normalizing constant via saddle-tilted Edgeworth (NRE) approximation.
%}
%}

%{
%{
 % @brief Normalizing constant via saddle-tilted Edgeworth (NRE) approximation.
 %
 % Evaluates the Norlund-Rice integral form of the limited-load-dependent
 % normalizing constant (Casale-Harrison-Ong, Perform. Eval. 152, 2021,
 % Thm. 6) by steepest descent instead of a Laplace approximation on the
 % untilted contour, as done by pfqn_nrl and pfqn_nrp. Two corrections are
 % applied over those methods:
 %
 %   1. The integrand is invariant under t -> t + c*1, since h is homogeneous
 %      of degree sum(N) in the class variables and that degree cancels
 %      against exp(-1i*N*t). The redundant direction is quotiented out, so
 %      the integral is (R-1)-dimensional, not R-dimensional. pfqn_nrl and
 %      pfqn_nrp integrate over R dimensions and let the logit/probit
 %      Jacobian supply curvature along the null direction, which is an
 %      artifact of the change of variables rather than of the integrand.
 %   2. The contour radii are tilted per class to the saddle point, i.e. to
 %      the X solving X_r*dlog(h)/dX_r = N_r, so that the origin is a
 %      stationary point of the phase. On the untilted contour X = 1 it is
 %      not, which is the leading source of bias in pfqn_nrl and pfqn_nrp.
 %
 % A second-order Edgeworth term built from the third and fourth cumulants
 % of the tilted distribution is then added, giving a relative error of
 % O(1/sum(N)^2) instead of O(1) for a heuristic Gaussian fit.
 %
 % All integrand evaluations are at real positive demands, so unlike
 % pfqn_nrl and pfqn_nrp this method needs no complex arithmetic and runs
 % entirely in the log domain through pfqn_lldsingle, whose cost is linear
 % rather than quadratic in the population wherever the rates settle.
 %
 % Cost is O(I*R^2 + R^4) evaluations of a single-class LLD normalizing
 % constant, hence polynomial in the number of classes.
 %
 % @fn pfqn_nre(L, N, Z, alpha, options)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR or DxR).
 % @param alpha Load-dependent rate matrix (Mx sum(N)).
 % @param options Solver options.
 % @param vfix Optional tilt to use instead of solving the saddle-point
 %        equation. Supplying the tilt obtained at a nearby population makes
 %        numerator and denominator of a ratio share one expansion point, which
 %        is the Tierney-Kadane arrangement; omit it for the standard estimator.
 % @return lG Logarithm of normalizing constant.
 % @return lGs Logarithm of the saddlepoint term alone, i.e. of lG with the
 %         Edgeworth factor omitted, so that lG-lGs is the correction.
 % @return vsad The tilt actually used, for reuse at a nearby population.
 % @return G Normalizing constant.
%}
%}
function [lG,G,lGs,vsad] = pfqn_nre(L,N,Z,alpha,options,vfix)
N = N(:)';
if sum(N)<0
    lG = -Inf; G = 0; lGs = lG; vsad = [];
    return
end
if sum(N)==0
    lG = 0; G = 1; lGs = lG; vsad = [];
    return
end
if nargin<3
    Z = zeros(1,length(N));
end
if nargin<5 || isempty(options)
    options = SolverNC.defaultOptions;
end
Nt = sum(N);
if size(alpha,2) < Nt
    line_error(mfilename,'the load-dependent rate matrix must have at least sum(N) columns.');
end
alpha = alpha(:,1:Nt); % trim so that every rate used downstream is positive
if ~isempty(Z) && sum(Z(:))>0
    L = [L; sum(Z,1)];
    alpha(end+1,1:Nt) = 1:Nt; % the delay is an infinite server station
end
[M,R] = size(L);
if M==1
    [~,lG] = pfqn_gld(L,N,alpha,options);
    G = exp(lG); lGs = lG; vsad = [];
    return
end

% scale demands in [0,1] per class, the residual factor is exact by homogeneity
Lmax = max(L,[],1);
Lmax(Lmax<=0) = 1;
L = L ./ repmat(Lmax,M,1);
lGscale = N*log(Lmax)';

if R==1
    % coefficient extraction is the identity in a single class
    lG = pfqn_lldsingle(L,Nt,alpha,options) + lGscale;
    G = exp(lG); lGs = lG; vsad = [];
    return
end

d = R - 1; % dimension of the quotient torus

hstep = 2e-2;    % finite-difference step, results are flat over [5e-3,5e-2]
% The stencil is memoised on the offset, which stays within [-4,4] along each
% dimension and leaves the origin in at most four dimensions at once, so the
% points it visits number O(d^4). Indexing them densely takes span^d entries,
% which is what used to stop the method at eight classes; a map costs space
% proportional to the work instead and leaves the arithmetic untouched.
mkcache = @() containers.Map('KeyType','char','ValueType','double');
cacheVal = mkcache();
vbase = zeros(1,d);
Nd = N(1:d);

% ---- saddle point: minimise the convex F(v) = K(v) - Nd*v ----
% A tilt supplied by the caller is used as given, so that a ratio of two
% constants can be expanded about one common point rather than two.
if nargin >= 6 && ~isempty(vfix)
    vbase = vfix(:)';
    converged = true;
else
converged = false;
for it = 1:100
    cacheVal = mkcache();
    grad = zeros(d,1);
    hess = zeros(d,d);
    for a = 1:d
        grad(a) = (cgf(unitoff(a,1)) - cgf(unitoff(a,-1)))/(2*hstep) - Nd(a);
    end
    for a = 1:d
        for b = 1:d
            hess(a,b) = (cgf(unitoff(a,1)+unitoff(b,1)) - cgf(unitoff(a,1)+unitoff(b,-1)) ...
                - cgf(unitoff(a,-1)+unitoff(b,1)) + cgf(unitoff(a,-1)+unitoff(b,-1)))/(4*hstep^2);
        end
    end
    step = -(hess\grad)';
    F0 = cgf(zeros(1,d)) - Nd*vbase';
    tau = 1;
    while tau > 1e-10
        vtry = vbase + tau*step;
        if cgfat(vtry) - Nd*vtry' <= F0
            break
        end
        tau = tau/2;
    end
    vbase = vbase + tau*step;
    % Newton converges to the root of the differenced gradient, whose own
    % O(hstep^2) bias puts any absolute gradient target out of reach
    if norm(tau*step) < 1e-10
        converged = true;
        break
    end
end
end
if ~converged
    line_warning(mfilename,'the saddle point search did not converge, the estimate may be inaccurate.');
end
vsad = vbase;

% ---- cumulants of the tilted distribution at the saddle ----
cacheVal = mkcache();
K0 = cgf(zeros(1,d));
Sigma = zeros(d,d);
for a = 1:d
    for b = 1:d
        Sigma(a,b) = (cgf(unitoff(a,1)+unitoff(b,1)) - cgf(unitoff(a,1)+unitoff(b,-1)) ...
            - cgf(unitoff(a,-1)+unitoff(b,1)) + cgf(unitoff(a,-1)+unitoff(b,-1)))/(4*hstep^2);
    end
end
Sigma = (Sigma + Sigma')/2;
if min(eig(Sigma)) <= 0
    line_error(mfilename,'the tilted covariance is singular, a class has no demand at any station.');
end

k3 = zeros(d,d,d);
for a = 1:d
    for b = 1:d
        for c = 1:d
            acc = 0;
            for s = 0:7
                sg = 1 - 2*[bitget(s,1),bitget(s,2),bitget(s,3)];
                acc = acc + prod(sg)*cgf(sg(1)*unitoff(a,1) + sg(2)*unitoff(b,1) + sg(3)*unitoff(c,1));
            end
            k3(a,b,c) = acc/(8*hstep^3);
        end
    end
end

k4 = zeros(d,d,d,d);
for a = 1:d
    for b = 1:d
        for c = 1:d
            for e = 1:d
                acc = 0;
                for s = 0:15
                    sg = 1 - 2*[bitget(s,1),bitget(s,2),bitget(s,3),bitget(s,4)];
                    acc = acc + prod(sg)*cgf(sg(1)*unitoff(a,1) + sg(2)*unitoff(b,1) ...
                        + sg(3)*unitoff(c,1) + sg(4)*unitoff(e,1));
                end
                k4(a,b,c,e) = acc/(16*hstep^4);
            end
        end
    end
end

% ---- second-order Edgeworth factor, see the derivation in the header ----
S = inv(Sigma);
rho4 = S(:)'*reshape(k4,[d*d,d*d])*S(:);
u = reshape(S(:)'*reshape(k3,[d*d,d]),1,d);
rhoA = u*S*u';
T = reshape(S*reshape(k3,d,d*d),[d d d]);
T = permute(reshape(S*reshape(permute(T,[2 1 3]),d,d*d),[d d d]),[2 1 3]);
T = permute(reshape(S*reshape(permute(T,[3 1 2]),d,d*d),[d d d]),[2 3 1]);
rhoB = sum(k3(:).*T(:));
corr = 1 + rho4/8 - (3*rhoA + 2*rhoB)/24;
if corr <= 0
    line_warning(mfilename,'the Edgeworth correction is non-positive, falling back on the saddlepoint term.');
    corr = 1;
end

lGs = K0 - Nd*vbase' - (d/2)*log(2*pi) - 0.5*log(det(Sigma)) + lGscale;
lG = lGs + log(corr);
G = exp(lG);

    function y = cgf(off)
        % cumulant generating function at vbase+off*hstep, memoised on the stencil
        key = char(off+5);   % offsets lie in [-4,4], so off+5 is a valid char code
        if isKey(cacheVal,key)
            y = cacheVal(key);
        else
            y = cgfat(vbase + off*hstep);
            cacheVal(key) = y;
        end
    end

    function y = cgfat(v)
        % log of the single-class LLD constant at the class tilt X=[exp(v),1]
        y = pfqn_lldsingle(L*[exp(v(:));1], Nt, alpha, options);
    end

    function off = unitoff(a,sgn)
        off = zeros(1,d);
        off(a) = sgn;
    end
end
