%{
%{
 % @file pfqn_mcmc.m
 % @brief Markov chain Monte Carlo (regularization) estimator of class
 %        throughputs and queue lengths in a closed multiclass product-form
 %        network.
%}
%}

%{
%{
 % @brief Chen-O'Cinneide REGULARIZATION: throughputs and queue lengths of a
 %        closed multiclass product-form network, obtained by simulating a
 %        regularized network that has the same steady-state distribution.
 % @fn pfqn_mcmc(L, N, Z, s, options)
 % @param L Service demand matrix (stations x classes).
 % @param N Closed population vector (1 x classes), finite and integer.
 % @param Z Think time (1 x classes, or a matrix summed over its rows).
 % @param s Server counts per station (stations x 1), Inf for an
 %        infinite-server station; [] means all stations single-server.
 % @param options Solver options; .samples is the number of simulated service
 %        completions, .seed the RNG seed, .config.mcmc_batches the batch count
 %        and .config.mcmc_burnin the discarded warm-up fraction.
 % @return X Throughput estimate per class, X(r) = G(N-e_r)/G(N).
 % @return Q Mean queue length estimate per station and class.
 % @return ci Batch-means standard errors and two-sigma intervals for X and Q.
%}
%}
function [X,Q,ci] = pfqn_mcmc(L,N,Z,s,options)
% [X,Q,CI] = PFQN_MCMC(L,N,Z,S,OPTIONS)
%
% Markov chain Monte Carlo estimator of the class throughputs
% X(r) = G(N-e_r)/G(N) and of the mean queue lengths Q(i,r) of a CLOSED
% multiclass product-form (BCMP, no type changes) network, by the
% REGULARIZATION algorithm of
%
%   W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm
%   for Closed Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.
%
% The three steps of the paper are:
%
%  I.  CONSTRUCT THE REGULARIZED NETWORK. Write rho(i,r) for the surrogate
%      traffic intensity of class r at station i -- here the service demand,
%      since rho = lambda/mu is a visit ratio over a service rate -- and
%      rho(r) = sum_i rho(i,r). The regularized network has the same stations,
%      classes and populations, UNIT service rates at every station, the
%      processor-sharing discipline, and a routing matrix that depends on the
%      destination only,
%
%          P*(i->m | class r) = rho(m,r)/rho(r).
%
%      By Theorem 2.1 it is a REVERSIBLE chain with the SAME steady-state
%      distribution as the original network, and its throughputs satisfy
%      Theta*(r) = rho(r)*Theta(r).
%
%  II. SIMULATE IT at service-completion epochs. With Y(i,r) the number of
%      class-r jobs at station i, Y(i) their total and Psi_i(k)=min(s_i,k) the
%      number of busy servers,
%
%          r(i,r) = Y(i,r)/Y(i) * Psi_i(Y(i)),   r(r) = sum_i r(i,r),
%          r      = sum_i Psi_i(Y(i)),
%
%      the next completion is of class r at station i with probability
%      r(i,r)/r, and the conditional expected time to it is 1/r. Equation (10)
%      of the paper is the holding-time weighted ratio estimator
%
%          Theta*(r) = sum_t r(r,t)/r(t)  /  sum_t 1/r(t),
%
%      and the same weights give the time-average queue lengths, which need no
%      transformation at all because the two networks share their steady state.
%
% III. TRANSFORM BACK: X(r) = Theta*(r)/rho(r).
%
% Because P* forgets the station of origin and every station serves at unit
% rate, the regularized chain has neither the slowly mixing routing chain nor
% the customer-trapping slow station that make the original chain converge
% slowly. The paper proves O(N^2*M^3) mixing in two special cases (Section 4)
% and reports the general behaviour experimentally (Section 5).
%
% Delay (infinite-server) demand enters as ONE extra station with s = Inf and
% demand Z. Aggregating infinite-server stations that way is exact in the
% product form, since their joint term is multinomial in the per-class totals.
%
% Confidence: the run is split into non-overlapping batches (Schmeiser 1982,
% 30 by default, the count used in the tables of the paper), the batch means of
% the ratio estimator give a standard error, and CI reports the paper's
% two-sigma interval. The estimator is a ratio of correlated averages, so it
% carries an O(1/samples) bias on top of the initialization bias; the paper
% ignores both, this implementation additionally discards a warm-up fraction
% (10% by default).
%
% Parameters:
%   L       - (M x R) per-class service demands at the M queueing stations.
%   N       - (1 x R) closed population vector; finite and integer.
%   Z       - (1 x R) aggregated think times, or a matrix summed over its
%             rows; [] or zeros if the model has no delay.
%   s       - (M x 1) number of servers at each queueing station, Inf for an
%             infinite server; [] (default) means all stations single-server.
%   options - solver options (optional). Fields used:
%               .samples              simulated service completions (1e5);
%               .seed                 RNG seed for reproducibility (optional);
%               .config.mcmc_batches  batch count for the CIs (30);
%               .config.mcmc_burnin   discarded warm-up fraction (0.1).
%
% Returns:
%   X  - (1 x R) throughput estimates G(N-e_r)/G(N).
%   Q  - (M x R) mean queue lengths at the queueing stations; the delay
%        aggregate is not returned, the caller recovers it as Z.*X.
%   ci - struct with fields Xse, Xlo, Xhi (1 x R), Qse, Qlo, Qhi (M x R),
%        batches, samples and burnin. The intervals are two-sigma, as in the
%        tables of the paper.
%
% Example (the three single-server stations plus IS station of Example 5.1):
%   mu = [0.2 0.5 0.8]; sets = {[1 2 3],[1 2],[1 3],[2 3]};
%   L = zeros(3,4); Z = zeros(1,4);
%   for c=1:4, L(sets{c},c) = 1./mu(sets{c})'; if c>1, Z(c) = 1/0.5; end, end
%   X = pfqn_mcmc(L, 3*ones(1,4), Z, [], struct('samples',1e5,'seed',23000))
%
% See also PFQN_NC, PFQN_MCI, PFQN_LS, PFQN_IS, PFQN_MVA.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3, Z = []; end
if nargin < 4, s = []; end
if nargin < 5, options = struct(); end
if isstruct(s) % called as pfqn_mcmc(L,N,Z,options)
    options = s;
    s = [];
end

L(~isfinite(L)) = 0;
[M,R] = size(L);
N = N(:)';
if isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
Z = Z(:)';
Z(~isfinite(Z)) = 0;

if any(~isfinite(N))
    line_error(mfilename,'PFQN_MCMC requires a closed model, but the population vector has an infinite entry.');
end
if any(abs(N - round(N)) > GlobalConstants.FineTol)
    % the chain lives on the integer lattice sum_i Y(i,r) = N(r), so a
    % fractional population has no state space at all; this is not a matter of
    % accuracy and must not be rounded away silently
    line_error(mfilename,'PFQN_MCMC simulates a state space of integer populations, but N = %s is fractional. Use an asymptotic method (''le'', ''ble'', ''kt'') on fractional populations.', mat2str(N));
end
N = round(N);

X = zeros(1,R);
Q = zeros(M,R);
ci = struct('Xse',zeros(1,R),'Xlo',zeros(1,R),'Xhi',zeros(1,R), ...
    'Qse',zeros(M,R),'Qlo',zeros(M,R),'Qhi',zeros(M,R), ...
    'batches',0,'samples',0,'burnin',0);
if sum(N) == 0
    return
end

if isempty(s)
    s = ones(M,1);
else
    s = s(:);
    if numel(s) ~= M
        line_error(mfilename,'PFQN_MCMC: the server count vector has %d entries but L has %d stations.', numel(s), M);
    end
end

%% Step I: the regularized network
% Only the surrogate traffic intensities rho(i,r) enter the product form, and
% scaling a whole class column by a constant leaves the steady-state
% distribution unchanged, so the demands are used as they are.
rho = L;
rho(rho<0) = 0;
svec = s(:);
if any(Z > 0)
    rho(end+1,:) = Z; %#ok<AGROW>
    % VERTICAL CONCATENATION, not svec(end+1): ones(M,1) is an unambiguous column
    % only for M >= 2. At M == 1 it is 1x1, and MATLAB grows a scalar along the
    % SECOND dimension, so svec(end+1) yields a 1x2 ROW. min(svec,Ytot) below then
    % implicitly expands 1x2 against 2x1 into a 2x2 instead of erroring, and the
    % find() over cumsum(Psi) returns a linear index past the rows of Y. The line
    % above is immune because (end+1,:) names the dimension; this one was not.
    svec = [svec; Inf]; %#ok<AGROW>
end
Mx = size(rho,1);
rhoTot = sum(rho,1);
if any(N > 0 & rhoTot <= 0)
    line_error(mfilename,'PFQN_MCMC: class %d has a positive population but no demand anywhere in the network.', find(N > 0 & rhoTot <= 0, 1));
end
% routing of the regularized network, P*(m|r) = rho(m,r)/rho(r), held as one
% column of cumulative probabilities per class
Pstar = rho ./ repmat(max(rhoTot,realmin),Mx,1);
cumP = cumsum(Pstar,1);
cumP(Mx,:) = 1; % guard the last bin against a floating-point shortfall

% run length, batching and warm-up
samples = 1e5;
if isfield(options,'samples') && ~isempty(options.samples) && isfinite(options.samples)
    samples = max(1,round(options.samples));
end
nbatches = 30; % Schmeiser (1982), the count used in the tables of the paper
if isfield(options,'config') && isstruct(options.config) ...
        && isfield(options.config,'mcmc_batches') && ~isempty(options.config.mcmc_batches)
    nbatches = max(1,round(options.config.mcmc_batches));
end
burninFrac = 0.1;
if isfield(options,'config') && isstruct(options.config) ...
        && isfield(options.config,'mcmc_burnin') && ~isempty(options.config.mcmc_burnin)
    burninFrac = min(0.9,max(0,options.config.mcmc_burnin));
end
batchLen = max(1,floor(samples/nbatches));
samples = batchLen*nbatches;
nburn = round(burninFrac*samples);
if isfield(options,'seed') && ~isempty(options.seed)
    rng(options.seed);
end

% Initial state: spread each class over the stations it can occupy in the
% proportions P*(.|r), by largest remainder. That is the marginal the
% regularized network would have with no queueing, so it costs nothing and
% starts the chain far closer to stationarity than a single-station state.
Y = zeros(Mx,R);
for r = 1:R
    if N(r) == 0
        continue
    end
    target = N(r)*Pstar(:,r);
    base = floor(target);
    short = N(r) - sum(base);
    if short > 0
        [~,ord] = sort(target - base,'descend');
        base(ord(1:short)) = base(ord(1:short)) + 1;
    end
    Y(:,r) = base;
end
Ytot = sum(Y,2);

%% Step II: simulate the regularized network at service-completion epochs
xnum = zeros(nbatches,R);    % sum_t r(r,t)/r(t) within the batch
qnum = zeros(Mx,R,nbatches); % sum_t Y(t)/r(t)   within the batch
den = zeros(nbatches,1);     % sum_t 1/r(t)      within the batch
stale = true;                % the rates below need recomputation
Psi = zeros(Mx,1); rvec = zeros(1,R); w = 0;
for t = 1:(nburn+samples)
    if stale
        % (8)-(9): busy servers, per-class completion rates, total rate
        Psi = min(svec,Ytot);
        rw = Psi./Ytot;
        rw(Ytot == 0) = 0;
        rvec = (rw')*Y;
        w = 1/sum(Psi);
        stale = false;
    end
    if t > nburn
        b = floor((t-nburn-1)/batchLen) + 1;
        den(b) = den(b) + w;
        xnum(b,:) = xnum(b,:) + w*rvec;
        qnum(:,:,b) = qnum(:,:,b) + w*Y;
    end
    % Pick the completing station with probability Psi(i)/r, then the
    % completing class within it with probability Y(i,r)/Y(i); the product is
    % the r(i,r)/r of the paper, since sum_r Y(i,r)/Y(i)*Psi(i) = Psi(i).
    cPsi = cumsum(Psi);
    i = find(cPsi >= rand*cPsi(Mx),1);
    if isempty(i)
        i = find(Psi > 0,1,'last');
    end
    cY = cumsum(Y(i,:));
    r = find(cY >= rand*cY(R),1);
    if isempty(r)
        r = find(Y(i,:) > 0,1,'last');
    end
    % Route it. A self-transition leaves the state, hence the rates and the
    % weight, unchanged: skipping the recomputation is the saving described at
    % the end of Section 2 of the paper.
    m = find(cumP(:,r) >= rand,1);
    if ~isempty(m) && m ~= i
        Y(i,r) = Y(i,r) - 1;
        Y(m,r) = Y(m,r) + 1;
        Ytot(i) = Ytot(i) - 1;
        Ytot(m) = Ytot(m) + 1;
        stale = true;
    end
end

%% Step III: back to the original network
% Theta(r) = Theta*(r)/rho(r) by (7) and (11); the queue lengths transfer
% unchanged, the two networks sharing their steady-state distribution.
Xb = zeros(nbatches,R);
Qb = zeros(Mx,R,nbatches);
for b = 1:nbatches
    Xb(b,:) = (xnum(b,:)./den(b))./rhoTot;
    Qb(:,:,b) = qnum(:,:,b)./den(b);
end
X = (sum(xnum,1)./sum(den))./rhoTot;
Qx = sum(qnum,3)./sum(den);
Q = Qx(1:M,:);

% batch-means standard error and the two-sigma interval of the paper
if nbatches > 1
    Xse = std(Xb,0,1)/sqrt(nbatches);
    Qse = std(Qb,0,3)/sqrt(nbatches);
else
    Xse = zeros(1,R);
    Qse = zeros(Mx,R);
end
ci.Xse = Xse;
ci.Xlo = X - 2*Xse;
ci.Xhi = X + 2*Xse;
ci.Qse = Qse(1:M,:);
ci.Qlo = Q - 2*ci.Qse;
ci.Qhi = Q + 2*ci.Qse;
ci.batches = nbatches;
ci.samples = samples;
ci.burnin = nburn;
end
