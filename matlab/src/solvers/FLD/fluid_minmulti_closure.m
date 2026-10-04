function [h, g, v] = fluid_minmulti_closure(mu, S, c)
% [H, G, V] = FLUID_MINMULTI_CLOSURE(MU, S, C)
%
% Min-normal moment closure of E[min(X_1,...,X_A,c)] for jointly normal X.
%
% FLUID_MIN_CLOSURE closes the TWO-argument min() that a queueing station
% needs, min(n_i,c_i). A Petri-net transition mode needs the many-argument
% one: its enabling degree is
%
%   e(m) = min_a ( m_a / w_a )
%
% over every input arc a, and the rate law then caps that at the mode's server
% count. There is no closed form for the expectation of a min of more than two
% correlated normals, so this function uses the recursion of Clark (1961): the
% running min is replaced at each step by the normal with its exact first two
% moments, and the next argument is folded in with the exact bivariate formulas.
% Each step is FLUID_MIN_CLOSURE's own expression, extended with the second
% moment and the cross-covariances Clark needs to carry the recursion forward:
%
%   th^2  = Var[Z] + Var[X] - 2*Cov[Z,X],  al = (E[Z]-E[X])/th,  p = Phi(-al)
%   E[W]  = E[Z]*p + E[X]*(1-p) - th*phi(al)
%   E[W^2]= (E[Z]^2+Var[Z])*p + (E[X]^2+Var[X])*(1-p) - (E[Z]+E[X])*th*phi(al)
%   Cov[W,Y] = Cov[Z,Y]*p + Cov[X,Y]*(1-p)
%
% with W = min(Z,X). The deterministic cap C is folded in last, by
% FLUID_MIN_CLOSURE itself, so a single-arc mode reduces EXACTLY to the closure
% the queueing methods already use and no second code path exists for it.
%
% THE RECURSION IS ORDER DEPENDENT, as Clark's approximation always is: only the
% first two moments of the running min are kept, so folding the arcs in a
% different order gives a slightly different answer. The order here is the
% caller's, i.e. increasing state coordinate, which is fixed by the layout and
% therefore reproducible.
%
% THE VARIANCES ARE HELD. The derivative G is taken with respect to the MEANS
% only, exactly as FLUID_DRIFT_JACOBIAN differentiates the queueing closure:
% the covariance is a separate unknown of the DAE, pinned by its own consistency
% row, not a function of the mean along the Newton step.
%
% Parameters:
%   mu - (A x 1) means of the arguments, already scaled by the arc weights
%   S  - (A x A) covariance of the arguments, symmetric positive semi-definite
%   c  - deterministic cap (the mode's server count); Inf for no cap
%
% Returns:
%   h - E[min(X_1,...,X_A,c)]
%   g - (A x 1) dH/dMU(a), the probability that arc a is the binding one
%   v - Var[min(X_1,...,X_A)] before the cap, carried for the caller's report
%
% -- Reference
% C. E. Clark, "The greatest of a finite set of random variables", Operations
% Research 9(2):145-162, 1961.
%
% See also FLUID_MIN_CLOSURE, FLUID_PETRI_TERMS, SOLVER_FLUID_PETRI.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(c)
    c = Inf;
end
mu = mu(:);
A = numel(mu);

% A mode with no input arc is enabled at degree one, which is the convention
% the exact engines use (SOLVER_SSA_NRM's SPNENDEGREE, STATE.AFTERGLOBALEVENT).
if A == 0
    h = min(1, c);
    g = zeros(0,1);
    v = 0;
    return
end
if isempty(S)
    S = zeros(A);
end

mz = mu(1);
vz = S(1,1);
covz = S(1,:).';   % Cov[Z, X_a] for every argument, updated as Z absorbs them
p = ones(A,1);     % p(k) = P(Z_{k-1} < X_k), the weight of the running min

for k = 2:A
    th2 = vz + S(k,k) - 2*covz(k);
    if th2 < 0
        th2 = 0; % a covariance beyond the Cauchy-Schwarz bound is not admissible
    end
    th = sqrt(th2);
    if th > 0
        al = (mz - mu(k))/th;
        Phi = 0.5*erfc(-al/sqrt(2));       % normcdf without the Statistics Toolbox
        phi = exp(-0.5*al^2)/sqrt(2*pi);
        pk = 1 - Phi;
        mw = mz*pk + mu(k)*Phi - th*phi;
        e2 = (mz^2 + vz)*pk + (mu(k)^2 + S(k,k))*Phi - (mz + mu(k))*th*phi;
        vw = e2 - mw^2;
        % Clark's moment match can leave a negative variance where the two
        % arguments are nearly identical; the min of two equal normals has the
        % variance of either, which is what the clamp restores.
        if vw < 0
            vw = 0;
        end
        cw = covz*pk + S(:,k)*Phi;
    else
        % Degenerate step: Z and X_k differ by a constant, so the min is the
        % smaller of the two exactly. The band matches FLUID_MIN_CLOSURE, where
        % it exists so that two codebases stopping either side of a kink read
        % the same indicator.
        if mu(k) - mz > GlobalConstants.FineTol*max(1, abs(mz))
            pk = 1;
        else
            pk = 0;
        end
        mw = min(mz, mu(k));
        vw = pk*vz + (1-pk)*S(k,k);
        cw = covz*pk + S(:,k)*(1-pk);
    end
    p(k) = pk;
    mz = mw; vz = vw; covz = cw;
end

v = vz;

% The cap, by the two-argument closure itself: C is deterministic, so its
% variance and its covariance with the running min are both zero.
[h, pcap] = fluid_min_closure(mz, c, vz, 0, 0);

% d/dMU(a): the running min carries weight p at every step it survives, and
% argument a enters with weight 1-p(a) at its own step. Argument 1 IS the
% running min at the start, so it enters with weight one.
g = zeros(A,1);
tail = pcap;
for a = A:-1:2
    g(a) = (1 - p(a))*tail;
    tail = tail*p(a);
end
g(1) = tail;
end
