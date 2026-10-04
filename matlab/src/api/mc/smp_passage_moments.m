function [Mall, m] = smp_passage_moments(P, hmom, pi0, target, nmax)
% [MALL, M] = SMP_PASSAGE_MOMENTS(P, HMOM, PI0, TARGET, NMAX)
%
% Moments of order 1..NMAX of the first passage time into the target state set
% for a semi-Markov chain with embedded transition matrix P and holding-time
% moments HMOM. MALL is (nstates x NMAX), row i for a passage started in state
% i, zero on the target. M is the PI0-weighted moment vector.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002, Sec. 3.2.
%
% HMOM selects which of the paper's two recursions runs, and they are NOT the
% same computation:
%
%   (nstates x NMAX) matrix    m_i(r), the holding time in i depends only on i.
%                              Eq. 7 with the u_i(r) recurrence of Eq. 8,
%                                  u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j),
%                                  u_i(0) = 1,
%                              which are the derivatives at the origin of
%                              1/h*_i(s). Cheaper: no per-pair moments.
%
%   (nstates x nstates) cell   HMOM{i,k} = [m_ik(1) ... m_ik(NMAX)], the r-th
%                              moment of the holding time in i WHEN THE NEXT
%                              STATE IS k. Eq. 6, the full Markov-renewal
%                              kernel. m_ik(0) = P(i,k) is implied and must not
%                              be supplied.
%
% Both solve one linear system of the size of the non-target block per order,
% with the lower-order terms already known -- the iteration the paper
% describes. Unlike the Markov case, the n-th moment needs every moment from
% 1 to n, so NMAX cannot be raised for free.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(nmax)
    nmax = 1;
end
n = size(P,1);
if size(P,2) ~= n
    line_error(mfilename, 'The embedded transition matrix must be square.');
end
if max(abs(sum(P,2) - 1)) > 1e-8
    line_error(mfilename, 'The embedded transition matrix rows must sum to one.');
end
target = unique(reshape(target,1,[]));
if isempty(target)
    line_error(mfilename, 'The target state set is empty: a first passage time into no state is undefined.');
end

isTarget = false(1,n);
isTarget(target) = true;
A = find(~isTarget);
nA = numel(A);
Mall = zeros(n, nmax);
if nA == 0
    m = zeros(1,nmax);
    return
end

PAA = P(A,A);
I = eye(nA);
Msub = zeros(nA, nmax);

if iscell(hmom)
    % --- Eq. 6, the full kernel -------------------------------------------
    if ~isequal(size(hmom), [n n])
        line_error(mfilename, 'A cell HMOM must be (nstates x nstates), entry (i,k) holding the moments of the sojourn in i when the next state is k.');
    end
    for q = 1:nmax
        b = zeros(nA,1);
        for a = 1:nA
            i = A(a);
            acc = 0;
            for c = 1:nA
                k = A(c);
                mik = local_kernel_moments(hmom, i, k, nmax);
                for r = 1:q
                    if r < q
                        acc = acc + nchoosek(q,r) * mik(r) * Msub(c, q-r);
                    else
                        acc = acc + nchoosek(q,r) * mik(r);   % M_k(0) = 1
                    end
                end
            end
            for k = target
                mik = local_kernel_moments(hmom, i, k, nmax);
                acc = acc + mik(q);
            end
            b(a) = acc;
        end
        Msub(:,q) = (I - PAA) \ b;
    end
else
    % --- Eq. 7 with the Eq. 8 recurrence ----------------------------------
    if size(hmom,1) ~= n
        line_error(mfilename, 'A matrix HMOM must carry one row per state.');
    end
    if size(hmom,2) < nmax
        line_error(mfilename, 'HMOM must carry at least NMAX holding-time moments per state.');
    end
    u = local_u(hmom(A,:), nmax);
    for q = 1:nmax
        b = zeros(nA,1);
        for r = 1:q
            if r < q
                b = b - nchoosek(q,r) * u(:,r) .* Msub(:,q-r);
            else
                b = b - nchoosek(q,r) * u(:,r);               % M_i(0) = 1
            end
        end
        Msub(:,q) = (I - PAA) \ b;
    end
end

Mall(A,:) = Msub;
if nargin < 3 || isempty(pi0)
    m = [];
else
    pi0 = reshape(pi0,1,[]);
    if numel(pi0) ~= n
        line_error(mfilename, 'PI0 must be a distribution over the state space, one entry per state.');
    end
    m = pi0 * Mall;
end
end

function u = local_u(mrows, nmax)
% Eq. 8: u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j), u_i(0) = 1. These are
% the derivatives at the origin of 1/h*_i(s), obtained from the moments of
% h*_i alone by differentiating h*(s) y(s) = 1 r times.
nA = size(mrows,1);
u = zeros(nA, nmax);
u0 = ones(nA,1);
for r = 1:nmax
    acc = zeros(nA,1);
    for j = 1:r
        if r-j == 0
            acc = acc + nchoosek(r,j) * mrows(:,j) .* u0;
        else
            acc = acc + nchoosek(r,j) * mrows(:,j) .* u(:,r-j);
        end
    end
    u(:,r) = -acc;
end
end

function mik = local_kernel_moments(hmom, i, k, nmax)
mik = hmom{i,k};
if isempty(mik)
    mik = zeros(1,nmax);
else
    mik = reshape(mik,1,[]);
end
end
