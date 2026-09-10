function L = smp_passage_lst(P, hlst, pi0, target, s)
% L = SMP_PASSAGE_LST(P, HLST, PI0, TARGET, S)
%
% Laplace-Stieltjes transform of the first passage time into the target state
% set for a semi-Markov chain, at the (possibly complex) points S.
%
% Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
% in Large Markov Chains", 2002, Eqs. 4-5:
%
%     L_i(s) = sum_{k not in B} r*_ik(s) L_k(s) + sum_{k in B} r*_ik(s)
%
% so (I - R*_AA(s)) L_A(s) = R*_AB(s) 1, one linear system per value of s.
%
% HLST is either
%   (nstates x 1) cell of handles h*_i(s)   the sojourn in i depends only on i,
%                                           so r*_ik(s) = P(i,k) h*_i(s) and
%                                           the complex numbers stay on the
%                                           DIAGONAL of the system (Eq. 5)
%   (nstates x nstates) cell of r*_ik(s)    the full Markov-renewal kernel; the
%                                           coefficients are then complex
%                                           functions of s throughout, which is
%                                           the harder case the paper flags
%
% Distribution objects supply their own transform: MARKOVIAN.evalLST gives the
% closed form pie (sI-D0)^{-1} (-D0) e for the phase-type family, so
% hlst{i} = @(s) dist.evalLST(s) is the intended way to build these.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = size(P,1);
target = unique(reshape(target,1,[]));
isTarget = false(1,n);
isTarget(target) = true;
A = find(~isTarget);
nA = numel(A);

if isempty(pi0)
    pi0 = ones(1,n)/n;
end
pi0 = reshape(pi0,1,[]);
if numel(pi0) ~= n
    line_error(mfilename, 'PI0 must be a distribution over the state space, one entry per state.');
end
atom = sum(pi0(target));

perState = iscell(hlst) && (size(hlst,2) == 1 || isvector(hlst));
L = zeros(size(s));
I = eye(nA);
for is = 1:numel(s)
    sv = s(is);
    if perState
        h = zeros(nA,1);
        for a = 1:nA
            h(a) = hlst{A(a)}(sv);
        end
        RAA = diag(h) * P(A,A);
        RAB = diag(h) * sum(P(A,target),2);
    else
        Rs = zeros(n,n);
        for i = 1:n
            for k = 1:n
                if ~isempty(hlst{i,k})
                    Rs(i,k) = hlst{i,k}(sv);
                end
            end
        end
        RAA = Rs(A,A);
        RAB = sum(Rs(A,target),2);
    end
    LA = (I - RAA) \ RAB;
    L(is) = pi0(A) * LA + atom;
end
end
