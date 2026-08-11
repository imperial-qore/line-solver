function d = dmap_dist_lag1(DMAPA,DMAPB,alA,alB)
% d=dmap_dist_lag1(DMAPA,DMAPB,alA,alB) - Lag-1 joint PMF L2 distance
% via Kronecker/Lyapunov formulation. Equivalent to dmap_dist(DMAPA,DMAPB,1).
%
% Reference (continuous-time formulation):
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
% Discrete-time extension by QORE Lab (https://qore.doc.ic.ac.uk/)
%
%  Input:
%  DMAPA: first D-MAP in the form of {D0,D1}
%  DMAPB: second D-MAP in the form of {D0,D1}
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: squared L2 distance of lag-1 joint PMFs

D0A=DMAPA{1}; D1A=DMAPA{2};
D0B=DMAPB{1}; D1B=DMAPB{2};
NA = size(D0A,1); NB = size(D0B,1);
if nargin<4, alB=dtmc_solve(inv(eye(NB)-D0B)*D1B); end
if nargin<3, alA=dtmc_solve(inv(eye(NA)-D0A)*D1A); end

a = sum(eye(NA)-D0A, 2);
b = sum(eye(NB)-D0B, 2);

Z_AB = dlyap(D0A', D0B, alA'*alB);
Z_AA = dlyap(D0A', D0A, alA'*alA);
Z_BB = dlyap(D0B', D0B, alB'*alB);

X_AB = dlyap(D0A, D0B', a*b');
X_AA = dlyap(D0A, D0A', a*a');
X_BB = dlyap(D0B, D0B', b*b');

vD1A = reshape(D1A, numel(D1A), 1);
vD1B = reshape(D1B, numel(D1B), 1);

d = vD1B'*kron(X_BB,Z_BB)*vD1B + vD1A'*kron(X_AA,Z_AA)*vD1A ...
  - 2*vD1A'*kron(X_AB,Z_AB)*vD1B;
end
