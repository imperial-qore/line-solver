function d = map_dist_lag1(MAPA,MAPB,alA,alB)
% d=map_dist_lag1(MAPA,MAPB,alA,alB) - Lag-1 joint density L2 distance
% via Kronecker/Lyapunov formulation. Equivalent to map_dist(MAPA,MAPB,1).
%
% Reference:
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
%  Input:
%  MAPA: first MAP in the form of {D0,D1}
%  MAPB: second MAP in the form of {D0,D1}
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: squared L2 distance of lag-1 joint densities

A0=MAPA{1}; A1=MAPA{2};
B0=MAPB{1}; B1=MAPB{2};
if nargin<4, alB=map_pie(MAPB); end
if nargin<3, alA=map_pie(MAPA); end

a = sum(-A0,2);
b = sum(-B0,2);

Z_AB = lyap(A0', B0, alA'*alB);
Z_AA = lyap(A0', A0, alA'*alA);
Z_BB = lyap(B0', B0, alB'*alB);

X_AB = lyap(A0, B0', a*b');
X_AA = lyap(A0, A0', a*a');
X_BB = lyap(B0, B0', b*b');

vA1 = reshape(A1, numel(A1), 1);
vB1 = reshape(B1, numel(B1), 1);

d = vB1'*kron(X_BB,Z_BB)*vB1 + vA1'*kron(X_AA,Z_AA)*vA1 ...
  - 2*vA1'*kron(X_AB,Z_AB)*vB1;
end
