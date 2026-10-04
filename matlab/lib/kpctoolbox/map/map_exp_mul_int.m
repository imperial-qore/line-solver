function d = map_exp_mul_int(MAPA,MAPB,L,alA,alB)
% d=map_exp_mul_int(MAPA,MAPB,L,alA,alB) - Joint density inner product
% of two MAPs via recursive Sylvester (Lyapunov) equations.
%
% Reference:
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
%  Input:
%  MAPA: first MAP in the form of {D0,D1}
%  MAPB: second MAP in the form of {D0,D1}
%  L: number of inter-arrival times in the joint density
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: inner product of joint densities

A0=MAPA{1}; A1=MAPA{2};
B0=MAPB{1}; B1=MAPB{2};
if nargin<5
    alB=map_pie(MAPB);
end
if nargin<4
    alA=map_pie(MAPA);
end
Z = lyap(B0', A0, alB'*alA);
for i=1:L-1
    Z = lyap(B0', A0, B1'*Z*A1);
end
d = sum(-B0,2)' * Z * sum(-A0,2);
end
