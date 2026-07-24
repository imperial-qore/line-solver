function d = map_geo_mul_sum(MAPA,MAPB,alA,alB)
% d=map_geo_mul_sum(MAPA,MAPB,alA,alB) - Geometric sum for
% autocorrelation distance computation.
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
%  d: geometric sum value

A0=MAPA{1}; A1=MAPA{2};
B0=MAPB{1}; B1=MAPB{2};
if nargin<4, alB=map_pie(MAPB); end
if nargin<3, alA=map_pie(MAPA); end

A0i = inv(-A0);
B0i = inv(-B0);
NA = size(A0,1);
NB = size(B0,1);

PAh = A0i*A1 - ones(NA,1)*alA;
PBh = B0i*B1 - ones(NB,1)*alB;

M = eye(NA*NB) - kron(PBh', PAh);
if rcond(M)<1e-10
    d = 1/rcond(M);
else
    X = dlyap(PAh, PBh, sum(A0i,2)*(alB*B0i));
    d = sum(alA*A0i*X*B0i);
end
end
