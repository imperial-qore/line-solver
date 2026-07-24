function d = map_dist(MAPA,MAPB,L,alA,alB)
% d=map_dist(MAPA,MAPB,L,alA,alB) - Squared L2 distance between the
% lag-L joint densities of two MAPs.
%
% Reference:
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
%  Input:
%  MAPA: first MAP in the form of {D0,D1}
%  MAPB: second MAP in the form of {D0,D1}
%  L: number of lags (L=1 for lag-1 joint density distance)
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: squared L2 distance

if nargin<5, alB=map_pie(MAPB); end
if nargin<4, alA=map_pie(MAPA); end
d = map_exp_mul_int(MAPA,MAPA,L+1,alA,alA) ...
  - 2*map_exp_mul_int(MAPA,MAPB,L+1,alA,alB) ...
  + map_exp_mul_int(MAPB,MAPB,L+1,alB,alB);
end
