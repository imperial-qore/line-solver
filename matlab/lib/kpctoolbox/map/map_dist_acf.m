function d = map_dist_acf(MAPA,MAPB,alA,alB)
% d=map_dist_acf(MAPA,MAPB,alA,alB) - Squared L2 distance between the
% autocorrelation functions of two MAPs.
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
%  d: squared L2 distance of autocorrelation functions

A0=MAPA{1}; B0=MAPB{1};
if nargin<4, alB=map_pie(MAPB); end
if nargin<3, alA=map_pie(MAPA); end

momA = map_moment(MAPA, [1 2]);
momB = map_moment(MAPB, [1 2]);
mA = momA(1); mB = momB(1);
varA = momA(2)-mA^2;
varB = momB(2)-mB^2;

d = (map_geo_mul_sum(MAPA,MAPA,alA,alA) - momA(2)^2/4) / varA^2 ...
  - 2*(map_geo_mul_sum(MAPA,MAPB,alA,alB) - momA(2)*momB(2)/4) / varA / varB ...
  + (map_geo_mul_sum(MAPB,MAPB,alB,alB) - momB(2)^2/4) / varB^2;
end
