function d = dmap_dist_acf(DMAPA,DMAPB,alA,alB)
% d=dmap_dist_acf(DMAPA,DMAPB,alA,alB) - Squared L2 distance between the
% autocorrelation functions of two discrete-time MAPs.
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
%  d: squared L2 distance of autocorrelation functions

D0A=DMAPA{1}; D0B=DMAPB{1};
NA = size(D0A,1); NB = size(D0B,1);
if nargin<4, alB=dtmc_solve(inv(eye(NB)-D0B)*DMAPB{2}); end
if nargin<3, alA=dtmc_solve(inv(eye(NA)-D0A)*DMAPA{2}); end

momA = dmap_moment(DMAPA, [1 2]);
momB = dmap_moment(DMAPB, [1 2]);
mA = momA(1); mB = momB(1);
varA = momA(2)-mA^2;
varB = momB(2)-mB^2;

% The constant removed from each geometric sum is that sum's own k=0 term,
% al*(I-D0)^-2*e, which in DISCRETE time is (m2+m1)/2 and NOT the continuous
% m2/2: (I-D0)^-1 replaces (-D0)^-1, so D0*D0i = D0i-I adds the extra -m1.
cA = (momA(2)+mA)/2;
cB = (momB(2)+mB)/2;

d = (dmap_geo_mul_sum(DMAPA,DMAPA,alA,alA) - cA^2) / varA^2 ...
  - 2*(dmap_geo_mul_sum(DMAPA,DMAPB,alA,alB) - cA*cB) / varA / varB ...
  + (dmap_geo_mul_sum(DMAPB,DMAPB,alB,alB) - cB^2) / varB^2;
end
