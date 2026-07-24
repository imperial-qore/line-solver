function d = dmap_dist(DMAPA,DMAPB,L,alA,alB)
% d=dmap_dist(DMAPA,DMAPB,L,alA,alB) - Squared L2 distance between the
% lag-L joint PMFs of two discrete-time MAPs.
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
%  L: number of lags (L=1 for lag-1 joint PMF distance)
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: squared L2 distance

NA = size(DMAPA{1},1); NB = size(DMAPB{1},1);
if nargin<5, alB=dtmc_solve(inv(eye(NB)-DMAPB{1})*DMAPB{2}); end
if nargin<4, alA=dtmc_solve(inv(eye(NA)-DMAPA{1})*DMAPA{2}); end
d = dmap_exp_mul_int(DMAPA,DMAPA,L+1,alA,alA) ...
  - 2*dmap_exp_mul_int(DMAPA,DMAPB,L+1,alA,alB) ...
  + dmap_exp_mul_int(DMAPB,DMAPB,L+1,alB,alB);
end
