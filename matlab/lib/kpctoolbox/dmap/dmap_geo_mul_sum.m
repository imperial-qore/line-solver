function d = dmap_geo_mul_sum(DMAPA,DMAPB,alA,alB)
% d=dmap_geo_mul_sum(DMAPA,DMAPB,alA,alB) - Geometric sum for
% autocorrelation distance computation of discrete-time MAPs.
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
%  d: geometric sum value

D0A=DMAPA{1}; D1A=DMAPA{2};
D0B=DMAPB{1}; D1B=DMAPB{2};
NA = size(D0A,1); NB = size(D0B,1);
if nargin<4, alB=dtmc_solve(inv(eye(NB)-D0B)*D1B); end
if nargin<3, alA=dtmc_solve(inv(eye(NA)-D0A)*D1A); end

IA = eye(NA); IB = eye(NB);
D0Ai = inv(IA-D0A);
D0Bi = inv(IB-D0B);

PAh = D0Ai*D1A - ones(NA,1)*alA;
PBh = D0Bi*D1B - ones(NB,1)*alB;

M = eye(NA*NB) - kron(PBh', PAh);
if rcond(M)<1e-10
    d = 1/rcond(M);
else
    X = dlyap(PAh, PBh, sum(D0Ai,2)*(alB*D0Bi));
    d = sum(alA*D0Ai*X*D0Bi);
end
end
