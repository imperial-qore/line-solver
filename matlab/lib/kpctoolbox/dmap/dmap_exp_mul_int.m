function d = dmap_exp_mul_int(DMAPA,DMAPB,L,alA,alB)
% d=dmap_exp_mul_int(DMAPA,DMAPB,L,alA,alB) - Joint PMF inner product
% of two discrete-time MAPs via recursive discrete Lyapunov equations.
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
%  L: number of inter-arrival times in the joint PMF
%  alA,alB: (optional) stationary vectors at arrival epochs
%
%  Output:
%  d: inner product of joint PMFs

D0A=DMAPA{1}; D1A=DMAPA{2};
D0B=DMAPB{1}; D1B=DMAPB{2};
NA = size(D0A,1); NB = size(D0B,1);
if nargin<5
    alB = dtmc_solve(inv(eye(NB)-D0B)*D1B);
end
if nargin<4
    alA = dtmc_solve(inv(eye(NA)-D0A)*D1A);
end
Z = dlyap(D0B', D0A, alB'*alA);
for i=1:L-1
    Z = dlyap(D0B', D0A, D1B'*Z*D1A);
end
dA = sum(eye(NA)-D0A, 2);
dB = sum(eye(NB)-D0B, 2);
d = dB' * Z * dA;
end
