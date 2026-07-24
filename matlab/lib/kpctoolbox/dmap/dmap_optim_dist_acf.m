function [D1B,d] = dmap_optim_dist_acf(DMAPA,alA,D0B,alB)
% [D1B,d]=dmap_optim_dist_acf(DMAPA,alA,D0B,alB) - Find D1B minimizing the
% autocorrelation function distance, given fixed D0B.
%
% Reference (continuous-time formulation):
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
% Discrete-time extension by QORE Lab (https://qore.doc.ic.ac.uk/)
%
%  Input:
%  DMAPA: reference D-MAP in the form of {D0,D1}
%  alA: stationary vector at arrivals for DMAPA
%  D0B: D0 matrix of the approximating D-MAP (fixed)
%  alB: stationary vector at arrivals for the approximating D-MAP
%
%  Output:
%  D1B: optimal D1 matrix
%  d: minimum distance achieved

NB = size(D0B,1);
IB = eye(NB);
b = sum(IB-D0B, 2);
Aeq = [kron(eye(NB),alB*inv(IB-D0B)) ; kron(ones(1,NB),eye(NB))];
beq = [alB';b];

    function di = myfun(x)
        D1Bx = reshape(x, NB, NB);
        di = dmap_dist_acf(DMAPA,{D0B,D1Bx},alA,alB);
    end

warning off;
options = optimset('Display','off');
[vD1B, d] = fmincon(@myfun, rand(NB*NB,1), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
D1B = reshape(vD1B, NB, NB);
warning on;
d = dmap_dist_acf(DMAPA,{D0B,D1B},alA,alB);
end
