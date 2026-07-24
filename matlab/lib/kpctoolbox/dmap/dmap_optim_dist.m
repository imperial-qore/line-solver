function [D1B,d] = dmap_optim_dist(DMAPA,alA,D0B,alB,L)
% [D1B,d]=dmap_optim_dist(DMAPA,alA,D0B,alB,L) - Find D1B minimizing the
% lag-L joint PMF distance, given fixed D0B.
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
%  L: number of lags
%
%  Output:
%  D1B: optimal D1 matrix
%  d: minimum distance achieved

D0A=DMAPA{1}; D1A=DMAPA{2};
NA = size(D0A,1); NB = size(D0B,1);
IA = eye(NA); IB = eye(NB);
a = sum(IA-D0A, 2);
b = sum(IB-D0B, 2);
Aeq = [kron(eye(NB),alB*inv(IB-D0B)) ; kron(ones(1,NB),eye(NB))];
beq = [alB';b];

    function di = myfun(x)
        D1Bx = reshape(x, NB, NB);
        di = dmap_dist(DMAPA,{D0B,D1Bx},L,alA,alB);
    end

warning off;
if L==1
    Z_AB = dlyap(D0A', D0B, alA'*alB);
    Z_AA = dlyap(D0A', D0A, alA'*alA);
    Z_BB = dlyap(D0B', D0B, alB'*alB);
    X_AB = dlyap(D0A, D0B', a*b');
    X_AA = dlyap(D0A, D0A', a*a');
    X_BB = dlyap(D0B, D0B', b*b');
    vD1A = reshape(D1A,numel(D1A),1);
    options = optimset('Display','off');
    vD1B = quadprog(kron(X_BB,Z_BB), -vD1A'*kron(X_AB,Z_AB), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
    D1B = reshape(vD1B, NB, NB);
    d = vD1B'*kron(X_BB,Z_BB)*vD1B + vD1A'*kron(X_AA,Z_AA)*vD1A - 2*vD1A'*kron(X_AB,Z_AB)*vD1B;
else
    options = optimset('Display','off');
    [vD1B, d] = fmincon(@myfun, rand(NB*NB,1), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
    D1B = reshape(vD1B, NB, NB);
end
warning on;
end
