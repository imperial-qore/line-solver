function [B1,d] = map_optim_dist(MAPA,alA,B0,alB,L)
% [B1,d]=map_optim_dist(MAPA,alA,B0,alB,L) - Find B1 minimizing the
% lag-L joint density distance, given fixed B0.
%
% Reference:
%   G. Horvath, "Measuring the distance between MAPs and some
%   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
%   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
%
%  Input:
%  MAPA: reference MAP in the form of {D0,D1}
%  alA: stationary vector at arrivals for MAPA
%  B0: D0 matrix of the approximating MAP (fixed)
%  alB: stationary vector at arrivals for the approximating MAP
%  L: number of lags
%
%  Output:
%  B1: optimal D1 matrix
%  d: minimum distance achieved

A0=MAPA{1}; A1=MAPA{2};
NB = size(B0,1);
a = sum(-A0,2);
b = sum(-B0,2);
Aeq = [kron(eye(NB),alB*inv(-B0)) ; kron(ones(1,NB),eye(NB))];
beq = [alB';b];

    function di = myfun(x)
        B1x = reshape(x, NB, NB);
        di = map_dist(MAPA,{B0,B1x},L,alA,alB);
    end

warning off;
if L==1
    Z_AB = lyap(A0', B0, alA'*alB);
    Z_AA = lyap(A0', A0, alA'*alA);
    Z_BB = lyap(B0', B0, alB'*alB);
    X_AB = lyap(A0, B0', a*b');
    X_AA = lyap(A0, A0', a*a');
    X_BB = lyap(B0, B0', b*b');
    vA1 = reshape(A1,numel(A1),1);
    options = optimset('Display','off');
    vB1 = quadprog(kron(X_BB,Z_BB), -vA1'*kron(X_AB,Z_AB), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
    B1 = reshape(vB1, NB, NB);
    d = vB1'*kron(X_BB,Z_BB)*vB1 + vA1'*kron(X_AA,Z_AA)*vA1 - 2*vA1'*kron(X_AB,Z_AB)*vB1;
else
    options = optimset('Display','off');
    [vB1, d] = fmincon(@myfun, rand(NB*NB,1), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
    B1 = reshape(vB1, NB, NB);
end
warning on;
end
