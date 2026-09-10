function [B1,d] = map_optim_dist_acf(MAPA,alA,B0,alB)
% [B1,d]=map_optim_dist_acf(MAPA,alA,B0,alB) - Find B1 minimizing the
% autocorrelation function distance, given fixed B0.
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
%
%  Output:
%  B1: optimal D1 matrix
%  d: minimum distance achieved

NB = size(B0,1);
b = sum(-B0,2);
Aeq = [kron(eye(NB),alB*inv(-B0)) ; kron(ones(1,NB),eye(NB))];
beq = [alB';b];

    function di = myfun(x)
        B1x = reshape(x, NB, NB);
        di = map_dist_acf(MAPA,{B0,B1x},alA,alB);
    end

warning off;
options = optimset('Display','off');
[vB1, d] = fmincon(@myfun, rand(NB*NB,1), [], [], Aeq, beq, 1e-6*ones(NB*NB,1), [], [], options);
B1 = reshape(vB1, NB, NB);
warning on;
d = map_dist_acf(MAPA,{B0,B1},alA,alB);
end
