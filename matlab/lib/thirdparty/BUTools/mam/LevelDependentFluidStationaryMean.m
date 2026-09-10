% res = LevelDependentFluidStationaryMean (masses, iniF, KF, cloF, iniB, KB, cloB, T)
% Calculates the stationary mean fluid level of first order and second order
% level dependent (multi-regime) fluid models, from the matrix-exponential
% building blocks returned by SecondOrderLevelDependentFluidSolve.
%
% * masses: list of point-mass vectors (K+1 vectors located at levels
%   T(1)=0, T(2), ..., T(K+1)=top threshold)
% * iniF, KF, cloF: initial vector, matrix exponent and closing matrix for
%   the forward direction for each regime (a list of size K)
% * iniB, KB, cloB: initial vector, matrix exponent and closing matrix for
%   the backward direction for each regime (a list of size K)
% * T: vector of regime thresholds (length K)
%
% res is the scalar mean fluid level E[X].
%
% Note: the level-dependent stationary density in regime k over the interval
% [T(k),T(k+1)] is
%   pi_k(x) = iniF{k}*expm(KF{k}*(x-T(k)))*cloF{k}
%           + iniB{k}*expm(KB{k}*(T(k+1)-x))*cloB{k}
% and E[X] is obtained in closed form by integrating x*pi_k(x) over each
% regime plus the contribution T(j)*sum(masses{j}) of the point masses.
function res = LevelDependentFluidStationaryMean (masses, iniF, KF, cloF, iniB, KB, cloB, T)

    K = length(T);
    T = [0 T];
    N = length(masses{1});
    h = ones(N,1);

    % integrals of a matrix exponential over [0,L] via nilpotent block
    % augmentation (robust even if M is singular, i.e. has a zero eigenvalue):
    %   J0 = int_0^L expm(M u) du,   J1 = int_0^L u*expm(M u) du
    function [J0, J1] = expIntMoments(M, L)
        n = size(M,1);
        Zn = zeros(n);
        In = eye(n);
        A = [M, In, Zn; Zn, Zn, In; Zn, Zn, Zn];
        W = expm(A*L);
        J0 = W(1:n, n+1:2*n);
        W13 = W(1:n, 2*n+1:3*n);
        J1 = L*J0 - W13;
    end

    res = 0;
    % contribution of the point masses (located at levels T(1..K+1))
    for j=1:K+1
        res = res + T(j) * sum(masses{j});
    end
    % contribution of the continuous density in each regime
    for k=1:K
        Tk = T(k+1) - T(k);
        [J0F, J1F] = expIntMoments(KF{k}, Tk);
        [J0B, J1B] = expIntMoments(KB{k}, Tk);
        res = res + iniF{k} * (T(k)*J0F + J1F)   * cloF{k} * h;
        res = res + iniB{k} * (T(k+1)*J0B - J1B) * cloB{k} * h;
    end
end
