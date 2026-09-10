function [a, ip, ier] = rodas_dec(n, a, ip)
% RODAS_DEC  DEC: LU by Gaussian elimination with partial pivoting
%
%   [A, IP, IER] = RODAS_DEC(N, A, IP)
%
%   PORTED THIRD-PARTY CODE (decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   Moler's ACM 423 factorisation, NOT MATLAB's backslash. The difference is
%   the point: LU here must produce the same factors in the same order as the
%   Fortran, because RODAS's step-size controller reads the residual of these
%   solves and a different pivot sequence moves the trajectory in the last
%   digits. LAPACK's dgetrf blocks the elimination and would not.
%
%   On return the upper triangle including the diagonal is U and the strict
%   lower triangle holds the NEGATED multipliers, which is what RODAS_SOL
%   expects. IP(K) is the K-th pivot row; IP(N) carries (-1)^(interchanges), or
%   0 when the matrix was found singular.
%
%   IER is 0 if A is nonsingular, else the stage K at which the pivot vanished.
%
%   See also RODAS_SOL, RODAS_DECB, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ier = 0;
ip(n) = 1;
if n ~= 1
    nm1 = n - 1;
    for k = 1:nm1
        kp1 = k + 1;
        % The Fortran scans i = k+1..n replacing on STRICTLY greater, so it
        % keeps the first index attaining the maximum, k included; MAX keeps
        % the first too.
        [~, mrel] = max(abs(a(k:n, k)));
        m = k + mrel - 1;
        ip(k) = m;
        t = a(m, k);
        if m ~= k
            ip(n) = -ip(n);
            a(m, k) = a(k, k);
            a(k, k) = t;
        end
        if t == 0
            ier = k; ip(n) = 0; return
        end
        t = 1/t;
        a(kp1:n, k) = -a(kp1:n, k)*t;
        for j = kp1:n
            t = a(m, j);
            a(m, j) = a(k, j);
            a(k, j) = t;
            if t ~= 0
                a(kp1:n, j) = a(kp1:n, j) + a(kp1:n, k)*t;
            end
        end
    end
end
if a(n, n) == 0
    ier = n; ip(n) = 0;
end
end
