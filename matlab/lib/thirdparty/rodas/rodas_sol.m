function b = rodas_sol(n, a, b, ip)
% RODAS_SOL  SOL: forward/back substitution against the factors of RODAS_DEC
%
%   B = RODAS_SOL(N, A, B, IP)
%
%   PORTED THIRD-PARTY CODE (decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   See also RODAS_DEC, RODAS_SOLB, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if n ~= 1
    nm1 = n - 1;
    for k = 1:nm1
        kp1 = k + 1;
        m = ip(k);
        t = b(m);
        b(m) = b(k);
        b(k) = t;
        b(kp1:n) = b(kp1:n) + a(kp1:n, k)*t;
    end
    for kb = 1:nm1
        km1 = n - kb;
        k = km1 + 1;
        b(k) = b(k)/a(k, k);
        t = -b(k);
        b(1:km1) = b(1:km1) + a(1:km1, k)*t;
    end
end
b(1) = b(1)/a(1, 1);
end
