function b = rodas_solb(n, a, ml, mu, b, ip)
% RODAS_SOLB  SOLB: substitution against the band factors of RODAS_DECB
%
%   B = RODAS_SOLB(N, A, ML, MU, B, IP)
%
%   PORTED THIRD-PARTY CODE (decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   See also RODAS_DECB, RODAS_SOL, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

md = ml + mu + 1;
md1 = md + 1;
mdm = md - 1;
nm1 = n - 1;
if ml ~= 0
    if n == 1
        b(1) = b(1)/a(md, 1);
        return
    end
    for k = 1:nm1
        m = ip(k);
        t = b(m);
        b(m) = b(k);
        b(k) = t;
        mdl = min(ml, n-k) + md;
        b(md1+k-md:mdl+k-md) = b(md1+k-md:mdl+k-md) + a(md1:mdl, k)*t;
    end
end
for kb = 1:nm1
    k = n + 1 - kb;
    b(k) = b(k)/a(md, k);
    t = -b(k);
    kmd = md - k;
    lm = max(1, kmd+1);
    if lm <= mdm
        b(lm-kmd:mdm-kmd) = b(lm-kmd:mdm-kmd) + a(lm:mdm, k)*t;
    end
end
b(1) = b(1)/a(md, 1);
end
