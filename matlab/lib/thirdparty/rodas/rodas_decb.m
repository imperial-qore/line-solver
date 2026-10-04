function [a, ip, ier] = rodas_decb(n, a, ml, mu, ip)
% RODAS_DECB  DECB: the banded counterpart of RODAS_DEC
%
%   [A, IP, IER] = RODAS_DECB(N, A, ML, MU, IP)
%
%   PORTED THIRD-PARTY CODE (decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   The matrix arrives in LINPACK band storage: its diagonals occupy rows
%   ML+1 through 2*ML+MU+1 of A, columns as columns.
%
%   See also RODAS_SOLB, RODAS_DEC, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ier = 0;
ip(n) = 1;
md = ml + mu + 1;
md1 = md + 1;
ju = 0;
if ml ~= 0 && n ~= 1
    if n >= mu + 2
        for j = (mu+2):n
            a(1:ml, j) = 0;
        end
    end
    nm1 = n - 1;
    for k = 1:nm1
        kp1 = k + 1;
        m = md;
        mdl = min(ml, n-k) + md;
        for i = md1:mdl
            if abs(a(i, k)) > abs(a(m, k)), m = i; end
        end
        ip(k) = m + k - md;
        t = a(m, k);
        if m ~= md
            ip(n) = -ip(n);
            a(m, k) = a(md, k);
            a(md, k) = t;
        end
        if t == 0
            ier = k; ip(n) = 0; return
        end
        t = 1/t;
        a(md1:mdl, k) = -a(md1:mdl, k)*t;
        ju = min(max(ju, mu + ip(k)), n);
        mm = md;
        if ju >= kp1
            for j = kp1:ju
                m = m - 1;
                mm = mm - 1;
                t = a(m, j);
                if m ~= mm
                    a(m, j) = a(mm, j);
                    a(mm, j) = t;
                end
                if t ~= 0
                    jk = j - k;
                    for i = md1:mdl
                        a(i-jk, j) = a(i-jk, j) + a(i, k)*t;
                    end
                end
            end
        end
    end
end
if a(md, n) == 0
    ier = n; ip(n) = 0;
end
end
