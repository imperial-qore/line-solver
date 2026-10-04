function [e1, ip1, ier] = rodas_decomr(n, fjac, fmas, mlmas, mumas, m1, m2, ...
    nm1, fac1, e1, ip1, ijob, lin)
% RODAS_DECOMR  DECOMR: assemble E1 = fac1*M - J for this IJOB and factor it
%
%   [E1, IP1, IER] = RODAS_DECOMR(N, FJAC, FMAS, MLMAS, MUMAS, M1, M2, NM1, ...
%                                 FAC1, E1, IP1, IJOB, LIN)
%
%   PORTED THIRD-PARTY CODE (dc_decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   IJOB names the storage combination:
%     1  M = I,      J full        2  M = I,      J banded
%     3  M banded,   J full        4  M banded,   J banded
%     5  M full,     J full
%   and 11..15 are the same five with the second-order block (M1 > 0) folded
%   in. 6 is "THIS OPTION IS NOT PROVIDED" upstream and 7..10 belong to
%   RADAU5's Hessenberg option, which RODAS never selects.
%
%   See also RODAS_SLVROD, RODAS_DEC, RODAS_DECB, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ier = 0;
switch ijob
    case 1
        % ---  B=IDENTITY, JACOBIAN A FULL MATRIX
        for j = 1:n
            e1(1:n, j) = -fjac(1:n, j);
            e1(j, j) = e1(j, j) + fac1;
        end
        [e1, ip1, ier] = rodas_dec(n, e1, ip1);
    case 11
        % ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
        for j = 1:nm1
            jm1 = j + m1;
            e1(1:nm1, j) = -fjac(1:nm1, jm1);
            e1(j, j) = e1(j, j) + fac1;
        end
        [e1, ip1, ier] = local_l45(fjac, m1, m2, nm1, fac1, e1, ip1);
    case 2
        % ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
        for j = 1:n
            e1(1+lin.mle:lin.mbjac+lin.mle, j) = -fjac(1:lin.mbjac, j);
            e1(lin.mdiag, j) = e1(lin.mdiag, j) + fac1;
        end
        [e1, ip1, ier] = rodas_decb(n, e1, lin.mle, lin.mue, ip1);
    case 12
        % ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
        for j = 1:nm1
            jm1 = j + m1;
            e1(1+lin.mle:lin.mbjac+lin.mle, j) = -fjac(1:lin.mbjac, jm1);
            e1(lin.mdiag, j) = e1(lin.mdiag, j) + fac1;
        end
        [e1, ip1, ier] = local_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin);
    case {3, 13}
        % ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
        if ijob == 3, nn = n; else, nn = nm1; end
        for j = 1:nn
            if ijob == 3, jm1 = j; else, jm1 = j + m1; end
            e1(1:nn, j) = -fjac(1:nn, jm1);
            for i = max(1, j-mumas):min(nn, j+mlmas)
                e1(i, j) = e1(i, j) + fac1*fmas(i-j+lin.mbdiag, j);
            end
        end
        if ijob == 3
            [e1, ip1, ier] = rodas_dec(n, e1, ip1);
        else
            [e1, ip1, ier] = local_l45(fjac, m1, m2, nm1, fac1, e1, ip1);
        end
    case {4, 14}
        % ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
        if ijob == 4, nn = n; else, nn = nm1; end
        for j = 1:nn
            if ijob == 4, jm1 = j; else, jm1 = j + m1; end
            e1(1+lin.mle:lin.mbjac+lin.mle, j) = -fjac(1:lin.mbjac, jm1);
            for i = 1:lin.mbb
                ib = i + lin.mdiff;
                e1(ib, j) = e1(ib, j) + fac1*fmas(i, j);
            end
        end
        if ijob == 4
            [e1, ip1, ier] = rodas_decb(n, e1, lin.mle, lin.mue, ip1);
        else
            [e1, ip1, ier] = local_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin);
        end
    case {5, 15}
        % ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
        if ijob == 5, nn = n; else, nn = nm1; end
        for j = 1:nn
            if ijob == 5, jm1 = j; else, jm1 = j + m1; end
            e1(1:nn, j) = fmas(1:nn, j)*fac1 - fjac(1:nn, jm1);
        end
        if ijob == 5
            [e1, ip1, ier] = rodas_dec(n, e1, ip1);
        else
            [e1, ip1, ier] = local_l45(fjac, m1, m2, nm1, fac1, e1, ip1);
        end
    otherwise
        % 6 is not provided upstream; 7..10 belong to RADAU5.
end
end

% ---------------------------------------------------------------------------
function [e1, ip1, ier] = local_l45(fjac, m1, m2, nm1, fac1, e1, ip1)
% DECOMR label 45: fold the second-order block into E1, then factor full.
mm = m1/m2;
for j = 1:m2
    for i = 1:nm1
        % A SEQUENTIAL SUM, as the Fortran accumulates it.
        s = 0;
        for k = 0:(mm-1)
            s = (s + fjac(i, j + k*m2))/fac1;
        end
        e1(i, j) = e1(i, j) - s;
    end
end
[e1, ip1, ier] = rodas_dec(nm1, e1, ip1);
end

% ---------------------------------------------------------------------------
function [e1, ip1, ier] = local_l46(fjac, m1, m2, nm1, fac1, e1, ip1, lin)
% DECOMR label 46: the same fold, banded.
mm = m1/m2;
for j = 1:m2
    for i = 1:lin.mbjac
        s = 0;
        for k = 0:(mm-1)
            s = (s + fjac(i, j + k*m2))/fac1;
        end
        e1(i+lin.mle, j) = e1(i+lin.mle, j) - s;
    end
end
[e1, ip1, ier] = rodas_decb(nm1, e1, lin.mle, lin.mue, ip1);
end
