function ak = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
    m2, nm1, fac1, e, ip, dy, ak, fx, ynew, hd, ijob, stage1, lin)
% RODAS_SLVROD  SLVROD: one Rosenbrock stage against the factored E
%
%   AK = RODAS_SLVROD(N, FJAC, MLJAC, MUJAC, FMAS, MLMAS, MUMAS, M1, M2, ...
%                     NM1, FAC1, E, IP, DY, AK, FX, YNEW, HD, IJOB, STAGE1, LIN)
%
%   PORTED THIRD-PARTY CODE (dc_decsol.f, Hairer & Wanner) -- see RODAS_CORE.
%
%   AK receives the solve. DY is the right hand side at the stage, FX the
%   derivative with respect to X (used only when the system is non-autonomous,
%   i.e. HD ~= 0), and YNEW the accumulated stage combination that STAGE1 adds
%   through the mass matrix. IJOB names the same storage combination as in
%   RODAS_DECOMR.
%
%   See also RODAS_DECOMR, RODAS_SOL, RODAS_SOLB, RODAS_CORE.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if hd == 0
    ak(1:n) = dy(1:n);
else
    ak(1:n) = dy(1:n) + hd*fx(1:n);
end

switch ijob
    case 1
        % ---  B=IDENTITY, JACOBIAN A FULL MATRIX
        if stage1, ak(1:n) = ak(1:n) + ynew(1:n); end
        ak = rodas_sol(n, e, ak, ip);
    case 11
        % ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
        if stage1, ak(1:n) = ak(1:n) + ynew(1:n); end
        ak = local_l48(fjac, m1, m2, nm1, fac1, e, ip, ak);
    case 2
        % ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
        if stage1, ak(1:n) = ak(1:n) + ynew(1:n); end
        ak = rodas_solb(n, e, lin.mle, lin.mue, ak, ip);
    case 12
        % ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
        if stage1, ak(1:n) = ak(1:n) + ynew(1:n); end
        ak = local_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin);
    case 3
        % ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
        if stage1
            for i = 1:n
                % A SEQUENTIAL SUM, as the Fortran accumulates it.
                s = 0;
                for j = max(1, i-mlmas):min(n, i+mumas)
                    s = s + fmas(i-j+lin.mbdiag, j)*ynew(j);
                end
                ak(i) = ak(i) + s;
            end
        end
        ak = rodas_sol(n, e, ak, ip);
    case {13, 14}
        % ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
        if stage1
            ak(1:m1) = ak(1:m1) + ynew(1:m1);
            for i = 1:nm1
                s = 0;
                for j = max(1, i-mlmas):min(nm1, i+mumas)
                    s = s + fmas(i-j+lin.mbdiag, j)*ynew(j+m1);
                end
                ak(i+m1) = ak(i+m1) + s;
            end
        end
        if ijob == 14
            ak = local_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin);
        else
            ak = local_l48(fjac, m1, m2, nm1, fac1, e, ip, ak);
        end
    case 4
        % ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
        if stage1
            for i = 1:n
                s = 0;
                for j = max(1, i-mlmas):min(n, i+mumas)
                    s = s + fmas(i-j+lin.mbdiag, j)*ynew(j);
                end
                ak(i) = ak(i) + s;
            end
        end
        ak = rodas_solb(n, e, lin.mle, lin.mue, ak, ip);
    case 5
        % ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
        if stage1
            for i = 1:n
                s = 0;
                for j = 1:n
                    s = s + fmas(i, j)*ynew(j);
                end
                ak(i) = ak(i) + s;
            end
        end
        ak = rodas_sol(n, e, ak, ip);
    case 15
        % ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
        if stage1
            ak(1:m1) = ak(1:m1) + ynew(1:m1);
            for i = 1:nm1
                s = 0;
                for j = 1:nm1
                    s = s + fmas(i, j)*ynew(j+m1);
                end
                ak(i+m1) = ak(i+m1) + s;
            end
        end
        ak = local_l48(fjac, m1, m2, nm1, fac1, e, ip, ak);
    case 6
        % ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
        % ---  THIS OPTION IS NOT PROVIDED. It solves only under STAGE1
        %      upstream; kept identical rather than tidied.
        if stage1
            for i = 1:n
                s = 0;
                for j = 1:n
                    s = s + fmas(i, j)*ynew(j);
                end
                ak(i) = ak(i) + s;
            end
            ak = rodas_solb(n, e, lin.mle, lin.mue, ak, ip);
        end
    otherwise
        % 7..10 belong to RADAU5.
end
end

% ---------------------------------------------------------------------------
function ak = local_l48(fjac, m1, m2, nm1, fac1, e, ip, ak)
% SLVROD label 48: the second-order elimination, full Jacobian.
mm = m1/m2;
for j = 1:m2
    s = 0;
    for k = (mm-1):-1:0
        jkm = j + k*m2;
        s = (ak(jkm) + s)/fac1;
        ak(1+m1:nm1+m1) = ak(1+m1:nm1+m1) + fjac(1:nm1, jkm)*s;
    end
end
sub = rodas_sol(nm1, e, ak(m1+1:m1+nm1), ip);
ak(m1+1:m1+nm1) = sub;
for i = m1:-1:1
    ak(i) = (ak(i) + ak(m2+i))/fac1;
end
end

% ---------------------------------------------------------------------------
function ak = local_l45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin)
% SLVROD label 45: the second-order elimination, banded Jacobian.
mm = m1/m2;
for j = 1:m2
    s = 0;
    for k = (mm-1):-1:0
        jkm = j + k*m2;
        s = (ak(jkm) + s)/fac1;
        for i = max(1, j-mujac):min(nm1, j+mljac)
            ak(i+m1) = ak(i+m1) + fjac(i+mujac+1-j, jkm)*s;
        end
    end
end
sub = rodas_solb(nm1, e, lin.mle, lin.mue, ak(m1+1:m1+nm1), ip);
ak(m1+1:m1+nm1) = sub;
for i = m1:-1:1
    ak(i) = (ak(i) + ak(m2+i))/fac1;
end
end
