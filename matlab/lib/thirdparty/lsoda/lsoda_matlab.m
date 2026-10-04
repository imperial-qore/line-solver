function [T, Y, stats] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords, opts)
% LSODA_MATLAB  Pure MATLAB LSODA integrator (no MEX, no compiler needed)
%
%   [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol)
%   [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep)
%   [T, Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol, mxstep, mxordn, mxords)
%   [T, Y] = lsoda_matlab(..., mxordn, mxords, opts)
%   [T, Y, stats] = lsoda_matlab(...)
%
%   Solves dy/dt = f(t,y) with automatic switching between the Adams
%   (nonstiff) and BDF (stiff) methods. The argument vector, the two output
%   modes and the failure behaviour replicate lsoda_mex, so this function is
%   a drop-in replacement for the MEX gateway on platforms without a
%   compiled binary.
%
%   Inputs:
%     odefun - function handle, dydt = odefun(t, y), y and dydt columns
%     tspan  - [t0 tf] (adaptive output, all internal steps returned) or
%              [t0 t1 ... tf] with more than 2 entries (output at those times)
%     y0     - initial condition vector
%     rtol   - relative tolerance, scalar or one per equation (default 1e-6)
%     atol   - absolute tolerance, scalar or one per equation (default 1e-9)
%     mxstep - max internal steps per output interval (default 500)
%     mxordn - max order of the Adams method, 1-12 (0 = library default 12)
%     mxords - max order of the BDF method, 1-5 (0 = library default 5)
%     opts   - (optional) struct with fields:
%              .hmax       - max step size (0 = unbounded), odeset MaxStep
%              .hmin       - min step size (default 0)
%              .h0         - initial step size (0 = LSODA's own heuristic)
%              .forceStiff - start on BDF and stay there, never switching to
%                            Adams. The Adams half loses its stability bound at
%                            a fixed point, where the corrector converges on the
%                            roundoff branch before `pdest` is ever formed, so
%                            `h` runs up to hmax and the trajectory wanders
%                            around the fixed point instead of settling on it.
%                            This is the same pin as LSODA.setForceStiff(true)
%                            in the JAR; see _kb/06-solver-catalog.md
%
%   Outputs:
%     T      - column vector of output times
%     Y      - solution, one row per output time
%     stats  - struct with fields nst, nfe, nje, nqu, meth
%
%   Port of liblsoda (C, Yu Feng and Hon Wah Tam) from the original LSODA of
%   Linda R. Petzold and Alan C. Hindmarsh (LLNL). The internal routines,
%   the coefficient tables and the step/order/method selection are a
%   line-for-line transcription of matlab/lib/thirdparty/lsoda/src.
%
%   The MIT License
%   Copyright (c) 2011 McWilliam Cosmology Center, Carnegie Mellon University.
%
%   See also: lsoda_solve, lsoda_fast, lsoda_accurate

%% ---- arguments ----
if nargin < 3
    error('lsoda:nargin', 'Usage: [T,Y] = lsoda_matlab(odefun, tspan, y0, rtol, atol [, mxstep [, mxordn, mxords]])');
end
if nargin < 4 || isempty(rtol), rtol = 1e-6; end
if nargin < 5 || isempty(atol), atol = 1e-9; end
if nargin < 6 || isempty(mxstep), mxstep = 500; end
if nargin < 7 || isempty(mxordn), mxordn = 0; end
if nargin < 8 || isempty(mxords), mxords = 0; end
if nargin < 9 || isempty(opts), opts = struct(); end
if ~isa(odefun, 'function_handle')
    error('lsoda:odefun', 'First argument must be a function handle.');
end

tspan = tspan(:).';
if numel(tspan) < 2
    error('lsoda:tspan', 'tspan must have at least 2 elements.');
end
y0 = y0(:).';
neq = numel(y0);
if neq < 1
    error('lsoda:neq', 'neq = %d is less than 1.', neq);
end

rtolv = expandTol(rtol, neq, 'rtol');
atolv = expandTol(atol, neq, 'atol');
if any(rtolv < 0)
    error('lsoda:rtol', 'rtol = %g is less than 0.', min(rtolv));
end
if any(atolv < 0)
    error('lsoda:atol', 'atol = %g is less than 0.', min(atolv));
end

mxstep = double(mxstep);
if mxstep <= 0, mxstep = 500; end
mxordn = double(mxordn);
if mxordn <= 0, mxordn = 100; end
mxordn = min(mxordn, 12);
mxords = double(mxords);
if mxords <= 0, mxords = 100; end
mxords = min(mxords, 5);

mxhnil = 10;
h0opt = optField(opts, 'h0', 0.0);
hmin = optField(opts, 'hmin', 0.0);
hmax = optField(opts, 'hmax', 0.0);
if hmax < 0.0 || hmin < 0.0
    error('lsoda:hbounds', 'hmax and hmin must be nonnegative.');
end
hmxi = 0.0;
if hmax > 0.0
    hmxi = 1.0 / hmax;
end
forceStiff = logical(optField(opts, 'forceStiff', false));
tcrit = 0.0;

%% ---- constants ----
ETA = 2.2204460492503131e-16;
SQRTETA = 1.4901161193847656e-08;
CCMAX = 0.3;
MAXCOR = 3;
MSBP = 20;
MXNCF = 10;
RATIO = 5.0;
% stability limits of the Adams method, C index 0..12, hence SM1(k+1)
SM1 = [0.0, 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025];
% method-switch constants, C index 1..12, leading dummy dropped
CM1 = [2.0, 5.999999999999999, 4.0, 1.5777777777777777, 0.44444444444444453, ...
       0.09721974866042459, 0.01755731922398589, 0.002677042079828849, ...
       0.00035127324808082747, 4.013484869498161e-05, 4.011251133543618e-06, ...
       3.541120754086548e-07];
CM2 = [2.0, 1.5, 0.6666666666666667, 0.20833333333333337, 0.04999999999999999, ...
       0.09721974866042459, 0.01755731922398589, 0.002677042079828849, ...
       0.00035127324808082747, 4.013484869498161e-05, 4.011251133543618e-06, ...
       3.541120754086548e-07];

%% ---- solver state, shared with the nested functions ----
lenyh = 1 + max(mxordn, mxords);
YH = zeros(lenyh, neq);
WM = zeros(neq, neq);
EWT = zeros(1, neq);
SAVF = zeros(1, neq);
ACOR = zeros(1, neq);
IPVT = zeros(1, neq);
EL = zeros(1, 14);
ELCO = zeros(12, 14);
TESCO = zeros(12, 3);

h = 0.0; hu = 0.0; rc = 0.0; tn = 0.0; tsw = 0.0; pdnorm = 0.0;
crate = 0.0; holdh = 0.0; rmax = 0.0; pdest = 0.0; pdlast = 0.0;
ialth = 0; ipup = 0; nslp = 0; icount = 0; irflag = 0;
imxer = 0; nhnil = 0; nslast = 0;
jcur = 0; meth = 0; mused = 0; nq = 0; nst = 0;
ncf = 0; nfe = 0; nje = 0; nqu = 0; miter = 0;
jstart = 0;
state = 1;
errmsg = '';
y = y0;
tcur = tspan(1);

%% ---- driver, mirroring the two output modes of lsoda_mex ----
ntout = numel(tspan);
if ntout == 2
    tf = tspan(2);
    T = tcur;
    Y = y;
    tdir = 1.0;
    if tf <= tspan(1), tdir = -1.0; end
    while (tf - tcur) * tdir > 0.0
        lsDrive(tf, 2);
        if state <= 0
            raiseFailure();
        end
        if (tf - tcur) * tdir >= 0.0
            T(end+1, 1) = tcur; %#ok<AGROW>
            Y(end+1, :) = y;    %#ok<AGROW>
        end
    end
    if T(end) ~= tf
        lsDrive(tf, 1);
        if state <= 0
            raiseFailure();
        end
        T(end+1, 1) = tcur;
        Y(end+1, :) = y;
    end
else
    T = zeros(ntout, 1);
    Y = zeros(ntout, neq);
    T(1) = tcur;
    Y(1, :) = y;
    for iout = 2:ntout
        lsDrive(tspan(iout), 1);
        if state <= 0
            raiseFailure();
        end
        T(iout) = tcur;
        Y(iout, :) = y;
    end
end

stats = struct('nst', nst, 'nfe', nfe, 'nje', nje, 'nqu', nqu, 'meth', meth);

%% ================= nested functions =================

    function raiseFailure()
        error('lsoda:failed', 'LSODA failed at t=%g (state=%d): %s', tcur, state, errmsg);
    end

    function ydot = lsF(t, yv)
        % right-hand side, user code sees columns
        dy = odefun(t, yv(:));
        ydot = dy(:).';
    end

    function vm = lsVmnorm(v)
        vm = max(abs(v) .* EWT);
    end

    function an = lsFnorm()
        an = max((sum(abs(WM) ./ EWT, 2)).' .* EWT);
    end

    function lsEwset(ycur)
        EWT = rtolv .* abs(ycur) + atolv;
        EWT = 1.0 ./ EWT;
    end

    function info = lsDgefa()
        % LINPACK dgefa, kji variant with row indexing, as in dgefa.c
        info = 0;
        for k = 1:(neq - 1)
            [~, jrel] = max(abs(WM(k, k:neq)));
            j = jrel + k - 1;
            IPVT(k) = j;
            if WM(k, j) == 0.0
                info = k;
                continue
            end
            if j ~= k
                tswap = WM(k, j);
                WM(k, j) = WM(k, k);
                WM(k, k) = tswap;
            end
            t = -1.0 / WM(k, k);
            WM(k, (k+1):neq) = t * WM(k, (k+1):neq);
            for i = (k+1):neq
                t = WM(i, j);
                if j ~= k
                    WM(i, j) = WM(i, k);
                    WM(i, k) = t;
                end
                WM(i, (k+1):neq) = WM(i, (k+1):neq) + t * WM(k, (k+1):neq);
            end
        end
        IPVT(neq) = neq;
        if WM(neq, neq) == 0.0
            info = neq;
        end
    end

    function b = lsDgesl(b)
        % solve a*x=b with the factors of lsDgefa, job = 0
        for k = 1:neq
            t = 0.0;
            if k > 1
                t = WM(k, 1:(k-1)) * b(1:(k-1)).';
            end
            b(k) = (b(k) - t) / WM(k, k);
        end
        for k = (neq - 1):-1:1
            b(k) = b(k) + WM(k, (k+1):neq) * b((k+1):neq).';
            j = IPVT(k);
            if j ~= k
                tswap = b(j);
                b(j) = b(k);
                b(k) = tswap;
            end
        end
    end

    function lsCfode(methval)
        if methval == 1
            ELCO(1, 1) = 1.0;
            ELCO(1, 2) = 1.0;
            TESCO(1, 1) = 0.0;
            TESCO(1, 2) = 2.0;
            TESCO(2, 1) = 1.0;
            TESCO(12, 3) = 0.0;
            pc = zeros(1, 14);
            pc(1) = 1.0;
            rqfac = 1.0;
            for nqi = 2:12
                rq1fac = rqfac;
                rqfac = rqfac / nqi;
                nqm1 = nqi - 1;
                fnqm1 = nqm1;
                nqp1 = nqi + 1;
                pc(nqi) = 0.0;
                for i = nqi:-1:2
                    pc(i) = pc(i-1) + fnqm1 * pc(i);
                end
                pc(1) = fnqm1 * pc(1);
                pint = pc(1);
                xpin = pc(1) / 2.0;
                tsign = 1.0;
                for i = 2:nqi
                    tsign = -tsign;
                    pint = pint + tsign * pc(i) / i;
                    xpin = xpin + tsign * pc(i) / (i + 1);
                end
                ELCO(nqi, 1) = pint * rq1fac;
                ELCO(nqi, 2) = 1.0;
                for i = 2:nqi
                    ELCO(nqi, i+1) = rq1fac * pc(i) / i;
                end
                agamq = rqfac * xpin;
                ragq = 1.0 / agamq;
                TESCO(nqi, 2) = ragq;
                if nqi < 12
                    TESCO(nqp1, 1) = ragq * rqfac / nqp1;
                end
                TESCO(nqm1, 3) = ragq;
            end
        else
            pc = zeros(1, 14);
            pc(1) = 1.0;
            rq1fac = 1.0;
            for nqi = 1:5
                fnq = nqi;
                nqp1 = nqi + 1;
                pc(nqp1) = 0.0;
                for i = (nqi+1):-1:2
                    pc(i) = pc(i-1) + fnq * pc(i);
                end
                pc(1) = pc(1) * fnq;
                for i = 1:nqp1
                    ELCO(nqi, i) = pc(i) / pc(2);
                end
                ELCO(nqi, 2) = 1.0;
                TESCO(nqi, 1) = rq1fac;
                TESCO(nqi, 2) = nqp1 / ELCO(nqi, 1);
                TESCO(nqi, 3) = (nqi + 2) / ELCO(nqi, 1);
                rq1fac = rq1fac / fnq;
            end
        end
    end

    function lsResetcoeff()
        el0 = EL(1);
        EL(1:(nq+1)) = ELCO(nq, 1:(nq+1));
        rc = rc * EL(1) / el0;
    end

    function lsScaleh(rh)
        rh = min(rh, rmax);
        rh = rh / max(1.0, abs(h) * hmxi * rh);
        if meth == 1
            irflag = 0;
            pdh = max(abs(h) * pdlast, 0.000001);
            if rh * pdh * 1.00001 >= SM1(nq+1)
                rh = SM1(nq+1) / pdh;
                irflag = 1;
            end
        end
        r = 1.0;
        for j = 2:(nq+1)
            r = r * rh;
            YH(j, :) = YH(j, :) * r;
        end
        h = h * rh;
        rc = rc * rh;
        ialth = nq + 1;
    end

    function [iflag, dky] = lsIntdy(t, k)
        dky = zeros(1, neq);
        iflag = 0;
        if k < 0 || k > nq
            iflag = -1;
            return
        end
        tp = tn - hu - 100.0 * ETA * (tn + hu);
        if (t - tp) * (t - tn) > 0.0
            iflag = -2;
            return
        end
        s = (t - tn) / h;
        ic = 1;
        for jj = (nq + 1 - k):nq
            ic = ic * jj;
        end
        co = ic;
        dky = co * YH(nq+1, :);
        for j = (nq-1):-1:k
            jp1 = j + 1;
            ic = 1;
            for jj = (jp1 - k):j
                ic = ic * jj;
            end
            co = ic;
            dky = co * YH(jp1, :) + s * dky;
        end
        if k ~= 0
            dky = dky * h^(-k);
        end
    end

    function ok = lsPrja()
        % P = I - h*el[1]*J by finite differences, then LU.
        %
        % UNDER forceStiff THE INCREMENT CARRIES NUMJAC'S FLOOR. LSODA sizes the
        % difference as max(sqrt(eps)*|y_j|, r0/ewt_j), and a component sitting
        % at EXACTLY zero under a tight atol takes the second branch with
        % r ~ 1e-19: that column is rounding noise divided by the increment and
        % the corrector converges to nonsense. The auto-switcher never meets it,
        % starting on Adams and taking a Jacobian only once the solution is
        % smooth; the pin takes one at t0, so it needs ode15s's own rule,
        % sqrt(eps)*max(|y_j|, atol_j/rtol_j). Applied ONLY under the pin, so
        % the C reference vectors the auto-switcher reproduces do not move.
        nje = nje + 1;
        hl0 = h * EL(1);
        if miter ~= 2
            ok = false;
            return
        end
        fac = lsVmnorm(SAVF);
        r0 = 1000.0 * abs(h) * ETA * neq * fac;
        if r0 == 0.0
            r0 = 1.0;
        end
        for j = 1:neq
            yj = y(j);
            r = max(SQRTETA * abs(yj), r0 / EWT(j));
            if forceStiff && rtolv(j) > 0.0
                r = max(r, SQRTETA * max(abs(yj), atolv(j) / rtolv(j)));
            end
            y(j) = y(j) + r;
            facv = -hl0 / r;
            ACOR = lsF(tn, y);
            WM(:, j) = ((ACOR - SAVF) * facv).';
            y(j) = yj;
        end
        nfe = nfe + neq;
        pdnorm = lsFnorm() / abs(hl0);
        WM(1:(neq+1):end) = WM(1:(neq+1):end) + 1.0;
        ier = lsDgefa();
        ok = (ier == 0);
    end

    function corflag = lsCorfailure(told)
        ncf = ncf + 1;
        rmax = 2.0;
        tn = told;
        for j = nq:-1:1
            for i1 = j:nq
                YH(i1, :) = YH(i1, :) - YH(i1+1, :);
            end
        end
        if abs(h) <= hmin * 1.00001 || ncf == MXNCF
            corflag = 2;
            return
        end
        ipup = miter;
        corflag = 1;
    end

    function [corflag, del, delp, m] = lsCorrection(pnorm, delp, told)
        m = 0;
        rate = 0.0;
        del = 0.0;
        y = YH(1, :);
        SAVF = lsF(tn, y);
        nfe = nfe + 1;
        while true
            if m == 0
                if ipup > 0
                    ierpj = lsPrja();
                    jcur = 1;
                    ipup = 0;
                    rc = 1.0;
                    nslp = nst;
                    crate = 0.7;
                    if ~ierpj
                        corflag = lsCorfailure(told);
                        return
                    end
                end
                ACOR = zeros(1, neq);
            end
            if miter == 0
                SAVF = h * SAVF - YH(2, :);
                y = SAVF - ACOR;
                del = lsVmnorm(y);
                y = YH(1, :) + EL(1) * SAVF;
                ACOR = SAVF;
            else
                y = h * SAVF - (YH(2, :) + ACOR);
                y = lsDgesl(y);
                del = lsVmnorm(y);
                ACOR = ACOR + y;
                y = YH(1, :) + EL(1) * ACOR;
            end
            if del <= 100.0 * pnorm * ETA
                break
            end
            if m ~= 0 || meth ~= 1
                if m ~= 0
                    rm = 1024.0;
                    if del <= 1024.0 * delp
                        rm = del / delp;
                    end
                    rate = max(rate, rm);
                    crate = max(0.2 * crate, rm);
                end
                conit = 0.5 / (nq + 2);
                dcon = del * min(1.0, 1.5 * crate) / (TESCO(nq, 2) * conit);
                if dcon <= 1.0
                    pdest = max(pdest, rate / abs(h * EL(1)));
                    if pdest ~= 0.0
                        pdlast = pdest;
                    end
                    break
                end
            end
            m = m + 1;
            if m == MAXCOR || (m >= 2 && del > 2.0 * delp)
                if miter == 0 || jcur == 1
                    corflag = lsCorfailure(told);
                    return
                end
                ipup = miter;
                m = 0;
                rate = 0.0;
                del = 0.0;
                y = YH(1, :);
                SAVF = lsF(tn, y);
                nfe = nfe + 1;
            else
                delp = del;
                SAVF = lsF(tn, y);
                nfe = nfe + 1;
            end
        end
        corflag = 0;
    end

    function rh = lsMethodswitch(dsm, pnorm, rh)
        if forceStiff
            return
        end
        if meth == 1
            if nq > 5
                return
            end
            if dsm <= 100.0 * pnorm * ETA || pdest == 0.0
                if irflag == 0
                    return
                end
                rh2 = 2.0;
                nqm2 = min(nq, mxords);
            else
                exsm = 1.0 / (nq + 1);
                rh1 = 1.0 / (1.2 * dsm^exsm + 0.0000012);
                rh1it = 2.0 * rh1;
                pdh = pdlast * abs(h);
                if pdh * rh1 > 0.00001
                    rh1it = SM1(nq+1) / pdh;
                end
                rh1 = min(rh1, rh1it);
                if nq > mxords
                    nqm2 = mxords;
                    lm2 = mxords + 1;
                    exm2 = 1.0 / lm2;
                    lm2p1 = lm2 + 1;
                    dm2 = lsVmnorm(YH(lm2p1, :)) / CM2(mxords);
                    rh2 = 1.0 / (1.2 * dm2^exm2 + 0.0000012);
                else
                    dm2 = dsm * (CM1(nq) / CM2(nq));
                    rh2 = 1.0 / (1.2 * dm2^exsm + 0.0000012);
                    nqm2 = nq;
                end
                if rh2 < RATIO * rh1
                    return
                end
            end
            rh = rh2;
            icount = 20;
            meth = 2;
            miter = 2;
            pdlast = 0.0;
            nq = nqm2;
            return
        end
        exsm = 1.0 / (nq + 1);
        if mxordn < nq
            nqm1 = mxordn;
            lm1 = mxordn + 1;
            exm1 = 1.0 / lm1;
            lm1p1 = lm1 + 1;
            dm1 = lsVmnorm(YH(lm1p1, :)) / CM1(mxordn);
            rh1 = 1.0 / (1.2 * dm1^exm1 + 0.0000012);
        else
            dm1 = dsm * (CM2(nq) / CM1(nq));
            rh1 = 1.0 / (1.2 * dm1^exsm + 0.0000012);
            nqm1 = nq;
            exm1 = exsm;
        end
        rh1it = 2.0 * rh1;
        pdh = pdnorm * abs(h);
        if pdh * rh1 > 0.00001
            rh1it = SM1(nqm1+1) / pdh;
        end
        rh1 = min(rh1, rh1it);
        rh2 = 1.0 / (1.2 * dsm^exsm + 0.0000012);
        if rh1 * RATIO < 5.0 * rh2
            return
        end
        alpha = max(0.001, rh1);
        dm1 = dm1 * alpha^exm1;
        if dm1 <= 1000.0 * ETA * pnorm
            return
        end
        rh = rh1;
        icount = 20;
        meth = 1;
        miter = 0;
        pdlast = 0.0;
        nq = nqm1;
    end

    function [orderflag, rh] = lsOrderswitch(rhup, dsm, kflag, maxord)
        exsm = 1.0 / (nq + 1);
        rhsm = 1.0 / (1.2 * dsm^exsm + 0.0000012);
        rhdn = 0.0;
        pdh = 0.0;
        if nq ~= 1
            ddn = lsVmnorm(YH(nq+1, :)) / TESCO(nq, 1);
            exdn = 1.0 / nq;
            rhdn = 1.0 / (1.3 * ddn^exdn + 0.0000013);
        end
        if meth == 1
            pdh = max(abs(h) * pdlast, 0.000001);
            if (nq + 1) < maxord + 1
                rhup = min(rhup, SM1(nq+2) / pdh);
            end
            rhsm = min(rhsm, SM1(nq+1) / pdh);
            if nq > 1
                rhdn = min(rhdn, SM1(nq) / pdh);
            end
            pdest = 0.0;
        end
        if rhsm >= rhup
            if rhsm >= rhdn
                newq = nq;
                rh = rhsm;
            else
                newq = nq - 1;
                rh = rhdn;
                if kflag < 0 && rh > 1.0
                    rh = 1.0;
                end
            end
        else
            if rhup <= rhdn
                newq = nq - 1;
                rh = rhdn;
                if kflag < 0 && rh > 1.0
                    rh = 1.0;
                end
            else
                rh = rhup;
                if rh >= 1.1
                    r = EL(nq+1) / (nq + 1);
                    nq = nq + 1;
                    YH(nq+1, :) = ACOR * r;
                    orderflag = 2;
                    return
                else
                    ialth = 3;
                    orderflag = 0;
                    return
                end
            end
        end
        if meth == 1
            if rh * pdh * 1.00001 < SM1(newq+1)
                if kflag == 0 && rh < 1.1
                    ialth = 3;
                    orderflag = 0;
                    return
                end
            end
        else
            if kflag == 0 && rh < 1.1
                ialth = 3;
                orderflag = 0;
                return
            end
        end
        if kflag <= -2
            rh = min(rh, 0.2);
        end
        if newq == nq
            orderflag = 1;
            return
        end
        nq = newq;
        orderflag = 2;
    end

    function kflag = lsStoda()
        kflag = 0;
        told = tn;
        ncf = 0;
        delp = 0.0;
        maxord = mxordn;
        if meth == 2
            maxord = mxords;
        end
        if jstart == 0
            nq = 1;
            ialth = 2;
            rmax = 10000.0;
            rc = 0.0;
            crate = 0.7;
            holdh = h;
            nslp = 0;
            ipup = miter;
            EL(1) = 1.0;
            icount = 20;
            irflag = 0;
            pdest = 0.0;
            pdlast = 0.0;
            % cfode(1) in the C, which assumes meth = 1 at the start; under
            % forceStiff the start is meth = 2 and the tables must match it
            lsCfode(meth);
            lsResetcoeff();
        end
        if jstart == -1
            ipup = miter;
            if ialth == 1
                ialth = 2;
            end
            if meth ~= mused
                lsCfode(meth);
                ialth = nq + 1;
                lsResetcoeff();
            end
            if h ~= holdh
                rh = h / holdh;
                h = holdh;
                lsScaleh(rh);
            end
        end
        if jstart == -2
            if h ~= holdh
                rh = h / holdh;
                h = holdh;
                lsScaleh(rh);
            end
        end
        dsm = 0.0;
        pnorm = 0.0;
        m = 0;
        del = 0.0;
        while true
            jcur = 0;
            while true
                if abs(rc - 1.0) > CCMAX
                    ipup = miter;
                end
                if nst >= nslp + MSBP
                    ipup = miter;
                end
                tn = tn + h;
                for j = nq:-1:1
                    for i1 = j:nq
                        YH(i1, :) = YH(i1, :) + YH(i1+1, :);
                    end
                end
                pnorm = lsVmnorm(YH(1, :));
                [corflag, del, delp, m] = lsCorrection(pnorm, delp, told);
                if corflag == 0
                    break
                end
                if corflag == 1
                    rh = max(0.25, hmin / abs(h));
                    lsScaleh(rh);
                    continue
                end
                if corflag == 2
                    kflag = -2;
                    holdh = h;
                    return
                end
            end
            if m == 0
                dsm = del / TESCO(nq, 2);
            end
            if m > 0
                dsm = lsVmnorm(ACOR) / TESCO(nq, 2);
            end
            if dsm <= 1.0
                kflag = 0;
                nst = nst + 1;
                hu = h;
                nqu = nq;
                mused = meth;
                for j = 1:(nq+1)
                    YH(j, :) = YH(j, :) + EL(j) * ACOR;
                end
                icount = icount - 1;
                if icount < 0
                    rh = lsMethodswitch(dsm, pnorm, 0.0);
                    if meth ~= mused
                        rh = max(rh, hmin / abs(h));
                        lsScaleh(rh);
                        rmax = 10.0;
                        lsEndstoda();
                        break
                    end
                end
                ialth = ialth - 1;
                if ialth == 0
                    rhup = 0.0;
                    if (nq + 1) ~= maxord + 1
                        SAVF = ACOR - YH(maxord+1, :);
                        dup = lsVmnorm(SAVF) / TESCO(nq, 3);
                        exup = 1.0 / (nq + 2);
                        rhup = 1.0 / (1.4 * dup^exup + 0.0000014);
                    end
                    [orderflag, rh] = lsOrderswitch(rhup, dsm, kflag, maxord);
                    if orderflag == 0
                        lsEndstoda();
                        break
                    end
                    if orderflag == 1
                        rh = max(rh, hmin / abs(h));
                        lsScaleh(rh);
                        rmax = 10.0;
                        lsEndstoda();
                        break
                    end
                    if orderflag == 2
                        lsResetcoeff();
                        rh = max(rh, hmin / abs(h));
                        lsScaleh(rh);
                        rmax = 10.0;
                        lsEndstoda();
                        break
                    end
                end
                if ialth > 1 || (nq + 1) == maxord + 1
                    lsEndstoda();
                    break
                end
                YH(maxord+1, :) = ACOR;
                lsEndstoda();
                break
            else
                kflag = kflag - 1;
                tn = told;
                for j = nq:-1:1
                    for i1 = j:nq
                        YH(i1, :) = YH(i1, :) - YH(i1+1, :);
                    end
                end
                rmax = 2.0;
                if abs(h) <= hmin * 1.00001
                    kflag = -1;
                    holdh = h;
                    break
                end
                if kflag > -3
                    [orderflag, rh] = lsOrderswitch(0.0, dsm, kflag, maxord);
                    if orderflag == 1 || orderflag == 0
                        if orderflag == 0
                            rh = min(rh, 0.2);
                        end
                        rh = max(rh, hmin / abs(h));
                        lsScaleh(rh);
                    end
                    if orderflag == 2
                        lsResetcoeff();
                        rh = max(rh, hmin / abs(h));
                        lsScaleh(rh);
                    end
                    continue
                else
                    if kflag == -10
                        kflag = -1;
                        holdh = h;
                        break
                    else
                        rh = 0.1;
                        rh = max(hmin / abs(h), rh);
                        h = h * rh;
                        y = YH(1, :);
                        SAVF = lsF(tn, y);
                        nfe = nfe + 1;
                        YH(2, :) = h * SAVF;
                        ipup = miter;
                        ialth = 5;
                        if nq == 1
                            continue
                        end
                        nq = 1;
                        lsResetcoeff();
                        continue
                    end
                end
            end
        end
    end

    function lsEndstoda()
        r = 1.0 / TESCO(nqu, 2);
        ACOR = ACOR * r;
        holdh = h;
    end

    function lsDrive(tout, itask)
        % one call of the C driver lsoda(), integrating from tcur towards tout
        h0 = h0opt;
        ihit = false;
        if state == 1
            if (tout - tcur) * h0 < 0.0
                error('lsoda:failed', '[lsoda] tout = %g behind t = %g, integration direction is given by %g', tout, tcur, h0);
            end
        end
        if state == 3
            jstart = -1;
        end
        if state == 1
            if forceStiff
                meth = 2;
                miter = 2;
            else
                meth = 1;
            end
            tn = tcur;
            tsw = tcur;
            if itask == 4 || itask == 5
                if (tcrit - tout) * (tout - tcur) < 0.0
                    error('lsoda:failed', '[lsoda] itask = 4 or 5 and tcrit behind tout');
                end
                if h0 ~= 0.0 && (tcur + h0 - tcrit) * h0 > 0.0
                    h0 = tcrit - tcur;
                end
            end
            jstart = 0;
            nq = 1;
            YH(2, :) = lsF(tcur, y);
            nfe = 1;
            YH(1, :) = y;
            lsEwset(y);
            if any(EWT <= 0.0)
                ibad = find(EWT <= 0.0, 1);
                error('lsoda:failed', '[lsoda] ewt[%d] = %g <= 0.', ibad, EWT(ibad));
            end
            if h0 == 0.0
                tdist = abs(tout - tcur);
                w0 = max(abs(tcur), abs(tout));
                if tdist < 2.0 * ETA * w0
                    error('lsoda:failed', '[lsoda] tout too close to t to start integration');
                end
                tol = max(rtolv);
                if tol <= 0.0
                    ay = abs(y);
                    nz = ay ~= 0.0;
                    if any(nz)
                        tol = max(tol, max(atolv(nz) ./ ay(nz)));
                    end
                end
                tol = max(tol, 100.0 * ETA);
                tol = min(tol, 0.001);
                sumv = lsVmnorm(YH(2, :));
                sumv = 1.0 / (tol * w0 * w0) + tol * sumv * sumv;
                h0 = 1.0 / sqrt(sumv);
                h0 = min(h0, tdist);
                if tout - tcur >= 0.0
                    h0 = h0 * 1.0;
                else
                    h0 = h0 * (-1.0);
                end
            end
            rh = abs(h0) * hmxi;
            if rh > 1.0
                h0 = h0 / rh;
            end
            h = h0;
            YH(2, :) = YH(2, :) * h0;
        end
        if state == 2 || state == 3
            % jstart is NOT reset to 1 here, though the C driver does reset it:
            % JSTART lives in the FORTRAN common block, dstoda leaves it at 1 and
            % the driver overrides it with -1 when a method switch must be
            % completed on the next step. Rebuilding it as 1 on every
            % continuation call drops that -1, which is invisible under itask=1
            % over sparse output times and fatal in STEPPING mode, where every
            % step returns: the switch is never completed and the elco tables
            % stay on the old method.
            nslast = nst;
            switch itask
                case 1
                    if (tn - tout) * h >= 0.0
                        doIntdyReturn(tout, itask);
                        return
                    end
                case 2
                    % nothing to check
                case 3
                    tp = tn - hu * (1.0 + 100.0 * ETA);
                    if (tp - tout) * h > 0.0
                        error('lsoda:failed', '[lsoda] itask = %d and tout behind tcur - hu', itask);
                    end
                    if (tn - tout) * h >= 0.0
                        doSuccessReturn(itask, false);
                        return
                    end
                case {4, 5}
                    % case 4 falls through into case 5 in the C driver
                    if itask == 4
                        if (tn - tcrit) * h > 0.0
                            error('lsoda:failed', '[lsoda] itask = 4 or 5 and tcrit behind tcur');
                        end
                        if (tcrit - tout) * h < 0.0
                            error('lsoda:failed', '[lsoda] itask = 4 or 5 and tcrit behind tout');
                        end
                        if (tn - tout) * h >= 0.0
                            doIntdyReturn(tout, itask);
                            return
                        end
                    else
                        if (tn - tcrit) * h > 0.0
                            error('lsoda:failed', '[lsoda] itask = 4 or 5 and tcrit behind tcur');
                        end
                    end
                    hmx = abs(tn) + abs(h);
                    ihit = abs(tn - tcrit) <= (100.0 * ETA * hmx);
                    if ihit
                        tcur = tcrit;
                        doSuccessReturn(itask, ihit);
                        return
                    end
                    tnext = tn + h * (1.0 + 4.0 * ETA);
                    if (tnext - tcrit) * h > 0.0
                        h = (tcrit - tn) * (1.0 - 4.0 * ETA);
                        if state == 2
                            jstart = -2;
                        end
                    end
                otherwise
                    error('lsoda:itask', '[lsoda] illegal itask = %d', itask);
            end
        end
        while true
            if state ~= 1 || nst ~= 0
                if (nst - nslast) >= mxstep
                    doSoftFailure(-1, sprintf('[lsoda] %d steps taken before reaching tout', mxstep));
                    return
                end
                lsEwset(YH(1, :));
                if any(EWT <= 0.0)
                    ibad = find(EWT <= 0.0, 1);
                    doSoftFailure(-6, sprintf('[lsoda] ewt[%d] = %g <= 0.', ibad, EWT(ibad)));
                    return
                end
            end
            tolsf = ETA * lsVmnorm(YH(1, :));
            if tolsf > 0.01
                tolsf = tolsf * 200.0;
                if nst == 0
                    error('lsoda:failed', ['[lsoda] at start of problem, too much accuracy requested ' ...
                        'for precision of machine, suggested scaling factor = %g'], tolsf);
                end
                doSoftFailure(-2, sprintf(['[lsoda] at t = %g, too much accuracy requested for precision ' ...
                    'of machine, suggested scaling factor = %g'], tcur, tolsf));
                return
            end
            if (tn + h) == tn
                nhnil = nhnil + 1;
                if nhnil <= mxhnil
                    warning('lsoda:hnil', '[lsoda] internal t = %g and h = %g are such that t + h = t on the next step', tn, h);
                end
            end
            kflag = lsStoda();
            if kflag == 0
                jstart = 1;
                if meth ~= mused
                    tsw = tn;
                    jstart = -1;
                end
                if itask == 1
                    if (tn - tout) * h < 0.0
                        continue
                    end
                    doIntdyReturn(tout, itask);
                    return
                end
                if itask == 2
                    doSuccessReturn(itask, ihit);
                    return
                end
                if itask == 3
                    if (tn - tout) * h >= 0.0
                        doSuccessReturn(itask, ihit);
                        return
                    end
                    continue
                end
                if itask == 4
                    if (tn - tout) * h >= 0.0
                        doIntdyReturn(tout, itask);
                        return
                    else
                        hmx = abs(tn) + abs(h);
                        ihit = abs(tn - tcrit) <= (100.0 * ETA * hmx);
                        if ihit
                            doSuccessReturn(itask, ihit);
                            return
                        end
                        tnext = tn + h * (1.0 + 4.0 * ETA);
                        if (tnext - tcrit) * h <= 0.0
                            continue
                        end
                        h = (tcrit - tn) * (1.0 - 4.0 * ETA);
                        jstart = -2;
                        continue
                    end
                end
                if itask == 5
                    hmx = abs(tn) + abs(h);
                    ihit = abs(tn - tcrit) <= (100.0 * ETA * hmx);
                    doSuccessReturn(itask, ihit);
                    return
                end
            end
            if kflag == -1 || kflag == -2
                big = 0.0;
                imxer = 1;
                for i = 1:neq
                    sz = abs(ACOR(i)) * EWT(i);
                    if big < sz
                        big = sz;
                        imxer = i;
                    end
                end
                if kflag == -1
                    doSoftFailure(-4, sprintf(['[lsoda] at t = %g and step size h = %g, the error test failed ' ...
                        'repeatedly or with abs(h) = hmin'], tn, h));
                    return
                end
                if kflag == -2
                    doSoftFailure(-5, sprintf(['[lsoda] at t = %g and step size h = %g, the corrector convergence ' ...
                        'failed repeatedly or with abs(h) = hmin'], tn, h));
                    return
                end
            end
        end
    end

    function doSuccessReturn(itask, ihit)
        y = YH(1, :);
        tcur = tn;
        if (itask == 4 || itask == 5) && ihit
            tcur = tcrit;
        end
        state = 2;
    end

    function doIntdyReturn(tout, itask)
        [iflag, dky] = lsIntdy(tout, 0);
        if iflag ~= 0
            warning('lsoda:intdy', '[lsoda] trouble from intdy, itask = %d, tout = %g', itask, tout);
            y = YH(1, :);
            tcur = tn;
        else
            y = dky;
        end
        tcur = tout;
        state = 2;
    end

    function doSoftFailure(code, msg)
        y = YH(1, :);
        tcur = tn;
        state = code;
        errmsg = msg;
    end

end

function v = optField(opts, name, defval)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = defval;
end
end

function v = expandTol(tolval, neq, name)
tolval = tolval(:).';
if numel(tolval) == 1
    v = repmat(tolval, 1, neq);
elseif numel(tolval) == neq
    v = tolval;
else
    error('lsoda:tol', '%s must be a scalar or have one entry per equation (%d).', name, neq);
end
end
