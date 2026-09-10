function [y, idid, stats, x] = rodas_core(n, fcn, x, y, xend, h, rtol, atol, itol, opt)
% RODAS_CORE  Rosenbrock method of order (3)4 for  M y' = f(x,y),  M singular
%
%   [Y, IDID, STATS, X] = RODAS_CORE(N, FCN, X, Y, XEND, H, RTOL, ATOL, ITOL, OPT)
%
%   PORTED THIRD-PARTY CODE -- do not edit to fix a LINE bug; fix the caller.
%
%   Source: E. Hairer and G. Wanner, rodas.f / dc_decsol.f / decsol.f, version
%   of October 28, 1996, as published with "Solving Ordinary Differential
%   Equations II. Stiff and Differential-Algebraic Problems", Springer Series
%   in Computational Mathematics 14.
%
%   WHY THIS EXISTS IN MATLAB, WHERE ODE15S ALREADY CARRIES A SINGULAR MASS
%   MATRIX. It is not that ode15s cannot do it -- it can, and RODAS.M leaves it
%   reachable. It is that the fluid 'dae' method now runs in four codebases,
%   and C++, native Python and the JAR have no ode15s: LSODA integrates y' = f
%   and cannot carry a mass matrix at all, and scipy's BDF/Radau take no M
%   either. RODAS is a FIXED SEQUENCE OF SIX LINEAR SOLVES against one real
%   matrix rather than an iteration with a convergence policy, so a port has no
%   iteration history to diverge on -- which is what lets all four agree in the
%   last digits on the DAE route. See _kb/06-solver-catalog.md.
%
%   THE TRANSLATION IS INDEX-FOR-INDEX. MATLAB indexes from 1 and so does
%   Fortran, so the loop bounds and array subscripts below are the original's
%   own, unshifted. That is deliberate: a re-based rewrite of 2000 lines of
%   Fortran is where an off-by-one hides, and it would be invisible until it
%   changed an answer nobody has a reference for.
%
%   WHERE A LOOP IS VECTORISED IT IS ELEMENTWISE, NEVER A REDUCTION. A
%   vectorised A(k+1:n,j) = A(k+1:n,j) + A(k+1:n,k)*t does the same multiply
%   and the same add on each element in the same order as the Fortran DO loop,
%   so it is bit identical. A vectorised SUM is not: MATLAB may pairwise-sum
%   where Fortran accumulates sequentially, and the two differ in the last
%   bits. Every reduction here is therefore still a loop -- the error norm, the
%   inner sums of DECOMR and SLVROD -- and that is not an oversight to be
%   optimised away.
%
%   Reachable IJOB values are 1..5 (and 11..15 when OPT.m1 > 0, the second
%   order form); RODAS itself never selects 6 or 7, which belong to RADAU5's
%   Hessenberg option and are absent here rather than transliterated dead.
%
%   Inputs
%     N     - dimension of the system
%     FCN   - @(x,y) -> f, the right hand side, length N
%     X     - initial abscissa
%     Y     - initial state, length N
%     XEND  - final abscissa (XEND-X may be positive or negative)
%     H     - initial step size guess (0 means 1e-6)
%     RTOL, ATOL, ITOL - tolerances; ITOL=0 scalar, ITOL=1 one per equation
%     OPT   - struct of the Fortran switches, all optional:
%       ifcn   0 f is autonomous (default), 1 f may depend on x
%       jac    @(x,y) -> dfy (LDJAC x N);  ijac 1 to use it, 0 for differences
%       mljac  lower bandwidth of the Jacobian, N for full (default N)
%       mujac  upper bandwidth of the Jacobian
%       dfx    @(x,y) -> df/dx;  idfx 1 to use it, 0 for differences
%       mas    @() -> am (LDMAS x NM1);  imas 1 to use it, 0 for M = I
%       mlmas, mumas  bandwidths of the mass matrix, N for full
%       solout @(nr,xold,x,y,dense) -> irtrn, called after every accepted step
%              when iout = 1; return a negative irtrn to stop. DENSE feeds
%              RODAS_CONTRO, the third-order interpolant valid over the step
%              just accepted.
%       iout   1 to call solout (default 0)
%       nmax   maximal number of steps (default 100000)
%       meth   coefficient set 1 (default), 2 or 3
%       pred   true for Gustafsson's predictive controller (default)
%       uround, hmax, fac1, fac2, safe, m1, m2 -- as in rodas.f
%
%   Outputs
%     Y     - the state at X
%     IDID  - 1 success, 2 stopped by SOLOUT, -2 more than NMAX steps,
%             -3 step size too small, -4 matrix repeatedly singular
%     STATS - struct: nfcn, njac, nstep, naccpt, nrejct, ndec, nsol
%     X     - the abscissa actually reached
%
%   See also RODAS, RODAS_CONTRO, ODE15S.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 10 || isempty(opt), opt = struct(); end
opt = local_defaults(opt, n);

y = y(:);
nm1 = n - opt.m1;
m1 = opt.m1;
m2 = opt.m2;
if m1 == 0, m2 = n; end
if m2 == 0, m2 = m1; end
if m1 < 0 || m2 < 0 || m1 + m2 > n
    error('rodas_core:input', 'CURIOUS INPUT FOR IWORK(9,10)= %d %d', m1, m2);
end
if opt.nmax <= 0
    error('rodas_core:input', 'WRONG INPUT IWORK(1)= %d', opt.nmax);
end
if opt.meth <= 0 || opt.meth >= 4
    error('rodas_core:input', 'CURIOUS INPUT IWORK(2)= %d', opt.meth);
end
uround = opt.uround;
if uround < 1e-16 || uround >= 1
    error('rodas_core:input', 'COEFFICIENTS HAVE 16 DIGITS, UROUND= %g', uround);
end
hmax = opt.hmax;
if isempty(hmax), hmax = xend - x; end
fac1 = opt.fac1; fac2 = opt.fac2;
if fac1 < 1 || fac2 > 1
    error('rodas_core:input', 'CURIOUS INPUT WORK(3,4)');
end
if opt.safe <= .001 || opt.safe >= 1
    error('rodas_core:input', 'CURIOUS INPUT FOR WORK(5)= %g', opt.safe);
end

if itol == 0
    rtolv = rtol(1); atolv = atol(1);
    if atolv <= 0 || rtolv <= uround*10
        error('rodas_core:input', 'TOLERANCES ARE TOO SMALL');
    end
    rtolv = repmat(rtolv, n, 1); atolv = repmat(atolv, n, 1);
else
    rtolv = rtol(:); atolv = atol(:);
    for i = 1:n
        if atolv(i) <= 0 || rtolv(i) <= uround*10
            error('rodas_core:input', 'TOLERANCES(%d) ARE TOO SMALL', i);
        end
    end
end

autnms = (opt.ifcn == 0);
implct = (opt.imas ~= 0);
mljac = opt.mljac; mujac = opt.mujac;
mlmas = opt.mlmas; mumas = opt.mumas;
jband = (mljac < nm1);

if jband
    ldjac = mljac + mujac + 1;
    lde = mljac + ldjac;
else
    mljac = nm1; mujac = nm1; ldjac = nm1; lde = nm1;
end
if implct
    if mlmas ~= nm1
        ldmas = mlmas + mumas + 1;
        if jband, ijob = 4; else, ijob = 3; end
    else
        ldmas = nm1;
        ijob = 5;
    end
    if mlmas > mljac || mumas > mujac
        error('rodas_core:input', ...
            'BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC"');
    end
else
    ldmas = 0;
    if jband, ijob = 2; else, ijob = 1; end
end
ldmas = max(1, ldmas);

% ---- the /LINAL/ COMMON block, as a struct ----
lin = struct('mle',0,'mue',0,'mbjac',0,'mbb',0,'mdiag',0,'mdiff',0,'mbdiag',0);
lin.mbdiag = mumas + 1;
if jband
    lin.mle = mljac;
    lin.mue = mujac;
    lin.mbjac = mljac + mujac + 1;
    lin.mbb = mlmas + mumas + 1;
    lin.mdiag = lin.mle + lin.mue + 1;
    lin.mdiff = lin.mle + lin.mue - mumas;
end

% ---- workspace ----
ynew = zeros(n,1); dy1 = zeros(n,1); dy = zeros(n,1);
ak1 = zeros(n,1); ak2 = zeros(n,1); ak3 = zeros(n,1);
ak4 = zeros(n,1); ak5 = zeros(n,1); ak6 = zeros(n,1);
fx = zeros(n,1); cont = zeros(4*n,1);
fjac = zeros(ldjac, n); e = zeros(lde, nm1); fmas = zeros(ldmas, nm1);
ip = zeros(nm1, 1);

if implct
    fmas = opt.mas();
end
c = rodas_coefficients(opt.meth);

if m1 > 0, ijob = ijob + 10; end

posneg = local_sign(1, xend - x);
hmaxn = min(abs(hmax), abs(xend - x));
if abs(h) <= uround*10, h = 1e-6; end
h = min(abs(h), hmaxn);
h = local_sign(h, posneg);
reject = false; last = false; nsing = 0;
hd1 = 0; hd2 = 0; hd3 = 0; hd4 = 0;
hacc = 0; erracc = 0; hopt = h;
nn2 = 2*n; nn3 = 3*n;

nfcn = 0; njac = 0; nstep = 0; naccpt = 0; nrejct = 0; ndec = 0; nsol = 0;
dense = struct('cont', cont, 'xold', x, 'h', h, 'n', n);
idid = 0;

if opt.iout ~= 0
    dense.xold = x; dense.h = h; dense.cont = cont;
    irtrn = opt.solout(naccpt+1, dense.xold, x, y, dense);
    if ~isempty(irtrn) && irtrn < 0
        idid = 2;
        stats = local_stats(nfcn, njac, nstep, naccpt, nrejct, ndec, nsol);
        return
    end
end

% ==== BASIC INTEGRATION STEP (label 1) ====
while true
    if nstep > opt.nmax
        idid = -2; break
    end
    if abs(h)*.1 <= abs(x)*uround
        idid = -3; break
    end
    if last
        h = hopt; idid = 1; break
    end
    hopt = h;
    if (x + h*1.0001 - xend)*posneg >= 0
        h = xend - x;
        last = true;
    end

    % ---- COMPUTATION OF THE JACOBIAN ----
    dy1 = fcn(x, y); dy1 = dy1(:);
    nfcn = nfcn + 1;
    njac = njac + 1;
    if opt.ijac == 0
        if jband
            mujacp = mujac + 1;
            md = min(lin.mbjac, n);
            for mm = 1:(m1/m2 + 1)
                for k = 1:md
                    j = k + (mm-1)*m2;
                    while true
                        ak2(j) = y(j);
                        ak3(j) = sqrt(uround*max(1e-5, abs(y(j))));
                        y(j) = y(j) + ak3(j);
                        j = j + md;
                        if j > mm*m2, break, end
                    end
                    ak1 = fcn(x, y); ak1 = ak1(:);
                    j = k + (mm-1)*m2;
                    j1 = k;
                    lbeg = max(1, j1 - mujac) + m1;
                    while true
                        lend = min(m2, j1 + mljac) + m1;
                        y(j) = ak2(j);
                        mujacj = mujacp - j1 - m1;
                        for l = lbeg:lend
                            fjac(l + mujacj, j) = (ak1(l) - dy1(l))/ak3(j);
                        end
                        j = j + md;
                        j1 = j1 + md;
                        lbeg = lend + 1;
                        if j > mm*m2, break, end
                    end
                end
            end
        else
            for i = 1:n
                ysafe = y(i);
                delt = sqrt(uround*max(1e-5, abs(ysafe)));
                y(i) = ysafe + delt;
                ak1 = fcn(x, y); ak1 = ak1(:);
                fjac(1:n-m1, i) = (ak1(m1+1:n) - dy1(m1+1:n))/delt;
                y(i) = ysafe;
            end
        end
    else
        fjac = opt.jac(x, y);
    end
    if ~autnms
        if opt.idfx == 0
            delt = sqrt(uround*max(1e-5, abs(x)));
            ak1 = fcn(x + delt, y); ak1 = ak1(:);
            fx = (ak1 - dy1)/delt;
        else
            fx = opt.dfx(x, y); fx = fx(:);
        end
    end

    % ==== COMPUTE THE STAGES (label 2) ====
    stepdone = false;
    while ~stepdone
        fac = 1/(h*c.gamma);
        [e, ip, ier] = rodas_decomr(n, fjac, fmas, mlmas, mumas, m1, m2, ...
            nm1, fac, e, ip, ijob, lin);
        if ier ~= 0
            % ---- SINGULAR MATRIX (label 80) ----
            nsing = nsing + 1;
            if nsing >= 5
                idid = -4; stepdone = true; break
            end
            h = h*.5;
            reject = true; last = false;
            continue
        end
        ndec = ndec + 1;

        hc21 = c.c21/h; hc31 = c.c31/h; hc32 = c.c32/h;
        hc41 = c.c41/h; hc42 = c.c42/h; hc43 = c.c43/h;
        hc51 = c.c51/h; hc52 = c.c52/h; hc53 = c.c53/h; hc54 = c.c54/h;
        hc61 = c.c61/h; hc62 = c.c62/h; hc63 = c.c63/h;
        hc64 = c.c64/h; hc65 = c.c65/h;
        if ~autnms
            hd1 = h*c.d1; hd2 = h*c.d2; hd3 = h*c.d3; hd4 = h*c.d4;
        end

        % ---- THE STAGES ----
        ak1 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy1, ak1, fx, ynew, hd1, ijob, false, lin);
        ynew = y + c.a21*ak1;
        dy = fcn(x + c.c2*h, ynew); dy = dy(:);
        ynew = hc21*ak1;
        ak2 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy, ak2, fx, ynew, hd2, ijob, true, lin);

        ynew = y + c.a31*ak1 + c.a32*ak2;
        dy = fcn(x + c.c3*h, ynew); dy = dy(:);
        ynew = hc31*ak1 + hc32*ak2;
        ak3 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy, ak3, fx, ynew, hd3, ijob, true, lin);

        ynew = y + c.a41*ak1 + c.a42*ak2 + c.a43*ak3;
        dy = fcn(x + c.c4*h, ynew); dy = dy(:);
        ynew = hc41*ak1 + hc42*ak2 + hc43*ak3;
        ak4 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy, ak4, fx, ynew, hd4, ijob, true, lin);

        ynew = y + c.a51*ak1 + c.a52*ak2 + c.a53*ak3 + c.a54*ak4;
        dy = fcn(x + h, ynew); dy = dy(:);
        ak6 = hc52*ak2 + hc54*ak4 + hc51*ak1 + hc53*ak3;
        ak5 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy, ak5, fx, ak6, 0, ijob, true, lin);

        % ---- EMBEDDED SOLUTION ----
        ynew = ynew + ak5;
        dy = fcn(x + h, ynew); dy = dy(:);
        cont(1:n) = hc61*ak1 + hc62*ak2 + hc65*ak5 + hc64*ak4 + hc63*ak3;
        ak6 = rodas_slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, ...
            m2, nm1, fac, e, ip, dy, ak6, fx, cont(1:n), 0, ijob, true, lin);

        % ---- NEW SOLUTION ----
        ynew = ynew + ak6;
        nsol = nsol + 6;
        nfcn = nfcn + 5;

        % ---- DENSE OUTPUT ----
        if opt.iout ~= 0
            cont(1:n) = y;
            cont(nn2+1:nn2+n) = c.d21*ak1 + c.d22*ak2 + c.d23*ak3 + ...
                c.d24*ak4 + c.d25*ak5;
            cont(nn3+1:nn3+n) = c.d31*ak1 + c.d32*ak2 + c.d33*ak3 + ...
                c.d34*ak4 + c.d35*ak5;
        end

        % ==== ERROR ESTIMATION ====
        nstep = nstep + 1;
        % A SEQUENTIAL SUM, as the Fortran accumulates it. A vectorised
        % sum((ak6./sk).^2) may pairwise-sum and differ in the last bits, which
        % is enough to move a step-size decision and desynchronise the run.
        err = 0;
        for i = 1:n
            sk = atolv(i) + rtolv(i)*max(abs(y(i)), abs(ynew(i)));
            q = ak6(i)/sk;
            err = err + q*q;
        end
        err = sqrt(err/n);

        % ---- COMPUTATION OF HNEW, .2 <= HNEW/H <= 6 ----
        fac = max(fac2, min(fac1, err^0.25/opt.safe));
        hnew = h/fac;

        if err <= 1
            % ---- STEP IS ACCEPTED ----
            naccpt = naccpt + 1;
            if opt.pred
                if naccpt > 1
                    facgus = (hacc/h)*(err*err/erracc)^0.25/opt.safe;
                    facgus = max(fac2, min(fac1, facgus));
                    fac = max(fac, facgus);
                    hnew = h/fac;
                end
                hacc = h;
                erracc = max(.01, err);
            end
            y = ynew;
            dense.xold = x;
            x = x + h;
            if opt.iout ~= 0
                cont(n+1:2*n) = y;
                dense.h = h; dense.cont = cont;
                irtrn = opt.solout(naccpt+1, dense.xold, x, y, dense);
                if ~isempty(irtrn) && irtrn < 0
                    idid = 2; stepdone = true; break
                end
            end
            if abs(hnew) > hmaxn, hnew = posneg*hmaxn; end
            if reject, hnew = posneg*min(abs(hnew), abs(h)); end
            reject = false;
            h = hnew;
            stepdone = true;   % back to label 1
        else
            % ---- STEP IS REJECTED ----
            reject = true;
            last = false;
            h = hnew;
            if naccpt >= 1, nrejct = nrejct + 1; end
            % back to label 2
        end
    end
    if idid ~= 0
        break
    end
end

stats = local_stats(nfcn, njac, nstep, naccpt, nrejct, ndec, nsol);
end

% ---------------------------------------------------------------------------
function s = local_stats(nfcn, njac, nstep, naccpt, nrejct, ndec, nsol)
s = struct('nfcn', nfcn, 'njac', njac, 'nstep', nstep, 'naccpt', naccpt, ...
    'nrejct', nrejct, 'ndec', ndec, 'nsol', nsol);
end

% ---------------------------------------------------------------------------
function v = local_sign(a, b)
% D_SIGN: |a| carrying the sign of b, with zero counted positive.
v = abs(a);
if b < 0, v = -v; end
end

% ---------------------------------------------------------------------------
function opt = local_defaults(opt, n)
d = struct('ifcn',0, 'jac',[], 'ijac',0, 'mljac',n, 'mujac',0, ...
    'dfx',[], 'idfx',0, 'mas',[], 'imas',0, 'mlmas',0, 'mumas',0, ...
    'solout',[], 'iout',0, 'nmax',100000, 'meth',1, 'pred',true, ...
    'uround',1e-16, 'hmax',[], 'fac1',5, 'fac2',.16666666666666666, ...
    'safe',0.9, 'm1',0, 'm2',0);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(opt, f{i}) || (isempty(opt.(f{i})) && ~isempty(d.(f{i})) ...
            && ~any(strcmp(f{i}, {'jac','dfx','mas','solout','hmax'})))
        opt.(f{i}) = d.(f{i});
    end
end
if ~isfield(opt,'hmax'), opt.hmax = []; end
end
