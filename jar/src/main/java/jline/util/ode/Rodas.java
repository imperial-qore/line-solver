package jline.util.ode;

/**
 * RODAS -- Rosenbrock method of order (3)4 for stiff and differential-algebraic
 * systems {@code M y' = f(x,y)}, including a SINGULAR mass matrix (index-1 DAE).
 *
 * <p>PORTED THIRD-PARTY CODE -- do not edit to fix a LINE bug; fix the caller.
 *
 * <p>Source: E. Hairer and G. Wanner, {@code rodas.f} / {@code dc_decsol.f} /
 * {@code decsol.f}, version of October 28, 1996, as published with "Solving
 * Ordinary Differential Equations II. Stiff and Differential-Algebraic
 * Problems", Springer Series in Computational Mathematics 14. Licence in
 * {@code matlab/lib/thirdparty/rodas/LICENSE}, which covers every port of it in
 * this repository.
 *
 * <p>WHY THIS EXISTS IN THE JAR. The fluid solver integrates with LSODA
 * ({@link jline.solvers.fluid.LSODAExt}), which solves {@code y' = f} and
 * cannot carry a mass matrix at all -- let alone a singular one. The {@code dae}
 * method writes population conservation as an ALGEBRAIC row, one per closed
 * chain, so its transient is an index-1 DAE with a singular M and has no LSODA
 * route. RODAS is a fixed sequence of six linear solves against one real matrix
 * rather than an iteration with a convergence policy, so a port has no
 * iteration history to diverge on -- which is what lets MATLAB, C++, native
 * Python and this agree in the last digits on that route. See
 * {@code _kb/06-solver-catalog.md}.
 *
 * <p>THE TRANSLATION IS INDEX-FOR-INDEX. Every array here is allocated one
 * longer than it needs to be and addressed from 1, exactly as the Fortran does,
 * and the loop bounds are the Fortran's own. That is deliberate: a 0-based
 * rewrite of 2000 lines of Fortran is where an off-by-one hides, and it would be
 * invisible until it changed an answer nobody has a reference for. Index 0 of
 * every internal array is unused. The CALLBACKS are the exception and are
 * 0-based, because they are the boundary a Java caller actually writes; the
 * copy at that boundary is O(n) against the O(n^3) factorisation it sits next
 * to.
 *
 * <p>EVERY REDUCTION IS A SEQUENTIAL LOOP, as the Fortran accumulates it -- the
 * error norm and the inner sums of DECOMR and SLVROD. Reassociating them would
 * change the last bits, which is enough to move a step-size decision and
 * desynchronise the whole trajectory from the other three codebases.
 *
 * <p>Reachable IJOB values are 1..5 (and 11..15 when {@code m1 > 0}, the second
 * order form); RODAS itself never selects 6 or 7, which belong to RADAU5's
 * Hessenberg option and are absent here rather than transliterated dead.
 *
 * <p>Java 8: no {@code var}, no {@code List.of}, no switch expressions.
 *
 * @see jline.solvers.fluid.analyzers.DaeAnalyzer
 */
public final class Rodas {

    private Rodas() {
    }

    /** IDID: the integration reached XEND. */
    public static final int IDID_SUCCESS = 1;
    /** IDID: SOLOUT asked for the integration to stop. */
    public static final int IDID_SOLOUT_STOP = 2;
    /** IDID: more than NMAX steps were needed. */
    public static final int IDID_NMAX = -2;
    /** IDID: the step size became too small. */
    public static final int IDID_STEP_TOO_SMALL = -3;
    /** IDID: the matrix was repeatedly singular. */
    public static final int IDID_SINGULAR = -4;

    /** The right hand side f(x,y). Arrays are 0-based and of length n. */
    public interface Fcn {
        void eval(double x, double[] y, double[] f);
    }

    /** df/dy. {@code dfy} is 0-based, (ldjac x n), row-major. */
    public interface Jac {
        void eval(double x, double[] y, double[][] dfy);
    }

    /** df/dx, for a non-autonomous system. */
    public interface Dfx {
        void eval(double x, double[] y, double[] fx);
    }

    /** The mass matrix. {@code am} is 0-based, (ldmas x nm1), row-major. */
    public interface Mas {
        void eval(double[][] am);
    }

    /** Called after each accepted step when {@code iout != 0}; negative stops. */
    public interface Solout {
        int eval(int nr, double xold, double x, double[] y, Dense dense);
    }

    /**
     * CONTRO: the third-order interpolant RODAS carries over the step just
     * accepted.
     *
     * <p>A Rosenbrock method chooses its step from the local error, so its
     * accepted points are wherever the stiffness put them and never the ones a
     * caller asked for. This is what makes an arbitrary output grid cost no
     * extra step. On the FIRST call -- made before any step, with {@code nr == 1}
     * -- {@code cont} holds no coefficients yet and the state there is {@code y}
     * itself, so callers must not interpolate at that point.
     */
    public static final class Dense {
        final double[] cont;
        final int n;
        double xold;
        double h;

        Dense(double[] cont, int n) {
            this.cont = cont;
            this.n = n;
        }

        /** Component {@code i} (0-BASED) of the solution at {@code x}. */
        public double value(int i, double x) {
            int ii = i + 1;
            double s = (x - xold) / h;
            return cont[ii] * (1 - s)
                    + s * (cont[ii + n] + (1 - s) * (cont[ii + 2 * n] + s * cont[ii + 3 * n]));
        }
    }

    /** The Fortran switches, with rodas.f's own defaults. */
    public static final class Options {
        /** 0: f is autonomous; 1: f may depend on x. */
        public int ifcn = 0;
        /** 1 to call {@link #jac}, 0 for finite differences. */
        public int ijac = 0;
        public Jac jac = null;
        /** Lower bandwidth of the Jacobian; n means full. -1 defaults to n. */
        public int mljac = -1;
        public int mujac = 0;
        /** 1 to call {@link #dfx}, 0 for finite differences. */
        public int idfx = 0;
        public Dfx dfx = null;
        /** 1 to call {@link #mas}, 0 for M = I. */
        public int imas = 0;
        public Mas mas = null;
        /** Lower bandwidth of the mass matrix; n means full. */
        public int mlmas = 0;
        public int mumas = 0;
        /** 1 to call {@link #solout} after each accepted step. */
        public int iout = 0;
        public Solout solout = null;
        public int nmax = 100000;
        /** Coefficient set: 1 (default), 2 or 3. */
        public int meth = 1;
        /** Gustafsson's predictive step-size controller. */
        public boolean pred = true;
        public double uround = 1e-16;
        /** Largest step; NaN defaults to xend - x. */
        public double hmax = Double.NaN;
        public double fac1 = 5.0;
        public double fac2 = .16666666666666666;
        public double safe = 0.9;
        /** Second-order form; 0 disables it, which is the usual case. */
        public int m1 = 0;
        public int m2 = 0;
    }

    /** Where the integration ended, and what it cost. */
    public static final class Result {
        public double[] y;
        public double x;
        public double h;
        public int idid;
        public int nfcn;
        public int njac;
        public int nstep;
        public int naccpt;
        public int nrejct;
        public int ndec;
        public int nsol;

        public boolean isSuccess() {
            return idid > 0;
        }

        @Override
        public String toString() {
            return "Rodas.Result(idid=" + idid + ", x=" + x + ", nfcn=" + nfcn
                    + ", njac=" + njac + ", nstep=" + nstep + ", naccpt=" + naccpt
                    + ", nrejct=" + nrejct + ")";
        }
    }

    /** Raised for an input rodas.f itself rejects before integrating. */
    public static class RodasInputException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public RodasInputException(String msg) {
            super(msg);
        }
    }

    /** The /LINAL/ COMMON block: band offsets shared by DECOMR and SLVROD. */
    private static final class Linal {
        int mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
    }

    // -----------------------------------------------------------------------
    // rodas.f: ROCOE, the method coefficients
    // -----------------------------------------------------------------------

    /** The Rosenbrock tableau. */
    static final class Coef {
        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
        double c21, c31, c32, c41, c42, c43, c51, c52, c53, c54;
        double c61, c62, c63, c64, c65;
        double gamma, c2, c3, c4, d1, d2, d3, d4;
        double d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
    }

    static Coef rocoe(int meth) {
        Coef c = new Coef();
        switch (meth) {
            case 1:
                c.c2 = .386; c.c3 = .21; c.c4 = .63;
                c.d1 = .25; c.d2 = -.1043; c.d3 = .1035; c.d4 = -.03620000000000023;
                c.a21 = 1.544;
                c.a31 = .9466785280815826; c.a32 = .2557011698983284;
                c.a41 = 3.314825187068521; c.a42 = 2.896124015972201;
                c.a43 = .9986419139977817;
                c.a51 = 1.221224509226641; c.a52 = 6.019134481288629;
                c.a53 = 12.53708332932087; c.a54 = -.687886036105895;
                c.c21 = -5.6688;
                c.c31 = -2.430093356833875; c.c32 = -.2063599157091915;
                c.c41 = -.1073529058151375; c.c42 = -9.594562251023355;
                c.c43 = -20.47028614809616;
                c.c51 = 7.496443313967647; c.c52 = -10.24680431464352;
                c.c53 = -33.99990352819905; c.c54 = 11.7089089320616;
                c.c61 = 8.083246795921522; c.c62 = -7.981132988064893;
                c.c63 = -31.52159432874371; c.c64 = 16.31930543123136;
                c.c65 = -6.058818238834054;
                c.gamma = .25;
                c.d21 = 10.12623508344586; c.d22 = -7.487995877610167;
                c.d23 = -34.80091861555747; c.d24 = -7.992771707568823;
                c.d25 = 1.025137723295662;
                c.d31 = -.6762803392801253; c.d32 = 6.087714651680015;
                c.d33 = 16.43084320892478; c.d34 = 24.76722511418386;
                c.d35 = -6.594389125716872;
                break;
            case 2:
                c.c2 = .3507221; c.c3 = .2557041; c.c4 = .681779;
                c.d1 = .25; c.d2 = -.06902209999999998;
                c.d3 = -9.671999999999459e-4; c.d4 = -.08797900000000025;
                c.a21 = 1.4028884;
                c.a31 = .6581212688557198; c.a32 = -1.320936088384301;
                c.a41 = 7.131197445744498; c.a42 = 16.02964143958207;
                c.a43 = -5.561572550509766;
                c.a51 = 22.73885722420363; c.a52 = 67.38147284535289;
                c.a53 = -31.2187749303856; c.a54 = .7285641833203814;
                c.c21 = -5.1043536;
                c.c31 = -2.899967805418783; c.c32 = 4.040399359702244;
                c.c41 = -32.64449927841361; c.c42 = -99.35311008728094;
                c.c43 = 49.99119122405989;
                c.c51 = -76.46023087151691; c.c52 = -278.5942120829058;
                c.c53 = 153.9294840910643; c.c54 = 10.97101866258358;
                c.c61 = -76.29701586804983; c.c62 = -294.2795630511232;
                c.c63 = 162.0029695867566; c.c64 = 23.6516690309527;
                c.c65 = -7.652977706771382;
                c.gamma = .25;
                c.d21 = -38.71940424117216; c.d22 = -135.8025833007622;
                c.d23 = 64.51068857505875; c.d24 = -4.192663174613162;
                c.d25 = -2.53193205033506;
                c.d31 = -14.99268484949843; c.d32 = -76.30242396627033;
                c.d33 = 58.65928432851416; c.d34 = 16.61359034616402;
                c.d35 = -.6758691794084156;
                break;
            case 3:
                // Coefficients for RODAS with order 4 for linear parabolic
                // problems, Gerd Steinebach (1993)
                c.gamma = .25;
                c.c2 = c.gamma * 3.; c.c3 = .21; c.c4 = .63;
                c.d1 = .25; c.d2 = -.5; c.d3 = -.023504; c.d4 = -.0362;
                c.a21 = 3.;
                c.a31 = 1.831036793486759; c.a32 = .4955183967433795;
                c.a41 = 2.304376582692669; c.a42 = -.05249275245743001;
                c.a43 = -1.176798761832782;
                c.a51 = -7.170454962423024; c.a52 = -4.741636671481785;
                c.a53 = -16.31002631330971; c.a54 = -1.062004044111401;
                c.c21 = -12.;
                c.c31 = -8.791795173947035; c.c32 = -2.207865586973518;
                c.c41 = 10.81793056857153; c.c42 = 6.780270611428266;
                c.c43 = 19.5348594464241;
                c.c51 = 34.19095006749676; c.c52 = 15.49671153725963;
                c.c53 = 54.7476087596413; c.c54 = 14.16005392148534;
                c.c61 = 34.62605830930532; c.c62 = 15.30084976114473;
                c.c63 = 56.99955578662667; c.c64 = 18.40807009793095;
                c.c65 = -5.714285714285717;
                c.d21 = 25.09876703708589; c.d22 = 11.62013104361867;
                c.d23 = 28.49148307714626; c.d24 = -5.664021568594133;
                c.d25 = 0.;
                c.d31 = 1.638054557396973; c.d32 = -.7373619806678748;
                c.d33 = 8.47791821923899; c.d34 = 15.9925314877952;
                c.d35 = -1.882352941176471;
                break;
            default:
                throw new RodasInputException("rodas: CURIOUS INPUT IWORK(2)=" + meth);
        }
        return c;
    }

    // -----------------------------------------------------------------------
    // decsol.f: the real factorisations RODAS uses. DECC/SOLC/DECHC/SOLBC
    // belong to RADAU5's complex pair and DECH/SOLH to its Hessenberg option;
    // RODAS reaches none of them, so they are absent rather than dead code.
    // -----------------------------------------------------------------------

    /**
     * DEC: LU by Gaussian elimination with partial pivoting, Moler's ACM 423.
     *
     * <p>Not a library factorisation, and that is the point: the factors must
     * come out in the same order as the Fortran's, because RODAS reads the
     * residual of these solves and a different pivot sequence moves the
     * trajectory in the last digits.
     *
     * @return 0 if nonsingular, else the stage at which the pivot vanished
     */
    static int dec(int n, double[][] a, int[] ip) {
        ip[n] = 1;
        if (n != 1) {
            int nm1 = n - 1;
            for (int k = 1; k <= nm1; k++) {
                int kp1 = k + 1;
                int m = k;
                for (int i = kp1; i <= n; i++) {
                    if (Math.abs(a[i][k]) > Math.abs(a[m][k])) {
                        m = i;
                    }
                }
                ip[k] = m;
                double t = a[m][k];
                if (m != k) {
                    ip[n] = -ip[n];
                    a[m][k] = a[k][k];
                    a[k][k] = t;
                }
                if (t == 0.0) {
                    ip[n] = 0;
                    return k;
                }
                t = 1.0 / t;
                for (int i = kp1; i <= n; i++) {
                    a[i][k] = -a[i][k] * t;
                }
                for (int j = kp1; j <= n; j++) {
                    t = a[m][j];
                    a[m][j] = a[k][j];
                    a[k][j] = t;
                    if (t != 0.0) {
                        for (int i = kp1; i <= n; i++) {
                            a[i][j] += a[i][k] * t;
                        }
                    }
                }
            }
        }
        if (a[n][n] == 0.0) {
            ip[n] = 0;
            return n;
        }
        return 0;
    }

    /** SOL: forward/back substitution against the factors DEC left in {@code a}. */
    static void sol(int n, double[][] a, double[] b, int[] ip) {
        if (n != 1) {
            int nm1 = n - 1;
            for (int k = 1; k <= nm1; k++) {
                int kp1 = k + 1;
                int m = ip[k];
                double t = b[m];
                b[m] = b[k];
                b[k] = t;
                for (int i = kp1; i <= n; i++) {
                    b[i] += a[i][k] * t;
                }
            }
            for (int kb = 1; kb <= nm1; kb++) {
                int km1 = n - kb;
                int k = km1 + 1;
                b[k] /= a[k][k];
                double t = -b[k];
                for (int i = 1; i <= km1; i++) {
                    b[i] += a[i][k] * t;
                }
            }
        }
        b[1] /= a[1][1];
    }

    /** DECB: the banded counterpart of DEC, LINPACK band storage. */
    static int decb(int n, double[][] a, int ml, int mu, int[] ip) {
        ip[n] = 1;
        int md = ml + mu + 1;
        int md1 = md + 1;
        int ju = 0;
        if (ml != 0 && n != 1) {
            if (n >= mu + 2) {
                for (int j = mu + 2; j <= n; j++) {
                    for (int i = 1; i <= ml; i++) {
                        a[i][j] = 0.0;
                    }
                }
            }
            int nm1 = n - 1;
            for (int k = 1; k <= nm1; k++) {
                int kp1 = k + 1;
                int m = md;
                int mdl = Math.min(ml, n - k) + md;
                for (int i = md1; i <= mdl; i++) {
                    if (Math.abs(a[i][k]) > Math.abs(a[m][k])) {
                        m = i;
                    }
                }
                ip[k] = m + k - md;
                double t = a[m][k];
                if (m != md) {
                    ip[n] = -ip[n];
                    a[m][k] = a[md][k];
                    a[md][k] = t;
                }
                if (t == 0.0) {
                    ip[n] = 0;
                    return k;
                }
                t = 1.0 / t;
                for (int i = md1; i <= mdl; i++) {
                    a[i][k] = -a[i][k] * t;
                }
                ju = Math.min(Math.max(ju, mu + ip[k]), n);
                int mm = md;
                if (ju >= kp1) {
                    for (int j = kp1; j <= ju; j++) {
                        m--;
                        mm--;
                        t = a[m][j];
                        if (m != mm) {
                            a[m][j] = a[mm][j];
                            a[mm][j] = t;
                        }
                        if (t != 0.0) {
                            int jk = j - k;
                            for (int i = md1; i <= mdl; i++) {
                                a[i - jk][j] += a[i][k] * t;
                            }
                        }
                    }
                }
            }
        }
        if (a[md][n] == 0.0) {
            ip[n] = 0;
            return n;
        }
        return 0;
    }

    /** SOLB: substitution against the band factors of DECB. */
    static void solb(int n, double[][] a, int ml, int mu, double[] b, int[] ip) {
        int md = ml + mu + 1;
        int md1 = md + 1;
        int mdm = md - 1;
        int nm1 = n - 1;
        if (ml != 0) {
            if (n == 1) {
                b[1] /= a[md][1];
                return;
            }
            for (int k = 1; k <= nm1; k++) {
                int m = ip[k];
                double t = b[m];
                b[m] = b[k];
                b[k] = t;
                int mdl = Math.min(ml, n - k) + md;
                for (int i = md1; i <= mdl; i++) {
                    b[i + k - md] += a[i][k] * t;
                }
            }
        }
        for (int kb = 1; kb <= nm1; kb++) {
            int k = n + 1 - kb;
            b[k] /= a[md][k];
            double t = -b[k];
            int kmd = md - k;
            int lm = Math.max(1, kmd + 1);
            for (int i = lm; i <= mdm; i++) {
                b[i - kmd] += a[i][k] * t;
            }
        }
        b[1] /= a[md][1];
    }

    // -----------------------------------------------------------------------
    // dc_decsol.f: build and factor E = fac1*M - J, and solve against it
    // -----------------------------------------------------------------------

    /** DECOMR: assemble E1 for this IJOB and factor it. */
    static int decomr(int n, double[][] fjac, double[][] fmas, int mlmas, int mumas,
                      int m1, int m2, int nm1, double fac1, double[][] e1, int[] ip1,
                      int ijob, Linal lin) {
        int nn;
        int jm1;
        switch (ijob) {
            case 1:
                // ---  B=IDENTITY, JACOBIAN A FULL MATRIX
                for (int j = 1; j <= n; j++) {
                    for (int i = 1; i <= n; i++) {
                        e1[i][j] = -fjac[i][j];
                    }
                    e1[j][j] += fac1;
                }
                return dec(n, e1, ip1);
            case 11:
                // ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
                for (int j = 1; j <= nm1; j++) {
                    jm1 = j + m1;
                    for (int i = 1; i <= nm1; i++) {
                        e1[i][j] = -fjac[i][jm1];
                    }
                    e1[j][j] += fac1;
                }
                return decomrL45(fjac, m1, m2, nm1, fac1, e1, ip1);
            case 2:
                // ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
                for (int j = 1; j <= n; j++) {
                    for (int i = 1; i <= lin.mbjac; i++) {
                        e1[i + lin.mle][j] = -fjac[i][j];
                    }
                    e1[lin.mdiag][j] += fac1;
                }
                return decb(n, e1, lin.mle, lin.mue, ip1);
            case 12:
                // ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
                for (int j = 1; j <= nm1; j++) {
                    jm1 = j + m1;
                    for (int i = 1; i <= lin.mbjac; i++) {
                        e1[i + lin.mle][j] = -fjac[i][jm1];
                    }
                    e1[lin.mdiag][j] += fac1;
                }
                return decomrL46(fjac, m1, m2, nm1, fac1, e1, ip1, lin);
            case 3:
            case 13:
                // ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
                nn = (ijob == 3) ? n : nm1;
                for (int j = 1; j <= nn; j++) {
                    jm1 = (ijob == 3) ? j : j + m1;
                    for (int i = 1; i <= nn; i++) {
                        e1[i][j] = -fjac[i][jm1];
                    }
                    int lo = Math.max(1, j - mumas);
                    int hi = Math.min(nn, j + mlmas);
                    for (int i = lo; i <= hi; i++) {
                        e1[i][j] += fac1 * fmas[i - j + lin.mbdiag][j];
                    }
                }
                if (ijob == 3) {
                    return dec(n, e1, ip1);
                }
                return decomrL45(fjac, m1, m2, nm1, fac1, e1, ip1);
            case 4:
            case 14:
                // ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
                nn = (ijob == 4) ? n : nm1;
                for (int j = 1; j <= nn; j++) {
                    jm1 = (ijob == 4) ? j : j + m1;
                    for (int i = 1; i <= lin.mbjac; i++) {
                        e1[i + lin.mle][j] = -fjac[i][jm1];
                    }
                    for (int i = 1; i <= lin.mbb; i++) {
                        int ib = i + lin.mdiff;
                        e1[ib][j] += fac1 * fmas[i][j];
                    }
                }
                if (ijob == 4) {
                    return decb(n, e1, lin.mle, lin.mue, ip1);
                }
                return decomrL46(fjac, m1, m2, nm1, fac1, e1, ip1, lin);
            case 5:
            case 15:
                // ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
                nn = (ijob == 5) ? n : nm1;
                for (int j = 1; j <= nn; j++) {
                    jm1 = (ijob == 5) ? j : j + m1;
                    for (int i = 1; i <= nn; i++) {
                        e1[i][j] = fmas[i][j] * fac1 - fjac[i][jm1];
                    }
                }
                if (ijob == 5) {
                    return dec(n, e1, ip1);
                }
                return decomrL45(fjac, m1, m2, nm1, fac1, e1, ip1);
            default:
                // 6 is "THIS OPTION IS NOT PROVIDED" upstream; 7..10 are RADAU5's.
                return 0;
        }
    }

    /** DECOMR label 45: fold the second-order block into E1, then factor full. */
    private static int decomrL45(double[][] fjac, int m1, int m2, int nm1, double fac1,
                                 double[][] e1, int[] ip1) {
        int mm = m1 / m2;
        for (int j = 1; j <= m2; j++) {
            for (int i = 1; i <= nm1; i++) {
                double sum = 0.0;
                for (int k = 0; k <= mm - 1; k++) {
                    sum = (sum + fjac[i][j + k * m2]) / fac1;
                }
                e1[i][j] -= sum;
            }
        }
        return dec(nm1, e1, ip1);
    }

    /** DECOMR label 46: the same fold, banded. */
    private static int decomrL46(double[][] fjac, int m1, int m2, int nm1, double fac1,
                                 double[][] e1, int[] ip1, Linal lin) {
        int mm = m1 / m2;
        for (int j = 1; j <= m2; j++) {
            for (int i = 1; i <= lin.mbjac; i++) {
                double sum = 0.0;
                for (int k = 0; k <= mm - 1; k++) {
                    sum = (sum + fjac[i][j + k * m2]) / fac1;
                }
                e1[i + lin.mle][j] -= sum;
            }
        }
        return decb(nm1, e1, lin.mle, lin.mue, ip1);
    }

    /** SLVROD: one Rosenbrock stage. {@code ak} receives the solve. */
    static void slvrod(int n, double[][] fjac, int mljac, int mujac, double[][] fmas,
                       int mlmas, int mumas, int m1, int m2, int nm1, double fac1,
                       double[][] e, int[] ip, double[] dy, double[] ak, double[] fx,
                       double[] ynew, double hd, int ijob, boolean stage1, Linal lin) {
        if (hd == 0.0) {
            for (int i = 1; i <= n; i++) {
                ak[i] = dy[i];
            }
        } else {
            for (int i = 1; i <= n; i++) {
                ak[i] = dy[i] + hd * fx[i];
            }
        }

        switch (ijob) {
            case 1:
                // ---  B=IDENTITY, JACOBIAN A FULL MATRIX
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        ak[i] += ynew[i];
                    }
                }
                sol(n, e, ak, ip);
                return;
            case 11:
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        ak[i] += ynew[i];
                    }
                }
                slvrodL48(fjac, m1, m2, nm1, fac1, e, ip, ak);
                return;
            case 2:
                // ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        ak[i] += ynew[i];
                    }
                }
                solb(n, e, lin.mle, lin.mue, ak, ip);
                return;
            case 12:
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        ak[i] += ynew[i];
                    }
                }
                slvrodL45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin);
                return;
            case 3:
                // ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        double sum = 0.0;
                        int lo = Math.max(1, i - mlmas);
                        int hi = Math.min(n, i + mumas);
                        for (int j = lo; j <= hi; j++) {
                            sum += fmas[i - j + lin.mbdiag][j] * ynew[j];
                        }
                        ak[i] += sum;
                    }
                }
                sol(n, e, ak, ip);
                return;
            case 13:
            case 14:
                // ---  B BANDED, JACOBIAN FULL, SECOND ORDER
                if (stage1) {
                    for (int i = 1; i <= m1; i++) {
                        ak[i] += ynew[i];
                    }
                    for (int i = 1; i <= nm1; i++) {
                        double sum = 0.0;
                        int lo = Math.max(1, i - mlmas);
                        int hi = Math.min(nm1, i + mumas);
                        for (int j = lo; j <= hi; j++) {
                            sum += fmas[i - j + lin.mbdiag][j] * ynew[j + m1];
                        }
                        ak[i + m1] += sum;
                    }
                }
                if (ijob == 14) {
                    slvrodL45(fjac, mljac, mujac, m1, m2, nm1, fac1, e, ip, ak, lin);
                } else {
                    slvrodL48(fjac, m1, m2, nm1, fac1, e, ip, ak);
                }
                return;
            case 4:
                // ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        double sum = 0.0;
                        int lo = Math.max(1, i - mlmas);
                        int hi = Math.min(n, i + mumas);
                        for (int j = lo; j <= hi; j++) {
                            sum += fmas[i - j + lin.mbdiag][j] * ynew[j];
                        }
                        ak[i] += sum;
                    }
                }
                solb(n, e, lin.mle, lin.mue, ak, ip);
                return;
            case 5:
                // ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        double sum = 0.0;
                        for (int j = 1; j <= n; j++) {
                            sum += fmas[i][j] * ynew[j];
                        }
                        ak[i] += sum;
                    }
                }
                sol(n, e, ak, ip);
                return;
            case 15:
                if (stage1) {
                    for (int i = 1; i <= m1; i++) {
                        ak[i] += ynew[i];
                    }
                    for (int i = 1; i <= nm1; i++) {
                        double sum = 0.0;
                        for (int j = 1; j <= nm1; j++) {
                            sum += fmas[i][j] * ynew[j + m1];
                        }
                        ak[i + m1] += sum;
                    }
                }
                slvrodL48(fjac, m1, m2, nm1, fac1, e, ip, ak);
                return;
            case 6:
                // ---  THIS OPTION IS NOT PROVIDED upstream, and it solves only
                //      under stage1 there. Kept identical rather than tidied.
                if (stage1) {
                    for (int i = 1; i <= n; i++) {
                        double sum = 0.0;
                        for (int j = 1; j <= n; j++) {
                            sum += fmas[i][j] * ynew[j];
                        }
                        ak[i] += sum;
                    }
                    solb(n, e, lin.mle, lin.mue, ak, ip);
                }
                return;
            default:
                // 7..10 belong to RADAU5.
        }
    }

    /** SLVROD label 48: the second-order elimination, full Jacobian. */
    private static void slvrodL48(double[][] fjac, int m1, int m2, int nm1, double fac1,
                                  double[][] e, int[] ip, double[] ak) {
        int mm = m1 / m2;
        for (int j = 1; j <= m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * m2;
                sum = (ak[jkm] + sum) / fac1;
                for (int i = 1; i <= nm1; i++) {
                    ak[i + m1] += fjac[i][jkm] * sum;
                }
            }
        }
        double[] sub = new double[nm1 + 1];
        System.arraycopy(ak, m1 + 1, sub, 1, nm1);
        sol(nm1, e, sub, ip);
        System.arraycopy(sub, 1, ak, m1 + 1, nm1);
        for (int i = m1; i >= 1; i--) {
            ak[i] = (ak[i] + ak[m2 + i]) / fac1;
        }
    }

    /** SLVROD label 45: the second-order elimination, banded Jacobian. */
    private static void slvrodL45(double[][] fjac, int mljac, int mujac, int m1, int m2,
                                  int nm1, double fac1, double[][] e, int[] ip,
                                  double[] ak, Linal lin) {
        int mm = m1 / m2;
        for (int j = 1; j <= m2; j++) {
            double sum = 0.0;
            for (int k = mm - 1; k >= 0; k--) {
                int jkm = j + k * m2;
                sum = (ak[jkm] + sum) / fac1;
                int lo = Math.max(1, j - mujac);
                int hi = Math.min(nm1, j + mljac);
                for (int i = lo; i <= hi; i++) {
                    ak[i + m1] += fjac[i + mujac + 1 - j][jkm] * sum;
                }
            }
        }
        double[] sub = new double[nm1 + 1];
        System.arraycopy(ak, m1 + 1, sub, 1, nm1);
        solb(nm1, e, lin.mle, lin.mue, sub, ip);
        System.arraycopy(sub, 1, ak, m1 + 1, nm1);
        for (int i = m1; i >= 1; i--) {
            ak[i] = (ak[i] + ak[m2 + i]) / fac1;
        }
    }

    // -----------------------------------------------------------------------
    // rodas.f: the driver and the core integrator
    // -----------------------------------------------------------------------

    /**
     * Integrate {@code M y' = f(x,y)} from {@code x} to {@code xend}.
     *
     * @param n     dimension of the system
     * @param fcn   the right hand side
     * @param x     initial abscissa
     * @param y     initial state, length n, 0-based; modified in place
     * @param xend  final abscissa
     * @param h     initial step size guess (0 means 1e-6)
     * @param rtol  relative tolerance, length 1 (itol 0) or n (itol 1)
     * @param atol  absolute tolerance, same shape as rtol
     * @param itol  0 scalar tolerances, 1 one per equation
     * @param opt   the Fortran switches; null for all defaults
     * @return where the integration ended and what it cost
     */
    public static Result integrate(int n, Fcn fcn, double x, double[] y, double xend,
                                   double h, double[] rtol, double[] atol, int itol,
                                   Options opt) {
        if (opt == null) {
            opt = new Options();
        }
        int mljac = (opt.mljac < 0) ? n : opt.mljac;
        int mujac = opt.mujac;
        int mlmas = opt.mlmas;
        int mumas = opt.mumas;
        int m1 = opt.m1;
        int m2 = opt.m2;
        int nm1 = n - m1;
        if (m1 == 0) {
            m2 = n;
        }
        if (m2 == 0) {
            m2 = m1;
        }
        if (opt.nmax <= 0) {
            throw new RodasInputException("rodas: WRONG INPUT IWORK(1)=" + opt.nmax);
        }
        if (opt.meth <= 0 || opt.meth >= 4) {
            throw new RodasInputException("rodas: CURIOUS INPUT IWORK(2)=" + opt.meth);
        }
        if (m1 < 0 || m2 < 0 || m1 + m2 > n) {
            throw new RodasInputException(
                    "rodas: CURIOUS INPUT FOR IWORK(9,10)=" + m1 + " " + m2);
        }
        double uround = opt.uround;
        if (uround < 1e-16 || uround >= 1.0) {
            throw new RodasInputException(
                    "rodas: COEFFICIENTS HAVE 16 DIGITS, UROUND=" + uround);
        }
        double hmax = Double.isNaN(opt.hmax) ? (xend - x) : opt.hmax;
        double fac1 = opt.fac1;
        double fac2 = opt.fac2;
        if (fac1 < 1.0 || fac2 > 1.0) {
            throw new RodasInputException("rodas: CURIOUS INPUT WORK(3,4)");
        }
        if (opt.safe <= .001 || opt.safe >= 1.0) {
            throw new RodasInputException("rodas: CURIOUS INPUT FOR WORK(5)=" + opt.safe);
        }

        double[] rtolv = new double[n + 1];
        double[] atolv = new double[n + 1];
        if (itol == 0) {
            if (atol[0] <= 0.0 || rtol[0] <= uround * 10.0) {
                throw new RodasInputException("rodas: TOLERANCES ARE TOO SMALL");
            }
            for (int i = 1; i <= n; i++) {
                rtolv[i] = rtol[0];
                atolv[i] = atol[0];
            }
        } else {
            for (int i = 1; i <= n; i++) {
                if (atol[i - 1] <= 0.0 || rtol[i - 1] <= uround * 10.0) {
                    throw new RodasInputException(
                            "rodas: TOLERANCES(" + i + ") ARE TOO SMALL");
                }
                rtolv[i] = rtol[i - 1];
                atolv[i] = atol[i - 1];
            }
        }

        boolean autnms = (opt.ifcn == 0);
        boolean implct = (opt.imas != 0);
        boolean jband = (mljac < nm1);

        int ldjac;
        int lde;
        if (jband) {
            ldjac = mljac + mujac + 1;
            lde = mljac + ldjac;
        } else {
            mljac = nm1;
            mujac = nm1;
            ldjac = nm1;
            lde = nm1;
        }
        int ldmas;
        int ijob;
        if (implct) {
            if (mlmas != nm1) {
                ldmas = mlmas + mumas + 1;
                ijob = jband ? 4 : 3;
            } else {
                ldmas = nm1;
                ijob = 5;
            }
            if (mlmas > mljac || mumas > mujac) {
                throw new RodasInputException(
                        "rodas: BANDWITH OF \"MAS\" NOT LARGER THAN BANDWITH OF \"JAC\"");
            }
        } else {
            ldmas = 0;
            ijob = jband ? 2 : 1;
        }
        ldmas = Math.max(1, ldmas);

        return roscor(n, fcn, x, y, xend, hmax, h, rtolv, atolv, itol, opt, mljac, mujac,
                mlmas, mumas, uround, ijob, fac1, fac2, autnms, implct, jband,
                ldjac, lde, ldmas, m1, m2, nm1);
    }

    private static Result roscor(int n, Fcn fcn, double x, double[] y0, double xend,
                                 double hmax, double h, double[] rtol, double[] atol,
                                 int itol, Options opt, int mljac, int mujac, int mlmas,
                                 int mumas, double uround, int ijob, double fac1,
                                 double fac2, boolean autnms, boolean implct,
                                 boolean jband, int ldjac, int lde, int ldmas,
                                 int m1, int m2, int nm1) {
        double[] y = new double[n + 1];
        System.arraycopy(y0, 0, y, 1, n);

        double[] ynew = new double[n + 1];
        double[] dy1 = new double[n + 1];
        double[] dy = new double[n + 1];
        double[] ak1 = new double[n + 1];
        double[] ak2 = new double[n + 1];
        double[] ak3 = new double[n + 1];
        double[] ak4 = new double[n + 1];
        double[] ak5 = new double[n + 1];
        double[] ak6 = new double[n + 1];
        double[] fx = new double[n + 1];
        double[] cont = new double[4 * n + 1];
        double[][] fjac = new double[ldjac + 1][n + 1];
        double[][] e = new double[lde + 1][nm1 + 1];
        double[][] fmas = new double[ldmas + 1][nm1 + 1];
        int[] ip = new int[nm1 + 1];

        // 0-based scratch for the callback boundary
        double[] yb = new double[n];
        double[] fb = new double[n];

        Linal lin = new Linal();
        Dense dense = new Dense(cont, n);
        int nn2 = 2 * n;
        int nn3 = 3 * n;

        if (implct) {
            double[][] amb = new double[ldmas][nm1];
            opt.mas.eval(amb);
            for (int i = 1; i <= ldmas; i++) {
                for (int j = 1; j <= nm1; j++) {
                    fmas[i][j] = amb[i - 1][j - 1];
                }
            }
        }
        Coef c = rocoe(opt.meth);
        if (m1 > 0) {
            ijob += 10;
        }

        double posneg = (xend - x >= 0) ? 1.0 : -1.0;
        double hmaxn = Math.min(Math.abs(hmax), Math.abs(xend - x));
        if (Math.abs(h) <= uround * 10.0) {
            h = 1e-6;
        }
        h = Math.min(Math.abs(h), hmaxn);
        h = (posneg >= 0) ? h : -h;
        boolean reject = false;
        boolean last = false;
        int nsing = 0;
        double hd1 = 0;
        double hd2 = 0;
        double hd3 = 0;
        double hd4 = 0;
        double hacc = 0;
        double erracc = 0;
        double hopt = h;

        lin.mbdiag = mumas + 1;
        if (jband) {
            lin.mle = mljac;
            lin.mue = mujac;
            lin.mbjac = mljac + mujac + 1;
            lin.mbb = mlmas + mumas + 1;
            lin.mdiag = lin.mle + lin.mue + 1;
            lin.mdiff = lin.mle + lin.mue - mumas;
        }

        Result res = new Result();
        res.idid = 0;

        if (opt.iout != 0) {
            dense.xold = x;
            dense.h = h;
            copyOut(y, yb, n);
            int irtrn = opt.solout.eval(res.naccpt + 1, dense.xold, x, yb, dense);
            if (irtrn < 0) {
                return finish(res, y, n, x, h, IDID_SOLOUT_STOP, y0);
            }
        }

        // ==== BASIC INTEGRATION STEP (label 1) ====
        while (true) {
            if (res.nstep > opt.nmax) {
                return finish(res, y, n, x, h, IDID_NMAX, y0);
            }
            if (Math.abs(h) * .1 <= Math.abs(x) * uround) {
                return finish(res, y, n, x, h, IDID_STEP_TOO_SMALL, y0);
            }
            if (last) {
                h = hopt;
                return finish(res, y, n, x, h, IDID_SUCCESS, y0);
            }
            hopt = h;
            if ((x + h * 1.0001 - xend) * posneg >= 0.0) {
                h = xend - x;
                last = true;
            }

            // ---- COMPUTATION OF THE JACOBIAN ----
            copyOut(y, yb, n);
            opt_eval(fcn, x, yb, fb);
            copyIn(fb, dy1, n);
            res.nfcn++;
            res.njac++;
            if (opt.ijac == 0) {
                if (jband) {
                    int mujacp = mujac + 1;
                    int md = Math.min(lin.mbjac, n);
                    for (int mm = 1; mm <= m1 / m2 + 1; mm++) {
                        for (int k = 1; k <= md; k++) {
                            int j = k + (mm - 1) * m2;
                            while (true) {
                                ak2[j] = y[j];
                                ak3[j] = Math.sqrt(uround
                                        * Math.max(1e-5, Math.abs(y[j])));
                                y[j] += ak3[j];
                                j += md;
                                if (j > mm * m2) {
                                    break;
                                }
                            }
                            copyOut(y, yb, n);
                            opt_eval(fcn, x, yb, fb);
                            copyIn(fb, ak1, n);
                            j = k + (mm - 1) * m2;
                            int j1 = k;
                            int lbeg = Math.max(1, j1 - mujac) + m1;
                            while (true) {
                                int lend = Math.min(m2, j1 + mljac) + m1;
                                y[j] = ak2[j];
                                int mujacj = mujacp - j1 - m1;
                                for (int l = lbeg; l <= lend; l++) {
                                    fjac[l + mujacj][j] = (ak1[l] - dy1[l]) / ak3[j];
                                }
                                j += md;
                                j1 += md;
                                lbeg = lend + 1;
                                if (j > mm * m2) {
                                    break;
                                }
                            }
                        }
                    }
                } else {
                    for (int i = 1; i <= n; i++) {
                        double ysafe = y[i];
                        double delt = Math.sqrt(uround * Math.max(1e-5, Math.abs(ysafe)));
                        y[i] = ysafe + delt;
                        copyOut(y, yb, n);
                        opt_eval(fcn, x, yb, fb);
                        copyIn(fb, ak1, n);
                        for (int j = m1 + 1; j <= n; j++) {
                            fjac[j - m1][i] = (ak1[j] - dy1[j]) / delt;
                        }
                        y[i] = ysafe;
                    }
                }
            } else {
                double[][] jb = new double[ldjac][n];
                copyOut(y, yb, n);
                opt.jac.eval(x, yb, jb);
                for (int i = 1; i <= ldjac; i++) {
                    for (int j = 1; j <= n; j++) {
                        fjac[i][j] = jb[i - 1][j - 1];
                    }
                }
            }
            if (!autnms) {
                if (opt.idfx == 0) {
                    double delt = Math.sqrt(uround * Math.max(1e-5, Math.abs(x)));
                    double xdelt = x + delt;
                    copyOut(y, yb, n);
                    opt_eval(fcn, xdelt, yb, fb);
                    copyIn(fb, ak1, n);
                    for (int j = 1; j <= n; j++) {
                        fx[j] = (ak1[j] - dy1[j]) / delt;
                    }
                } else {
                    copyOut(y, yb, n);
                    opt.dfx.eval(x, yb, fb);
                    copyIn(fb, fx, n);
                }
            }

            // ==== COMPUTE THE STAGES (label 2) ====
            boolean stepdone = false;
            while (!stepdone) {
                double fac = 1.0 / (h * c.gamma);
                int ier = decomr(n, fjac, fmas, mlmas, mumas, m1, m2, nm1, fac, e, ip,
                        ijob, lin);
                if (ier != 0) {
                    // ---- SINGULAR MATRIX (label 80) ----
                    nsing++;
                    if (nsing >= 5) {
                        return finish(res, y, n, x, h, IDID_SINGULAR, y0);
                    }
                    h *= .5;
                    reject = true;
                    last = false;
                    continue;
                }
                res.ndec++;

                double hc21 = c.c21 / h;
                double hc31 = c.c31 / h;
                double hc32 = c.c32 / h;
                double hc41 = c.c41 / h;
                double hc42 = c.c42 / h;
                double hc43 = c.c43 / h;
                double hc51 = c.c51 / h;
                double hc52 = c.c52 / h;
                double hc53 = c.c53 / h;
                double hc54 = c.c54 / h;
                double hc61 = c.c61 / h;
                double hc62 = c.c62 / h;
                double hc63 = c.c63 / h;
                double hc64 = c.c64 / h;
                double hc65 = c.c65 / h;
                if (!autnms) {
                    hd1 = h * c.d1;
                    hd2 = h * c.d2;
                    hd3 = h * c.d3;
                    hd4 = h * c.d4;
                }

                // ---- THE STAGES ----
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy1, ak1, fx, ynew, hd1, ijob, false, lin);
                for (int i = 1; i <= n; i++) {
                    ynew[i] = y[i] + c.a21 * ak1[i];
                }
                copyOut(ynew, yb, n);
                opt_eval(fcn, x + c.c2 * h, yb, fb);
                copyIn(fb, dy, n);
                for (int i = 1; i <= n; i++) {
                    ynew[i] = hc21 * ak1[i];
                }
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy, ak2, fx, ynew, hd2, ijob, true, lin);

                for (int i = 1; i <= n; i++) {
                    ynew[i] = y[i] + c.a31 * ak1[i] + c.a32 * ak2[i];
                }
                copyOut(ynew, yb, n);
                opt_eval(fcn, x + c.c3 * h, yb, fb);
                copyIn(fb, dy, n);
                for (int i = 1; i <= n; i++) {
                    ynew[i] = hc31 * ak1[i] + hc32 * ak2[i];
                }
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy, ak3, fx, ynew, hd3, ijob, true, lin);

                for (int i = 1; i <= n; i++) {
                    ynew[i] = y[i] + c.a41 * ak1[i] + c.a42 * ak2[i] + c.a43 * ak3[i];
                }
                copyOut(ynew, yb, n);
                opt_eval(fcn, x + c.c4 * h, yb, fb);
                copyIn(fb, dy, n);
                for (int i = 1; i <= n; i++) {
                    ynew[i] = hc41 * ak1[i] + hc42 * ak2[i] + hc43 * ak3[i];
                }
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy, ak4, fx, ynew, hd4, ijob, true, lin);

                for (int i = 1; i <= n; i++) {
                    ynew[i] = y[i] + c.a51 * ak1[i] + c.a52 * ak2[i] + c.a53 * ak3[i]
                            + c.a54 * ak4[i];
                }
                copyOut(ynew, yb, n);
                opt_eval(fcn, x + h, yb, fb);
                copyIn(fb, dy, n);
                for (int i = 1; i <= n; i++) {
                    ak6[i] = hc52 * ak2[i] + hc54 * ak4[i] + hc51 * ak1[i]
                            + hc53 * ak3[i];
                }
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy, ak5, fx, ak6, 0.0, ijob, true, lin);

                // ---- EMBEDDED SOLUTION ----
                for (int i = 1; i <= n; i++) {
                    ynew[i] += ak5[i];
                }
                copyOut(ynew, yb, n);
                opt_eval(fcn, x + h, yb, fb);
                copyIn(fb, dy, n);
                for (int i = 1; i <= n; i++) {
                    cont[i] = hc61 * ak1[i] + hc62 * ak2[i] + hc65 * ak5[i]
                            + hc64 * ak4[i] + hc63 * ak3[i];
                }
                slvrod(n, fjac, mljac, mujac, fmas, mlmas, mumas, m1, m2, nm1, fac, e,
                        ip, dy, ak6, fx, cont, 0.0, ijob, true, lin);

                // ---- NEW SOLUTION ----
                for (int i = 1; i <= n; i++) {
                    ynew[i] += ak6[i];
                }
                res.nsol += 6;
                res.nfcn += 5;

                // ---- DENSE OUTPUT ----
                if (opt.iout != 0) {
                    for (int i = 1; i <= n; i++) {
                        cont[i] = y[i];
                        cont[i + nn2] = c.d21 * ak1[i] + c.d22 * ak2[i] + c.d23 * ak3[i]
                                + c.d24 * ak4[i] + c.d25 * ak5[i];
                        cont[i + nn3] = c.d31 * ak1[i] + c.d32 * ak2[i] + c.d33 * ak3[i]
                                + c.d34 * ak4[i] + c.d35 * ak5[i];
                    }
                }

                // ==== ERROR ESTIMATION ====
                res.nstep++;
                // A SEQUENTIAL SUM, as the Fortran accumulates it.
                double err = 0.0;
                for (int i = 1; i <= n; i++) {
                    double sk = atol[i] + rtol[i]
                            * Math.max(Math.abs(y[i]), Math.abs(ynew[i]));
                    double q = ak6[i] / sk;
                    err += q * q;
                }
                err = Math.sqrt(err / n);

                // ---- COMPUTATION OF HNEW, .2 <= HNEW/H <= 6 ----
                fac = Math.max(fac2, Math.min(fac1, Math.pow(err, 0.25) / opt.safe));
                double hnew = h / fac;

                if (err <= 1.0) {
                    // ---- STEP IS ACCEPTED ----
                    res.naccpt++;
                    if (opt.pred) {
                        if (res.naccpt > 1) {
                            double facgus = (hacc / h)
                                    * Math.pow(err * err / erracc, 0.25) / opt.safe;
                            facgus = Math.max(fac2, Math.min(fac1, facgus));
                            fac = Math.max(fac, facgus);
                            hnew = h / fac;
                        }
                        hacc = h;
                        erracc = Math.max(.01, err);
                    }
                    for (int i = 1; i <= n; i++) {
                        y[i] = ynew[i];
                    }
                    dense.xold = x;
                    x += h;
                    if (opt.iout != 0) {
                        for (int i = 1; i <= n; i++) {
                            cont[n + i] = y[i];
                        }
                        dense.h = h;
                        copyOut(y, yb, n);
                        int irtrn = opt.solout.eval(res.naccpt + 1, dense.xold, x, yb,
                                dense);
                        if (irtrn < 0) {
                            return finish(res, y, n, x, h, IDID_SOLOUT_STOP, y0);
                        }
                    }
                    if (Math.abs(hnew) > hmaxn) {
                        hnew = posneg * hmaxn;
                    }
                    if (reject) {
                        hnew = posneg * Math.min(Math.abs(hnew), Math.abs(h));
                    }
                    reject = false;
                    h = hnew;
                    stepdone = true;   // back to label 1
                } else {
                    // ---- STEP IS REJECTED ----
                    reject = true;
                    last = false;
                    h = hnew;
                    if (res.naccpt >= 1) {
                        res.nrejct++;
                    }
                    // back to label 2
                }
            }
        }
    }

    private static void opt_eval(Fcn fcn, double x, double[] yb, double[] fb) {
        fcn.eval(x, yb, fb);
    }

    private static void copyOut(double[] src1, double[] dst0, int n) {
        System.arraycopy(src1, 1, dst0, 0, n);
    }

    private static void copyIn(double[] src0, double[] dst1, int n) {
        System.arraycopy(src0, 0, dst1, 1, n);
    }

    private static Result finish(Result res, double[] y, int n, double x, double h,
                                 int idid, double[] y0) {
        res.x = x;
        res.h = h;
        res.idid = idid;
        res.y = new double[n];
        System.arraycopy(y, 1, res.y, 0, n);
        System.arraycopy(y, 1, y0, 0, n);
        return res;
    }
}
