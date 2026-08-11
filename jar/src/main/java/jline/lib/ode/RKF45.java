/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Runge-Kutta-Fehlberg 4(5) adaptive ODE solver.
 *
 * Ported from the C version by John Burkardt, which was itself a port of the
 * original FORTRAN77 code by Herman Watts and Lawrence Shampine (Sandia).
 */
package jline.lib.ode;

/**
 * Runge-Kutta-Fehlberg 4(5) adaptive ODE integrator.
 */
public final class RKF45 {

    private static final int MAXNFE = 3000;
    private static final double REMIN = 1.0e-12;

    private static double sign(double x) { return x < 0.0 ? -1.0 : 1.0; }

    private static double epsilon() {
        double eps = 1.0;
        while (1.0 + eps > 1.0) {
            eps /= 2.0;
        }
        return 2.0 * eps;
    }

    private final int neqn;

    private double h = -1.0;
    private int init = -1000;
    private int kflag = -1000;
    private int kop = -1;
    private int nfe = -1;
    private int flagSave = -1000;
    private double relerrSave = -1.0;
    private double abserrSave = -1.0;

    public RKF45(int neqn) {
        this.neqn = neqn;
    }

    /** Reset internal state so the next call is treated as a fresh start. */
    public void reset() {
        h = -1.0;
        init = -1000;
        kflag = -1000;
        kop = -1;
        nfe = -1;
        flagSave = -1000;
        relerrSave = -1.0;
        abserrSave = -1.0;
    }

    public RKF45Result integrate(
            OdeFunction f,
            double[] yIn,
            double[] ypIn,
            double tIn,
            double tout,
            double relerr,
            double abserr,
            int flag) {
        double[] y = yIn.clone();
        double[] yp = ypIn.clone();
        double t = tIn;
        double relErr = relerr;
        int flagReturn = flag;

        double eps = epsilon();

        if (neqn < 1 || relErr < 0.0 || abserr < 0.0
                || flagReturn == 0 || flagReturn < -2 || flagReturn > 8) {
            return new RKF45Result(y, yp, t, 8, relErr);
        }

        int mflag = Math.abs(flagReturn);

        if (mflag != 1) {
            if (t == tout && kflag != 3) {
                return new RKF45Result(y, yp, t, 8, relErr);
            }

            if (mflag == 2) {
                if (kflag == 3) {
                    flagReturn = flagSave;
                    mflag = Math.abs(flagReturn);
                } else if (init == 0) {
                    flagReturn = flagSave;
                } else if (kflag == 4) {
                    nfe = 0;
                } else if (kflag == 5 && abserr == 0.0) {
                    throw new IllegalStateException("RKF45: KFLAG = 5 and ABSERR = 0.0");
                } else if (kflag == 6 && relErr <= relerrSave && abserr <= abserrSave) {
                    throw new IllegalStateException(
                            "RKF45: KFLAG = 6 and RELERR <= RELERR_SAVE and ABSERR <= ABSERR_SAVE");
                }
            } else {
                if (flagReturn == 3) {
                    flagReturn = flagSave;
                    if (kflag == 3) mflag = Math.abs(flagReturn);
                } else if (flagReturn == 4) {
                    nfe = 0;
                    flagReturn = flagSave;
                    if (kflag == 3) mflag = Math.abs(flagReturn);
                } else if (flagReturn == 5 && abserr > 0.0) {
                    flagReturn = flagSave;
                    if (kflag == 3) mflag = Math.abs(flagReturn);
                } else {
                    throw new IllegalStateException(
                            "RKF45: Integration cannot be continued. " +
                                    "User did not respond to FLAG = " + flagReturn);
                }
            }
        }

        flagSave = flagReturn;
        kflag = 0;
        relerrSave = relErr;
        abserrSave = abserr;

        double relerrMin = 2.0 * eps + REMIN;
        if (relErr < relerrMin) {
            relErr = relerrMin;
            kflag = 3;
            return new RKF45Result(y, yp, t, 3, relErr);
        }

        double dt = tout - t;

        double[] f1 = new double[neqn];
        double[] f2 = new double[neqn];
        double[] f3 = new double[neqn];
        double[] f4 = new double[neqn];
        double[] f5 = new double[neqn];

        if (mflag == 1) {
            init = 0;
            kop = 0;
            f.evaluate(t, y, yp);
            nfe = 1;

            if (t == tout) {
                return new RKF45Result(y, yp, t, 2, relErr);
            }
        }

        if (init == 0) {
            init = 1;
            h = Math.abs(dt);
            double toln = 0.0;

            for (int k = 0; k < neqn; k++) {
                double tol = relErr * Math.abs(y[k]) + abserr;
                if (tol > 0.0) {
                    toln = tol;
                    double ypk = Math.abs(yp[k]);
                    if (tol < ypk * Math.pow(h, 5)) {
                        h = Math.pow(tol / ypk, 0.2);
                    }
                }
            }

            if (toln <= 0.0) {
                h = 0.0;
            }

            h = Math.max(h, 26.0 * eps * Math.max(Math.abs(t), Math.abs(dt)));
            flagSave = (flagReturn < 0) ? -2 : 2;
        }

        h = sign(dt) * Math.abs(h);

        if (2.0 * Math.abs(dt) <= Math.abs(h)) {
            kop++;
        }

        if (kop == 100) {
            kop = 0;
            return new RKF45Result(y, yp, t, 7, relErr);
        }

        if (Math.abs(dt) <= 26.0 * eps * Math.abs(t)) {
            t = tout;
            for (int i = 0; i < neqn; i++) {
                y[i] = y[i] + dt * yp[i];
            }
            f.evaluate(t, y, yp);
            nfe++;
            return new RKF45Result(y, yp, t, 2, relErr);
        }

        boolean output = false;
        double scale = 2.0 / relErr;
        double ae = scale * abserr;

        while (true) {
            boolean hfaild = false;
            double hmin = 26.0 * eps * Math.abs(t);

            dt = tout - t;

            if (2.0 * Math.abs(h) > Math.abs(dt)) {
                if (Math.abs(dt) <= Math.abs(h)) {
                    output = true;
                    h = dt;
                } else {
                    h = 0.5 * dt;
                }
            }

            double esttol = 0.0;

            while (true) {
                if (nfe > MAXNFE) {
                    kflag = 4;
                    return new RKF45Result(y, yp, t, 4, relErr);
                }

                fehl(f, y, t, h, yp, f1, f2, f3, f4, f5, f1);
                nfe += 5;

                double eeoet = 0.0;
                for (int k = 0; k < neqn; k++) {
                    double et = Math.abs(y[k]) + Math.abs(f1[k]) + ae;
                    if (et <= 0.0) {
                        return new RKF45Result(y, yp, t, 5, relErr);
                    }
                    double ee = Math.abs(
                            (-2090.0 * yp[k] + (21970.0 * f3[k] - 15048.0 * f4[k]))
                                    + (22528.0 * f2[k] - 27360.0 * f5[k]));
                    eeoet = Math.max(eeoet, ee / et);
                }

                esttol = Math.abs(h) * eeoet * scale / 752400.0;
                if (esttol <= 1.0) break;

                hfaild = true;
                output = false;

                double s = (esttol < 59049.0) ? 0.9 / Math.pow(esttol, 0.2) : 0.1;
                h *= s;

                if (Math.abs(h) < hmin) {
                    kflag = 6;
                    return new RKF45Result(y, yp, t, 6, relErr);
                }
            }

            t += h;
            for (int i = 0; i < neqn; i++) {
                y[i] = f1[i];
            }
            f.evaluate(t, y, yp);
            nfe++;

            double s2 = (esttol > 0.0001889568) ? 0.9 / Math.pow(esttol, 0.2) : 5.0;
            if (hfaild) s2 = Math.min(s2, 1.0);
            h = sign(h) * Math.max(s2 * Math.abs(h), hmin);

            if (output) {
                t = tout;
                return new RKF45Result(y, yp, t, 2, relErr);
            }

            if (flagReturn <= 0) {
                return new RKF45Result(y, yp, t, -2, relErr);
            }
        }
    }

    private void fehl(
            OdeFunction f,
            double[] y,
            double t,
            double h,
            double[] yp,
            double[] f1,
            double[] f2,
            double[] f3,
            double[] f4,
            double[] f5,
            double[] s) {
        double ch;

        ch = h / 4.0;
        for (int i = 0; i < neqn; i++) {
            f5[i] = y[i] + ch * yp[i];
        }
        f.evaluate(t + ch, f5, f1);

        ch = 3.0 * h / 32.0;
        for (int i = 0; i < neqn; i++) {
            f5[i] = y[i] + ch * (yp[i] + 3.0 * f1[i]);
        }
        f.evaluate(t + 3.0 * h / 8.0, f5, f2);

        ch = h / 2197.0;
        for (int i = 0; i < neqn; i++) {
            f5[i] = y[i] + ch * (1932.0 * yp[i] + (7296.0 * f2[i] - 7200.0 * f1[i]));
        }
        f.evaluate(t + 12.0 * h / 13.0, f5, f3);

        ch = h / 4104.0;
        for (int i = 0; i < neqn; i++) {
            f5[i] = y[i] + ch * (
                    (8341.0 * yp[i] - 845.0 * f3[i])
                            + (29440.0 * f2[i] - 32832.0 * f1[i]));
        }
        f.evaluate(t + h, f5, f4);

        ch = h / 20520.0;
        for (int i = 0; i < neqn; i++) {
            f1[i] = y[i] + ch * (
                    (-6080.0 * yp[i] + (9295.0 * f3[i] - 5643.0 * f4[i]))
                            + (41040.0 * f1[i] - 28352.0 * f2[i]));
        }
        f.evaluate(t + h / 2.0, f1, f5);

        ch = h / 7618050.0;
        for (int i = 0; i < neqn; i++) {
            s[i] = y[i] + ch * (
                    (902880.0 * yp[i] + (3855735.0 * f3[i] - 1371249.0 * f4[i]))
                            + (3953664.0 * f2[i] + 277020.0 * f5[i]));
        }
    }
}
