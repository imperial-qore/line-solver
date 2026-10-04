package jline.api.qsys;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.complex.Complex;

/**
 * Sojourn time distribution of the M/G/1 processor-sharing queue.
 *
 * <p>Jobs arrive in a Poisson stream of rate lambda at a single egalitarian
 * processor-sharing server whose service requirement has Laplace-Stieltjes
 * transform bhat(tau) and mean m1. Writing V(x) for the sojourn time of a tagged
 * job of service requirement x and rho = lambda*m1 &lt; 1, Ott (1984) and Yashkov
 * (1983) express the conditional transform as</p>
 *
 * <pre>   E[exp(-s V(x))] = (1-rho) / D(s,x),</pre>
 *
 * <p>where D(s,x) is the inverse Laplace transform, evaluated at x, of</p>
 *
 * <pre>   f(tau;s) = [ (1-rho)*tau^2 - (1-rho)*lambda*(1-bhat(tau))*tau
 *                + s*rho*tau - s*lambda*(1-bhat(tau)) ]
 *              / [ tau^2 * (tau - s - lambda*(1-bhat(tau))) ].</pre>
 *
 * <p>The transform is exact but implicit, since f must be inverted in tau. For
 * phase-type service f(tau;s) is a proper rational function of tau, the double
 * pole at tau = 0 cancels, and D(s,x) is obtained in closed form as a finite sum
 * of residues (or, for repeated poles, from the matrix exponential of the
 * companion realization). This makes the M/PH/1-PS queue, hence every service law
 * that can be fitted with a phase-type distribution, exactly solvable. For a
 * service transform supplied as a callback, f is inverted in tau numerically on a
 * Bromwich contour placed to the right of the dominant singularity tau*(s), the
 * unique root of tau = s + lambda*(1-bhat(tau)) in the right half plane, which the
 * fixed-point iteration of that equation reaches at geometric rate rho.</p>
 *
 * <p>The conditional sojourn time is atomic on the lattice t = (k+1)*x, for
 * k = 0,1,2,...: processor sharing gives every job in the system the same amount of
 * work, so if the k jobs present on arrival all outlive the tagged job and no
 * arrival intervenes, the sojourn is exactly (k+1)*x. For exponential service the
 * masses are A_k = (1-rho)*rho^k*exp(-k*mu*x)*exp(-lambda*(k+1)*x), the k = 0 term
 * being the probability of finding the system empty and sharing it with nobody,
 * which is the only one that stays exact for general service. The k = 0 atom is
 * removed before inverting in s; the remaining atoms make cdfCond jump and leave no
 * density, so pdfCond is NaN on the lattice.</p>
 *
 * <p>References: T. J. Ott, "The sojourn-time distribution in the M/G/1 queue with
 * processor sharing", J. Appl. Prob. 21(2), 1984, pp. 360-378; S. F. Yashkov, "A
 * derivation of response time distribution for an M/G/1 processor-sharing queue",
 * Probl. Contr. Inform. Theory 12, 1983, pp. 133-148; Q. Zhen, C. Knessl,
 * "Asymptotic expansions for the sojourn time distribution in the M/G/1-PS queue",
 * Math. Meth. Oper. Res. 74, 2011, equations (2.2)-(2.5).</p>
 *
 * <p>Port of MATLAB qsys_mg1_ps.m.</p>
 */
public class Qsys_mg1_ps {

    private Qsys_mg1_ps() {
    }

    /** Service Laplace-Stieltjes transform, which must accept complex arguments. */
    public interface ServiceLST {
        Complex value(Complex tau);
    }

    /** Service density. */
    public interface ServicePDF {
        double value(double y);
    }

    private static final int NPANEL = 8;
    private static final int NGAUSS = 32;

    /**
     * Sojourn time distribution with phase-type service PH(alpha,T).
     *
     * @param lambda Poisson arrival rate
     * @param alpha  phase-type initial probability vector
     * @param T      phase-type subgenerator
     * @param x      service requirements to condition on, may be null
     * @param s      transform arguments to tabulate, may be null
     * @param t      times at which to evaluate the distribution, may be null
     * @param nterms function evaluations per numerical Laplace inversion, odd
     * @return the tabulated result, still callable for other arguments
     */
    public static QsysMg1PsResult qsys_mg1_ps(double lambda, double[] alpha, final double[][] T,
                                              double[] x, double[] s, double[] t, int nterms) {
        checkCommon(lambda, nterms);
        final int n = alpha.length;
        if (T.length != n) {
            throw new RuntimeException("qsys_mg1_ps: T must be " + n + " x " + n + " to match alpha");
        }
        double asum = 0;
        for (int i = 0; i < n; i++) {
            if (T[i].length != n) {
                throw new RuntimeException("qsys_mg1_ps: T must be " + n + " x " + n + " to match alpha");
            }
            if (alpha[i] < -1e-12) {
                throw new RuntimeException("qsys_mg1_ps: alpha must be a probability vector");
            }
            asum += alpha[i];
        }
        if (Math.abs(asum - 1) > 1e-8) {
            throw new RuntimeException("qsys_mg1_ps: alpha must be a probability vector");
        }
        final double[] exitrate = new double[n];
        for (int i = 0; i < n; i++) {
            double rowsum = 0;
            for (int j = 0; j < n; j++) {
                rowsum += T[i][j];
            }
            exitrate[i] = -rowsum;
            if (exitrate[i] < -1e-10 || T[i][i] >= 0) {
                throw new RuntimeException("qsys_mg1_ps: T must be a proper phase-type subgenerator");
            }
        }
        Matrix mT = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                mT.set(i, j, T[i][j]);
            }
        }
        Matrix negTinv = mT.scale(-1.0).inv();
        double m1 = 0;
        double m2 = 0;
        double[] v1 = new double[n];
        double[] v2 = new double[n];
        for (int i = 0; i < n; i++) {
            double r1 = 0;
            for (int j = 0; j < n; j++) {
                r1 += negTinv.get(i, j);
            }
            v1[i] = r1;
        }
        for (int i = 0; i < n; i++) {
            double r2 = 0;
            for (int j = 0; j < n; j++) {
                r2 += negTinv.get(i, j) * v1[j];
            }
            v2[i] = r2;
        }
        for (int i = 0; i < n; i++) {
            m1 += alpha[i] * v1[i];
            m2 += 2 * alpha[i] * v2[i];
        }

        final double[] alphaCopy = alpha.clone();
        ServiceLST bhat = new ServiceLST() {
            public Complex value(Complex tau) {
                return phLst(alphaCopy, T, exitrate, tau);
            }
        };
        final Matrix mTf = mT;
        ServicePDF bpdf = new ServicePDF() {
            public double value(double y) {
                Matrix E = mTf.scale(y).expm();
                double acc = 0;
                for (int i = 0; i < n; i++) {
                    double r = 0;
                    for (int j = 0; j < n; j++) {
                        r += E.get(i, j) * exitrate[j];
                    }
                    acc += alphaCopy[i] * r;
                }
                return acc;
            }
        };

        // Faddeev-LeVerrier gives det(tau*I-T) and the adjugate in one sweep, so
        // bhat(tau) = nb(tau)/db(tau) as polynomials of degree n-1 and n
        double[] db = new double[n + 1];
        double[] nb = new double[n];
        db[0] = 1.0;
        double[][] Mk = eye(n);
        for (int k = 1; k <= n; k++) {
            double acc = 0;
            for (int i = 0; i < n; i++) {
                double r = 0;
                for (int j = 0; j < n; j++) {
                    r += Mk[i][j] * exitrate[j];
                }
                acc += alpha[i] * r;
            }
            nb[k - 1] = acc;
            double[][] TM = matmul(T, Mk);
            double tr = 0;
            for (int i = 0; i < n; i++) {
                tr += TM[i][i];
            }
            db[k] = -tr / k;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Mk[i][j] = TM[i][j] + (i == j ? db[k] : 0.0);
                }
            }
        }

        double minrate = Double.POSITIVE_INFINITY;
        for (Complex ev : mT.eig()) {
            minrate = Math.min(minrate, -ev.getReal());
        }
        double ymax = Math.max(40 * m1, 40 / minrate);
        return build(lambda, m1, m2, true, db, nb, bhat, bpdf, ymax, x, s, t, nterms);
    }

    /**
     * Sojourn time distribution with the service law given by its transform.
     *
     * <p>Only the transforms and their moments are available on this path: the time
     * grid is rejected, because inverting a numerically inverted transform is
     * unstable in double precision.</p>
     *
     * @param lambda Poisson arrival rate
     * @param bhat   service Laplace-Stieltjes transform
     * @param m1     mean service requirement
     * @param x      service requirements to condition on, may be null
     * @param s      transform arguments to tabulate, may be null
     * @param nterms function evaluations per numerical Laplace inversion, odd
     * @param pdf    service density, needed to remove the conditioning, may be null
     * @return the tabulated result, still callable for other arguments
     */
    public static QsysMg1PsResult qsys_mg1_ps(double lambda, ServiceLST bhat, double m1,
                                              double[] x, double[] s, int nterms, ServicePDF pdf) {
        checkCommon(lambda, nterms);
        if (!Double.isFinite(m1) || m1 <= 0) {
            throw new RuntimeException("qsys_mg1_ps: the mean service time must be a finite positive scalar");
        }
        return build(lambda, m1, Double.NaN, false, null, null, bhat, pdf, 60 * m1, x, s, null, nterms);
    }

    private static void checkCommon(double lambda, int nterms) {
        if (!Double.isFinite(lambda) || lambda <= 0) {
            throw new RuntimeException("qsys_mg1_ps: lambda must be a finite positive scalar");
        }
        if (nterms % 2 == 0 || nterms < 11) {
            throw new RuntimeException("qsys_mg1_ps: nterms must be an odd integer of at least 11");
        }
    }

    private static QsysMg1PsResult build(double lambda, double m1, double m2, boolean isPH,
                                         double[] db, double[] nb, ServiceLST bhat, ServicePDF bpdf,
                                         double ymax, double[] x, double[] s, double[] t, int nterms) {
        double rho = lambda * m1;
        if (rho >= 1) {
            throw new RuntimeException(String.format("qsys_mg1_ps: system is unstable: utilization %.6f >= 1", rho));
        }
        QsysMg1PsResult res = new QsysMg1PsResult();
        res.lambda = lambda;
        res.rho = rho;
        res.m1 = m1;
        res.m2 = m2;
        res.isPH = isPH;
        res.db = db;
        res.nb = nb;
        res.bhat = bhat;
        res.bpdf = bpdf;
        res.nterms = nterms;
        res.shift = 0.25 * (1 - rho) / m1;
        res.meanUncond = m1 / (1 - rho);
        res.atomUncond = (1 - rho) * bhat.value(new Complex(lambda, 0)).getReal();
        res.x = x == null ? new double[0] : x.clone();
        res.s = s == null ? new double[0] : s.clone();
        res.t = t == null ? new double[0] : t.clone();

        // fixed panelled Gauss-Legendre rule for removing the conditioning: the
        // nodes do not move with s, so the service density is sampled only once
        double[][] rule = quadNodes(ymax);
        res.yq = rule[0];
        res.wq = rule[1];
        boolean hasPdf = bpdf != null;
        if (hasPdf) {
            res.bq = new double[res.yq.length];
            for (int k = 0; k < res.yq.length; k++) {
                res.bq[k] = bpdf.value(res.yq[k]);
            }
        }

        int nx = res.x.length;
        int ns = res.s.length;
        int nt = res.t.length;
        res.lstCondVal = new double[nx][ns];
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ns; j++) {
                res.lstCondVal[i][j] = lstCond(res, new Complex(res.s[j], 0), res.x[i], false).getReal();
            }
        }
        res.lstUncondVal = new double[ns];
        for (int j = 0; j < ns; j++) {
            res.lstUncondVal[j] = hasPdf ? lstUncond(res, new Complex(res.s[j], 0)).getReal() : Double.NaN;
        }

        res.atomCond = new double[nx];
        res.meanCond = new double[nx];
        res.m2Cond = new double[nx];
        res.varCond = new double[nx];
        for (int i = 0; i < nx; i++) {
            res.atomCond[i] = (1 - rho) * Math.exp(-lambda * res.x[i]);
            res.meanCond[i] = res.x[i] / (1 - rho);
            res.m2Cond[i] = secondMoment(res, res.x[i], res.meanCond[i]);
            res.varCond[i] = res.m2Cond[i] - res.meanCond[i] * res.meanCond[i];
        }
        if (hasPdf) {
            // integrate the conditional second moment, which is far better conditioned
            // than differentiating the quadrature that removes the conditioning
            double acc = 0;
            for (int k = 0; k < res.yq.length; k++) {
                acc += res.wq[k] * res.bq[k] * secondMoment(res, res.yq[k], res.yq[k] / (1 - rho));
            }
            res.m2Uncond = acc;
            res.varUncond = acc - res.meanUncond * res.meanUncond;
        }

        if (nt > 0 && !isPH) {
            throw new RuntimeException("qsys_mg1_ps: the sojourn time distribution needs phase-type "
                    + "service, since inverting a numerically inverted transform is unstable in double "
                    + "precision; with a transform handle only the LST and its moments are available");
        }
        res.pdfCond = new double[nx][nt];
        res.cdfCond = new double[nx][nt];
        for (int i = 0; i < nx; i++) {
            final double atom = res.atomCond[i];
            final double xi = res.x[i];
            final QsysMg1PsResult fres = res;
            // V(x) >= x with an atom at x, so invert the excess V(x)-x net of its atom
            LaplaceFun gpdf = new LaplaceFun() {
                public Complex value(Complex u) {
                    return lstCond(fres, u, xi, true).subtract(atom);
                }
            };
            LaplaceFun gcdf = new LaplaceFun() {
                public Complex value(Complex u) {
                    return lstCond(fres, u, xi, true).subtract(atom).divide(u);
                }
            };
            for (int j = 0; j < nt; j++) {
                if (res.t[j] < xi) {
                    continue;
                }
                if (res.t[j] == xi) {
                    res.cdfCond[i][j] = atom;
                    continue;
                }
                res.cdfCond[i][j] = ilt(gcdf, res.t[j] - xi, nterms).getReal() + atom;
                // V(x) is atomic on the lattice (k+1)*x, where no density exists
                double ratio = res.t[j] / xi;
                if (Math.abs(ratio - Math.round(ratio)) < 1e-9) {
                    res.pdfCond[i][j] = Double.NaN;
                } else {
                    res.pdfCond[i][j] = ilt(gpdf, res.t[j] - xi, nterms).getReal();
                }
            }
        }
        res.pdfUncond = new double[nt];
        res.cdfUncond = new double[nt];
        final QsysMg1PsResult fres = res;
        LaplaceFun un = new LaplaceFun() {
            public Complex value(Complex u) {
                return lstUncond(fres, u);
            }
        };
        LaplaceFun uncdf = new LaplaceFun() {
            public Complex value(Complex u) {
                return lstUncond(fres, u).divide(u);
            }
        };
        for (int j = 0; j < nt; j++) {
            if (!hasPdf) {
                res.pdfUncond[j] = Double.NaN;
                res.cdfUncond[j] = Double.NaN;
            } else if (res.t[j] <= 0) {
                res.pdfUncond[j] = 0;
                res.cdfUncond[j] = 0;
            } else {
                res.pdfUncond[j] = ilt(un, res.t[j], nterms).getReal();
                res.cdfUncond[j] = ilt(uncdf, res.t[j], nterms).getReal();
            }
        }
        return res;
    }

    /** Laplace-domain callback used by the inversion. */
    interface LaplaceFun {
        Complex value(Complex s);
    }

    static Complex lstCond(QsysMg1PsResult r, Complex s, double x, boolean excess) {
        if (x == 0) {
            return Complex.ONE;
        }
        Complex[] vs = r.isPH ? denomPh(r, s, new double[]{x}) : denomGen(r, s, new double[]{x});
        Complex val = vs[0];
        double scale = vs[1].getReal();
        Complex e = (excess ? s : Complex.ZERO).subtract(scale).multiply(x).exp();
        return e.multiply(1 - r.rho).divide(val);
    }

    static Complex lstUncond(QsysMg1PsResult r, Complex s) {
        if (r.bq == null) {
            throw new RuntimeException("qsys_mg1_ps: the service density is required to remove the "
                    + "conditioning, pass it as the pdf argument");
        }
        Complex[] vs = r.isPH ? denomPh(r, s, r.yq) : denomGen(r, s, r.yq);
        double scale = vs[vs.length - 1].getReal();
        Complex acc = Complex.ZERO;
        for (int k = 0; k < r.yq.length; k++) {
            Complex e = new Complex(-scale * r.yq[k], 0).exp();
            acc = acc.add(e.multiply((1 - r.rho) * r.wq[k] * r.bq[k]).divide(vs[k]));
        }
        return acc;
    }

    static Complex dominantRoot(QsysMg1PsResult r, Complex s) {
        // unique root of tau = s + lambda*(1-bhat(tau)) in the right half plane
        Complex tau = s;
        int maxit = Math.max(200, (int) Math.ceil(3 * Math.log(1e-15) / Math.log(Math.max(r.rho, 1e-3))));
        for (int it = 0; it < maxit; it++) {
            Complex taunew = s.add(Complex.ONE.subtract(r.bhat.value(tau)).multiply(r.lambda));
            if (taunew.subtract(tau).abs() <= 1e-14 * Math.max(1, taunew.abs())) {
                return taunew;
            }
            tau = taunew;
        }
        throw new RuntimeException("qsys_mg1_ps: the dominant root iteration did not converge");
    }

    /**
     * D(s,x) = exp(scale*x)*val for phase-type service, exactly, from the residues of
     * f(tau;s), whose double pole at the origin cancels. The returned array holds one
     * val per requested x followed by the scale.
     */
    private static Complex[] denomPh(QsysMg1PsResult r, Complex s, double[] x) {
        int n = r.nb.length;
        Complex[] dbc = toComplex(r.db);
        Complex[] nbc = toComplex(r.nb);
        Complex[] dm = new Complex[n + 1];
        dm[0] = dbc[0];
        for (int i = 1; i <= n; i++) {
            dm[i] = dbc[i].subtract(nbc[i - 1]);
        }
        Complex[] P = conv(new Complex[]{Complex.ONE, s.add(r.lambda).negate()}, dbc);
        for (int i = 0; i < n; i++) {
            P[i + 2] = P[i + 2].add(nbc[i].multiply(r.lambda));
        }
        Complex[] A = conv(new Complex[]{new Complex(1 - r.rho, 0), s.multiply(r.rho), Complex.ZERO}, dbc);
        Complex[] B = conv(new Complex[]{new Complex(1 - r.rho, 0), s}, scale(dm, r.lambda));
        for (int i = 0; i < B.length; i++) {
            A[i + 1] = A[i + 1].subtract(B[i]);
        }
        double tail = A[A.length - 1].abs() + A[A.length - 2].abs();
        double anorm = 0;
        for (int i = 0; i < A.length; i++) {
            anorm = Math.max(anorm, A[i].abs());
        }
        if (tail > 1e-6 * Math.max(1, anorm)) {
            throw new RuntimeException("qsys_mg1_ps: the double pole at the origin did not cancel");
        }
        Complex[] Ahat = new Complex[n + 1];
        System.arraycopy(A, 0, Ahat, 0, n + 1);
        Complex[] rts = roots(P);
        double sc = Double.NEGATIVE_INFINITY;
        double rmax = 0;
        for (int i = 0; i < rts.length; i++) {
            sc = Math.max(sc, rts[i].getReal());
            rmax = Math.max(rmax, rts[i].abs());
        }
        double sep = Double.POSITIVE_INFINITY;
        for (int i = 0; i < rts.length; i++) {
            for (int j = i + 1; j < rts.length; j++) {
                sep = Math.min(sep, rts[i].subtract(rts[j]).abs());
            }
        }
        Complex[] out = new Complex[x.length + 1];
        out[x.length] = new Complex(sc, 0);
        if (sep > 1e-7 * Math.max(1, rmax)) {
            Complex[] dP = polyder(P);
            Complex[] coef = new Complex[rts.length];
            for (int i = 0; i < rts.length; i++) {
                coef[i] = polyval(Ahat, rts[i]).divide(polyval(dP, rts[i]));
            }
            for (int k = 0; k < x.length; k++) {
                Complex acc = Complex.ZERO;
                for (int i = 0; i < rts.length; i++) {
                    acc = acc.add(coef[i].multiply(rts[i].subtract(sc).multiply(x[k]).exp()));
                }
                out[k] = acc;
            }
        } else {
            // repeated poles: use the companion realization of Ahat/P instead
            int m = n + 1;
            Complex[][] Ac = new Complex[m][m];
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    Ac[i][j] = Complex.ZERO;
                }
            }
            for (int j = 0; j < m; j++) {
                Ac[0][j] = P[j + 1].divide(P[0]).negate();
            }
            for (int i = 1; i < m; i++) {
                Ac[i][i - 1] = Complex.ONE;
            }
            for (int i = 0; i < m; i++) {
                Ac[i][i] = Ac[i][i].subtract(sc);
            }
            for (int k = 0; k < x.length; k++) {
                Complex[][] E = expmComplex(scale2(Ac, x[k]));
                Complex acc = Complex.ZERO;
                for (int i = 0; i < m; i++) {
                    acc = acc.add(Ahat[i].divide(P[0]).multiply(E[i][0]));
                }
                out[k] = acc;
            }
        }
        return out;
    }

    /**
     * D(s,x) for a general service transform, by inverting f(tau;s) in tau on a
     * contour placed just to the right of the dominant singularity.
     */
    private static Complex[] denomGen(final QsysMg1PsResult r, final Complex s, double[] x) {
        final double sc = dominantRoot(r, s).getReal() + r.shift;
        LaplaceFun f = new LaplaceFun() {
            public Complex value(Complex u) {
                Complex tau = u.add(sc);
                Complex bh = r.bhat.value(tau);
                Complex one = Complex.ONE.subtract(bh);
                Complex num = tau.multiply(tau).multiply(1 - r.rho)
                        .subtract(one.multiply(tau).multiply((1 - r.rho) * r.lambda))
                        .add(tau.multiply(s).multiply(r.rho))
                        .subtract(one.multiply(s).multiply(r.lambda));
                Complex den = tau.multiply(tau).multiply(tau.subtract(s).subtract(one.multiply(r.lambda)));
                return num.divide(den);
            }
        };
        Complex[] out = new Complex[x.length + 1];
        for (int k = 0; k < x.length; k++) {
            out[k] = ilt(f, x[k], r.nterms);
        }
        out[x.length] = new Complex(sc, 0);
        return out;
    }

    /**
     * Second moment from the transform curvature at the origin. The stencil is
     * one-sided so that the transform is never sampled at negative arguments, where
     * it need not converge, the step is scaled by the conditional mean but capped by
     * the unconditional one so that it stays finite as x -&gt; 0, and Richardson
     * extrapolation over h and h/2 removes the leading truncation.
     */
    private static double secondMoment(QsysMg1PsResult r, double x, double meanref) {
        if (meanref <= 0) {
            return 0;
        }
        double h = Math.min(1e-2 / meanref, 1 / r.meanUncond);
        return (16 * d2Forward(r, x, h / 2) - d2Forward(r, x, h)) / 15;
    }

    private static double d2Forward(QsysMg1PsResult r, double x, double h) {
        double[] f = new double[6];
        for (int j = 0; j < 6; j++) {
            f[j] = lstCond(r, new Complex(j * h, 0), x, false).getReal();
        }
        return (45 * f[0] - 154 * f[1] + 214 * f[2] - 156 * f[3] + 61 * f[4] - 10 * f[5]) / (12 * h * h);
    }

    /**
     * Abate-Whitt Euler inversion, symmetrized so that complex-valued time functions
     * are handled as well as real-valued ones.
     */
    static Complex ilt(LaplaceFun fun, double t, int nterms) {
        int ne = (nterms - 1) / 2;
        double[] eta = new double[2 * ne + 1];
        eta[0] = 0.5;
        for (int i = 1; i <= ne; i++) {
            eta[i] = 1.0;
        }
        eta[2 * ne] = Math.pow(2.0, -ne);
        for (int k = 1; k < ne; k++) {
            eta[2 * ne - k] = eta[2 * ne - k + 1]
                    + Math.exp(lgamma(ne + 1) - ne * Math.log(2) - lgamma(k + 1) - lgamma(ne - k + 1));
        }
        double pre = Math.pow(10.0, ne / 3.0);
        Complex g = Complex.ZERO;
        for (int k = 0; k <= 2 * ne; k++) {
            double e = pre * (1 - (k % 2) * 2) * eta[k];
            Complex bj = new Complex(ne * Math.log(10) / 3 / t, Math.PI * k / t);
            g = g.add(fun.value(bj).add(fun.value(bj.conjugate())).multiply(0.5 * e));
        }
        return g.divide(t);
    }

    /**
     * Panelled Gauss-Legendre rule on [0,ymax], with the panels growing geometrically
     * so that both ends of an exponentially decaying density are resolved. The rule is
     * fixed, so it is identical in every codebase.
     */
    static double[][] quadNodes(double ymax) {
        double[][] gl = legendre(NGAUSS);
        double[] edges = new double[NPANEL + 2];
        edges[0] = 0;
        for (int k = 0; k <= NPANEL; k++) {
            edges[k + 1] = ymax * Math.pow(2.0, k - NPANEL);
        }
        int np = edges.length - 1;
        double[] y = new double[np * NGAUSS];
        double[] w = new double[np * NGAUSS];
        for (int k = 0; k < np; k++) {
            double a = edges[k];
            double b = edges[k + 1];
            for (int i = 0; i < NGAUSS; i++) {
                y[k * NGAUSS + i] = 0.5 * (a + b) + 0.5 * (b - a) * gl[0][i];
                w[k * NGAUSS + i] = 0.5 * (b - a) * gl[1][i];
            }
        }
        return new double[][]{y, w};
    }

    /** Golub-Welsch nodes and weights of the ng-point Legendre rule on [-1,1]. */
    static double[][] legendre(int ng) {
        Matrix J = new Matrix(ng, ng);
        for (int k = 1; k < ng; k++) {
            double bk = k / Math.sqrt(4.0 * k * k - 1);
            J.set(k - 1, k, bk);
            J.set(k, k - 1, bk);
        }
        jline.io.Ret.Eigs ed = J.eigvec();
        Matrix V = ed.vectors;
        Matrix D = ed.values;
        double[] xs = new double[ng];
        double[] ws = new double[ng];
        Integer[] idx = new Integer[ng];
        final double[] diag = new double[ng];
        for (int i = 0; i < ng; i++) {
            diag[i] = D.get(i);
            idx[i] = i;
        }
        java.util.Arrays.sort(idx, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(diag[a], diag[b]);
            }
        });
        for (int i = 0; i < ng; i++) {
            xs[i] = diag[idx[i]];
            double v0 = V.get(0, idx[i]);
            double nrm = 0;
            for (int j = 0; j < ng; j++) {
                nrm += V.get(j, idx[i]) * V.get(j, idx[i]);
            }
            ws[i] = 2 * v0 * v0 / nrm;
        }
        return new double[][]{xs, ws};
    }

    private static Complex phLst(double[] alpha, double[][] T, double[] exitrate, Complex tau) {
        int n = alpha.length;
        Complex[][] M = new Complex[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M[i][j] = new Complex(-T[i][j], 0);
            }
            M[i][i] = M[i][i].add(tau);
        }
        Complex[] b = new Complex[n];
        for (int i = 0; i < n; i++) {
            b[i] = new Complex(exitrate[i], 0);
        }
        Complex[] z = solve(M, b);
        Complex acc = Complex.ZERO;
        for (int i = 0; i < n; i++) {
            acc = acc.add(z[i].multiply(alpha[i]));
        }
        return acc;
    }

    /** Gaussian elimination with partial pivoting on a complex system. */
    static Complex[] solve(Complex[][] A, Complex[] b) {
        int n = b.length;
        Complex[][] M = new Complex[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int i = c + 1; i < n; i++) {
                if (M[i][c].abs() > M[piv][c].abs()) {
                    piv = i;
                }
            }
            Complex[] tmp = M[c];
            M[c] = M[piv];
            M[piv] = tmp;
            if (M[c][c].abs() == 0) {
                throw new RuntimeException("qsys_mg1_ps: singular system in the phase-type transform");
            }
            for (int i = c + 1; i < n; i++) {
                Complex fac = M[i][c].divide(M[c][c]);
                for (int j = c; j <= n; j++) {
                    M[i][j] = M[i][j].subtract(fac.multiply(M[c][j]));
                }
            }
        }
        Complex[] z = new Complex[n];
        for (int i = n - 1; i >= 0; i--) {
            Complex acc = M[i][n];
            for (int j = i + 1; j < n; j++) {
                acc = acc.subtract(M[i][j].multiply(z[j]));
            }
            z[i] = acc.divide(M[i][i]);
        }
        return z;
    }

    /** Durand-Kerner iteration, which handles complex coefficients. */
    static Complex[] roots(Complex[] p) {
        int n = p.length - 1;
        Complex[] a = new Complex[p.length];
        for (int i = 0; i < p.length; i++) {
            a[i] = p[i].divide(p[0]);
        }
        Complex[] z = new Complex[n];
        Complex seed = new Complex(0.4, 0.9);
        z[0] = Complex.ONE;
        for (int i = 1; i < n; i++) {
            z[i] = z[i - 1].multiply(seed);
        }
        for (int it = 0; it < 500; it++) {
            double move = 0;
            for (int i = 0; i < n; i++) {
                Complex den = Complex.ONE;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        den = den.multiply(z[i].subtract(z[j]));
                    }
                }
                if (den.abs() == 0) {
                    continue;
                }
                Complex d = polyval(a, z[i]).divide(den);
                z[i] = z[i].subtract(d);
                move = Math.max(move, d.abs());
            }
            if (move < 1e-14) {
                break;
            }
        }
        return z;
    }

    static Complex[] conv(Complex[] a, Complex[] b) {
        Complex[] c = new Complex[a.length + b.length - 1];
        for (int i = 0; i < c.length; i++) {
            c[i] = Complex.ZERO;
        }
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                c[i + j] = c[i + j].add(a[i].multiply(b[j]));
            }
        }
        return c;
    }

    static Complex polyval(Complex[] p, Complex z) {
        Complex acc = Complex.ZERO;
        for (int i = 0; i < p.length; i++) {
            acc = acc.multiply(z).add(p[i]);
        }
        return acc;
    }

    static Complex[] polyder(Complex[] p) {
        int n = p.length - 1;
        if (n == 0) {
            return new Complex[]{Complex.ZERO};
        }
        Complex[] d = new Complex[n];
        for (int i = 0; i < n; i++) {
            d[i] = p[i].multiply(n - i);
        }
        return d;
    }

    private static Complex[] toComplex(double[] v) {
        Complex[] c = new Complex[v.length];
        for (int i = 0; i < v.length; i++) {
            c[i] = new Complex(v[i], 0);
        }
        return c;
    }

    private static Complex[] scale(Complex[] v, double f) {
        Complex[] c = new Complex[v.length];
        for (int i = 0; i < v.length; i++) {
            c[i] = v[i].multiply(f);
        }
        return c;
    }

    private static Complex[][] scale2(Complex[][] A, double f) {
        int n = A.length;
        Complex[][] B = new Complex[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                B[i][j] = A[i][j].multiply(f);
            }
        }
        return B;
    }

    /**
     * Matrix exponential of a complex matrix, via the real embedding
     * [[Re A, -Im A],[Im A, Re A]] so that the tested real routine does the work.
     */
    static Complex[][] expmComplex(Complex[][] A) {
        int n = A.length;
        Matrix M = new Matrix(2 * n, 2 * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M.set(i, j, A[i][j].getReal());
                M.set(i, j + n, -A[i][j].getImaginary());
                M.set(i + n, j, A[i][j].getImaginary());
                M.set(i + n, j + n, A[i][j].getReal());
            }
        }
        Matrix E = M.expm();
        Complex[][] out = new Complex[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                out[i][j] = new Complex(E.get(i, j), E.get(i + n, j));
            }
        }
        return out;
    }

    private static double[][] eye(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++) {
            I[i][i] = 1.0;
        }
        return I;
    }

    private static double[][] matmul(double[][] A, double[][] B) {
        int n = A.length;
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                double aik = A[i][k];
                if (aik == 0) {
                    continue;
                }
                for (int j = 0; j < n; j++) {
                    C[i][j] += aik * B[k][j];
                }
            }
        }
        return C;
    }

    private static double lgamma(int k) {
        double acc = 0;
        for (int i = 2; i < k; i++) {
            acc += Math.log(i);
        }
        return acc;
    }
}
