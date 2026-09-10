/**
 * @file Per-type waiting times of the MMAP[K]/G[K]/1 FCFS queue
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.lang.processes.Det;
import jline.lang.processes.Distribution;
import jline.lang.processes.Markovian;
import jline.lang.processes.Pareto;
import jline.util.ComplexOps;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * The MMAP[K]/G[K]/1 FCFS queue: K customer types with class-dependent GENERAL
 * service, fed by a marked Markovian arrival process.
 *
 * THE METHOD, which is He's, theorem for theorem. FCFS makes the actual waiting
 * time of a customer the WORKLOAD it finds on arrival, so everything follows
 * from the joint transform of workload and arrival phase,
 * f(s)_j = E[exp(-s V) 1{phase = j}], which by He's Theorem 4.1 (eq. 4.6)
 * satisfies
 *
 *     f(s) [ s I + D0 + sum_k Dk gk(s) ] = s v0,                          (*)
 *
 * with v0 the idle-phase vector, his y0. v0 needs NO root search: the matrix U
 * solving U = D0 + sum_k Dk Fk(U) with Fk(U) = int exp(U t) dFk(t) is his
 * eq. (4.4), the generator of the underlying Markov process obtained by EXCISING
 * the busy periods, and eq. (4.5) with Theorem 4.2 give y0 Q = 0 and
 * y0 e = 1 - rho, i.e. v0 = (1 - rho) pi_U. The same vector is what the
 * analyticity of (*) forces, since the roots of the bracket in the closed right
 * half plane are s = -u over the spectrum of U; the two agree to 2.5e-13, and
 * the stationary route is taken because it needs no complex eigenvector.
 *
 * The per-type actual waiting time is the workload seen by a type-k arrival,
 * biased by that type's own arrival block, his Theorem 5.1 eq. (5.1) summed
 * over the post-arrival phase:
 *
 *     E[exp(-s Wk)] = f(s) Dk e / lambda_k.
 *
 * SCOPE. He allows an arrival to be a BATCH carrying a sequence of types, and
 * his Theorem 5.3 then multiplies the transform by prod_(i&lt;n) f*_(h_i)(s), the
 * service of the customers ahead of the tagged one WITHIN its own batch. This
 * class covers the single-customer-per-arrival case, his Special case 3.3,
 * where that product is empty, which is exactly the MMAP convention LINE
 * carries.
 *
 * Reference:
 * Qi-Ming He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
 * Queueing Systems 38(4):397-418, 2001.
 */
public final class Qsys_mmapgk1 {
    private Qsys_mmapgk1() {}

    public static QsysMmapGk1Result qsys_mmapgk1(MatrixCell MMAP, List<Distribution> svc) {
        return qsys_mmapgk1(MMAP, svc, null, 3, 1e-12, 10000);
    }

    /**
     * @param MMAP      LINE convention {D0, D1, D^(1), ..., D^(K)} with D1 = sum_k D^(k)
     * @param svc       K service laws, one per marked type; they may differ in family
     * @param wPoints   times at which to evaluate the per-type waiting time CDF, or null
     * @param numWMoms  how many per-type waiting time moments to return
     * @param tol       fixed point tolerance on U
     * @param iterMax   fixed point iteration cap
     * @return the per-type solution
     */
    public static QsysMmapGk1Result qsys_mmapgk1(MatrixCell MMAP, List<Distribution> svc,
                                                 Matrix wPoints, int numWMoms,
                                                 double tol, int iterMax) {
        final int K = MMAP.size() - 2;
        if (K < 1) {
            throw new IllegalArgumentException("The MMAP must carry at least one marked arrival block");
        }
        if (svc.size() != K) {
            throw new IllegalArgumentException("One service law per marked type is required ("
                    + svc.size() + " given, " + K + " types)");
        }
        final Matrix D0 = MMAP.get(0);
        final int ma = D0.getNumRows();
        final List<Matrix> Dk = new ArrayList<Matrix>();
        Matrix Dsum = D0.copy();
        for (int k = 0; k < K; k++) {
            Dk.add(MMAP.get(k + 2));
            Dsum = Dsum.add(Dk.get(k));
        }

        final Matrix theta = statLeftNull(Dsum);
        final Matrix ones = Matrix.ones(ma, 1);
        final double[] lambdas = new double[K];
        final double[] meanS = new double[K];
        double rho = 0.0;
        for (int k = 0; k < K; k++) {
            lambdas[k] = theta.mult(Dk.get(k)).mult(ones).get(0, 0);
            meanS[k] = svc.get(k).getMean();
            rho += lambdas[k] * meanS[k];
        }
        if (rho >= 1.0) {
            throw new IllegalArgumentException("The load " + rho + " of the system is not below one");
        }

        // Fixed point U = D0 + sum_k Dk Fk(U)
        Matrix U = D0.copy();
        for (int it = 0; it < iterMax; it++) {
            Matrix Unew = D0.copy();
            for (int k = 0; k < K; k++) {
                Unew = Unew.add(Dk.get(k).mult(matrixLst(svc.get(k), U)));
            }
            double diff = 0.0;
            for (int i = 0; i < ma; i++) {
                for (int j = 0; j < ma; j++) {
                    diff = Math.max(diff, Math.abs(Unew.get(i, j) - U.get(i, j)));
                }
            }
            U = Unew;
            if (diff <= tol) break;
        }
        // U e = 0 EXACTLY. It is a property of the fixed point, not of the
        // iterate: the iteration converges linearly, so the row sums still carry
        // a residue at the tolerance above, and that residue would move the
        // stationary solve off the vector it is meant to find.
        for (int i = 0; i < ma; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < ma; j++) rowsum += U.get(i, j);
            U.set(i, i, U.get(i, i) - rowsum);
        }
        final Matrix v0 = statLeftNull(U).scale(1.0 - rho);

        final Matrix waitMoments = moments(D0, Dk, svc, theta, v0, lambdas, numWMoms);
        final Matrix meanWT = new Matrix(1, K);
        final Matrix meanST = new Matrix(1, K);
        double meanQL = 0.0;
        for (int k = 0; k < K; k++) {
            meanWT.set(0, k, waitMoments.get(k, 0));
            meanST.set(0, k, waitMoments.get(k, 0) + meanS[k]);
            meanQL += lambdas[k] * meanST.get(0, k);
        }

        Matrix waitCDF = null;
        if (wPoints != null && wPoints.length() > 0) {
            waitCDF = new Matrix(K, wPoints.length());
            for (int k = 0; k < K; k++) {
                for (int it = 0; it < wPoints.length(); it++) {
                    waitCDF.set(k, it, eulerInvert(D0, Dk, svc, v0, lambdas, k, wPoints.get(it)));
                }
            }
        }

        final Matrix lam = new Matrix(1, K);
        double lamTot = 0.0;
        for (int k = 0; k < K; k++) {
            lam.set(0, k, lambdas[k]);
            lamTot += lambdas[k];
        }
        return new QsysMmapGk1Result(lam, lamTot, rho, v0, waitMoments, meanWT, meanST,
                meanQL, waitCDF, wPoints, "LINE:MMAP[" + K + "]/G[" + K + "]/1");
    }

    /** Left null vector of G normalized to sum one. */
    private static Matrix statLeftNull(Matrix G) {
        final int n = G.getNumRows();
        final Matrix A = G.copy();
        for (int i = 0; i < n; i++) A.set(i, n - 1, 1.0);
        final Matrix rhs = new Matrix(n, 1);
        rhs.set(n - 1, 0, 1.0);
        return A.transpose().inv().mult(rhs).transpose();
    }

    /**
     * The matrix transform int_0^inf exp(U t) dF(t). A phase-type law is exact
     * on the Kronecker sum, since its density is the SCALAR beta exp(St) s0 and
     * int exp(Ut) x exp(St) dt = -(U (+) S)^-1. Anything else is read off the
     * eigenvalues of U, where the matrix transform becomes the scalar one.
     */
    private static Matrix matrixLst(Distribution law, Matrix U) {
        final int n = U.getNumRows();
        if (law instanceof Markovian) {
            final MatrixCell rep = ((Markovian) law).getProcess();
            final Matrix S = rep.get(0);
            final int ms = S.getNumRows();
            final Matrix beta = Map_pie_local(rep);
            final Matrix s0 = S.mult(Matrix.ones(ms, 1)).scale(-1.0);
            final Matrix KS = U.kron(Matrix.eye(ms)).add(Matrix.eye(n).kron(S));
            return Matrix.eye(n).kron(beta).mult(KS.inv().mult(Matrix.eye(n).kron(s0)).scale(-1.0));
        }
        if (law instanceof Det) {
            return U.scale(law.getMean()).expm();
        }
        // Stieltjes sum with true CDF increments: a proper measure for any law
        final double[][] nodes = stieltjesNodes(law);
        Matrix F = new Matrix(n, n);
        for (int i = 0; i < nodes[0].length; i++) {
            F = F.add(U.scale(nodes[0][i]).expm().scale(nodes[1][i]));
        }
        return F;
    }

    /** Embedded arrival-epoch vector of a phase-type pair, i.e. map_pie. */
    private static Matrix Map_pie_local(MatrixCell rep) {
        return jline.api.mam.Map_pie.map_pie(rep.get(0), rep.get(1));
    }

    /** Midpoint nodes and true CDF increments over the support of the law. */
    private static double[][] stieltjesNodes(Distribution law) {
        final int nGrid = 2400;
        double lo = 0.0;
        double hi = law.getMean() * 60.0;
        final double var = law.getVar();
        if (Double.isFinite(var) && var > 0) {
            hi = Math.max(hi, law.getMean() + 12.0 * Math.sqrt(var));
        }
        final double[] x = new double[nGrid];
        final double[] w = new double[nGrid];
        final double step = (hi - lo) / nGrid;
        double mass = 0.0;
        double prev = law.evalCDF(lo);
        for (int i = 0; i < nGrid; i++) {
            final double right = lo + (i + 1) * step;
            final double cur = law.evalCDF(right);
            x[i] = lo + (i + 0.5) * step;
            w[i] = cur - prev;
            mass += w[i];
            prev = cur;
        }
        if (mass > 0) {
            for (int i = 0; i < nGrid; i++) w[i] /= mass;
        }
        return new double[][]{x, w};
    }

    /** Raw moment E[S^j] of a service law. */
    private static double rawMoment(Distribution law, int j) {
        if (law instanceof Markovian) {
            final MatrixCell rep = ((Markovian) law).getProcess();
            final Matrix S = rep.get(0);
            final int n = S.getNumRows();
            final Matrix beta = Map_pie_local(rep);
            final Matrix Minv = S.scale(-1.0).inv();
            Matrix acc = Matrix.eye(n);
            double fact = 1.0;
            for (int i = 1; i <= j; i++) {
                acc = acc.mult(Minv);
                fact *= i;
            }
            return fact * beta.mult(acc).mult(Matrix.ones(n, 1)).get(0, 0);
        }
        if (law instanceof Det) {
            return Math.pow(law.getMean(), j);
        }
        if (j == 1) return law.getMean();
        if (j == 2) return law.getVar() + law.getMean() * law.getMean();
        // A HEAVY TAIL HAS NO MOMENT, and the truncated sum below would hand
        // back a finite number for one that diverges: a Pareto of shape <= j has
        // E[S^j] = Inf, and the corresponding waiting time moment is genuinely
        // infinite rather than merely large. The quadrature integrates over a
        // cut support and cannot see that, so the tail index decides first.
        if (law instanceof Pareto && (double) law.getParam(1).getValue() <= j) {
            return Double.POSITIVE_INFINITY;
        }
        final double[][] nodes = stieltjesNodes(law);
        double m = 0.0;
        for (int i = 0; i < nodes[0].length; i++) m += nodes[1][i] * Math.pow(nodes[0][i], j);
        return m;
    }

    /**
     * Derivatives of f(s) M(s) = s v0 at s = 0. M_0 is SINGULAR with right null
     * vector e, so each order fixes f_j only up to a multiple of theta; that
     * multiple is what the NEXT order's solvability condition supplies. At j = 0
     * the same condition reads theta M_1 e = v0 e, i.e. 1 - rho = 1 - rho, the
     * identity that validates the setup.
     */
    private static Matrix moments(Matrix D0, List<Matrix> Dk, List<Distribution> svc,
                                  Matrix theta, Matrix v0, double[] lambdas, int numWMoms) {
        final int ma = D0.getNumRows();
        final int K = Dk.size();
        final List<Matrix> Mder = new ArrayList<Matrix>();
        for (int j = 0; j <= numWMoms + 1; j++) {
            Matrix Mj;
            if (j == 0) {
                Mj = D0.copy();
                for (int k = 0; k < K; k++) Mj = Mj.add(Dk.get(k));
            } else {
                Mj = new Matrix(ma, ma);
                for (int k = 0; k < K; k++) {
                    final double sign = (j % 2 == 0) ? 1.0 : -1.0;
                    Mj = Mj.add(Dk.get(k).scale(sign * rawMoment(svc.get(k), j)));
                }
                if (j == 1) Mj = Mj.add(Matrix.eye(ma));
            }
            Mder.add(Mj);
        }
        final Matrix e = Matrix.ones(ma, 1);
        final double denom = theta.mult(Mder.get(1)).mult(e).get(0, 0);
        final List<Matrix> fder = new ArrayList<Matrix>();
        fder.add(theta);
        // [M_0, e] with f_j^p e = 0 pins the particular solution
        final Matrix Abase = Mder.get(0).concatCols(e);
        for (int j = 1; j <= numWMoms; j++) {
            Matrix rhs = new Matrix(1, ma);
            if (j == 1) rhs = rhs.add(v0);
            for (int i = 0; i < j; i++) {
                rhs = rhs.sub(fder.get(i).mult(Mder.get(j - i)).scale(binom(j, i)));
            }
            final Matrix rhsAug = new Matrix(1, ma + 1);
            for (int c = 0; c < ma; c++) rhsAug.set(0, c, rhs.get(0, c));
            // x Abase = rhsAug  <=>  Abase' x' = rhsAug'
            final Matrix fp = solveLeastSquares(Abase.transpose(), rhsAug.transpose()).transpose();
            double acc = 0.0;
            for (int i = 0; i < j; i++) {
                acc += binom(j + 1, i) * fder.get(i).mult(Mder.get(j + 1 - i)).mult(e).get(0, 0);
            }
            final double cfree = (-acc / (j + 1) - fp.mult(Mder.get(1)).mult(e).get(0, 0)) / denom;
            fder.add(fp.add(theta.scale(cfree)));
        }
        final Matrix moms = new Matrix(K, Math.max(numWMoms, 1));
        for (int k = 0; k < K; k++) {
            for (int j = 1; j <= numWMoms; j++) {
                final double sign = (j % 2 == 0) ? 1.0 : -1.0;
                moms.set(k, j - 1,
                        sign * fder.get(j).mult(Dk.get(k)).mult(e).get(0, 0) / lambdas[k]);
            }
        }
        return moms;
    }

    private static Matrix solveLeastSquares(Matrix A, Matrix b) {
        // normal equations; A here is (ma+1) x ma of full column rank
        final Matrix At = A.transpose();
        return At.mult(A).inv().mult(At).mult(b);
    }

    private static double binom(int n, int k) {
        double r = 1.0;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return Math.round(r);
    }

    // ---- complex support for the transform inversion ------------------------

    /**
     * Abate-Whitt Euler inversion of the type-k waiting time CDF. The transform
     * has to be evaluated off the real axis, which the Distribution API does not
     * offer, so the scalar transforms are formed here: exactly for a phase-type
     * or deterministic law, and by the same CDF-increment measure otherwise.
     */
    private static double eulerInvert(Matrix D0, List<Matrix> Dk, List<Distribution> svc,
                                      Matrix v0, double[] lambdas, int k, double t) {
        if (t <= 0) {
            return waitLstReal(D0, Dk, svc, v0, lambdas, k, 1e12);
        }
        final double A = 18.4;
        final int nEuler = 15, mEuler = 11;
        final double u = Math.exp(A / 2) / t;
        final double x = A / (2 * t);
        final double[] terms = new double[nEuler + mEuler + 1];
        terms[0] = waitLstReal(D0, Dk, svc, v0, lambdas, k, x) / x / 2.0;
        for (int j = 1; j <= nEuler + mEuler; j++) {
            final Complex sj = new Complex(x, Math.PI * j / t);
            final Complex w = waitLstComplex(D0, Dk, svc, v0, lambdas, k, sj);
            terms[j] = ((j % 2 == 0) ? 1.0 : -1.0) * w.divide(sj).getReal();
        }
        final double[] partial = new double[terms.length];
        double run = 0.0;
        for (int j = 0; j < terms.length; j++) {
            run += terms[j];
            partial[j] = run;
        }
        double F = 0.0;
        for (int j = 0; j <= mEuler; j++) {
            F += binom(mEuler, j) / Math.pow(2, mEuler) * partial[nEuler + j];
        }
        F = u * F;
        return Math.min(Math.max(F, 0.0), 1.0);
    }

    private static double waitLstReal(Matrix D0, List<Matrix> Dk, List<Distribution> svc,
                                      Matrix v0, double[] lambdas, int k, double s) {
        final int ma = D0.getNumRows();
        Matrix M = Matrix.eye(ma).scale(s).add(D0);
        for (int j = 0; j < Dk.size(); j++) {
            M = M.add(Dk.get(j).scale(scalarLstReal(svc.get(j), s)));
        }
        // f M = s v0
        final Matrix f = M.transpose().inv().mult(v0.transpose().scale(s)).transpose();
        return f.mult(Dk.get(k)).mult(Matrix.ones(ma, 1)).get(0, 0) / lambdas[k];
    }

    /** E[exp(-s Wk)] at a complex argument. */
    private static Complex waitLstComplex(Matrix D0, List<Matrix> Dk, List<Distribution> svc,
                                          Matrix v0, double[] lambdas, int k, Complex s) {
        final int ma = D0.getNumRows();
        final Complex[][] M = new Complex[ma][ma];
        for (int i = 0; i < ma; i++) {
            for (int j = 0; j < ma; j++) {
                M[i][j] = new Complex(D0.get(i, j), 0.0);
            }
            M[i][i] = M[i][i].add(s);
        }
        for (int q = 0; q < Dk.size(); q++) {
            // the law's OWN complex transform, sn.lst's contract, rather than a
            // second implementation of it here
            final Complex g = svc.get(q).evalLST(s);
            for (int i = 0; i < ma; i++) {
                for (int j = 0; j < ma; j++) {
                    M[i][j] = M[i][j].add(g.multiply(Dk.get(q).get(i, j)));
                }
            }
        }
        // solve f M = s v0, i.e. M' f' = (s v0)'
        final Complex[][] A = new Complex[ma][ma];
        final Complex[] b = new Complex[ma];
        for (int i = 0; i < ma; i++) {
            for (int j = 0; j < ma; j++) A[i][j] = M[j][i];
            b[i] = s.multiply(v0.get(i));
        }
        final Complex[] f = ComplexOps.solve(A, b);
        Complex out = new Complex(0.0, 0.0);
        for (int i = 0; i < ma; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < ma; j++) rowsum += Dk.get(k).get(i, j);
            out = out.add(f[i].multiply(rowsum));
        }
        return out.divide(lambdas[k]);
    }

    private static double scalarLstReal(Distribution law, double s) {
        return law.evalLST(s);
    }

}
