/**
 * @file Equilibrium analysis of a general QBD process with Rational Arrival
 * Process components
 *
 * <p>The process is specified directly by its level-independent blocks
 * (A0,A1,A2) and its boundary blocks (B0,B1), where A0 drives level increases,
 * A2 drives level decreases and A1 the within-level evolution. Unlike a
 * Markovian QBD the blocks need not be nonnegative: they are only required to
 * be conservative, (A0+A1+A2)*e = 0, and to define a genuine RAP through the
 * prediction-process interpretation. This makes qbd_rap strictly more general
 * than {@link Qbd_raprap1}, which builds a product-space QBD from two
 * INDEPENDENT RAPs; here the arrival process and the sequence of service times
 * may be driven from a shared phase space and therefore be cross-correlated.</p>
 *
 * <p>References: N. G. Bean and B. F. Nielsen, "Quasi-Birth-and-Death Processes
 * with Rational Arrival Process Components", Stochastic Models, 26(3), 2010,
 * pp. 309-334 (DTU technical report IMM-2007-20). The equilibrium construction
 * below is their Theorem 7 and the stability test is their Corollary 8. The
 * argument rests on the prediction-process interpretation of a RAP due to
 * Asmussen and Bladt, which is what allows a QBD argument to be carried over to
 * matrices that are not nonnegative; the same prediction process underlies the
 * conditional-vector RAP sampler in {@link Rap_sample}.</p>
 *
 * <p>Algorithm (Theorem 7): solve A0*G^2 + A1*G + A2 = 0 for G, set
 * U = A1 + A0*G and R = A0*inv(-U); find the row vector pihat0 with
 * pihat0*(B1 + R*A2) = 0 normalised to pihat0*e = 1; set pi0 = K*pihat0 with K
 * such that pi0*inv(I-R)*e = 1; then pi_n = pi0*R^n and the marginal level
 * probability is pi_n*e.</p>
 *
 * <p>Computation of G: the blocks are not nonnegative, so the probabilistic
 * iterations used for Markovian QBDs (logarithmic reduction, cyclic reduction)
 * carry no convergence guarantee here, and the paper explicitly leaves the
 * general case open ("The issue of justifying algorithms for the evaluation of
 * the matrix G for such processes has not been undertaken", Section 6). Two
 * paths are therefore taken. If A2 has rank one, A2 = u*v, then G = e*v/(v*e)
 * solves the equation exactly: G is idempotent and, by conservativity,
 * (A0+A1)*e = -A2*e = -u*(v*e), so (A0+A1)*e*v/(v*e) = -u*v = -A2. This is the
 * case deliberately chosen in the paper's example. Otherwise the quadratic
 * matrix equation is solved numerically by the natural functional iteration
 * G &lt;- inv(-A1)*(A2 + A0*G^2) used as a warm start, followed by Newton's
 * method on the Sylvester-form Jacobian, (A0*G+A1)*H + A0*H*G = -F(G), solved
 * through its Kronecker expansion. If the residual does not reach roundoff
 * level, or the iterate does not satisfy G*e = e, an exception is raised rather
 * than an unconverged G being returned.</p>
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qbd_rap {
    private Qbd_rap() {}

    /** Default number of levels reported when the caller does not choose one. */
    public static final int DEFAULT_NUM_LEVELS = 20;

    /**
     * Equilibrium analysis of a QBD with RAP components.
     *
     * @param A0        level-up block (m x m)
     * @param A1        local block (m x m)
     * @param A2        level-down block (m x m)
     * @param B0        boundary level-up block, null for A0
     * @param B1        boundary local block, null for A1
     * @param numLevels highest level reported
     * @return the level distribution, mean queue length, R, G, U, Sp(R) and pi_0
     */
    public static QbdRapResult qbd_rap(Matrix A0, Matrix A1, Matrix A2, Matrix B0, Matrix B1, int numLevels) {
        if (A0 == null || A1 == null || A2 == null) {
            throw new IllegalArgumentException("qbd_rap requires at least the blocks A0, A1 and A2.");
        }
        Matrix bb0 = (B0 == null) ? A0.copy() : B0;
        Matrix bb1 = (B1 == null) ? A1.copy() : B1;
        int m = A1.getNumRows();
        if (A1.getNumCols() != m
                || A0.getNumRows() != m || A0.getNumCols() != m
                || A2.getNumRows() != m || A2.getNumCols() != m
                || bb0.getNumRows() != m || bb0.getNumCols() != m
                || bb1.getNumRows() != m || bb1.getNumCols() != m) {
            throw new IllegalArgumentException("All QBD blocks must be square and of the same order.");
        }
        if (numLevels < 0) {
            throw new IllegalArgumentException("numLevels must be a nonnegative integer.");
        }

        Matrix e = Matrix.ones(m, 1);
        Matrix I = Matrix.eye(m);
        double blockScale = Math.max(1.0, Math.max(A0.normFrobenius(),
                Math.max(A1.normFrobenius(), A2.normFrobenius())));

        // conservativity check: see _kb/03-api-layer.md ("Qbd_rap")
        double consA = normInf(A0.add(A1).add(A2).mult(e));
        if (consA > 1e-8 * blockScale) {
            throw new RuntimeException("The repeating blocks are not conservative: ||(A0+A1+A2)*e||_inf = "
                    + consA + ". A QBD with RAP components requires (A0+A1+A2)*e = 0.");
        }
        double consB = normInf(bb0.add(bb1).mult(e));
        if (consB > 1e-8 * blockScale) {
            throw new RuntimeException("The boundary blocks are not conservative: ||(B0+B1)*e||_inf = "
                    + consB + ". A QBD with RAP components requires (B0+B1)*e = 0 at level 0.");
        }

        // Step 1: matrix G.
        Matrix G = qbd_rap_g(A0, A1, A2, blockScale);

        // Steps 2 and 3: U and R.
        Matrix U = A1.add(A0.mult(G));
        Matrix R;
        try {
            R = A0.mult(U.scale(-1.0).inv());
        } catch (RuntimeException ex) {
            throw new RuntimeException("The matrix U = A1 + A0*G is singular, R = A0*inv(-U) does not exist.");
        }

        // Corollary 8(i): positive recurrence.
        double spr = spectralRadius(R);
        // 1e-12 margin at the null-recurrent boundary: see _kb/03-api-layer.md ("Qbd_rap")
        if (spr >= 1.0 - 1e-12) {
            throw new RuntimeException("The process is not positive recurrent: Sp(R) = " + spr
                    + " >= 1 (Corollary 8 of Bean and Nielsen, 2010).");
        }

        // Step 4: boundary vector via SVD null space (relative tolerance): see _kb/03-api-layer.md ("Qbd_rap")
        Matrix V = bb1.add(R.mult(A2));
        Ret.SVD svdV = V.transpose().svd();
        Matrix sv = svdV.s;
        Matrix W = svdV.v;
        double svMax = sv.get(0, 0);
        double nullTol = 1e-8 * Math.max(svMax, 1.0);
        double svMin = sv.get(sv.getNumRows() - 1, 0);
        if (svMin > nullTol) {
            throw new RuntimeException("The boundary equation x*(B1 + R*A2) = 0 has no nontrivial solution "
                    + "(smallest singular value " + svMin + " against tolerance " + nullTol
                    + "), so the process is not positive recurrent (Corollary 8(ii) of Bean and Nielsen, 2010).");
        }
        if (m > 1 && sv.get(sv.getNumRows() - 2, 0) <= nullTol) {
            throw new RuntimeException("The boundary equation x*(B1 + R*A2) = 0 has a solution space of "
                    + "dimension greater than one, the equilibrium vector is not unique.");
        }
        Matrix pihat0 = new Matrix(1, m);
        for (int j = 0; j < m; j++) {
            pihat0.set(0, j, W.get(j, W.getNumCols() - 1));
        }
        double den = pihat0.mult(e).get(0, 0);
        if (Math.abs(den) < 1e-12 * normInf(pihat0)) {
            throw new RuntimeException("The boundary vector cannot be normalised, x*e = 0.");
        }
        pihat0 = pihat0.scale(1.0 / den);

        // Step 5: level-0 vector.
        Matrix ImR = I.sub(R);
        Matrix ImRinv = ImR.inv();
        double K = 1.0 / pihat0.mult(ImRinv).mult(e).get(0, 0);
        Matrix pi0 = pihat0.scale(K);

        // boundary consistency check (Theorem 7): see _kb/03-api-layer.md ("Qbd_rap")
        Matrix bal = pi0.mult(bb0).add(pi0.mult(R).mult(A1)).add(pi0.mult(R).mult(R).mult(A2));
        double balNorm = normInf(bal);
        if (balNorm > 1e-8 * blockScale * Math.max(normInf(pi0), 1.0)) {
            throw new RuntimeException("The boundary block B0 is inconsistent with the repeating blocks: "
                    + "||pi0*B0 + pi1*A1 + pi2*A2||_inf = " + balNorm + ". The level-0 balance equation of "
                    + "Theorem 7 requires pi0*(B0-A0) = 0.");
        }

        // Step 6: level vectors and marginal level distribution.
        Matrix pqueue = new Matrix(numLevels + 1, m);
        Matrix levelProb = new Matrix(1, numLevels + 1);
        Matrix pin = pi0.copy();
        for (int n = 0; n <= numLevels; n++) {
            double lp = 0.0;
            for (int j = 0; j < m; j++) {
                pqueue.set(n, j, pin.get(0, j));
                lp += pin.get(0, j);
            }
            levelProb.set(0, n, lp);
            pin = pin.mult(R);
        }

        // exact mean queue length in closed form: see _kb/03-api-layer.md ("Qbd_rap")
        double QN = pi0.mult(R).mult(ImRinv).mult(ImRinv).mult(e).get(0, 0);

        return new QbdRapResult(levelProb, QN, R, G, U, spr, pqueue, pi0);
    }

    public static QbdRapResult qbd_rap(Matrix A0, Matrix A1, Matrix A2, Matrix B0, Matrix B1) {
        return qbd_rap(A0, A1, A2, B0, B1, DEFAULT_NUM_LEVELS);
    }

    public static QbdRapResult qbd_rap(Matrix A0, Matrix A1, Matrix A2) {
        return qbd_rap(A0, A1, A2, null, null, DEFAULT_NUM_LEVELS);
    }

    /**
     * Solves A0*G^2 + A1*G + A2 = 0 for the matrix G, using the exact rank-one
     * closed form when A2 has rank one and otherwise functional iteration
     * followed by Newton's method. Never returns an unconverged iterate.
     */
    public static Matrix qbd_rap_g(Matrix A0, Matrix A1, Matrix A2, double blockScale) {
        int m = A1.getNumRows();
        Matrix e = Matrix.ones(m, 1);
        double resTol = 1e-10 * blockScale;

        // rank-one closed form for G: see _kb/03-api-layer.md ("Qbd_rap")
        Ret.SVD svdA2 = A2.svd();
        Matrix s2 = svdA2.s;
        if (s2.getNumRows() > 1 && s2.get(0, 0) > 0 && s2.get(1, 0) <= 1e-10 * s2.get(0, 0)) {
            Matrix W = svdA2.v;
            Matrix v = new Matrix(1, m);
            for (int j = 0; j < m; j++) {
                v.set(0, j, W.get(j, 0));
            }
            double ve = v.mult(e).get(0, 0);
            if (Math.abs(ve) < 1e-12 * normInf(v)) {
                throw new RuntimeException("A2 has rank one but its right factor v satisfies v*e = 0, "
                        + "so the closed form G = e*v/(v*e) is undefined.");
            }
            Matrix G = e.mult(v.scale(1.0 / ve));
            double res = residual(A0, A1, A2, G);
            if (res > resTol) {
                throw new RuntimeException("The rank-one closed form for G leaves a residual "
                        + "||A0*G^2 + A1*G + A2||_F = " + res + ", which is above the roundoff level "
                        + resTol + ".");
            }
            return G;
        }

        // functional iteration warm start for G: see _kb/03-api-layer.md ("Qbd_rap")
        Matrix mA1inv;
        try {
            mA1inv = A1.scale(-1.0).inv();
        } catch (RuntimeException ex) {
            throw new RuntimeException("The local block A1 is singular, the iteration for G cannot be started.");
        }
        Matrix G = new Matrix(m, m);
        for (int it = 0; it < 200; it++) {
            Matrix Gnew = mA1inv.mult(A2.add(A0.mult(G).mult(G)));
            if (!isFinite(Gnew)) {
                break;
            }
            double step = Gnew.sub(G).normFrobenius();
            G = Gnew;
            if (step <= 1e-14 * Math.max(G.normFrobenius(), 1.0)) {
                break;
            }
        }
        if (!isFinite(G)) {
            G = new Matrix(m, m);
        }

        // Newton's method via Kronecker expansion: see _kb/03-api-layer.md ("Qbd_rap")
        Matrix Im = Matrix.eye(m);
        for (int it = 0; it < 100; it++) {
            Matrix res = A0.mult(G).mult(G).add(A1.mult(G)).add(A2);
            if (res.normFrobenius() <= resTol) {
                break;
            }
            Matrix J = Im.kron(A0.mult(G).add(A1)).add(G.transpose().kron(A0)).toDense();
            Matrix rhs = new Matrix(m * m, 1);
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < m; i++) {
                    rhs.set(j * m + i, 0, -res.get(i, j));
                }
            }
            Matrix sol = new Matrix(m * m, 1);
            try {
                Matrix.solve(J, rhs.toDense(), sol);
            } catch (RuntimeException ex) {
                break;
            }
            Matrix H = new Matrix(m, m);
            for (int j = 0; j < m; j++) {
                for (int i = 0; i < m; i++) {
                    H.set(i, j, sol.get(j * m + i, 0));
                }
            }
            G = G.add(H);
            if (!isFinite(G)) {
                break;
            }
        }

        double res = Double.POSITIVE_INFINITY;
        double geErr = Double.POSITIVE_INFINITY;
        if (isFinite(G)) {
            res = residual(A0, A1, A2, G);
            geErr = normInf(G.mult(e).sub(e));
        }
        if (!(res <= resTol) || geErr > 1e-8) {
            throw new RuntimeException("Could not compute the matrix G for this QBD with RAP components: "
                    + "residual ||A0*G^2 + A1*G + A2||_F = " + res + " against a tolerance of " + resTol
                    + ", and ||G*e-e||_inf = " + geErr + ". The blocks are not nonnegative, so neither the "
                    + "functional iteration nor Newton's method is guaranteed to converge, and the "
                    + "justification of algorithms for G in this setting is left as an open problem in "
                    + "Section 6 of N. G. Bean and B. F. Nielsen, \"Quasi-Birth-and-Death Processes with "
                    + "Rational Arrival Process Components\", Stochastic Models, 26(3), 2010, pp. 309-334. "
                    + "Supply a model with a rank-one A2, for which G is available in closed form.");
        }
        return G;
    }

    private static double residual(Matrix A0, Matrix A1, Matrix A2, Matrix G) {
        return A0.mult(G).mult(G).add(A1.mult(G)).add(A2).normFrobenius();
    }

    private static boolean isFinite(Matrix A) {
        return !A.hasNaN() && !A.hasInfinite();
    }

    private static double normInf(Matrix A) {
        double maxAbs = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < A.getNumCols(); j++) {
                double v = Math.abs(A.get(i, j));
                if (v > maxAbs) {
                    maxAbs = v;
                }
            }
        }
        return maxAbs;
    }

    /**
     * Spectral radius (maximum absolute eigenvalue) of a matrix.
     */
    private static double spectralRadius(Matrix A) {
        int n = A.getNumRows();
        DMatrixRMaj dm = new DMatrixRMaj(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dm.set(i, j, A.get(i, j));
            }
        }
        EigenDecomposition_F64<DMatrixRMaj> evd = DecompositionFactory_DDRM.eig(n, false);
        evd.decompose(dm);
        double maxAbs = 0.0;
        for (int i = 0; i < evd.getNumberOfEigenvalues(); i++) {
            Complex_F64 ev = evd.getEigenvalue(i);
            double absVal = Math.sqrt(ev.real * ev.real + ev.imaginary * ev.imaginary);
            if (absVal > maxAbs) {
                maxAbs = absVal;
            }
        }
        return maxAbs;
    }
}
