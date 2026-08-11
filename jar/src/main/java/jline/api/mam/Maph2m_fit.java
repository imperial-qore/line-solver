/**
 * @file Second-order MAPH[m] fitting (class probabilities and backward moments)
 *
 * Port of m3a maph2m_fit.m / maph2m_fit_multiclass.m. The underlying APH(2) is
 * supplied by aph2_fit; the class marking probabilities q(j,c) are then affine
 * in the fitted backward moments, so the fit is a small quadratic program with
 * the same objective, constraints and bounds as the MATLAB reference.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.quantego.josqp.OSQP;
import com.quantego.josqp.Model;

import jline.io.Ret;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Maph2m_fit {
    private Maph2m_fit() {}

    private static final double DEGENTOL = 1e-6;
    private static final double FEASTOL = 1e-8;

    /**
     * Computes the second-order MAPH[m] fitting the first three moments, the class
     * probabilities (always fitted exactly) and the first-order backward moments.
     *
     * @param M1 first moment of the inter-arrival times
     * @param M2 second moment of the inter-arrival times
     * @param M3 third moment of the inter-arrival times
     * @param P  class probabilities
     * @param B  first-order backward moments
     * @return the fitted MAPH(2,m) as {D0, D1, D11, ..., D1m}
     */
    public static MatrixCell maph2m_fit(double M1, double M2, double M3, Matrix P, Matrix B) {
        Matrix backward = (B.getNumRows() == 1) ? B.transpose() : B;

        Ret.mamAPH2Fit aphResult = Aph2_fit.aph2_fit(M1, M2, M3);
        Collection<MatrixCell> aphsCol = aphResult.APHS.values();
        List<MatrixCell> aphs = new ArrayList<MatrixCell>(aphsCol);

        List<MatrixCell> maphs = new ArrayList<MatrixCell>();
        List<Double> errors = new ArrayList<Double>();

        for (MatrixCell aph : aphs) {
            try {
                Pair<MatrixCell, Matrix> res = maph2m_fit_multiclass(aph, P, backward, null);
                maphs.add(res.getLeft());
                Matrix fB = res.getRight();
                double error = 0.0;
                for (int c = 0; c < backward.getNumRows(); c++) {
                    double rel = fB.get(c, 0) / backward.get(c, 0) - 1.0;
                    error += rel * rel;
                }
                errors.add(Double.valueOf(error));
            } catch (RuntimeException e) {
                // this candidate APH(2) form admits no feasible marking
                maphs.add(null);
                errors.add(Double.valueOf(Double.MAX_VALUE));
            }
        }

        int best = -1;
        double bestErr = Double.MAX_VALUE;
        for (int i = 0; i < errors.size(); i++) {
            if (maphs.get(i) != null && errors.get(i).doubleValue() <= bestErr) {
                bestErr = errors.get(i).doubleValue();
                best = i;
            }
        }
        if (best < 0) {
            throw new RuntimeException("Fitting MAPH(2,m): feasibility could not be restored");
        }
        return maphs.get(best);
    }

    /** Convenience overload with unit class weights. */
    public static Pair<MatrixCell, Matrix> maph2m_fit_multiclass(MatrixCell aph, Matrix P, Matrix B) {
        return maph2m_fit_multiclass(aph, P, B, null);
    }

    /**
     * Marks a canonical acyclic APH(2) so that the class probabilities are matched
     * exactly and the backward moments as closely as possible.
     *
     * @param aph          underlying APH(2) in canonical acyclic form
     * @param P            class probabilities
     * @param B            target first-order backward moments
     * @param classWeights per-class weights in the objective (null for unit weights)
     * @return the fitted MAPH and the backward moments it realizes
     */
    public static Pair<MatrixCell, Matrix> maph2m_fit_multiclass(MatrixCell aph, Matrix P, Matrix B,
                                                                 double[] classWeights) {
        Matrix D0 = aph.get(0);
        Matrix D1 = aph.get(1);
        if (D0.getNumRows() != 2) {
            throw new IllegalArgumentException("Underlying APH must be of second-order");
        }
        if (D0.get(1, 0) != 0.0) {
            throw new IllegalArgumentException("Underlying APH must be acyclic");
        }
        if (D1.get(0, 1) != 0.0 || D1.get(1, 1) != 0.0) {
            throw new IllegalArgumentException("Underlying APH must be in canonical acyclic form");
        }

        int k = (P.getNumRows() == 1) ? P.getNumCols() : P.getNumRows();
        double[] p = new double[k];
        for (int c = 0; c < k; c++) {
            p[c] = (P.getNumRows() == 1) ? P.get(0, c) : P.get(c, 0);
        }
        Matrix backward = (B.getNumRows() == 1) ? B.transpose() : B;

        double[] w = new double[k];
        if (classWeights == null) {
            for (int c = 0; c < k; c++) w[c] = 1.0;
        } else {
            System.arraycopy(classWeights, 0, w, 0, k);
        }

        double h1 = -1.0 / D0.get(0, 0);
        double h2 = -1.0 / D0.get(1, 1);
        double r1 = D0.get(0, 1) * h1;

        double[][] q = new double[2][k];
        Matrix fB = new Matrix(k, 1);

        if (Math.abs(1.0 - r1) < DEGENTOL) {
            // degenerate form: one degree of freedom, match the class probabilities
            for (int c = 0; c < k; c++) {
                q[0][c] = p[c];
                q[1][c] = p[c];
            }
        } else {
            // q(j,c) = fB(c) * q_b(j,c) + q_0(j,c)
            double[][] qb = new double[2][k];
            double[][] q0 = new double[2][k];
            for (int c = 0; c < k; c++) {
                qb[0][c] = p[c] * (1.0 / (h2 * (r1 - 1.0)));
                q0[0][c] = p[c] * (-(h1 + h2) / (h2 * (r1 - 1.0)));
                qb[1][c] = p[c] * (1.0 / (h2 * r1));
                q0[1][c] = p[c] * (-h1 / (h2 * r1));
            }

            // inequality constraints: 0 <= q(j,c) <= 1
            double[][] A = new double[4 * k][k];
            double[] b = new double[4 * k];
            for (int c = 0; c < k; c++) {
                for (int j = 0; j < 2; j++) {
                    int row = c * 4 + j * 2;
                    A[row][c] = qb[j][c];
                    b[row] = 1.0 - q0[j][c];
                    A[row + 1][c] = -qb[j][c];
                    b[row + 1] = q0[j][c];
                }
            }

            // equality constraints: the marking probabilities sum to one in each phase
            double[][] Aeq = new double[2][k];
            double[] beq = new double[]{1.0, 1.0};
            for (int c = 0; c < k; c++) {
                for (int j = 0; j < 2; j++) {
                    Aeq[j][c] = qb[j][c];
                    beq[j] -= q0[j][c];
                }
            }

            // objective: sum_c w_c (x_c/B_c - 1)^2, up to an additive constant
            double[][] H = new double[k][k];
            double[] h = new double[k];
            for (int c = 0; c < k; c++) {
                double bc = backward.get(c, 0);
                H[c][c] = 2.0 / (bc * bc) * w[c];
                h[c] = -2.0 / bc * w[c];
            }

            double[] x = solveQP(H, h, A, b, Aeq, beq, 1e-6, 1e6);
            for (int c = 0; c < k; c++) {
                fB.set(c, 0, x[c]);
                for (int j = 0; j < 2; j++) {
                    q[j][c] = x[c] * qb[j][c] + q0[j][c];
                }
            }
        }

        for (int j = 0; j < 2; j++) {
            if (!isFeasible(q[j])) {
                throw new RuntimeException("Fitting MAPH(2,m): Feasibility could not be restored");
            }
            q[j] = clipAndNormalize(q[j]);
        }

        MatrixCell maph = new MatrixCell(2 + k);
        maph.set(0, D0.copy());
        maph.set(1, D1.copy());
        for (int c = 0; c < k; c++) {
            // maph{2+c} = D1 .* [q(1,c) 0; q(2,c) 0]
            Matrix Dc = Matrix.zeros(2, 2);
            Dc.set(0, 0, D1.get(0, 0) * q[0][c]);
            Dc.set(1, 0, D1.get(1, 0) * q[1][c]);
            maph.set(2 + c, Dc);
        }

        if (Math.abs(1.0 - r1) < DEGENTOL) {
            // realized backward moments of the degenerate form
            MatrixCell cell = new MatrixCell(maph.size());
            for (int i = 0; i < maph.size(); i++) cell.set(i, maph.get(i));
            Matrix realized = Mmap_backward_moment.mmap_backward_moment(cell, Matrix.ones(1, 1));
            for (int c = 0; c < k; c++) {
                fB.set(c, 0, realized.get(c, 0));
            }
        }

        return new Pair<MatrixCell, Matrix>(maph, fB);
    }

    private static boolean isFeasible(double[] qj) {
        double sum = 0.0;
        for (int i = 0; i < qj.length; i++) {
            if (qj[i] < -FEASTOL) return false;
            sum += qj[i];
        }
        return sum <= 1.0 + FEASTOL;
    }

    private static double[] clipAndNormalize(double[] qj) {
        double[] out = new double[qj.length];
        double sum = 0.0;
        for (int i = 0; i < qj.length; i++) {
            out[i] = Math.max(qj[i], 0.0);
            sum += out[i];
        }
        for (int i = 0; i < out.length; i++) {
            out[i] /= sum;
        }
        return out;
    }

    /** Convex QP with inequality, equality and box constraints, as in MATLAB QP(). */
    private static double[] solveQP(double[][] H, double[] h, double[][] A, double[] b,
                                    double[][] Aeq, double[] beq, double lb, double ub) {
        int n = h.length;
        Model.Builder builder = Model.getBuilder();
        Model.Variable[] vars = new Model.Variable[n];
        for (int i = 0; i < n; i++) {
            vars[i] = builder.addVariable().lb(lb).ub(ub);
        }
        Model.Objective obj = builder.setObjective();
        for (int i = 0; i < n; i++) {
            obj.add(h[i], vars[i]);
            for (int j = i; j < n; j++) {
                double coeff = (i == j) ? H[i][j] : H[i][j] + H[j][i];
                if (Math.abs(coeff) > 1e-12) {
                    obj.add(coeff / 2.0, vars[i], vars[j]);
                }
            }
        }
        obj.minimize();
        for (int i = 0; i < A.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            boolean any = false;
            for (int j = 0; j < n; j++) {
                if (Math.abs(A[i][j]) > 1e-12) {
                    ctr.add(A[i][j], vars[j]);
                    any = true;
                }
            }
            if (any) ctr.leq(b[i]);
        }
        for (int i = 0; i < Aeq.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            boolean any = false;
            for (int j = 0; j < n; j++) {
                if (Math.abs(Aeq[i][j]) > 1e-12) {
                    ctr.add(Aeq[i][j], vars[j]);
                    any = true;
                }
            }
            if (any) ctr.eq(beq[i]);
        }
        Model model = builder.build();
        model.getParam().verbose = GlobalConstants.getVerbose() != VerboseLevel.SILENT;
        OSQP.Status status = model.solve();
        if (status != OSQP.Status.SOLVED) {
            throw new RuntimeException("Fitting MAPH(2,m): quadratic programming solver failed: " + status);
        }
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = model.getSolution(vars[i].getIndex());
        }
        return x;
    }
}
