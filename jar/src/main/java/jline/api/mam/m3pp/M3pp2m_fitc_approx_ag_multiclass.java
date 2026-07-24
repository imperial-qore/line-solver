/**
 * @file M3PP(2,m) auto-gamma multiclass fitting with variance constraints
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import com.quantego.josqp.Model;
import com.quantego.josqp.OSQP;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.*;
import jline.io.InputOutput;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;

public final class M3pp2m_fitc_approx_ag_multiclass {

    private M3pp2m_fitc_approx_ag_multiclass() {}

    public static Matrix[] m3pp2m_fitc_approx_ag_multiclass(
            Matrix[] mmpp, double[] ai, double[] gt3, double t3) {

        int m = ai.length;
        Matrix[] fit = mmpp;
        double sumAi = 0.0;
        for (double v : ai) sumAi += v;

        double a = Map_count_mean.map_count_mean(new MatrixCell(mmpp[0], mmpp[1]), 1.0);

        if (Math.abs(a - sumAi) > 1e-8) {
            throw new IllegalArgumentException("Inconsistent per-class arrival rates.");
        }

        if (fit[0].getNumRows() == 1) {
            Matrix D0 = fit[0];
            Matrix D1 = fit[1];
            Matrix[] result = new Matrix[2 + m];
            result[0] = D0;
            result[1] = D1;
            for (int i = 0; i < m; i++) {
                double pi = ai[i] / a;
                result[2 + i] = D1.scale(pi);
            }
            return result;
        }

        if (m == 1) {
            return new Matrix[]{fit[0], fit[1], fit[1]};
        }

        double l1 = fit[1].get(0, 0);
        double l2 = fit[1].get(1, 1);
        double r1 = fit[0].get(0, 1);
        double r2 = fit[0].get(1, 0);
        double t = t3;

        double f1num = l1 * r2 * (2 * l2 * r1 - 2 * l1 * r1 + Math.pow(r1, 3) * t + Math.pow(r2, 3) * t +
                2 * l1 * Math.pow(r1, 2) * t - 2 * l2 * Math.pow(r1, 2) * t + 3 * r1 * Math.pow(r2, 2) * t +
                3 * Math.pow(r1, 2) * r2 * t + 2 * l1 * r1 * Math.exp(-r1 * t - r2 * t) -
                2 * l2 * r1 * Math.exp(-r1 * t - r2 * t) + 2 * l1 * r1 * r2 * t - 2 * l2 * r1 * r2 * t);
        double f1 = f1num / Math.pow(r1 + r2, 4);

        double f2num = l2 * r1 * (2 * l1 * r2 - 2 * l2 * r2 + Math.pow(r1, 3) * t + Math.pow(r2, 3) * t -
                2 * l1 * Math.pow(r2, 2) * t + 2 * l2 * Math.pow(r2, 2) * t + 3 * r1 * Math.pow(r2, 2) * t +
                3 * Math.pow(r1, 2) * r2 * t - 2 * l1 * r2 * Math.exp(-r1 * t - r2 * t) +
                2 * l2 * r2 * Math.exp(-r1 * t - r2 * t) - 2 * l1 * r1 * r2 * t + 2 * l2 * r1 * r2 * t);
        double f2 = f2num / Math.pow(r1 + r2, 4);

        double tmp = f1 * l2 * r1 - f2 * l1 * r2;

        double q1i_ai = -(f2 * (r1 + r2)) / tmp;
        double q1i_gi = (l2 * r1) / tmp;
        double q1i_const = 0.0;

        double q2i_ai = (f1 * (r1 + r2)) / tmp;
        double q2i_gi = -(l1 * r2) / tmp;
        double q2i_const = 0.0;

        int n = m;

        // m3pp2m_fitc_approx_ag_multiclass.m optimizes over G, with q(1,i) = q1i_ai*ai(i)
        // + q1i_gi*Gi + q1i_const. The slope q1i_gi is routinely O(1e5), so an equality
        // residual that any iterative QP solver would call converged in G translates into
        // a gross violation of sum_i q(j,i) = 1, which is what makes the D1 partition
        // hold. The problem is therefore posed in z(i) = q(1,i), an affine and invertible
        // change of variables that leaves the optimum unchanged but gives the partition
        // unit coefficients. G is recovered from z afterwards.
        double[] beta = new double[n];   // z_i at Gi = 0
        for (int i = 0; i < n; i++) beta[i] = q1i_ai * ai[i] + q1i_const;
        double ratio = q2i_gi / q1i_gi;  // dq(2,i)/dq(1,i)
        if (q1i_gi == 0.0 || !Double.isFinite(ratio)) {
            // q(1,i) does not depend on G: the fitting problem carries no information on
            // how to split the classes, so fall back to the uniform split, which is the
            // only choice that keeps the class matrices a partition of D1.
            Matrix D1deg = fit[1];
            Matrix[] degenerate = new Matrix[2 + m];
            degenerate[0] = fit[0];
            degenerate[1] = D1deg;
            for (int i = 0; i < m; i++) degenerate[2 + i] = D1deg.scale(1.0 / m);
            return degenerate;
        }
        double[] alpha = new double[n];  // q(2,i) at z_i = 0
        for (int i = 0; i < n; i++) alpha[i] = q2i_ai * ai[i] + q2i_const - ratio * beta[i];

        // Objective sum_i (Gi/gt3(i) - 1)^2 = sum_i (c_i*z_i + d_i)^2 after substitution.
        // Non-finite entries are zeroed, as the MATLAB code does; substituting a surrogate
        // gt3 would silently fit a different target.
        double[] c = new double[n];
        double[] d = new double[n];
        for (int i = 0; i < n; i++) {
            c[i] = 1.0 / (q1i_gi * gt3[i]);
            d[i] = -beta[i] * c[i] - 1.0;
        }
        double[][] H = new double[n][n];
        double[] h = new double[n];
        for (int i = 0; i < n; i++) {
            H[i][i] = 2.0 * c[i] * c[i];
            h[i] = 2.0 * c[i] * d[i];
        }
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(h[i])) h[i] = 0.0;
            for (int j = 0; j < n; j++) {
                if (!Double.isFinite(H[i][j])) H[i][j] = 0.0;
            }
        }
        // Rescale the objective so its curvature is O(1); a positive factor leaves the
        // argmin untouched but keeps OSQP's regularization from swamping the quadratic term.
        double hMax = 0.0;
        for (int i = 0; i < n; i++) hMax = Math.max(hMax, H[i][i]);
        if (hMax > 0.0 && Double.isFinite(hMax)) {
            for (int i = 0; i < n; i++) {
                H[i][i] /= hMax;
                h[i] /= hMax;
            }
        }

        // Inequalities: q(1,i) >= 0 and q(2,i) >= 0.
        List<double[]> constraintsList = new ArrayList<double[]>();
        List<Double> boundsList = new ArrayList<Double>();
        for (int i = 0; i < m; i++) {
            double[] c1 = new double[n]; c1[i] = -1.0;
            constraintsList.add(c1);
            boundsList.add(0.0);

            double[] c2 = new double[n]; c2[i] = -ratio;
            constraintsList.add(c2);
            boundsList.add(alpha[i]);
        }

        // Equalities: sum_i q(1,i) = 1 and sum_i q(2,i) = 1.
        List<double[]> eqConstraintsList = new ArrayList<double[]>();
        List<Double> eqBoundsList = new ArrayList<Double>();

        double[] eq1 = new double[n];
        for (int i = 0; i < n; i++) eq1[i] = 1.0;
        eqConstraintsList.add(eq1);
        eqBoundsList.add(1.0);

        double[] eq2 = new double[n];
        for (int i = 0; i < n; i++) eq2[i] = ratio;
        double alphaSum = 0.0;
        for (int i = 0; i < n; i++) alphaSum += alpha[i];
        eqConstraintsList.add(eq2);
        eqBoundsList.add(1.0 - alphaSum);

        for (int i = 0; i < boundsList.size(); i++) {
            if (!Double.isFinite(boundsList.get(i))) boundsList.set(i, 0.0);
        }
        for (int i = 0; i < eqBoundsList.size(); i++) {
            if (!Double.isFinite(eqBoundsList.get(i))) eqBoundsList.set(i, 1.0);
        }

        double[][] A = constraintsList.toArray(new double[0][]);
        double[] b = new double[boundsList.size()];
        for (int i = 0; i < b.length; i++) b[i] = boundsList.get(i);
        // Both equality rows are multiples of the all-ones vector, so the pair is rank
        // deficient. quadprog factorizes Aeq and works on its row space; OSQP does not,
        // and a singular KKT system stalls it at max_iter. Drop rows that duplicate an
        // earlier one up to a scalar factor: an inconsistent duplicate cannot be met by
        // any z anyway and is reported downstream by mmap_isfeasible.
        dropDependentRows(eqConstraintsList, eqBoundsList);

        double[][] Aeq = eqConstraintsList.toArray(new double[0][]);
        double[] beq = new double[eqBoundsList.size()];
        for (int i = 0; i < beq.length; i++) beq[i] = eqBoundsList.get(i);

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                if (!Double.isFinite(A[i][j])) A[i][j] = 0.0;
            }
        }
        for (int i = 0; i < Aeq.length; i++) {
            for (int j = 0; j < Aeq[i].length; j++) {
                if (!Double.isFinite(Aeq[i][j])) Aeq[i][j] = 0.0;
            }
        }

        // Fallback if the QP does not solve: the uniform split, the one point that is
        // always on the equality manifold, so the class matrices still partition D1.
        double[] uniform = new double[n];
        for (int i = 0; i < n; i++) uniform[i] = 1.0 / n;
        double[] z = solveQP(H, h, A, b, Aeq, beq, uniform);

        // Rate fitted exactly, z from the optimization. The equality constraints make
        // sum_i q(j,i) = 1, so the class matrices partition D1; clamping q here would
        // break that partition and is deliberately not done.
        double[][] q = new double[2][m];
        for (int i = 0; i < m; i++) {
            q[0][i] = z[i];
            q[1][i] = alpha[i] + ratio * z[i];
        }

        Matrix D0 = fit[0];
        Matrix D1 = fit[1];
        Matrix[] result = new Matrix[2 + m];
        result[0] = D0;
        result[1] = D1;
        for (int i = 0; i < m; i++) {
            Matrix diag = new Matrix(2, 2);
            diag.set(0, 0, q[0][i]);
            diag.set(1, 1, q[1][i]);
            result[2 + i] = D1.elementMult(diag);
        }
        if (!Mmap_isfeasible.mmap_isfeasible(new MatrixCell(result))) {
            InputOutput.line_warning("m3pp2m_fitc_approx_ag_multiclass", "Infeasible fitted M3PP");
        }
        return result;
    }

    /**
     * Solves min 0.5 x' H x + h' x subject to A x &lt;= b and Aeq x = beq.
     *
     * <p>The constraint rows carry the coefficients q*i_gi, which are routinely O(1e5)
     * while the variables are O(1e2), so the unscaled system is far outside the range
     * where OSQP's relative termination criterion is meaningful: it declares SOLVED with
     * an equality residual large enough to destroy sum_i q(j,i) = 1, and tightening the
     * tolerance alone only turns that into MAX_ITER_REACHED. Each row is therefore
     * normalized by its infinity norm before the solve, which leaves the feasible set
     * unchanged and brings all data to O(1).
     */
    private static double[] solveQP(double[][] H, double[] h, double[][] A, double[] b,
                                    double[][] Aeq, double[] beq, double[] initialGuess) {
        int n = h.length;
        double[] bScaled = b.clone();
        double[] beqScaled = beq.clone();
        A = normalizeRows(A, bScaled);
        Aeq = normalizeRows(Aeq, beqScaled);
        b = bScaled;
        beq = beqScaled;
        Model.Builder builder = Model.getBuilder();
        Model.Variable[] vars = new Model.Variable[n];
        // No variable bounds: m3pp2m_fitc_approx_ag_multiclass.m calls quadprog with
        // empty lb/ub, and the sign of Gi is already governed by the constraints above.
        for (int i = 0; i < n; i++) vars[i] = builder.addVariable();

        Model.Objective obj = builder.setObjective();
        for (int i = 0; i < n; i++) obj.add(h[i], vars[i]);
        for (int i = 0; i < n; i++) {
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
            for (int j = 0; j < n; j++) {
                if (Math.abs(A[i][j]) > 1e-12) ctr.add(A[i][j], vars[j]);
            }
            ctr.leq(b[i]);
        }
        for (int i = 0; i < Aeq.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            for (int j = 0; j < n; j++) {
                if (Math.abs(Aeq[i][j]) > 1e-12) ctr.add(Aeq[i][j], vars[j]);
            }
            ctr.eq(beq[i]);
        }

        Model model = builder.build();
        model.getParam().verbose = GlobalConstants.getVerbose() != VerboseLevel.SILENT;
        // q(j,i) is an affine map of G with slope q*i_gi, which is routinely O(1e5) here,
        // so OSQP's default 1e-3 termination tolerance leaves an equality residual large
        // enough to break sum_i q(j,i) = 1. Solve to machine accuracy and polish the
        // active set, which quadprog's interior-point-convex algorithm does by default.
        model.getParam().eps_abs = 1e-12;
        model.getParam().eps_rel = 1e-12;
        model.getParam().max_iter = 100000;
        model.getParam().polish = true;
        model.getParam().polish_refine_iter = 10;
        OSQP.Status status = model.solve();
        if (status == OSQP.Status.SOLVED) {
            double[] result = new double[n];
            for (int i = 0; i < n; i++) result[i] = model.getSolution(vars[i].getIndex());
            return result;
        }
        return initialGuess;
    }

    /**
     * Removes equality rows that are a scalar multiple of a row kept earlier, together with
     * their right-hand sides. Rows are compared after normalization by their infinity norm.
     */
    private static void dropDependentRows(List<double[]> rows, List<Double> rhs) {
        List<double[]> kept = new ArrayList<double[]>();
        List<Double> keptRhs = new ArrayList<Double>();
        for (int i = 0; i < rows.size(); i++) {
            double[] row = rows.get(i);
            double norm = 0.0;
            for (int j = 0; j < row.length; j++) norm = Math.max(norm, Math.abs(row[j]));
            if (norm == 0.0 || !Double.isFinite(norm)) continue;
            boolean dependent = false;
            for (int k = 0; k < kept.size() && !dependent; k++) {
                double[] ref = kept.get(k);
                double refNorm = 0.0;
                for (int j = 0; j < ref.length; j++) refNorm = Math.max(refNorm, Math.abs(ref[j]));
                // Sign of the scalar factor, taken from the entry that sets the norm.
                double sign = 0.0;
                for (int j = 0; j < row.length && sign == 0.0; j++) {
                    if (Math.abs(Math.abs(ref[j]) - refNorm) <= 0.0 && ref[j] != 0.0) {
                        sign = Math.signum(row[j] / norm) * Math.signum(ref[j] / refNorm);
                    }
                }
                if (sign == 0.0) sign = 1.0;
                dependent = true;
                for (int j = 0; j < row.length; j++) {
                    if (Math.abs(row[j] / norm - sign * ref[j] / refNorm) > 1e-10) {
                        dependent = false;
                        break;
                    }
                }
            }
            if (!dependent) {
                kept.add(row);
                keptRhs.add(rhs.get(i));
            }
        }
        rows.clear();
        rows.addAll(kept);
        rhs.clear();
        rhs.addAll(keptRhs);
    }

    /**
     * Returns a copy of the constraint matrix with every row divided by its infinity norm,
     * scaling the corresponding right-hand side in place. Rows that are entirely zero are
     * left untouched.
     */
    private static double[][] normalizeRows(double[][] M, double[] rhs) {
        double[][] out = new double[M.length][];
        for (int i = 0; i < M.length; i++) {
            double norm = 0.0;
            for (int j = 0; j < M[i].length; j++) {
                norm = Math.max(norm, Math.abs(M[i][j]));
            }
            out[i] = M[i].clone();
            if (norm > 0.0 && Double.isFinite(norm)) {
                for (int j = 0; j < out[i].length; j++) out[i][j] /= norm;
                rhs[i] /= norm;
            }
        }
        return out;
    }
}
