/**
 * @file M3PP(2,m) approximate count-based fitting with optimization
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import com.quantego.josqp.Model;
import com.quantego.josqp.OSQP;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.*;
import jline.api.mam.Mmpp2_fitc_approx;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;

public final class M3pp2m_fitc_approx {

    private M3pp2m_fitc_approx() {}

    public static Matrix[] m3pp2m_fitc_approx(double a, double bt1, double bt2, double binf,
                                              double m3t2, double t1, double t2,
                                              double[] ai, double[] dvt3, double t3) {
        double sumAi = 0.0;
        for (double v : ai) sumAi += v;
        if (Math.abs(a - sumAi) > 1e-8) {
            throw new IllegalArgumentException("Inconsistent per-class arrival rates.");
        }

        int m = ai.length;
        Matrix[] mmppFit = Mmpp2_fitc_approx.mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);

        if (mmppFit[0].getNumRows() == 1) {
            Matrix D0 = mmppFit[0];
            Matrix D1 = mmppFit[1];
            Matrix[] fit = new Matrix[2 + m];
            fit[0] = D0;
            fit[1] = D1;
            for (int i = 0; i < m; i++) {
                fit[2 + i] = D1.scale(ai[i] / a);
            }
            return fit;
        }

        if (m == 1) {
            return new Matrix[]{mmppFit[0], mmppFit[1], mmppFit[1]};
        }

        double l1 = mmppFit[1].get(0, 0);
        double l2 = mmppFit[1].get(1, 1);
        double r1 = mmppFit[0].get(0, 1);
        double r2 = mmppFit[0].get(1, 0);
        double t = t3;

        double sinhTerm = Math.sinh((r1 * t) / 2 + (r2 * t) / 2);
        double expTerm = Math.exp(-(r1 * t) / 2 - (r2 * t) / 2);
        double sinhExp = sinhTerm * expTerm;

        double q1i_ai = (Math.pow(r1, 4) * t / 2 + Math.pow(r2, 4) * t / 2 - l1 * Math.pow(r2, 3) * t +
                l2 * Math.pow(r2, 3) * t + 2 * r1 * Math.pow(r2, 3) * t + 2 * Math.pow(r1, 3) * r2 * t +
                3 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 2 * l1 * Math.pow(r2, 2) * sinhExp -
                2 * l2 * Math.pow(r2, 2) * sinhExp - 2 * l1 * r1 * Math.pow(r2, 2) * t -
                l1 * Math.pow(r1, 2) * r2 * t + 2 * l2 * r1 * Math.pow(r2, 2) * t +
                l2 * Math.pow(r1, 2) * r2 * t + 2 * l1 * r1 * r2 * sinhExp -
                2 * l2 * r1 * r2 * sinhExp) /
                (l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp -
                        l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t));

        double q1i_dvi = -(Math.pow(r1, 4) + 4 * Math.pow(r1, 3) * r2 + 6 * Math.pow(r1, 2) * Math.pow(r2, 2) +
                4 * r1 * Math.pow(r2, 3) + Math.pow(r2, 4)) /
                (4 * l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp -
                        l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t));

        double q1i_const = -(l1 * Math.pow(r2, 4) * t + l2 * Math.pow(r1, 4) * t + 3 * l1 * r1 * Math.pow(r2, 3) * t +
                l1 * Math.pow(r1, 3) * r2 * t + l2 * r1 * Math.pow(r2, 3) * t + 3 * l2 * Math.pow(r1, 3) * r2 * t +
                3 * l1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 2 * Math.pow(l1, 2) * r1 * Math.pow(r2, 2) * t +
                2 * Math.pow(l1, 2) * Math.pow(r1, 2) * r2 * t + 3 * l2 * Math.pow(r1, 2) * Math.pow(r2, 2) * t +
                2 * Math.pow(l2, 2) * r1 * Math.pow(r2, 2) * t + 2 * Math.pow(l2, 2) * Math.pow(r1, 2) * r2 * t -
                4 * Math.pow(l1, 2) * r1 * r2 * sinhExp - 4 * Math.pow(l2, 2) * r1 * r2 * sinhExp -
                4 * l1 * l2 * r1 * Math.pow(r2, 2) * t - 4 * l1 * l2 * Math.pow(r1, 2) * r2 * t +
                8 * l1 * l2 * r1 * r2 * sinhExp) /
                (4 * l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp - 2 * l2 * sinhExp -
                        l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t));

        double denom2 = (r1 + r2) * (Math.pow(l2, 2) * Math.pow(r1, 2) * t - 2 * Math.pow(l2, 2) * r1 * sinhExp -
                l1 * l2 * Math.pow(r1, 2) * t + Math.pow(l2, 2) * r1 * r2 * t +
                2 * l1 * l2 * r1 * sinhExp - l1 * l2 * r1 * r2 * t);

        double q2i_ai = -(Math.pow(r1, 4) * t / 2 + Math.pow(r2, 4) * t / 2 + l1 * Math.pow(r1, 3) * t -
                l2 * Math.pow(r1, 3) * t + 2 * r1 * Math.pow(r2, 3) * t + 2 * Math.pow(r1, 3) * r2 * t +
                3 * Math.pow(r1, 2) * Math.pow(r2, 2) * t - 2 * l1 * Math.pow(r1, 2) * sinhExp +
                2 * l2 * Math.pow(r1, 2) * sinhExp + l1 * r1 * Math.pow(r2, 2) * t +
                2 * l1 * Math.pow(r1, 2) * r2 * t - l2 * r1 * Math.pow(r2, 2) * t -
                2 * l2 * Math.pow(r1, 2) * r2 * t - 2 * l1 * r1 * r2 * sinhExp +
                2 * l2 * r1 * r2 * sinhExp) / denom2;

        double q2i_dvi = (Math.pow(r1, 4) + 4 * Math.pow(r1, 3) * r2 + 6 * Math.pow(r1, 2) * Math.pow(r2, 2) +
                4 * r1 * Math.pow(r2, 3) + Math.pow(r2, 4)) / (4 * denom2);

        double q2i_const = (l1 * Math.pow(r2, 4) * t + l2 * Math.pow(r1, 4) * t + 3 * l1 * r1 * Math.pow(r2, 3) * t +
                l1 * Math.pow(r1, 3) * r2 * t + l2 * r1 * Math.pow(r2, 3) * t + 3 * l2 * Math.pow(r1, 3) * r2 * t +
                3 * l1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 2 * Math.pow(l1, 2) * r1 * Math.pow(r2, 2) * t +
                2 * Math.pow(l1, 2) * Math.pow(r1, 2) * r2 * t + 3 * l2 * Math.pow(r1, 2) * Math.pow(r2, 2) * t +
                2 * Math.pow(l2, 2) * r1 * Math.pow(r2, 2) * t + 2 * Math.pow(l2, 2) * Math.pow(r1, 2) * r2 * t -
                4 * Math.pow(l1, 2) * r1 * r2 * sinhExp - 4 * Math.pow(l2, 2) * r1 * r2 * sinhExp -
                4 * l1 * l2 * r1 * Math.pow(r2, 2) * t - 4 * l1 * l2 * Math.pow(r1, 2) * r2 * t +
                8 * l1 * l2 * r1 * r2 * sinhExp) / (4 * denom2);

        int n = m;
        double[][] H = new double[n][n];
        for (int i = 0; i < n; i++) H[i][i] = 2.0 / (dvt3[i] * dvt3[i]);
        double[] h = new double[n];
        for (int i = 0; i < n; i++) h[i] = -2.0 / dvt3[i];

        List<double[]> constraintsList = new ArrayList<double[]>();
        List<Double> boundsList = new ArrayList<Double>();
        for (int i = 0; i < m; i++) {
            double[] c1 = new double[n]; c1[i] = -q1i_dvi;
            constraintsList.add(c1); boundsList.add(q1i_ai * ai[i] + q1i_const);
            double[] c2 = new double[n]; c2[i] = q1i_dvi;
            constraintsList.add(c2); boundsList.add(1.0 - q1i_ai * ai[i] - q1i_const);
            double[] c3 = new double[n]; c3[i] = -q2i_dvi;
            constraintsList.add(c3); boundsList.add(q2i_ai * ai[i] + q2i_const);
            double[] c4 = new double[n]; c4[i] = q2i_dvi;
            constraintsList.add(c4); boundsList.add(1.0 - q2i_ai * ai[i] - q2i_const);
        }

        List<double[]> eqConstraintsList = new ArrayList<double[]>();
        List<Double> eqBoundsList = new ArrayList<Double>();
        double[] eq1 = new double[n];
        for (int i = 0; i < n; i++) eq1[i] = q1i_dvi;
        eqConstraintsList.add(eq1); eqBoundsList.add(1.0 - m * q1i_const - q1i_ai * a);
        double[] eq2 = new double[n];
        for (int i = 0; i < n; i++) eq2[i] = q2i_dvi;
        eqConstraintsList.add(eq2); eqBoundsList.add(1.0 - m * q2i_const - q2i_ai * a);

        double[][] A = constraintsList.toArray(new double[0][]);
        double[] b = new double[boundsList.size()];
        for (int i = 0; i < b.length; i++) b[i] = boundsList.get(i);
        double[][] Aeq = eqConstraintsList.toArray(new double[0][]);
        double[] beq = new double[eqBoundsList.size()];
        for (int i = 0; i < beq.length; i++) beq[i] = eqBoundsList.get(i);

        double[] optimalDv = solveQP(H, h, A, b, Aeq, beq, dvt3);

        double[][] q = new double[2][m];
        for (int i = 0; i < m; i++) {
            double Ai = ai[i];
            double Dvi = optimalDv[i];
            q[0][i] = q1i_ai * Ai + q1i_dvi * Dvi + q1i_const;
            q[1][i] = q2i_ai * Ai + q2i_dvi * Dvi + q2i_const;
        }

        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < m; i++) q[j][i] = Math.max(0.0, q[j][i]);
            double sum = 0.0;
            for (int i = 0; i < m; i++) sum += q[j][i];
            if (sum > 1.0) for (int i = 0; i < m; i++) q[j][i] /= sum;
        }

        Matrix D0 = mmppFit[0];
        Matrix D1 = mmppFit[1];
        Matrix[] fit = new Matrix[2 + m];
        fit[0] = D0;
        fit[1] = D1;
        for (int i = 0; i < m; i++) {
            Matrix diag = new Matrix(2, 2);
            diag.set(0, 0, q[0][i]);
            diag.set(1, 1, q[1][i]);
            fit[2 + i] = D1.elementMult(diag);
        }
        if (!Mmap_isfeasible.mmap_isfeasible(new MatrixCell(fit))) {
            throw new IllegalStateException("Infeasible fitted M3PP");
        }
        return fit;
    }

    private static double[] solveQP(double[][] H, double[] h, double[][] A, double[] b,
                                    double[][] Aeq, double[] beq, double[] initialGuess) {
        int n = h.length;
        Model.Builder builder = Model.getBuilder();
        Model.Variable[] vars = new Model.Variable[n];
        for (int i = 0; i < n; i++) vars[i] = builder.addVariable().lb(1e-6);

        Model.Objective obj = builder.setObjective();
        for (int i = 0; i < n; i++) obj.add(h[i], vars[i]);
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double coeff = (i == j) ? H[i][j] : H[i][j] + H[j][i];
                if (Math.abs(coeff) > 1e-12) obj.add(coeff / 2.0, vars[i], vars[j]);
            }
        }
        obj.minimize();
        for (int i = 0; i < A.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            for (int j = 0; j < n; j++) if (Math.abs(A[i][j]) > 1e-12) ctr.add(A[i][j], vars[j]);
            ctr.leq(b[i]);
        }
        for (int i = 0; i < Aeq.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            for (int j = 0; j < n; j++) if (Math.abs(Aeq[i][j]) > 1e-12) ctr.add(Aeq[i][j], vars[j]);
            ctr.eq(beq[i]);
        }
        Model model = builder.build();
        model.getParam().verbose = GlobalConstants.getVerbose() != VerboseLevel.SILENT;
        OSQP.Status status = model.solve();
        if (status == OSQP.Status.SOLVED) {
            double[] result = new double[n];
            for (int i = 0; i < n; i++) result[i] = model.getSolution(vars[i].getIndex());
            return result;
        }
        return initialGuess;
    }
}
