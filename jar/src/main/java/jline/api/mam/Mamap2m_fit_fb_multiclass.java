/**
 * @file MAMAP forward-backward multiclass fitting.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import com.quantego.josqp.Model;
import com.quantego.josqp.OSQP;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Triple of (fitted MMAP, feasible forward moments, feasible backward moments).
 */
public final class Mamap2m_fit_fb_multiclass {
    private Mamap2m_fit_fb_multiclass() {}

    public static final class FitResult {
        public final MatrixCell mmap;
        public final double[] feasibleForwardMoments;
        public final double[] feasibleBackwardMoments;

        public FitResult(MatrixCell mmap, double[] fF, double[] fB) {
            this.mmap = mmap;
            this.feasibleForwardMoments = fF;
            this.feasibleBackwardMoments = fB;
        }
    }

    public static FitResult mamap2m_fit_fb_multiclass(
            MatrixCell map,
            double[] p,
            double[] F,
            double[] B,
            double[] classWeights,
            double[] fbWeights) {

        if (map.get(0).getNumRows() != 2) {
            throw new IllegalArgumentException("Underlying MAP must be of second-order");
        }
        if (map.get(0).get(1, 0) != 0.0) {
            throw new IllegalArgumentException("Underlying MAP must be acyclic");
        }

        // Canonical form, as in mamap2m_fit_fb_multiclass.m: D1(1,2)==0 is form 1,
        // D1(1,1)==0 is form 2. Every branch below is written against those two
        // labels, so swapping them here inverts the whole algorithm.
        int form;
        if (map.get(1).get(0, 1) == 0.0) {
            form = 1;
        } else if (map.get(1).get(0, 0) == 0.0) {
            form = 2;
        } else {
            throw new IllegalArgumentException("Underlying MAP must be in canonical acyclic form");
        }

        int k = p.length;

        double[] actualClassWeights;
        if (classWeights != null) {
            actualClassWeights = classWeights;
        } else {
            actualClassWeights = new double[k];
            for (int i = 0; i < k; i++) actualClassWeights[i] = 1.0;
        }
        double[] actualFbWeights = (fbWeights != null) ? fbWeights : new double[]{1.0, 1.0};

        MatrixCell mmap = new MatrixCell(2 + k);
        mmap.set(0, map.get(0).copy());
        mmap.set(1, map.get(1).copy());

        double h1 = -1.0 / map.get(0).get(0, 0);
        double h2 = -1.0 / map.get(0).get(1, 1);
        double r1 = map.get(0).get(0, 1) * h1;
        double r2 = map.get(1).get(1, 1) * h2;

        double degentol = 1e-8;

        boolean isPoissonCase =
                (form == 1 && (r1 < degentol || r2 > 1 - degentol
                        || Math.abs(h2 - h1 * r2) < degentol
                        || Math.abs(h1 - h2 + h2 * r1) < degentol))
                || (form == 2 && (r2 > 1 - degentol
                        || Math.abs(h1 - h2 + h2 * r1) < degentol
                        || Math.abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol));

        if (isPoissonCase) {
            return fitPoissonProcess(mmap, p, k);
        }
        if (form == 2 && r2 < degentol && Math.abs(1 - r1) < degentol) {
            return fitDegeneratePhaseType(mmap, p, k);
        }
        if (form == 1 && r2 < degentol) {
            return fitCanonicalPhaseType(map, p, B, actualClassWeights, k);
        }
        if ((form == 1 && Math.abs(1 - r1) < degentol) || (form == 2 && Math.abs(1 - r1) < degentol)) {
            return fitNonCanonicalPhaseType(mmap, p, F, B, actualClassWeights, actualFbWeights, h1, h2, r1, r2, k);
        }
        if (form == 2 && r2 < degentol) {
            return fitDegenerateGamma(mmap, p, F, B, actualClassWeights, actualFbWeights, h1, h2, r1, r2, k);
        }
        return fitGeneralCase(mmap, p, F, B, actualClassWeights, actualFbWeights, h1, h2, r1, r2, form, k);
    }

    public static FitResult mamap2m_fit_fb_multiclass(MatrixCell map, double[] p, double[] F, double[] B) {
        return mamap2m_fit_fb_multiclass(map, p, F, B, null, null);
    }

    private static FitResult fitPoissonProcess(MatrixCell mmap, double[] p, int k) {
        double h = Map_mean.map_mean(mmap);

        MatrixCell result = new MatrixCell(2 + k);
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -1.0 / h);
        result.set(0, d0);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, 1.0 / h);
        result.set(1, d1);

        for (int c = 0; c < k; c++) {
            result.set(2 + c, result.get(1).scale(p[c]));
        }

        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fF = Mmap_forward_moment.mmap_forward_moment(result, moments);
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(result, moments);

        return new FitResult(result, fF.toArray1D(), fB.toArray1D());
    }

    private static FitResult fitDegeneratePhaseType(MatrixCell mmap, double[] p, int k) {
        for (int c = 0; c < k; c++) {
            Matrix D1c = new Matrix(2, 2);
            D1c.set(0, 0, p[c]);
            D1c.set(0, 1, p[c]);
            D1c.set(1, 0, p[c]);
            D1c.set(1, 1, p[c]);
            mmap.set(2 + c, mmap.get(1).elementMult(D1c));
        }

        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fF = Mmap_forward_moment.mmap_forward_moment(mmap, moments);
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(mmap, moments);
        return new FitResult(mmap, fF.toArray1D(), fB.toArray1D());
    }

    private static FitResult fitCanonicalPhaseType(MatrixCell map, double[] p, double[] B, double[] classWeights, int k) {
        MatrixCell aph = new MatrixCell(2);
        aph.set(0, map.get(0).copy());
        aph.set(1, map.get(1).copy());
        aph.get(1).set(1, 1, 0.0);
        MatrixCell normalizedAph = Map_normalize.map_normalize(aph);

        MatrixCell mmap = new MatrixCell(2 + k);
        mmap.set(0, normalizedAph.get(0));
        mmap.set(1, normalizedAph.get(1));

        for (int c = 0; c < k; c++) {
            mmap.set(2 + c, mmap.get(1).scale(p[c]));
        }

        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fF = Mmap_forward_moment.mmap_forward_moment(mmap, moments);
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(mmap, moments);
        return new FitResult(mmap, fF.toArray1D(), fB.toArray1D());
    }

    private static FitResult fitNonCanonicalPhaseType(
            MatrixCell mmap, double[] p, double[] F, double[] B,
            double[] classWeights, double[] fbWeights,
            double h1, double h2, double r1, double r2, int k) {
        double[][] qf = new double[2][k];
        double[][] q0 = new double[2][k];

        for (int c = 0; c < k; c++) {
            qf[0][c] = p[c] * (-1.0 / ((h1 + h2 * (r1 - 1)) * (r2 - 1) * (r1 + r2 - r1 * r2)));
            q0[0][c] = p[c] * (h2 / ((r2 - 1) * (r1 + r2 - r1 * r2) * (h1 - h2 + h2 * r1)));
            qf[1][c] = p[c] * (-1.0 / (r2 * (h1 + h2 * (r1 - 1)) * (r1 + r2 - r1 * r2)));
            q0[1][c] = p[c] * ((h1 + h2 * r1) / (r2 * (r1 + r2 - r1 * r2) * (h1 - h2 + h2 * r1)));
        }

        double[] fF = solveForwardOptimization(F, qf, q0, classWeights, fbWeights, k);

        double[][] q = new double[3][k];
        for (int c = 0; c < k; c++) {
            q[0][c] = 1.0 / k;
            q[1][c] = fF[c] * qf[0][c] + q0[0][c];
            q[2][c] = fF[c] * qf[1][c] + q0[1][c];
        }

        for (int c = 0; c < k; c++) {
            Matrix D1c = new Matrix(2, 2);
            D1c.set(0, 0, q[0][c]);
            D1c.set(0, 1, 0.0);
            D1c.set(1, 0, q[1][c]);
            D1c.set(1, 1, q[2][c]);
            mmap.set(2 + c, mmap.get(1).elementMult(D1c));
        }

        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(mmap, moments);
        return new FitResult(mmap, fF, fB.toArray1D());
    }

    private static FitResult fitDegenerateGamma(
            MatrixCell mmap, double[] p, double[] F, double[] B,
            double[] classWeights, double[] fbWeights,
            double h1, double h2, double r1, double r2, int k) {
        if (fbWeights[0] >= fbWeights[1]) {
            return fitForwardMoments(mmap, p, F, classWeights, fbWeights, h1, h2, r1, r2, k);
        } else {
            return fitBackwardMoments(mmap, p, B, classWeights, fbWeights, h1, h2, r1, r2, k);
        }
    }

    private static FitResult fitGeneralCase(
            MatrixCell mmap, double[] p, double[] F, double[] B,
            double[] classWeights, double[] fbWeights,
            double h1, double h2, double r1, double r2, int form, int k) {
        double[] fF = F.clone();
        double[] fB = B.clone();
        for (int c = 0; c < k; c++) {
            mmap.set(2 + c, mmap.get(1).scale(p[c]));
        }
        return new FitResult(mmap, fF, fB);
    }

    private static FitResult fitForwardMoments(
            MatrixCell mmap, double[] p, double[] F,
            double[] classWeights, double[] fbWeights,
            double h1, double h2, double r1, double r2, int k) {
        double[] fF = F.clone();
        for (int c = 0; c < k; c++) {
            Matrix D1c = new Matrix(2, 2);
            D1c.set(0, 0, p[c]);
            D1c.set(0, 1, p[c]);
            D1c.set(1, 0, p[c]);
            D1c.set(1, 1, 1.0 / k);
            mmap.set(2 + c, mmap.get(1).elementMult(D1c));
        }
        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(mmap, moments);
        return new FitResult(mmap, fF, fB.toArray1D());
    }

    private static FitResult fitBackwardMoments(
            MatrixCell mmap, double[] p, double[] B,
            double[] classWeights, double[] fbWeights,
            double h1, double h2, double r1, double r2, int k) {
        double[] fB = B.clone();
        for (int c = 0; c < k; c++) {
            Matrix D1c = new Matrix(2, 2);
            D1c.set(0, 0, p[c]);
            D1c.set(0, 1, p[c]);
            D1c.set(1, 0, p[c]);
            D1c.set(1, 1, 1.0 / k);
            mmap.set(2 + c, mmap.get(1).elementMult(D1c));
        }
        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix fF = Mmap_forward_moment.mmap_forward_moment(mmap, moments);
        return new FitResult(mmap, fF.toArray1D(), fB);
    }

    private static double[] solveForwardOptimization(
            double[] F, double[][] qf, double[][] q0,
            double[] classWeights, double[] fbWeights, int k) {
        int n = k;
        double[][] H = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                H[i][j] = (i == j) ? 2.0 * classWeights[i] * fbWeights[0] : 0.0;
            }
        }

        double[] h = new double[n];
        for (int i = 0; i < n; i++) {
            h[i] = -2.0 * F[i] * classWeights[i] * fbWeights[0];
        }

        double[][] A = new double[0][n];
        double[] b = new double[0];
        double[][] Aeq = new double[0][n];
        double[] beq = new double[0];

        QPResult qpResult = solveQP(H, h, A, b, Aeq, beq);
        if (qpResult.success) {
            return qpResult.solution;
        }
        return F.clone();
    }

    private static QPResult solveQP(
            double[][] H, double[] h,
            double[][] A, double[] b,
            double[][] Aeq, double[] beq) {
        int n = h.length;

        Model.Builder builder = Model.getBuilder();
        Model.Variable[] vars = new Model.Variable[n];
        for (int i = 0; i < n; i++) {
            vars[i] = builder.addVariable().lb(0.0);
        }

        Model.Objective obj = builder.setObjective();
        for (int i = 0; i < n; i++) {
            obj.add(h[i], vars[i]);
        }
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                double coeff = (i == j) ? H[i][j] : H[i][j] + H[j][i];
                if (Math.abs(coeff) > 1e-12) {
                    if (i == j) obj.add(coeff / 2.0, vars[i], vars[j]);
                    else obj.add(coeff / 2.0, vars[i], vars[j]);
                }
            }
        }
        obj.minimize();

        for (int i = 0; i < A.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            for (int j = 0; j < n; j++) {
                if (Math.abs(A[i][j]) > 1e-12) {
                    ctr.add(A[i][j], vars[j]);
                }
            }
            ctr.leq(b[i]);
        }

        for (int i = 0; i < Aeq.length; i++) {
            Model.Constraint ctr = builder.addConstraint();
            for (int j = 0; j < n; j++) {
                if (Math.abs(Aeq[i][j]) > 1e-12) {
                    ctr.add(Aeq[i][j], vars[j]);
                }
            }
            ctr.eq(beq[i]);
        }

        Model model = builder.build();
        model.getParam().verbose = GlobalConstants.getVerbose() != VerboseLevel.SILENT;
        OSQP.Status status = model.solve();

        if (status == OSQP.Status.SOLVED) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = model.getSolution(vars[i].getIndex());
            }
            double objVal = 0.0;
            for (int i = 0; i < n; i++) {
                objVal += h[i] * x[i];
                for (int j = 0; j < n; j++) {
                    objVal += 0.5 * x[i] * H[i][j] * x[j];
                }
            }
            return new QPResult(x, objVal, true);
        } else {
            double[] fallback = new double[n];
            for (int i = 0; i < n; i++) fallback[i] = 1.0 / n;
            return new QPResult(fallback, Double.MAX_VALUE, false);
        }
    }
}
