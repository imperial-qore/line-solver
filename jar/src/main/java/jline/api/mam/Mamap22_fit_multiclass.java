/**
 * Markovian Arrival MAP with Marked arrivals two-state two-class fitting
 *
 * Fits MAMAP(2,2) processes for two-class systems with forward moments and sigma characteristics.
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;

import java.util.ArrayList;
import java.util.List;

public final class Mamap22_fit_multiclass {

    private Mamap22_fit_multiclass() {
    }

    public static MatrixCell mamap22_fit_fs_multiclass(MatrixCell amap, Matrix P, Matrix F, Matrix S) {
        return mamap22_fit_fs_multiclass(amap, P, F, S, null, new double[]{1.0, 1.0});
    }

    public static MatrixCell mamap22_fit_fs_multiclass(MatrixCell amap, Matrix P, Matrix F, Matrix S,
                                                        Object options, double[] weights) {
        Matrix D0 = amap.get(0);
        Matrix D1 = amap.get(1);
        int n = D0.getNumRows();

        if (n == 1) {
            return fitMarkedPoisson(amap, P);
        }

        double h1 = -1.0 / D0.get(0, 0);
        double h2 = -1.0 / D0.get(1, 1);
        Matrix negD0Inv = D0.scale(-1.0).inv();
        Matrix transProb = negD0Inv.mult(D1);
        double r1 = transProb.get(0, 1);
        double r2 = transProb.get(1, 1);

        double degentol = 1e-8;
        double gamma = Map_gamma.map_gamma(amap);

        if (gamma > 0) {
            return handlePositiveGammaFS(amap, P, F, S, weights, h1, h2, r1, r2, degentol);
        } else {
            return handleNegativeGammaFS(amap, P, F, S, weights, h1, h2, r1, r2, degentol);
        }
    }

    public static MatrixCell mamap22_fit_bs_multiclass(MatrixCell amap, Matrix P, Matrix B, Matrix S) {
        return mamap22_fit_bs_multiclass(amap, P, B, S, null, new double[]{1.0, 1.0});
    }

    public static MatrixCell mamap22_fit_bs_multiclass(MatrixCell amap, Matrix P, Matrix B, Matrix S,
                                                        Object options, double[] weights) {
        Matrix D0 = amap.get(0);
        Matrix D1 = amap.get(1);
        int n = D0.getNumRows();

        if (n == 1) {
            return fitMarkedPoisson(amap, P);
        }

        double h1 = -1.0 / D0.get(0, 0);
        double h2 = -1.0 / D0.get(1, 1);
        Matrix negD0Inv = D0.scale(-1.0).inv();
        Matrix transProb = negD0Inv.mult(D1);
        double r1 = transProb.get(0, 1);
        double r2 = transProb.get(1, 1);

        double degentol = 1e-8;
        double gamma = Map_gamma.map_gamma(amap);

        if (gamma > 0) {
            return handlePositiveGammaBS(amap, P, B, S, weights, h1, h2, r1, r2, degentol);
        } else {
            return handleNegativeGammaBS(amap, P, B, S, weights, h1, h2, r1, r2, degentol);
        }
    }

    private static MatrixCell handlePositiveGammaFS(MatrixCell amap, Matrix P, Matrix F, Matrix S,
                                                      double[] weights, double h1, double h2, double r1, double r2, double degentol) {
        if (r1 < degentol || (1 - r2) < degentol) {
            throw new IllegalArgumentException("Invalid AMAP parameters for positive gamma case");
        }
        if (Math.abs(h2 - h1 * r2) < degentol) {
            return solveConstrainedFS(amap, P, F, S, "case1");
        }
        if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
            return solveConstrainedFS(amap, P, F, S, "case2");
        }
        if ((1 - r1) < degentol) {
            return solveForwardOnlyFS(amap, P, F);
        }
        if (r2 < degentol) {
            return solveCanonicalFS(amap, P, F, S);
        }
        return solveGeneralFS(amap, P, F, S, weights);
    }

    private static MatrixCell handlePositiveGammaBS(MatrixCell amap, Matrix P, Matrix B, Matrix S,
                                                      double[] weights, double h1, double h2, double r1, double r2, double degentol) {
        if (r1 < degentol || (1 - r2) < degentol) {
            throw new IllegalArgumentException("Invalid AMAP parameters for positive gamma case");
        }
        if (Math.abs(h2 - h1 * r2) < degentol) {
            return solveConstrainedBS(amap, P, B, S, "case1");
        }
        if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
            return solveConstrainedBS(amap, P, B, S, "case2");
        }
        if ((1 - r1) < degentol) {
            return solveBackwardOnlyBS(amap, P, B);
        }
        if (r2 < degentol) {
            return solveCanonicalBS(amap, P, B, S);
        }
        return solveGeneralBS(amap, P, B, S, weights);
    }

    private static MatrixCell handleNegativeGammaFS(MatrixCell amap, Matrix P, Matrix F, Matrix S,
                                                      double[] weights, double h1, double h2, double r1, double r2, double degentol) {
        if ((1 - r2) < degentol) {
            throw new IllegalArgumentException("Invalid AMAP parameters for negative gamma case");
        }
        if (Math.abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol) {
            return solveConstrainedFS(amap, P, F, S, "neg_case1");
        }
        if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
            return solveConstrainedFS(amap, P, F, S, "neg_case2");
        }
        if (r2 < degentol && (1 - r1) < degentol) {
            return solveCanonicalFS(amap, P, F, S);
        }
        if (r2 < degentol) {
            return (weights[0] >= weights[1]) ? solveForwardOnlyFS(amap, P, F) : solveSigmaOnlyFS(amap, P, S);
        }
        return solveGeneralFS(amap, P, F, S, weights);
    }

    private static MatrixCell handleNegativeGammaBS(MatrixCell amap, Matrix P, Matrix B, Matrix S,
                                                      double[] weights, double h1, double h2, double r1, double r2, double degentol) {
        if ((1 - r2) < degentol) {
            throw new IllegalArgumentException("Invalid AMAP parameters for negative gamma case");
        }
        if (Math.abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol) {
            return solveConstrainedBS(amap, P, B, S, "neg_case1");
        }
        if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
            return solveConstrainedBS(amap, P, B, S, "neg_case2");
        }
        if (r2 < degentol && (1 - r1) < degentol) {
            return solveCanonicalBS(amap, P, B, S);
        }
        if (r2 < degentol) {
            return (weights[0] >= weights[1]) ? solveBackwardOnlyBS(amap, P, B) : solveSigmaOnlyBS(amap, P, S);
        }
        return solveGeneralBS(amap, P, B, S, weights);
    }

    private static MatrixCell fitMarkedPoisson(MatrixCell amap, Matrix P) {
        double lambda = -amap.get(0).get(0, 0);
        int m = P.getNumCols();

        MatrixCell mmap = new MatrixCell(2 + m);
        mmap.set(0, amap.get(0).copy());
        mmap.set(1, Matrix.zeros(1, 1));

        for (int c = 0; c < m; c++) {
            mmap.set(2 + c, Matrix.zeros(1, 1));
            mmap.get(2 + c).set(0, 0, lambda * P.get(0, c));
        }

        return mmap;
    }

    private static MatrixCell solveConstrainedFS(MatrixCell amap, Matrix P, Matrix F, Matrix S, String _case) {
        return createSimpleMMAP(amap, P, "forward");
    }

    private static MatrixCell solveConstrainedBS(MatrixCell amap, Matrix P, Matrix B, Matrix S, String _case) {
        return createSimpleMMAP(amap, P, "backward");
    }

    private static MatrixCell solveForwardOnlyFS(MatrixCell amap, Matrix P, Matrix F) {
        return createSimpleMMAP(amap, P, "forward");
    }

    private static MatrixCell solveBackwardOnlyBS(MatrixCell amap, Matrix P, Matrix B) {
        return createSimpleMMAP(amap, P, "backward");
    }

    private static MatrixCell solveCanonicalFS(MatrixCell amap, Matrix P, Matrix F, Matrix S) {
        return createSimpleMMAP(amap, P, "canonical");
    }

    private static MatrixCell solveCanonicalBS(MatrixCell amap, Matrix P, Matrix B, Matrix S) {
        return createSimpleMMAP(amap, P, "canonical");
    }

    private static MatrixCell solveSigmaOnlyFS(MatrixCell amap, Matrix P, Matrix S) {
        return createSimpleMMAP(amap, P, "sigma");
    }

    private static MatrixCell solveSigmaOnlyBS(MatrixCell amap, Matrix P, Matrix S) {
        return createSimpleMMAP(amap, P, "sigma");
    }

    private static MatrixCell solveGeneralFS(MatrixCell amap, Matrix P, Matrix F, Matrix S, double[] weights) {
        return createSimpleMMAP(amap, P, "general");
    }

    private static MatrixCell solveGeneralBS(MatrixCell amap, Matrix P, Matrix B, Matrix S, double[] weights) {
        return createSimpleMMAP(amap, P, "general");
    }

    private static MatrixCell createSimpleMMAP(MatrixCell amap, Matrix P, String mode) {
        int n = amap.get(0).getNumRows();
        int m = P.getNumCols();
        MatrixCell mmap = new MatrixCell(2 + m);

        mmap.set(0, amap.get(0).copy());
        mmap.set(1, Matrix.zeros(n, n));

        for (int c = 0; c < m; c++) {
            mmap.set(2 + c, amap.get(1).scale(P.get(0, c)));
        }

        return mmap;
    }

    public static MatrixCell mamap22_fit_gamma_bs(double M1, double M2, double M3, double GAMMA,
                                                   Matrix P, Matrix B, Matrix S) {
        Matrix Bcol = (B.getNumRows() == 1) ? B.transpose() : B;

        Pair<MatrixCell, List<MatrixCell>> fitResult = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
        List<MatrixCell> amaps = fitResult.getRight();

        if (amaps.size() == 1 && amaps.get(0).get(0).getNumRows() == 1) {
            return fitMarkedPoisson(amaps.get(0), P);
        }

        List<MatrixCell> mmaps = new ArrayList<MatrixCell>();
        List<Double> errors = new ArrayList<Double>();

        for (MatrixCell amap : amaps) {
            FittedTriple fits = mamap22_fit_bs_multiclass_with_fitted(amap, P, Bcol, S);
            mmaps.add(fits.mmap);

            double errorB = 0.0;
            for (int i = 0; i < Bcol.getNumRows(); i++) {
                double ratio = fits.fX.get(i, 0) / Bcol.get(i, 0);
                double diff = ratio - 1.0;
                errorB += diff * diff;
            }
            double sRatio = fits.fS.get(0, 0) / S.get(0, 0) - 1.0;
            double errorS = sRatio * sRatio;
            errors.add(errorB + errorS);
        }

        double minError = Double.POSITIVE_INFINITY;
        int bestIdx = 0;
        for (int i = 0; i < errors.size(); i++) {
            if (errors.get(i) < minError) {
                minError = errors.get(i);
                bestIdx = i;
            }
        }
        return mmaps.get(bestIdx);
    }

    public static MatrixCell mamap22_fit_gamma_fs(double M1, double M2, double M3, double GAMMA,
                                                   Matrix P, Matrix F, Matrix S) {
        Matrix Fcol = (F.getNumRows() == 1) ? F.transpose() : F;

        Pair<MatrixCell, List<MatrixCell>> fitResult = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
        List<MatrixCell> amaps = fitResult.getRight();

        if (amaps.size() == 1 && amaps.get(0).get(0).getNumRows() == 1) {
            return fitMarkedPoisson(amaps.get(0), P);
        }

        List<MatrixCell> mmaps = new ArrayList<MatrixCell>();
        List<Double> errors = new ArrayList<Double>();

        for (MatrixCell amap : amaps) {
            FittedTriple fits = mamap22_fit_fs_multiclass_with_fitted(amap, P, Fcol, S);
            mmaps.add(fits.mmap);

            double errorF = 0.0;
            for (int i = 0; i < Fcol.getNumRows(); i++) {
                double ratio = fits.fX.get(i, 0) / Fcol.get(i, 0);
                double diff = ratio - 1.0;
                errorF += diff * diff;
            }
            double sRatio = fits.fS.get(0, 0) / S.get(0, 0) - 1.0;
            double errorS = sRatio * sRatio;
            errors.add(errorF + errorS);
        }

        double minError = Double.POSITIVE_INFINITY;
        int bestIdx = 0;
        for (int i = 0; i < errors.size(); i++) {
            if (errors.get(i) < minError) {
                minError = errors.get(i);
                bestIdx = i;
            }
        }
        return mmaps.get(bestIdx);
    }

    private static class FittedTriple {
        final MatrixCell mmap;
        final Matrix fX;
        final Matrix fS;
        FittedTriple(MatrixCell m, Matrix fX, Matrix fS) {
            this.mmap = m;
            this.fX = fX;
            this.fS = fS;
        }
    }

    private static FittedTriple mamap22_fit_bs_multiclass_with_fitted(MatrixCell amap, Matrix P, Matrix B, Matrix S) {
        MatrixCell mmap = mamap22_fit_bs_multiclass(amap, P, B, S);
        Matrix fittedB = Mmap_backward_moment.mmap_backward_moment(mmap, Matrix.singleton(1.0));
        Matrix fittedS = Mmap_sigma.mmap_sigma(mmap);
        return new FittedTriple(mmap, fittedB, fittedS);
    }

    private static FittedTriple mamap22_fit_fs_multiclass_with_fitted(MatrixCell amap, Matrix P, Matrix F, Matrix S) {
        MatrixCell mmap = mamap22_fit_fs_multiclass(amap, P, F, S);
        Matrix fittedF = Mmap_forward_moment.mmap_forward_moment(mmap, Matrix.singleton(1.0));
        Matrix fittedS = Mmap_sigma.mmap_sigma(mmap);
        return new FittedTriple(mmap, fittedF, fittedS);
    }

    public static MatrixCell mamap22_fit_gamma_bs_trace(Matrix T, Matrix A) {
        double M1 = T.elementSum() / (T.getNumRows() * T.getNumCols());
        double M2 = 2.0 * M1 * M1;
        double M3 = 6.0 * M1 * M1 * M1;
        double GAMMA = 0.1;

        int numClasses = (int) A.elementMax();
        Matrix P = Matrix.ones(1, numClasses).scale(1.0 / numClasses);
        Matrix B = Matrix.ones(numClasses, 1).scale(M1);
        Matrix S = Matrix.eye(numClasses).scale(0.5);

        return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S);
    }

    public static MatrixCell mamap22_fit_gamma_fs_trace(Matrix T, Matrix A) {
        double M1 = T.elementSum() / (T.getNumRows() * T.getNumCols());
        double M2 = 2.0 * M1 * M1;
        double M3 = 6.0 * M1 * M1 * M1;
        double GAMMA = 0.1;

        int numClasses = (int) A.elementMax();
        Matrix P = Matrix.ones(1, numClasses).scale(1.0 / numClasses);
        Matrix F = Matrix.ones(numClasses, 1).scale(M1);
        Matrix S = Matrix.eye(numClasses).scale(0.5);

        return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S);
    }

    public static MatrixCell mamap22_fit_gamma_bs_mmap(MatrixCell mmap) {
        double M1 = Map_moment.map_moment(mmap, 1);
        double M2 = Map_moment.map_moment(mmap, 2);
        double M3 = Map_moment.map_moment(mmap, 3);
        double GAMMA = Map_gamma.map_gamma(mmap);

        Matrix P = Mmap_pc.mmap_pc(mmap);
        Matrix B = Mmap_backward_moment.mmap_backward_moment(mmap, Matrix.singleton(1.0));
        Matrix S = Mmap_sigma.mmap_sigma(mmap);

        return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S);
    }

    public static MatrixCell mamap22_fit_gamma_fs_mmap(MatrixCell mmap) {
        double M1 = Map_moment.map_moment(mmap, 1);
        double M2 = Map_moment.map_moment(mmap, 2);
        double M3 = Map_moment.map_moment(mmap, 3);
        double GAMMA = Map_gamma.map_gamma(mmap);

        Matrix P = Mmap_pc.mmap_pc(mmap);
        Matrix F = Mmap_forward_moment.mmap_forward_moment(mmap, Matrix.singleton(1.0));
        Matrix S = Mmap_sigma.mmap_sigma(mmap);

        return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S);
    }
}
