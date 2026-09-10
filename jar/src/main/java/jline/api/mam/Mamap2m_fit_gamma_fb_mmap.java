/**
 * @file Markovian Arrival MAP with Marked arrivals gamma forward-backward MMAP fitting
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import jline.util.Triple;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mamap2m_fit_gamma_fb_mmap {
    private Mamap2m_fit_gamma_fb_mmap() {}

    /**
     * Fits a second-order acyclic MMAP[m] to match the characteristics of the input MMAP.
     */
    public static MatrixCell mamap2m_fit_gamma_fb_mmap(MatrixCell mmap) {
        double M1 = Map_moment.map_moment(mmap.get(0), mmap.get(1), 1);
        double M2 = Map_moment.map_moment(mmap.get(0), mmap.get(1), 2);
        double M3 = Map_moment.map_moment(mmap.get(0), mmap.get(1), 3);
        double GAMMA = Map_gamma.map_gamma(mmap);

        Matrix P = Mmap_pc.mmap_pc(mmap);
        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix F = Mmap_forward_moment.mmap_forward_moment(mmap, moments);
        Matrix B = Mmap_backward_moment.mmap_backward_moment(mmap, moments);

        return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P.toArray1D(), F.toArray1D(), B.toArray1D());
    }

    /**
     * Computes the second-order MAMAP[m] fitting the given moments.
     */
    public static MatrixCell mamap2m_fit_gamma_fb(double M1, double M2, double M3, double GAMMA,
                                                  double[] P, double[] F, double[] B) {
        double[] Fcol = new double[F.length];
        System.arraycopy(F, 0, Fcol, 0, F.length);
        double[] Bcol = new double[B.length];
        System.arraycopy(B, 0, Bcol, 0, B.length);

        Pair<MatrixCell, List<MatrixCell>> amapResult = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
        List<MatrixCell> MAPS = amapResult.getSecond();

        if (MAPS.size() == 1 && MAPS.get(0).get(0).getNumRows() == 1) {
            MatrixCell MAP = MAPS.get(0);
            int m = P.length;

            MatrixCell result = new MatrixCell(2 + m);
            result.set(0, MAP.get(0).copy());
            result.set(1, MAP.get(1).copy());

            for (int c = 0; c < m; c++) {
                result.set(2 + c, result.get(1).scale(P[c]));
            }

            return result;
        }

        if (MAPS.isEmpty()) {
            return createFallbackMMAP(M1, P);
        }

        List<MatrixCell> MMAPS = new ArrayList<MatrixCell>();
        List<Double> ERRORS = new ArrayList<Double>();

        for (int j = 0; j < MAPS.size(); j++) {
            try {
                Mamap2m_fit_fb_multiclass.FitResult fitResult =
                        Mamap2m_fit_fb_multiclass.mamap2m_fit_fb_multiclass(MAPS.get(j), P, Fcol, Bcol);
                MatrixCell fittedMMAP = fitResult.mmap;
                double[] fF = fitResult.feasibleForwardMoments;
                double[] fB = fitResult.feasibleBackwardMoments;
                MMAPS.add(fittedMMAP);

                double forwardError = computeRelativeError(fF, Fcol);
                double backwardError = computeRelativeError(fB, Bcol);
                ERRORS.add(Double.valueOf(forwardError + backwardError));
            } catch (Exception e) {
                MMAPS.add(createFallbackMMAP(M1, P));
                ERRORS.add(Double.valueOf(Double.MAX_VALUE));
            }
        }

        int bestIndex = 0;
        double bestErr = ERRORS.get(0).doubleValue();
        for (int i = 1; i < ERRORS.size(); i++) {
            if (ERRORS.get(i).doubleValue() < bestErr) {
                bestErr = ERRORS.get(i).doubleValue();
                bestIndex = i;
            }
        }
        return MMAPS.get(bestIndex);
    }

    /**
     * Creates a fallback MMAP when fitting fails.
     */
    private static MatrixCell createFallbackMMAP(double M1, double[] P) {
        int m = P.length;
        MatrixCell result = new MatrixCell(2 + m);

        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -1.0 / M1);
        result.set(0, d0);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, 1.0 / M1);
        result.set(1, d1);

        for (int c = 0; c < m; c++) {
            result.set(2 + c, result.get(1).scale(P[c]));
        }

        return result;
    }

    /**
     * Computes relative error between fitted and target moments.
     */
    private static double computeRelativeError(double[] fitted, double[] target) {
        double error = 0.0;
        for (int i = 0; i < fitted.length; i++) {
            if (target[i] != 0.0) {
                double relError = (fitted[i] / target[i] - 1.0);
                error += relError * relError;
            }
        }
        return error;
    }
}
