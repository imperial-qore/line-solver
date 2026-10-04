/**
 * @file MAMAP(2,m) fitting from moments, autocorrelation decay and class statistics
 *
 * Port of m3a mamap2m_fit.m. The dispatcher below reproduces the MATLAB one
 * exactly: underlying AMAP(2) from amap2_fit_gamma, then per-candidate marking
 * through the degenerate-form branches (F+S, B+S, MAPH) or the general F+B fit,
 * with the candidate of least weighted error returned.
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mamap2m_fit {
    private Mamap2m_fit() {}

    private static final double GAMMATOL = 1e-4;
    private static final double DEGENTOL = 1e-8;

    public static Matrix[] mamap2m_fit(double M1, double M2, double M3, double GAMMA,
                                       double[] P, double[] F, double[] B, Matrix S) {
        return mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S, new double[]{1.0, 1.0, 1.0});
    }

    /**
     * Fits a MAPH(2,m) or MAMAP(2,m) matching the inter-arrival moments and decay
     * rate, the class probabilities (always exactly) and, as far as the form
     * allows, the forward moments, backward moments and class transition
     * probabilities.
     *
     * @param fbsWeights weights of forward moments, backward moments and sigma
     */
    public static Matrix[] mamap2m_fit(double M1, double M2, double M3, double GAMMA,
                                       double[] P, double[] F, double[] B, Matrix S,
                                       double[] fbsWeights) {
        double[] fbWeights = {fbsWeights[0], fbsWeights[1]};
        double[] fsWeights = {fbsWeights[0], fbsWeights[2]};
        double[] bsWeights = {fbsWeights[1], fbsWeights[2]};
        int m = P.length;

        if (m > 2) {
            return toArray(mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B));
        }

        Matrix Pm = column(P);
        Matrix Fm = column(F);
        Matrix Bm = column(B);

        if (Math.abs(GAMMA) < GAMMATOL) {
            return toArray(Maph2m_fit.maph2m_fit(M1, M2, M3, Pm, Bm));
        }

        List<MatrixCell> maps = Amap2_fit_gamma.amap2_fitall_gamma(M1, M2, M3, GAMMA);
        if (maps.isEmpty()) {
            Pair<MatrixCell, List<MatrixCell>> fit = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
            maps = new ArrayList<MatrixCell>();
            maps.add(fit.getLeft());
        }

        if (maps.size() == 1 && maps.get(0).get(0).getNumRows() == 1) {
            // the underlying process is Poisson: perturb just above the exponential
            // to recover a second-order form, else return a marked Poisson process
            double M2a = M2 * (1 + 1e-4);
            double M3a = M3 * Math.pow(M2a / M2, 1.5);
            List<MatrixCell> maps2 = Amap2_fit_gamma.amap2_fitall_gamma(M1, M2a, M3a, GAMMA);
            if (!maps2.isEmpty()) {
                List<MatrixCell> normalized = new ArrayList<MatrixCell>();
                for (MatrixCell map : maps2) {
                    normalized.add(Map_normalize.map_normalize(map));
                }
                maps = normalized;
            } else {
                return markedPoisson(maps.get(0), P);
            }
        }

        List<MatrixCell> mmaps = new ArrayList<MatrixCell>();
        List<Double> errors = new ArrayList<Double>();

        for (MatrixCell map : maps) {
            Matrix D0 = map.get(0);
            Matrix D1 = map.get(1);
            double h1 = -1.0 / D0.get(0, 0);
            double h2 = -1.0 / D0.get(1, 1);
            double r1 = h1 * D0.get(0, 1);
            double r2 = h2 * D1.get(1, 1);

            MatrixCell fitted = null;
            boolean degen = true;
            if (GAMMA > 0) {
                if (r1 < DEGENTOL || (1 - r2) < DEGENTOL) {
                    throw new RuntimeException("Fitting MAMAP(2,m): should not happen");
                } else if (Math.abs(h2 - h1 * r2) < DEGENTOL) {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(map, Pm, Fm, S, null, fsWeights);
                } else if (Math.abs(h1 - h2 + h2 * r1) < DEGENTOL) {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(map, Pm, Bm, S, null, bsWeights);
                } else if ((1 - r1) < DEGENTOL) {
                    // non-canonical APH(2): only the forward moments can be fitted
                    fitted = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(map, Pm, Fm, S, null, fsWeights);
                } else if (r2 < DEGENTOL) {
                    // canonical APH(2)
                    fitted = Maph2m_fit.maph2m_fit_multiclass(map, Pm, Bm).getLeft();
                } else {
                    degen = false;
                }
            } else {
                if ((1 - r2) < DEGENTOL) {
                    throw new RuntimeException("Fitting MAMAP(2,m): should not happen");
                } else if (Math.abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < DEGENTOL) {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(map, Pm, Fm, S, null, fsWeights);
                } else if (Math.abs(h1 - h2 + h2 * r1) < DEGENTOL) {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(map, Pm, Bm, S, null, bsWeights);
                } else if (r2 < DEGENTOL && (1 - r1) < DEGENTOL) {
                    fitted = Maph2m_fit.maph2m_fit_multiclass(map, Pm, Bm).getLeft();
                } else if (r2 < DEGENTOL) {
                    if (fbsWeights[0] >= fbsWeights[1]) {
                        fitted = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(map, Pm, Fm, S, null, fsWeights);
                    } else {
                        fitted = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(map, Pm, Bm, S, null, bsWeights);
                    }
                } else {
                    degen = false;
                }
            }

            if (!degen) {
                if (fbsWeights[0] >= fbsWeights[2] && fbsWeights[1] >= fbsWeights[2]) {
                    Mamap2m_fit_fb_multiclass.FitResult r =
                        Mamap2m_fit_fb_multiclass.mamap2m_fit_fb_multiclass(map, P, F, B, null, fbWeights);
                    fitted = r.mmap;
                } else if (fbsWeights[0] >= fbsWeights[1]) {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(map, Pm, Fm, S, null, fsWeights);
                } else {
                    fitted = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(map, Pm, Bm, S, null, bsWeights);
                }
            }

            mmaps.add(fitted);
            errors.add(Double.valueOf(fittingError(fitted, F, B, S, fbsWeights)));
        }

        int best = 0;
        double bestErr = errors.get(0).doubleValue();
        for (int i = 1; i < errors.size(); i++) {
            if (errors.get(i).doubleValue() < bestErr) {
                bestErr = errors.get(i).doubleValue();
                best = i;
            }
        }
        return toArray(mmaps.get(best));
    }

    /**
     * MAMAP(2,m) fitted on the forward and backward moments only, used when the
     * number of classes exceeds two (mamap2m_fit_gamma_fb.m).
     */
    public static MatrixCell mamap2m_fit_gamma_fb(double M1, double M2, double M3, double GAMMA,
                                                  double[] P, double[] F, double[] B) {
        List<MatrixCell> maps = Amap2_fit_gamma.amap2_fitall_gamma(M1, M2, M3, GAMMA);
        if (maps.isEmpty()) {
            Pair<MatrixCell, List<MatrixCell>> fit = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
            maps = new ArrayList<MatrixCell>();
            maps.add(fit.getLeft());
        }
        if (maps.size() == 1 && maps.get(0).get(0).getNumRows() == 1) {
            Matrix[] poisson = markedPoisson(maps.get(0), P);
            MatrixCell out = new MatrixCell(poisson.length);
            for (int i = 0; i < poisson.length; i++) out.set(i, poisson[i]);
            return out;
        }

        MatrixCell best = null;
        double bestErr = Double.MAX_VALUE;
        for (MatrixCell map : maps) {
            Mamap2m_fit_fb_multiclass.FitResult r =
                Mamap2m_fit_fb_multiclass.mamap2m_fit_fb_multiclass(map, P, F, B);
            double err = 0.0;
            for (int c = 0; c < F.length; c++) {
                double e1 = r.feasibleForwardMoments[c] / F[c] - 1.0;
                double e2 = r.feasibleBackwardMoments[c] / B[c] - 1.0;
                err += e1 * e1 + e2 * e2;
            }
            if (err < bestErr) {
                bestErr = err;
                best = r.mmap;
            }
        }
        return best;
    }

    private static double fittingError(MatrixCell mmap, double[] F, double[] B, Matrix S,
                                       double[] fbsWeights) {
        Matrix fF = Mmap_forward_moment.mmap_forward_moment(mmap, Matrix.ones(1, 1));
        Matrix fB = Mmap_backward_moment.mmap_backward_moment(mmap, Matrix.ones(1, 1));
        Matrix fS = Mmap_sigma.mmap_sigma(mmap);
        double eF = F[0] / fF.get(0, 0) - 1.0;
        double eB = B[0] / fB.get(0, 0) - 1.0;
        double eS = S.get(0, 0) / fS.get(0, 0) - 1.0;
        return fbsWeights[0] * eF * eF + fbsWeights[1] * eB * eB + fbsWeights[2] * eS * eS;
    }

    private static Matrix[] markedPoisson(MatrixCell map, double[] P) {
        int m = P.length;
        Matrix[] out = new Matrix[2 + m];
        out[0] = map.get(0).copy();
        out[1] = map.get(1).copy();
        for (int c = 0; c < m; c++) {
            out[2 + c] = map.get(1).scale(P[c]);
        }
        return out;
    }

    private static Matrix column(double[] v) {
        Matrix M = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            M.set(i, 0, v[i]);
        }
        return M;
    }

    private static Matrix[] toArray(MatrixCell cell) {
        Matrix[] out = new Matrix[cell.size()];
        for (int i = 0; i < cell.size(); i++) {
            out[i] = cell.get(i);
        }
        return out;
    }
}
