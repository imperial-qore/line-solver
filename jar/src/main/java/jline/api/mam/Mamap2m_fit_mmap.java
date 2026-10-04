package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Markovian Arrival MAP with Marked arrivals MMAP-based fitting.
 */
public final class Mamap2m_fit_mmap {
    private Mamap2m_fit_mmap() {}

    public static MatrixCell mamap2m_fit_mmap(MatrixCell MMAP) {
        return mamap2m_fit_mmap(MMAP, new double[]{1.0, 1.0, 1.0});
    }

    public static MatrixCell mamap2m_fit_mmap(MatrixCell MMAP, double[] fbsWeights) {
        double M1 = Map_moment.map_moment(MMAP.get(0), MMAP.get(1), 1);
        double M2 = Map_moment.map_moment(MMAP.get(0), MMAP.get(1), 2);
        double M3 = Map_moment.map_moment(MMAP.get(0), MMAP.get(1), 3);
        double GAMMA = Map_gamma.map_gamma(MMAP);
        Matrix P = Mmap_pc.mmap_pc(MMAP);
        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix F = Mmap_forward_moment.mmap_forward_moment(MMAP, moments);
        Matrix B = Mmap_backward_moment.mmap_backward_moment(MMAP, moments);
        Matrix S = Mmap_sigma.mmap_sigma(MMAP);
        return mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S, fbsWeights);
    }

    private static MatrixCell mamap2m_fit(double M1, double M2, double M3, double GAMMA,
                                          Matrix P, Matrix F, Matrix B, Matrix S, double[] fbsWeights) {
        double gammatol = 1e-4;
        double degentol = 1e-8;
        int m = P.getNumCols();
        if (m > 2) {
            System.out.println("Fitting MAMAP(2,m): fitting F+B because m > 2");
            return Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P.toArray1D(), F.toArray1D(), B.toArray1D());
        }
        double cv2 = M2 / (M1 * M1) - 1.0;
        if (Math.abs(cv2) < degentol && Math.abs(GAMMA) < gammatol) {
            System.out.println("Fitting MAMAP(2,m): converting Poisson to second-order");
            double lambda = 1.0 / M1;
            double switchRate = lambda * 0.1;
            Matrix D0 = new Matrix(2, 2);
            D0.set(0, 0, -(lambda + switchRate));
            D0.set(0, 1, switchRate);
            D0.set(1, 0, switchRate);
            D0.set(1, 1, -(lambda + switchRate));
            Matrix D1 = new Matrix(2, 2);
            D1.set(0, 0, lambda * 0.5); D1.set(0, 1, lambda * 0.5);
            D1.set(1, 0, lambda * 0.5); D1.set(1, 1, lambda * 0.5);
            MatrixCell poissonMap = new MatrixCell(2);
            poissonMap.set(0, D0);
            poissonMap.set(1, D1);
            java.util.List<MatrixCell> amaps = java.util.Collections.singletonList(poissonMap);
            return handleMultipleAmaps(amaps, P, F, B, S, fbsWeights, GAMMA, degentol, m);
        }
        if (Math.abs(GAMMA) < gammatol) {
            System.out.println("Fitting MAMAP(2,m): fitting MAPH because gamma = " + GAMMA);
            if (fbsWeights[0] > fbsWeights[1]) {
                System.out.println("Converting PH-renewal to non-canonical form for forward moment fitting");
                return maph2m_fit_noncanonical(M1, M2, M3, P, F, fbsWeights);
            }
            return Maph2m_fit.maph2m_fit(M1, M2, M3, P, B);
        }
        jline.util.Pair<MatrixCell, java.util.List<MatrixCell>> amapsRet =
                Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, GAMMA);
        java.util.List<MatrixCell> amaps = amapsRet.getRight();
        if (amaps.size() == 1 && amaps.get(0).get(0).getNumRows() == 1) {
            System.out.println("Fitting MAMAP(2,m): fitting marked Poisson because the underlying process has one state");
            MatrixCell map = amaps.get(0);
            MatrixCell mmap = new MatrixCell(2 + m);
            mmap.set(0, map.get(0).copy());
            mmap.set(1, map.get(1).copy());
            for (int c = 0; c < m; c++) mmap.set(2 + c, mmap.get(1).scale(P.get(0, c)));
            return mmap;
        }
        return handleMultipleAmaps(amaps, P, F, B, S, fbsWeights, GAMMA, degentol, m);
    }

    private static MatrixCell handleMultipleAmaps(java.util.List<MatrixCell> amaps,
                                                  Matrix P, Matrix F, Matrix B, Matrix S, double[] fbsWeights,
                                                  double GAMMA, double degentol, int m) {
        double[] fbWeights = new double[]{fbsWeights[0], fbsWeights[1]};
        double[] fsWeights = new double[]{fbsWeights[0], fbsWeights[2]};
        double[] bsWeights = new double[]{fbsWeights[1], fbsWeights[2]};

        java.util.List<MatrixCell> mmaps = new java.util.ArrayList<MatrixCell>();
        java.util.List<Double> errors = new java.util.ArrayList<Double>();

        for (MatrixCell amap : amaps) {
            double h1 = -1.0 / amap.get(0).get(0, 0);
            double h2 = -1.0 / amap.get(0).get(1, 1);
            Matrix negD0Inv = amap.get(0).scale(-1.0).inv();
            Matrix transProb = negD0Inv.mult(amap.get(1));
            double r1 = h1 * amap.get(0).get(0, 1);
            double r2 = h2 * amap.get(1).get(1, 1);

            boolean degen = true;
            MatrixCell mmap;

            if (GAMMA > 0) {
                if (r1 < degentol || Math.abs(1 - r2) < degentol) {
                    throw new IllegalArgumentException("Should not happen for positive gamma");
                } else if (Math.abs(h2 - h1 * r2) < degentol) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(amap, P, F, S, null, fsWeights);
                } else if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(amap, P, B, S, null, bsWeights);
                } else if (Math.abs(1 - r1) < degentol) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(amap, P, F, S, null, fsWeights);
                } else if (r2 < degentol) {
                    mmap = maph2m_fit_multiclass_simple(amap, P, B);
                } else {
                    degen = false;
                    mmap = new MatrixCell(2 + m);
                }
            } else {
                if (Math.abs(1 - r2) < degentol) {
                    throw new IllegalArgumentException("Should not happen for negative gamma");
                } else if (Math.abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < degentol) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(amap, P, F, S, null, fsWeights);
                } else if (Math.abs(h1 - h2 + h2 * r1) < degentol) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(amap, P, B, S, null, bsWeights);
                } else if (r2 < degentol && Math.abs(1 - r1) < degentol) {
                    mmap = maph2m_fit_multiclass_simple(amap, P, B);
                } else if (r2 < degentol) {
                    if (fbsWeights[0] >= fbsWeights[1]) {
                        mmap = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(amap, P, F, S, null, fsWeights);
                    } else {
                        mmap = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(amap, P, B, S, null, bsWeights);
                    }
                } else {
                    degen = false;
                    mmap = new MatrixCell(2 + m);
                }
            }

            if (!degen) {
                if (fbsWeights[0] >= fbsWeights[2] && fbsWeights[1] >= fbsWeights[2]) {
                    double[] sliced = new double[]{fbWeights[0], fbWeights[1]};
                    Mamap2m_fit_fb_multiclass.FitResult retObj = Mamap2m_fit_fb_multiclass.mamap2m_fit_fb_multiclass(amap,
                            P.toArray1D(), F.toArray1D(), B.toArray1D(), null, sliced);
                    mmap = retObj.mmap;
                } else if (fbsWeights[0] >= fbsWeights[1]) {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_fs_multiclass(amap, P, F, S, null, fsWeights);
                } else {
                    mmap = Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(amap, P, B, S, null, bsWeights);
                }
            }

            mmaps.add(mmap);
            try {
                Matrix fittedF = Mmap_forward_moment.mmap_forward_moment(mmap, Matrix.ones(1, 1));
                Matrix fittedB = Mmap_backward_moment.mmap_backward_moment(mmap, Matrix.ones(1, 1));
                Matrix fittedS = Mmap_sigma.mmap_sigma(mmap);
                double error = 0.0;
                for (int c = 0; c < m; c++) {
                    double fError = (F.get(0, c) / fittedF.get(c, 0)) - 1.0;
                    double bError = (B.get(0, c) / fittedB.get(c, 0)) - 1.0;
                    error += fbsWeights[0] * fError * fError;
                    error += fbsWeights[1] * bError * bError;
                    if (c < S.getNumRows() && c < S.getNumCols() && c < fittedS.getNumRows() && c < fittedS.getNumCols()) {
                        double sError = (S.get(c, c) / fittedS.get(c, c)) - 1.0;
                        error += fbsWeights[2] * sError * sError;
                    }
                }
                errors.add(error);
            } catch (Exception e) {
                errors.add(Double.MAX_VALUE);
            }
        }

        int bestIndex = 0;
        double minError = errors.get(0);
        for (int i = 1; i < errors.size(); i++) {
            if (errors.get(i) < minError) { minError = errors.get(i); bestIndex = i; }
        }
        return mmaps.get(bestIndex);
    }

    private static MatrixCell maph2m_fit_multiclass_simple(MatrixCell aph, Matrix P, Matrix B) {
        int n = aph.get(0).getNumRows();
        int m = P.getNumCols();
        MatrixCell mmap = new MatrixCell(2 + m);
        mmap.set(0, aph.get(0).copy());
        mmap.set(1, Matrix.zeros(n, n));
        for (int c = 0; c < m; c++) mmap.set(2 + c, aph.get(1).scale(P.get(0, c)));
        return mmap;
    }

    private static MatrixCell maph2m_fit_noncanonical(double M1, double M2, double M3,
                                                      Matrix P, Matrix F, double[] fbsWeights) {
        int m = P.getNumCols();
        double scv = M2 / (M1 * M1) - 1.0;
        if (scv <= 1.0) {
            if (Math.abs(scv) < 1e-8) {
                double lambda = 1.0 / M1;
                double p = 0.7;
                Matrix D0 = new Matrix(2, 2);
                D0.set(0, 0, -lambda / p); D0.set(0, 1, lambda * (1 - p) / p);
                D0.set(1, 0, 0.0); D0.set(1, 1, -lambda);
                Matrix D1 = new Matrix(2, 2);
                D1.set(0, 0, lambda / p); D1.set(0, 1, 0.0);
                D1.set(1, 0, 0.0); D1.set(1, 1, lambda);
                Matrix T = new Matrix(2, 2);
                T.set(0, 0, 1.0); T.set(0, 1, 0.2);
                T.set(1, 0, 0.0); T.set(1, 1, 1.0);
                Matrix Tinv = T.inv();
                Matrix D0trans = Tinv.mult(D0).mult(T);
                Matrix D1trans = Tinv.mult(D1).mult(T);
                MatrixCell mmap = new MatrixCell(2 + m);
                mmap.set(0, D0trans);
                mmap.set(1, Matrix.zeros(2, 2));
                for (int c = 0; c < m; c++) mmap.set(2 + c, D1trans.scale(P.get(0, c)));
                return mmap;
            }
        }
        return Maph2m_fit.maph2m_fit(M1, M2, M3, P, F);
    }

    /** MAMAP 2m fit mmap algorithms. */
    public static final class Mamap2mFitMmapAlgo {}
}
