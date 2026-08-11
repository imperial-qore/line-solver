package jline.api.sn;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class SnGetDemandsChain {
    private SnGetDemandsChain() {}

    /**
     * Calculate new queueing network parameters after aggregating classes into chains.
     */
    public static Ret.snGetDemands snGetDemandsChain(NetworkStruct sn_in) {
        // see _kb/03-api-layer.md for rationale
        NetworkStruct sn = sn_in;
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix N = sn.njobs;

        Matrix scv = sn.scv.copy();
        scv.apply(Double.NaN, 1.0, "equal");

        Matrix ST = new Matrix(0, 0);
        sn.rates.divide(1.0, ST, false);
        ST.removeNaN();

        Matrix alpha = new Matrix(M, K);
        Matrix Vchain = new Matrix(M, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            if (sn.refclass.get(0, c) > -1) {
                for (int i = 0; i < M; i++) {
                    Matrix visits = sn.visits.get(Integer.valueOf(c));
                    double res = 0.0;
                    int iIdx = (int) sn.stationToStateful.get(i);
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        res += visits.get(iIdx, (int) inchain.get(0, col));
                    }
                    Vchain.set(i, c, res / visits.get((int) sn.stationToStateful.get((int) sn.refstat.get((int) inchain.value(), 0)),
                            (int) sn.refclass.get(0, c)));
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int k = (int) inchain.get(0, col);
                        alpha.set(i, k, alpha.get(i, k) + visits.get(iIdx, k) / res);
                    }
                }
            } else {
                for (int i = 0; i < M; i++) {
                    Matrix visits = sn.visits.get(Integer.valueOf(c));
                    double res1 = 0.0;
                    double res2 = 0.0;
                    int refIdx = (int) sn.stationToStateful.get((int) sn.refstat.get((int) inchain.value(), 0));
                    int iIdx = (int) sn.stationToStateful.get(i);
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        res1 += visits.get(iIdx, idx);
                        res2 += visits.get(refIdx, idx);
                    }
                    Vchain.set(i, c, res1 / res2);
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int k = (int) inchain.get(0, col);
                        alpha.set(i, k, alpha.get(i, k) + visits.get(iIdx, k) / res1);
                    }
                }
            }
        }

        Vchain.apply(GlobalConstants.Inf, 0.0, "equal");
        Vchain.apply(Double.NaN, 0.0, "equal");
        for (int c = 0; c < C; c++) {
            double vchainRef = Vchain.get((int) sn.refstat.get((int) sn.inchain.get(Integer.valueOf(c)).value(), 0), c);
            int[] colIdx = Vchain.getColIndexes();
            double[] nzVals = Vchain.getNonZeroValues();
            for (int i = colIdx[c]; i < colIdx[c + 1]; i++) {
                nzVals[i] /= vchainRef;
            }
        }
        alpha.apply(GlobalConstants.Inf, 0.0, "equal");
        alpha.apply(Double.NaN, 0.0, "equal");
        alpha.apply(GlobalConstants.Zero, 0.0, "less");

        Matrix Lchain = new Matrix(M, C);
        Matrix STchain = new Matrix(M, C);
        Matrix SCVchain = new Matrix(M, C);
        Matrix Nchain = new Matrix(1, C);
        Matrix refstatchain = new Matrix(C, 1);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(Integer.valueOf(c));
            boolean isOpenChain = false;
            double sum = 0.0;
            for (int col = 0; col < inchain.getNumCols(); col++) {
                sum += N.get((int) inchain.get(0, col));
                if (Utils.isInf(sum)) {
                    isOpenChain = true;
                    break;
                }
            }
            Nchain.set(0, c, sum);

            for (int i = 0; i < M; i++) {
                sum = 0.0;
                if (isOpenChain && (double) i == sn.refstat.get((int) inchain.value(), 0)) {
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        double rateValue = sn.rates.get(i, (int) inchain.get(0, col));
                        if (Double.isFinite(rateValue)) sum += rateValue;
                    }
                    STchain.set(i, c, 1.0 / sum);
                } else {
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        sum += ST.get(i, idx) * alpha.get(i, idx);
                    }
                    STchain.set(i, c, sum);
                }
                Lchain.set(i, c, Vchain.get(i, c) * STchain.get(i, c));
                double alphachain = 0.0;
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int idx = (int) inchain.get(0, col);
                    double scvValue = scv.get(i, idx);
                    if (Double.isFinite(scvValue)) alphachain += alpha.get(i, idx);
                }
                if (alphachain > GlobalConstants.Zero) {
                    sum = 0.0;
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        sum += scv.get(i, idx) * alpha.get(i, idx);
                    }
                    SCVchain.set(i, c, sum / alphachain);
                }
            }
            refstatchain.set(c, 0, sn.refstat.get((int) inchain.value(), 0));
            for (int col = 1; col < inchain.getNumCols(); col++) {
                int classIdx = (int) inchain.get(0, col);
                if (sn.refstat.get(classIdx, 0) != refstatchain.get(c, 0)) {
                    throw new RuntimeException("Class have different reference station");
                }
            }
        }
        Lchain.apply(GlobalConstants.Inf, 0.0, "equal");
        Lchain.apply(Double.NaN, 0.0, "equal");
        STchain.apply(GlobalConstants.Inf, 0.0, "equal");
        STchain.apply(Double.NaN, 0.0, "equal");
        return new Ret.snGetDemands(Lchain, STchain, Vchain, alpha, Nchain, SCVchain, refstatchain);
    }
}
