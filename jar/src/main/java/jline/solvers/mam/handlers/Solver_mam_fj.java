/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import jline.solvers.mam.MAMFJResult;

import java.util.ArrayList;
import java.util.List;

import jline.api.fj.FJConvert;
import jline.api.fj.FJInfo;
import jline.api.fj.FJValidation;
import jline.lang.NetworkStruct;
import jline.lib.fjcodes.FJArrival;
import jline.lib.fjcodes.MainFJ.FJPercentileResult;
import jline.lib.fjcodes.FJService;
import jline.lib.fjcodes.MainFJ;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Solver_mam_fj {
    private Solver_mam_fj() {}

    public static MAMFJResult solver_mam_fj(NetworkStruct sn) {
        return solver_mam_fj(sn, new double[]{0.50, 0.90, 0.95, 0.99}, 100, "NARE");
    }

    public static MAMFJResult solver_mam_fj(NetworkStruct sn, double[] percentiles) {
        return solver_mam_fj(sn, percentiles, 100, "NARE");
    }

    public static MAMFJResult solver_mam_fj(NetworkStruct sn, double[] percentiles, int C) {
        return solver_mam_fj(sn, percentiles, C, "NARE");
    }

    /**
     * Solve Fork-Join network using FJ_codes algorithm.
     */
    public static MAMFJResult solver_mam_fj(NetworkStruct sn, double[] percentiles, int C, String tMode) {
        Pair<Boolean, FJInfo> homogCheck = FJValidation.isHomogeneous(sn);
        boolean isValid = homogCheck.getLeft();
        FJInfo fjInfo = homogCheck.getRight();
        if (!isValid || fjInfo == null) {
            throw new IllegalArgumentException(
                    "Network is not a valid Fork-Join topology. "
                            + "Required: Source -> Fork -> K Queues -> Join -> Sink");
        }

        Pair<java.util.List<FJArrival>, java.util.List<FJService>> params = FJConvert.extractFJParams(sn, fjInfo);
        FJArrival[] arrivals = params.getLeft().toArray(new FJArrival[0]);
        FJService[] services = params.getRight().toArray(new FJService[0]);

        List<FJPercentileResult> percentileResults = new ArrayList<FJPercentileResult>();

        // Highest requested percentile (0-1 probability scale), used to size the
        // QBD truncation.
        double pMax = 0.0;
        for (double p : percentiles) {
            if (p > pMax) pMax = p;
        }

        for (int r = 0; r < sn.nclasses; r++) {
            try {
                // see _kb/06-solver-catalog.md for rationale
                int Cr = C;
                if (Cr <= 0) {
                    double rho = arrivals[r].lambda / services[r].mu;
                    if (!(rho > 0.0) || rho >= 1.0) {
                        Cr = 100;  // saturated/unknown: fall back to conservative
                    } else {
                        double tail = Math.max(1e-6, (1.0 - pMax) / 10.0);
                        int est = (int) Math.ceil(Math.log(tail) / Math.log(rho));
                        Cr = Math.max(20, Math.min(100, est));
                    }
                }
                FJPercentileResult result = MainFJ.mainFJ(arrivals[r], services[r], percentiles, fjInfo.K, Cr, tMode);
                percentileResults.add(result);
            } catch (Exception e) {
                throw new RuntimeException("FJ_codes failed for class " + r + ": " + e.getMessage(), e);
            }
        }

        double[] meanRT = new double[sn.nclasses];
        for (int r = 0; r < sn.nclasses; r++) {
            FJPercentileResult result = percentileResults.get(r);
            double sum = 0.0;
            for (int i = 0; i < result.RTp.length - 1; i++) {
                double dp = result.percentiles[i + 1] - result.percentiles[i];
                sum += (result.RTp[i] + result.RTp[i + 1]) / 2.0 * dp / 100.0;
            }
            meanRT[r] = sum;
        }

        Matrix QN = new Matrix(sn.nstations, sn.nclasses);
        Matrix RN = new Matrix(sn.nstations, sn.nclasses);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix UN = new Matrix(sn.nstations, sn.nclasses);
        Matrix XN = new Matrix(1, sn.nclasses);

        for (int r = 0; r < sn.nclasses; r++) {
            double lambda = arrivals[r].lambda;
            double mu = services[r].mu;

            XN.set(0, r, lambda);

            for (int qIdx : fjInfo.queueIndices) {
                int station = (int) sn.nodeToStation.get(qIdx);

                RN.set(station, r, meanRT[r] / fjInfo.K);
                TN.set(station, r, lambda / fjInfo.K);
                QN.set(station, r, lambda * meanRT[r] / fjInfo.K);
                UN.set(station, r, lambda / (fjInfo.K * mu));
            }
        }

        return new MAMFJResult(QN, UN, RN, TN, XN, percentileResults);
    }

    /**
     * Interpolate percentile from stored results.
     */
    public static double interpolatePercentile(double[] storedPercentiles, double[] storedValues, double targetPercentile) {
        for (int i = 0; i < storedPercentiles.length - 1; i++) {
            if (targetPercentile >= storedPercentiles[i] && targetPercentile <= storedPercentiles[i + 1]) {
                double t = (targetPercentile - storedPercentiles[i])
                        / (storedPercentiles[i + 1] - storedPercentiles[i]);
                return storedValues[i] + t * (storedValues[i + 1] - storedValues[i]);
            }
        }

        if (targetPercentile < storedPercentiles[0]) {
            return storedValues[0];
        } else {
            return storedValues[storedValues.length - 1];
        }
    }
}
