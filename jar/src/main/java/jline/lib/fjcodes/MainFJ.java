/**
 * Main orchestrator for FJ_codes Fork-Join percentile analysis
 *
 * Implements the approximation method from "Beyond the Mean in Fork-Join Queues:
 * Efficient Approximation for Response-Time Tails" (IFIP Performance 2015).
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import jline.util.matrix.Matrix;

public final class MainFJ {
    private MainFJ() {}

    /**
     * Result for a specific K value
     */
    public static final class FJPercentileResult {
        public final int K;                     // Number of parallel queues
        public final double[] percentiles;      // Percentile levels (0-100 scale)
        public final double[] RTp;              // Response time percentile values

        public FJPercentileResult(int K, double[] percentiles, double[] RTp) {
            this.K = K;
            this.percentiles = percentiles;
            this.RTp = RTp;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof FJPercentileResult)) return false;
            FJPercentileResult that = (FJPercentileResult) o;
            return K == that.K && java.util.Arrays.equals(percentiles, that.percentiles)
                    && java.util.Arrays.equals(RTp, that.RTp);
        }

        @Override
        public int hashCode() {
            int r = Objects.hash(K);
            r = 31 * r + java.util.Arrays.hashCode(percentiles);
            r = 31 * r + java.util.Arrays.hashCode(RTp);
            return r;
        }
    }

    public static List<FJPercentileResult> mainFJ(FJArrival arrival, FJService service, double[] pers, int[] K) {
        return mainFJ(arrival, service, pers, K, new int[]{100}, "NARE");
    }

    public static List<FJPercentileResult> mainFJ(FJArrival arrival, FJService service, double[] pers, int[] K, int[] Cs) {
        return mainFJ(arrival, service, pers, K, Cs, "NARE");
    }

    /**
     * Compute response time percentiles for K-node Fork-Join queue.
     */
    public static List<FJPercentileResult> mainFJ(FJArrival arrival, FJService service, double[] pers,
                                                  int[] K, int[] Cs, String tMode) {
        // Check stability condition
        double load = arrival.lambda / service.mu;
        if (load >= 1.0) {
            throw new IllegalArgumentException("System not stable: mean arrival rate " + arrival.lambda
                    + " >= mean service rate " + service.mu + " (load = " + load + ")");
        }

        // Use the last C value from Cs array
        int C = Cs[Cs.length - 1];

        // Compute percentiles for K=1 (exact)
        Matrix percentileRT_1 = ReturnRT1.returnRT1(arrival, service, pers);

        // Compute percentiles for K=2 (approximation)
        Matrix percentileRT_2 = ReturnRT2.returnRT2(arrival, service, pers, C, tMode);

        List<FJPercentileResult> results = new ArrayList<FJPercentileResult>();

        for (int k : K) {
            double[] percentilesScaled = new double[pers.length];
            for (int i = 0; i < pers.length; i++) {
                percentilesScaled[i] = pers[i] * 100.0;
            }
            double[] RTp = new double[pers.length];

            for (int p = 0; p < pers.length; p++) {
                if (k == 1) {
                    RTp[p] = percentileRT_1.get(p, 1);
                } else if (k == 2) {
                    RTp[p] = percentileRT_2.get(p, 1);
                } else {
                    double rt1 = percentileRT_1.get(p, 1);
                    double rt2 = percentileRT_2.get(p, 1);
                    RTp[p] = rt1 + (rt2 - rt1) * Math.log((double) k) / Math.log(2.0);
                }
            }

            results.add(new FJPercentileResult(k, percentilesScaled, RTp));
        }

        return results;
    }

    /**
     * Convenience overload for single K value
     */
    public static FJPercentileResult mainFJ(FJArrival arrival, FJService service, double[] pers, int K) {
        return mainFJ(arrival, service, pers, K, 100, "NARE");
    }

    public static FJPercentileResult mainFJ(FJArrival arrival, FJService service, double[] pers, int K, int C) {
        return mainFJ(arrival, service, pers, K, C, "NARE");
    }

    public static FJPercentileResult mainFJ(FJArrival arrival, FJService service, double[] pers,
                                            int K, int C, String tMode) {
        return mainFJ(arrival, service, pers, new int[]{K}, new int[]{C}, tMode).get(0);
    }
}
