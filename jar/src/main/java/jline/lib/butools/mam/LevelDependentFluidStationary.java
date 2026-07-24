/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 *
 * Stationary distribution and mean of first/second-order level-dependent fluid
 * queues, from the building blocks produced by SecondOrderLevelDependentFluidSolve.
 */
package jline.lib.butools.mam;

import jline.util.matrix.Matrix;

public final class LevelDependentFluidStationary {
    private LevelDependentFluidStationary() {}

    /**
     * Stationary distribution ('pdf','pdfd','cdf','cdfm') at the requested
     * points. Returns a (points.length x N) matrix of per-state values.
     */
    public static Matrix stationaryDistr(LevelDependentFluidSolution sol, double[] T, String what, double[] points) {
        int K = T.length;
        double[] Tarr = new double[K + 1];
        for (int k = 0; k < K; k++) Tarr[k + 1] = T[k];
        int N = sol.masses.get(0).getNumCols();
        boolean cumm = what.equals("cdf") || what.equals("cdfm");
        boolean cdfm = what.equals("cdfm");

        Matrix res = Matrix.zeros(points.length, N);
        for (int pi = 0; pi < points.length; pi++) {
            double p = points[pi];
            Matrix pres = Matrix.zeros(1, N);
            int k = 0;
            while (k < K && p >= Tarr[k]) {
                if (cumm) {
                    if (k > 0) {
                        Matrix[] ie = FluidTools.integExp2(sol.KF.get(k - 1), sol.KB.get(k - 1), Tarr[k] - Tarr[k - 1]);
                        pres = pres.add(sol.iniF.get(k - 1).mult(ie[0]).mult(sol.cloF.get(k - 1)))
                                .add(sol.iniB.get(k - 1).mult(ie[1]).mult(sol.cloB.get(k - 1)));
                    }
                    if (p > Tarr[k] || cdfm) pres = pres.add(sol.masses.get(k));
                }
                k++;
            }
            if (k == K && p == Tarr[k] && cdfm) pres = pres.add(sol.masses.get(k));
            int ki = k - 1;
            double prem = p - Tarr[ki];
            double Tk = Tarr[ki + 1] - Tarr[ki];
            if (what.equals("pdf")) {
                pres = sol.iniF.get(ki).mult(sol.KF.get(ki).scale(prem).expm()).mult(sol.cloF.get(ki))
                        .add(sol.iniB.get(ki).mult(sol.KB.get(ki).scale(Tk - prem).expm()).mult(sol.cloB.get(ki)));
            } else if (what.equals("pdfd")) {
                pres = sol.iniF.get(ki).mult(sol.KF.get(ki)).mult(sol.KF.get(ki).scale(prem).expm()).mult(sol.cloF.get(ki))
                        .sub(sol.iniB.get(ki).mult(sol.KB.get(ki)).mult(sol.KB.get(ki).scale(Tk - prem).expm()).mult(sol.cloB.get(ki)));
            } else {
                Matrix[] ie = FluidTools.integExp2(sol.KF.get(ki), sol.KB.get(ki), prem);
                pres = pres.add(sol.iniF.get(ki).mult(ie[0]).mult(sol.cloF.get(ki)))
                        .add(sol.iniB.get(ki).mult(sol.KB.get(ki).scale(Tk - prem).expm()).mult(ie[1]).mult(sol.cloB.get(ki)));
            }
            for (int j = 0; j < N; j++) res.set(pi, j, pres.get(0, j));
        }
        return res;
    }

    /** Stationary mean fluid level E[X] (closed form). */
    public static double stationaryMean(LevelDependentFluidSolution sol, double[] T) {
        int K = T.length;
        double[] Tarr = new double[K + 1];
        for (int k = 0; k < K; k++) Tarr[k + 1] = T[k];
        int N = sol.masses.get(0).getNumCols();
        Matrix h = Matrix.ones(N, 1);

        double res = 0.0;
        for (int j = 0; j <= K; j++) res += Tarr[j] * sol.masses.get(j).elementSum();
        for (int k = 0; k < K; k++) {
            double Tk = Tarr[k + 1] - Tarr[k];
            Matrix[] fF = FluidTools.expIntMoments(sol.KF.get(k), Tk);
            Matrix[] fB = FluidTools.expIntMoments(sol.KB.get(k), Tk);
            res += sol.iniF.get(k).mult(fF[0].scale(Tarr[k]).add(fF[1])).mult(sol.cloF.get(k)).mult(h).get(0, 0);
            res += sol.iniB.get(k).mult(fB[0].scale(Tarr[k + 1]).sub(fB[1])).mult(sol.cloB.get(k)).mult(h).get(0, 0);
        }
        return res;
    }
}
