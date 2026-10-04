/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

/**
 * Modified-offered-load and pointwise-stationary approximations for a
 * time-varying multiserver system.
 *
 * <p>THE ONE IDEA. A stationary loss system with offered load a blocks with
 * probability B(s,a). In a time-varying system the question is WHICH LOAD goes
 * into that formula. PSA uses the instantaneous one, lambda(t)E[S]. MOL uses the
 * offered load of the corresponding INFINITE-SERVER system,
 * m(t) = int_0^Inf lambda(t-x)P(S&gt;x)dx, which is EXACT there and therefore
 * carries the time lag and the smoothing the finite-server system also has. The
 * difference between the two is precisely the lag: PSA peaks when the arrival
 * rate peaks, MOL peaks later, and the real system peaks later too.
 *
 * <p>WHAT TO EXPECT. Against the exact time-varying birth-death chain on a
 * sinusoidal rate, MOL cuts the mean RELATIVE error roughly threefold (0.13
 * against 0.44 at s = 100); it does not always win on ABSOLUTE error, which is
 * dominated by the peak of the cycle. Under constant input MOL is exact.
 *
 * <p>Port of MATLAB qsys_mtgs0_mol.m.
 *
 * <p>Reference: W. A. Massey, W. Whitt (1994). An analysis of the modified
 * offered load approximation for the nonstationary Erlang loss model. Annals of
 * Applied Probability 4(4), 1145-1160; W. Whitt (1991). Management Science
 * 37(3), 307-314.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mtgs0_mol {

    private Qsys_mtgs0_mol() {
    }

    /**
     * Erlang B by the recursion B_j = a B_(j-1)/(j + a B_(j-1)), which never
     * forms a^s/s! and so never overflows.
     *
     * @param s number of servers
     * @param a offered load in erlangs
     * @return the probability that all servers are busy
     */
    public static double erlangB(int s, double a) {
        double b = 1.0;
        for (int j = 1; j <= s; j++) {
            b = a * b / (j + a * b);
        }
        return b;
    }

    /**
     * Erlang C from the same recursion; 1 when the load saturates the servers.
     *
     * @param s number of servers
     * @param a offered load in erlangs
     * @return the probability that an arrival waits
     */
    public static double erlangC(int s, double a) {
        if (a >= s) {
            return 1.0;
        }
        double b = erlangB(s, a);
        double rho = a / s;
        return b / (1.0 - rho * (1.0 - b));
    }

    /**
     * @param lambdaFun   the arrival rate
     * @param serviceCcdf G^c(x) = P(S &gt; x)
     * @param ES          the mean service time
     * @param s           number of servers
     * @param tvals       times at which to evaluate
     * @param startTime   time the system started empty; -Inf assumes an infinite past
     * @param delay       use Erlang C rather than Erlang B
     * @return map with times, offeredLoad, instantLoad, probBlockMOL,
     *         probBlockPSA, meanBusyMOL and arrivalRate, all double[]
     */
    public static Map<String, double[]> qsys_mtgs0_mol(DoubleUnaryOperator lambdaFun,
                                                        DoubleUnaryOperator serviceCcdf, double ES,
                                                        int s, double[] tvals, double startTime,
                                                        boolean delay) {
        if (s < 1) {
            throw new RuntimeException("qsys_mtgs0_mol: the number of servers s must be at least 1");
        }
        QsysMtginfResult inf = Qsys_mtginf.qsys_mtginf(lambdaFun, serviceCcdf, ES, tvals, startTime,
                Double.NaN, null, Qsys_mtginf.DEFAULT_TOL, Qsys_mtginf.DEFAULT_PANELS, 1e12);
        int n = inf.times.length;
        double[] mol = new double[n];
        double[] psa = new double[n];
        double[] busy = new double[n];
        for (int i = 0; i < n; i++) {
            mol[i] = delay ? erlangC(s, inf.meanNumber[i]) : erlangB(s, inf.meanNumber[i]);
            psa[i] = delay ? erlangC(s, inf.offeredLoadPSA[i]) : erlangB(s, inf.offeredLoadPSA[i]);
            busy[i] = delay ? Math.min(inf.meanNumber[i], s) : inf.meanNumber[i] * (1.0 - mol[i]);
        }
        Map<String, double[]> res = new HashMap<String, double[]>();
        res.put("times", inf.times);
        res.put("offeredLoad", inf.meanNumber);
        res.put("instantLoad", inf.offeredLoadPSA);
        res.put("probBlockMOL", mol);
        res.put("probBlockPSA", psa);
        res.put("meanBusyMOL", busy);
        res.put("arrivalRate", inf.arrivalRate);
        return res;
    }
}
