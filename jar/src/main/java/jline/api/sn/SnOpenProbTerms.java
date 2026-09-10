package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

import static jline.util.Maths.factln;

/**
 * Open-class contribution to an aggregate state probability at one station.
 *
 * <p>The mixed branch shared by {@code getProbAggr} and {@code getProbSysAggr}:
 * a Delay carries an independent Poisson per open class, a queueing station
 * carries the multinomial-geometric product form in its utilizations, and a
 * Source (EXT) carries neither, its population belonging to the environment.
 * The system probability is the sum of these terms over the stations, plus the
 * closed classes' multinomial, so both getters call this one method.
 *
 * <p>Mirrors {@code @SolverMVA/getProbAggr.m} and {@code @SolverFLD/getProbAggr.m}
 * and the native Python {@code _open_prob_aggr_terms}.
 *
 * @since LINE 3.0
 */
public final class SnOpenProbTerms {
    private SnOpenProbTerms() {}

    /**
     * Log-probability contributed by the open classes at one station.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @param Q mean queue lengths, stations by classes
     * @param U utilizations, stations by classes
     * @param nir per-class job counts at station ist, as a row
     * @param ist station index (0-based)
     * @return the log contribution, or negative infinity for an impossible state
     */
    public static double snOpenProbTerms(NetworkStruct sn, Matrix Q, Matrix U, Matrix nir, int ist) {
        SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
        if (sched == SchedStrategy.EXT) {
            return 0.0;
        }
        double total = 0.0;
        if (sched == SchedStrategy.INF) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (!Double.isInfinite(sn.njobs.get(r))) {
                    continue;
                }
                double nr = nir.get(0, r);
                if (Q.get(ist, r) > 0) {
                    total += nr * Math.log(Q.get(ist, r)) - Q.get(ist, r) - factln((int) nr);
                } else if (nr > 0) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
            return total;
        }
        double rhoTotal = 0.0;
        double nTotal = 0.0;
        for (int r = 0; r < sn.nclasses; r++) {
            if (!Double.isInfinite(sn.njobs.get(r))) {
                continue;
            }
            rhoTotal += U.get(ist, r);
            nTotal += nir.get(0, r);
        }
        if (!(rhoTotal < 1.0)) {
            return Double.NEGATIVE_INFINITY;
        }
        total += Math.log(1.0 - rhoTotal) + factln((int) nTotal);
        for (int r = 0; r < sn.nclasses; r++) {
            if (!Double.isInfinite(sn.njobs.get(r))) {
                continue;
            }
            double nr = nir.get(0, r);
            if (nr <= 0) {
                continue;
            }
            if (U.get(ist, r) <= 0) {
                return Double.NEGATIVE_INFINITY;
            }
            total += nr * Math.log(U.get(ist, r)) - factln((int) nr);
        }
        return total;
    }
}
