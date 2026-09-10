package jline.api.sn;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class SnHasProductFormNotHetFCFS {
    private SnHasProductFormNotHetFCFS() {}

    /**
     * Checks if the network satisfies product-form assumptions (does not have heterogeneous FCFS)
     *
     * @param sn NetworkStruct object for the queueing network model
     * @return true if the network satisfies the product-form assumptions
     */
    public static boolean snHasProductFormNotHetFCFS(NetworkStruct sn) {
        return snHasProductFormNotHetFCFS(sn, true);
    }

    /**
     * Checks if the network satisfies product-form assumptions (does not have heterogeneous FCFS)
     *
     * @param sn NetworkStruct object for the queueing network model
     * @param checkMeans also demand class-independent FCFS service means. Pass false only for
     *                   an algorithm that models class-dependent FCFS itself (ab, schmidt,
     *                   schmidt-ext), for which the exclusion is the whole point.
     * @return true if the network satisfies the product-form assumptions
     */
    public static boolean snHasProductFormNotHetFCFS(NetworkStruct sn, boolean checkMeans) {
        boolean ret = true;
        for (int i = 0; i < sn.sched.size(); i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            ret = ret && (s == SchedStrategy.INF
                    || s == SchedStrategy.PS
                    || s == SchedStrategy.FCFS
                    || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.EXT);
        }
        ret = ret && !SnHasPriorities.snHasPriorities(sn);
        ret = ret && !SnHasForkJoin.snHasForkJoin(sn);
        ret = ret && !SnHasSDRouting.snHasSDRouting(sn);

        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                for (int j = 0; j < sn.scv.getNumCols(); j++) {
                    double scv = sn.scv.get(i, j);
                    if (Double.isFinite(scv) && scv > 0) {
                        ret = ret && (scv > 1 - GlobalConstants.FineTol) && (scv < 1 + GlobalConstants.FineTol);
                    }
                }
                // BCMP type 1 asks the FCFS service to be exponential AND
                // class-independent, so the means must agree too: with unequal
                // means the product-form solve returns a wait proportional to
                // each class's own demand where FCFS makes every class wait
                // behind the same queue. The comparison is between CHAIN
                // service times (visit-weighted over the classes that actually
                // visit the station): a class that never visits cannot break
                // product form, and within-chain heterogeneity is invisible to
                // both the product-form and the qd branch, which deaggregate a
                // chain result proportionally to each class's own demand, so
                // only between-chain heterogeneity warrants the divert. LN
                // layers carry seeded rates for classes with zero visits, which
                // a raw per-class comparison mistakes for heterogeneity.
                double stmin = 0.0, stmax = 0.0;
                boolean anyserved = false;
                int iIdx = (int) sn.stationToStateful.get(i);
                for (int c = 0; checkMeans && c < sn.nchains; c++) {
                    Matrix visitsC = sn.visits.get(Integer.valueOf(c));
                    double num = 0.0, den = 0.0;
                    for (int j = 0; j < sn.rates.getNumCols(); j++) {
                        if (sn.chains.get(c, j) == 0) {
                            continue;
                        }
                        double w = visitsC.get(iIdx, j);
                        double rate = sn.rates.get(i, j);
                        if (w > GlobalConstants.Zero && Double.isFinite(rate) && rate > 0) {
                            num += w / rate;
                            den += w;
                        }
                    }
                    if (den > 0) {
                        double st = num / den;
                        if (!anyserved) {
                            stmin = st;
                            stmax = st;
                            anyserved = true;
                        } else {
                            stmin = Math.min(stmin, st);
                            stmax = Math.max(stmax, st);
                        }
                    }
                }
                if (anyserved && stmax - stmin > GlobalConstants.CoarseTol * stmax) {
                    ret = false;
                }
            }
        }

        return ret;
    }
}
