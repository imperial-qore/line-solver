/**
 * @file Convert Place throughputs from firing events to tokens
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

/**
 * Place throughput, arrival rate and response time in tokens.
 *
 * A Place is a station and a token is the job it holds, so a firing that
 * consumes two tokens is two departures, not one. The CTMC and SSA analyzers
 * count firing events instead, which for unit arc multiplicities is the same
 * number and for weighted arcs is not: the reported throughput is then not a
 * token rate, and QLen over it is not a sojourn time.
 *
 * This function rescales the Place rows to tokens:
 *
 * <ul>
 * <li>TN(p,k) tokens consumed from the Place per unit time</li>
 * <li>AN(p,k) tokens produced into the Place per unit time</li>
 * <li>RN(p,k) QN(p,k) / TN(p,k), Little's law over the Place</li>
 * </ul>
 *
 * Rows that do not belong to a Place are left untouched, so a mixed
 * Queue/Place model keeps its queueing metrics. When the firing rates cannot be
 * recovered from the throughputs the inputs are left unchanged rather than
 * replaced by a guess.
 *
 * Port of matlab/src/api/sn/sn_pn_avg_rates.m
 */
public final class SnPnAvgRates {
    private SnPnAvgRates() {}

    /**
     * Rescales the Place rows of TN, AN and RN to tokens, in place.
     *
     * @param sn network structure
     * @param QN average queue lengths, i.e. mean token counts at the Places
     * @param TN average throughputs at stations, counting firing events
     * @param AN average arrival rates at stations, or null when not yet computed
     * @param RN average response times at stations, or null when not needed
     */
    public static void snPnAvgRates(NetworkStruct sn, Matrix QN, Matrix TN, Matrix AN, Matrix RN) {
        if (TN == null || TN.isEmpty()) {
            return;
        }
        boolean hasPlace = false;
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Place) {
                hasPlace = true;
                break;
            }
        }
        if (!hasPlace) {
            return;
        }

        // The analyzers hand over event counts, hence the false.
        SnPnFiringRates.Ret ret = SnPnFiringRates.snPnFiringRates(sn, TN, false);
        if (ret.rates == null) {
            return;
        }

        int R = sn.nclasses;
        for (int pp = 0; pp < ret.placeNodes.size(); pp++) {
            int ist = (int) sn.nodeToStation.get(ret.placeNodes.get(pp));
            if (ist < 0) {
                continue;
            }
            for (int k = 0; k < R; k++) {
                double tput = 0.0;
                double arvr = 0.0;
                for (int mm = 0; mm < ret.rates.getNumRows(); mm++) {
                    double x = ret.rates.get(mm, 0);
                    tput += ret.consumed[mm][pp][k] * x;
                    arvr += ret.produced[mm][pp][k] * x;
                }
                TN.set(ist, k, tput);
                if (AN != null && !AN.isEmpty()) {
                    AN.set(ist, k, arvr);
                }
                if (RN != null && !RN.isEmpty()) {
                    if (tput > 0) {
                        RN.set(ist, k, QN.get(ist, k) / tput);
                    } else {
                        RN.set(ist, k, 0.0);
                    }
                }
            }
        }
    }
}
