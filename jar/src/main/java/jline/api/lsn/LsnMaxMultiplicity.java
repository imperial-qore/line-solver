/**
 * @file Load Sharing Network Maximum Multiplicity Computation
 *
 * Computes maximum multiplicity constraints for load sharing network (LSN) analysis.
 *
 * @since LINE 3.0
 */
package jline.api.lsn;

import jline.GlobalConstants;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.util.graph.DirectedGraph;
import jline.util.matrix.Matrix;

public final class LsnMaxMultiplicity {
    private LsnMaxMultiplicity() {}

    /**
     * Computes the maximum multiplicity that can be sustained by each task.
     *
     * @param lsn The layered network structure
     * @return Matrix of maximum multiplicities for each task
     */
    public static Matrix lsnMaxMultiplicity(LayeredNetworkStruct lsn) {
        Matrix ag = new Matrix(lsn.dag);
        int n = lsn.dag.getNumRows();
        for (int ist = 0; ist < n; ist++) {
            for (int jst = 0; jst < n; jst++) {
                if (lsn.dag.get(ist, jst) > 0) {
                    ag.set(ist, jst, 1.0);
                }
            }
        }
        Matrix mult = new Matrix(lsn.mult);
        Matrix type = lsn.type;
        Matrix isref = lsn.isref;

        Matrix order = DirectedGraph.kahn(ag);

        Matrix inflow = new Matrix(n, 1);
        for (int ist = 0; ist < n; ist++) {
            if (type.get(ist) == (double) LayeredNetworkElement.TASK && isref.get(ist) != 0.0) {
                inflow.set(ist, 0, mult.get(ist));
            } else if (type.get(ist) == (double) LayeredNetworkElement.ENTRY
                    && lsn.arrival != null && lsn.arrival.get(ist) != null) {
                // an open arrival needs one thread of the entry's parent task, and is
                // the only inflow source when no reference task exists
                // (lsn_max_multiplicity.m:47-55)
                inflow.set(ist, 0, 1.0);
            }
        }

        Matrix outflow = new Matrix(n, 1);

        if (mult.length() < n) {
            int cur_end = mult.length();
            mult.expandMatrix(1, n, n);
            for (int ist = cur_end; ist < n; ist++) {
                mult.set(ist, GlobalConstants.Inf);
            }
        }

        for (int k = 0; k < n; k++) {
            int ist = (int) order.get(k);
            double inflowVal = inflow.get(ist);
            double multVal = mult.get(ist);
            boolean isSetupTask = lsn.hassetup != null && lsn.hassetup.length() > ist
                    && lsn.hassetup.get(0, ist) == 1.0;
            if (isSetupTask && inflowVal > 0) {
                // see _kb/03-api-layer.md for rationale (fj/ and lsn/ additions)
                outflow.set(ist, 0, multVal);
            } else {
                outflow.set(ist, 0, Math.min(inflowVal, multVal));
            }

            for (int jst = 0; jst < n; jst++) {
                if (jst != ist && ag.get(ist, jst) != 0.0) {
                    inflow.set(jst, inflow.get(jst) + outflow.get(ist));
                }
            }
        }

        for (int ist = 0; ist < n; ist++) {
            if (type.get(ist) == (double) LayeredNetworkElement.TASK
                    && mult.get(ist) == GlobalConstants.Inf
                    && isref.get(ist) == 0.0) {
                outflow.set(ist, GlobalConstants.Inf);
            }
        }

        return outflow;
    }
}
