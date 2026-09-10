/**
 * @file NPFQN Deterministic (round-robin) Traffic Split Degrees
 *
 * Computes, for every station-class, the degree of the deterministic
 * round-robin split its departure stream undergoes, which the QNA/MNA traffic
 * equations use in place of the Markovian splitting formula.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.NodeParam;
import jline.GlobalConstants;
import jline.lang.constant.RoutingStrategy;
import jline.util.matrix.Matrix;

public final class Npfqn_traffic_split_rr {
    private Npfqn_traffic_split_rr() {}

    /**
     * Deterministic (round-robin) split degree of the departure stream of each
     * station-class. An entry k &gt; 1 means that the class-r departures of that
     * station are dispatched one-in-k by a round-robin node, so that a
     * downstream flow carrying a fraction p of them is the k-fold convolution
     * thinned with probability q = k*p and has SCV 1+p*(d2-k), against the
     * Markovian 1+p*(d2-1). An entry of 1 marks an ordinary probabilistic split.
     *
     * Only two topologies admit the deterministic rule: the station dispatches
     * round-robin itself, or it feeds with probability one a router that does
     * and whose pointer no other flow advances. Anything else falls back to 1.
     *
     * @param sn network structure
     * @return (M,K) matrix of split degrees
     */
    public static Matrix npfqn_traffic_split_rr(NetworkStruct sn) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix kRR = new Matrix(M, K);
        kRR.fill(1.0);
        if (sn.routing == null || sn.routing.isEmpty()) {
            return kRR;
        }
        boolean anyRR = false;
        for (Node nd : sn.routing.keySet()) {
            for (RoutingStrategy rs : sn.routing.get(nd).values()) {
                if (rs == RoutingStrategy.RROBIN) {
                    anyRR = true;
                    break;
                }
            }
            if (anyRR) {
                break;
            }
        }
        if (!anyRR) {
            return kRR;
        }
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int r = 0; r < K; r++) {
                if (routingOf(sn, ind, r) == RoutingStrategy.RROBIN) {
                    kRR.set(ist, r, rrDegree(sn, ind, r));
                    continue;
                }
                if (sn.rtnodes == null || sn.rtnodes.isEmpty()) {
                    continue;
                }
                // a sure transition into a router that dispatches round-robin
                int dest = -1;
                int ndest = 0;
                for (int col = 0; col < sn.rtnodes.getNumCols(); col++) {
                    if (sn.rtnodes.get(ind * K + r, col) > 0) {
                        dest = col;
                        ndest++;
                    }
                }
                if (ndest != 1) {
                    continue;
                }
                int jnd = dest / K;
                int s = dest - jnd * K;
                if (sn.isstation.get(jnd) > 0 || routingOf(sn, jnd, s) != RoutingStrategy.RROBIN) {
                    continue;
                }
                if (sn.rtnodes.get(ind * K + r, dest) < 1 - GlobalConstants.FineTol) {
                    continue;
                }
                // the round-robin pointer must be advanced by this stream alone
                int nfeed = 0;
                for (int row = 0; row < sn.rtnodes.getNumRows(); row++) {
                    if (sn.rtnodes.get(row, dest) > 0) {
                        nfeed++;
                    }
                }
                if (nfeed != 1) {
                    continue;
                }
                kRR.set(ist, r, rrDegree(sn, jnd, s));
            }
        }
        return kRR;
    }

    private static RoutingStrategy routingOf(NetworkStruct sn, int ind, int r) {
        Node nd = sn.nodes.get(ind);
        JobClass jc = sn.jobclasses.get(r);
        if (sn.routing.containsKey(nd) && sn.routing.get(nd).containsKey(jc)) {
            return sn.routing.get(nd).get(jc);
        }
        return null;
    }

    /** Number of destinations the round-robin pointer of (ind,r) cycles through. */
    private static double rrDegree(NetworkStruct sn, int ind, int r) {
        Node nd = sn.nodes.get(ind);
        JobClass jc = sn.jobclasses.get(r);
        if (sn.nodeparam != null && sn.nodeparam.containsKey(nd)) {
            NodeParam np = sn.nodeparam.get(nd);
            if (np != null && np.outlinks != null && np.outlinks.containsKey(jc)) {
                Matrix ol = np.outlinks.get(jc);
                if (ol != null && ol.length() > 0) {
                    return atLeastOne(ol.length());
                }
            }
        }
        if (sn.connmatrix != null && ind < sn.connmatrix.getNumRows()) {
            int k = 0;
            for (int j = 0; j < sn.connmatrix.getNumCols(); j++) {
                if (sn.connmatrix.get(ind, j) > 0) {
                    k++;
                }
            }
            return atLeastOne(k);
        }
        return 1.0;
    }

    private static double atLeastOne(int k) {
        return k < 1 ? 1.0 : (double) k;
    }
}
