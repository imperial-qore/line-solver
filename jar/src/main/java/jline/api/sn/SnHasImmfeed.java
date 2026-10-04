package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

public final class SnHasImmfeed {
    private SnHasImmfeed() {}

    /**
     * Whether immediate feedback is EFFECTIVE anywhere in the model.
     *
     * A declaration alone is not enough. Queue.setImmediateFeedback marks a station
     * and JobClass.setImmediateFeedback marks a class -- and the class spelling marks
     * EVERY station, since sn.immfeed is the OR of the two -- so a model with no
     * self-loop at all can carry a full sn.immfeed matrix while the feature changes
     * nothing. Reading the raw matrix made every solver that consults it warn, or
     * refuse, on a plain M/M/1 that merely mentioned the flag.
     *
     * Immediate feedback is effective at (station i, class r) when sn.immfeed(i,r)
     * holds AND the routing table has a self-loop INTO (i,r) from some class s at the
     * same station, which is the only way a job can come back to the server it just
     * left. A class switch on the way round is folded into sn.rt by refreshRouting, so
     * the incoming class s need not be r.
     *
     * Solvers that handle immediate feedback look at the SYNCHRONIZATION instead; see
     * State.immfeedSelfLoop, which applies the same test per sync.
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasImmfeed(NetworkStruct sn) {
        if (sn == null || sn.immfeed == null || sn.immfeed.isEmpty()) {
            return false;
        }
        boolean declared = false;
        for (int i = 0; i < sn.immfeed.getNumRows() && !declared; i++) {
            for (int r = 0; r < sn.immfeed.getNumCols(); r++) {
                if (sn.immfeed.get(i, r) > 0.0) {
                    declared = true;
                    break;
                }
            }
        }
        if (!declared) {
            return false;
        }
        final Matrix rt = sn.rt;
        // Without a routing table there is nothing to qualify the declaration with, so
        // report it as declared rather than silently dropping it.
        if (rt == null || rt.isEmpty()) {
            return true;
        }
        final int R = sn.nclasses;
        for (int ist = 0; ist < Math.min(sn.nstations, sn.immfeed.getNumRows()); ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            if (ind < 0 || ind >= sn.nnodes || sn.isstateful.get(ind, 0) != 1.0) {
                continue;
            }
            int isf = (int) sn.nodeToStateful.get(ind);
            for (int r = 0; r < Math.min(R, sn.immfeed.getNumCols()); r++) {
                if (sn.immfeed.get(ist, r) <= 0.0) {
                    continue;
                }
                int col = isf * R + r;
                if (col >= rt.getNumCols()) {
                    continue;
                }
                for (int s = 0; s < R; s++) {
                    int row = isf * R + s;
                    if (row < rt.getNumRows() && rt.get(row, col) > 0.0) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
