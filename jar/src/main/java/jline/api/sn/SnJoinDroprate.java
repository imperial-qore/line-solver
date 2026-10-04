package jline.api.sn;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * Rate at which sibling tasks are discarded at each Join.
 *
 * <p>Port of matlab/src/api/fj/sn_join_droprate.m.</p>
 */
public final class SnJoinDroprate {
    private SnJoinDroprate() {}

    /**
     * A Join is the one station where the loss identity ArvR - Tput does NOT hold, because the two
     * rates are in different units: AN counts the SIBLINGS offered to the join (N per parent job)
     * while TN counts the PARENT jobs released by it (one per synchronisation). Reading
     * ArvR - Tput there reports (N-1)/N of the offered traffic as lost at every join, standard
     * joins included, when a standard join loses nothing at all.
     *
     * <p>The siblings a join actually consumes are K per synchronisation, where K is the quorum
     * (K = N on a standard join), so DropRateJoin = max(0, AN - K*TN), which is 0 for a standard
     * join and (N-K)*TN for a quorum.</p>
     *
     * <p>This is the DERIVED value, exact given TN and AN. A solver that MEASURES the discards on
     * its own sample path (SolverLDES) reports its own.</p>
     *
     * @param sn the network structure
     * @param TN station throughputs
     * @param AN station arrival rates
     * @return an (nstations x nclasses) matrix, zero away from the Join rows
     */
    public static Matrix snJoinDroprate(NetworkStruct sn, Matrix TN, Matrix AN) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix out = new Matrix(M, K);
        out.fill(0.0);
        if (TN == null || AN == null || sn.fj == null || sn.fj.getNumRows() == 0
                || sn.nodes == null) {
            return out;
        }
        for (int ind = 0; ind < sn.nodes.size(); ind++) {
            if (sn.nodetype == null || ind >= sn.nodetype.size()
                    || sn.nodetype.get(ind) != NodeType.Join) {
                continue;
            }
            Node joinNode = sn.nodes.get(ind);
            // sn.nodeToStation is a 1 x nnodes ROW vector (Network.java:7872), so it
            // is read with the linear accessor, as the rest of the JAR reads it.
            if (sn.nodeToStation == null || ind >= sn.nodeToStation.getNumElements()) {
                continue;
            }
            int ist = (int) sn.nodeToStation.get(ind);
            if (ist < 0 || ist >= M) {
                continue;
            }
            for (int r = 0; r < K && r < sn.jobclasses.size(); r++) {
                double a = ist < AN.getNumRows() && r < AN.getNumCols() ? AN.get(ist, r) : 0.0;
                double t = ist < TN.getNumRows() && r < TN.getNumCols() ? TN.get(ist, r) : 0.0;
                if (!Double.isFinite(a) || !Double.isFinite(t) || a <= 0) {
                    continue;
                }
                // PER CLASS: a variable forking level makes the sibling count differ
                // between classes, so it cannot be hoisted out of this loop.
                int nsib = SnJoinSiblings.snJoinSiblings(sn, joinNode, r);
                if (nsib <= 0) {
                    continue;
                }
                JobClass jobClass = sn.jobclasses.get(r);
                int kreq = SnJoinQuorum.snJoinQuorum(sn, joinNode, jobClass, nsib);
                out.set(ist, r, Math.max(0.0, a - kreq * t));
            }
        }
        return out;
    }
}
