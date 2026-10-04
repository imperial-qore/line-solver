package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * Number of sibling tasks forked per parent job on a fork-join pair.
 *
 * <p>Port of matlab/src/api/fj/sn_join_siblings.m.</p>
 */
public final class SnJoinSiblings {
    private SnJoinSiblings() {}

    /** Class-blind form: the widest fork over the classes. */
    public static int snJoinSiblings(NetworkStruct sn, Node joinNode) {
        return snJoinSiblings(sn, joinNode, -1);
    }

    /**
     * Siblings are counted at the FORK, as the simulation engines count them.
     *
     * <p>THE COUNT IS PER LINK, not the out-degree times a node-wide scalar. A fork carries a
     * VARIABLE FORKING LEVEL: setTasksPerLink(class, n [, dest]) sets one link of one class,
     * setTasksPerLinkDistribution makes the degree a draw, and setBranchProb makes a link taken
     * with probability below one. All three land in {@code ForkNodeParam.fanOutLink} and
     * {@code .fanOutProb}, both (nnodes x nclasses) and indexed by DESTINATION NODE, with the
     * DISTRIBUTION case storing its mean, so</p>
     *
     * <pre>    N = sum_d fanOutProb(d,r) * fanOutLink(d,r)</pre>
     *
     * <p>which is the EXPECTED sibling count and reduces to out-degree times the node-wide
     * {@code fanOut} on a model that sets none of the three. Falls back to that older product when
     * {@code fanOutLink} is absent, and to the Join's in-degree when the matched fork cannot be
     * identified.</p>
     *
     * <p>The result is what a quorum is measured against: the join fires on the k-th of N siblings
     * and the remaining N-k are discarded when they arrive.</p>
     *
     * @param sn       the network structure
     * @param joinNode the Join node
     * @param r        the class index, or -1 for the widest fork over the classes
     * @return the number of siblings forked per parent job, 0 when it cannot be determined
     */
    public static int snJoinSiblings(NetworkStruct sn, Node joinNode, int r) {
        if (sn == null || joinNode == null || sn.connmatrix == null) {
            return 0;
        }
        int j = joinNode.getNodeIndex();
        if (j < 0 || j >= sn.connmatrix.getNumCols()) {
            return 0;
        }
        int n = 0;
        for (int a = 0; a < sn.connmatrix.getNumRows(); a++) {
            if (sn.connmatrix.get(a, j) > 0) {
                n++;
            }
        }
        if (sn.fj == null || j >= sn.fj.getNumCols()) {
            return n;
        }
        int f = -1;
        for (int a = 0; a < sn.fj.getNumRows(); a++) {
            if (sn.fj.get(a, j) > 0) {
                f = a;
                break;
            }
        }
        if (f < 0 || f >= sn.connmatrix.getNumRows()) {
            return n;
        }
        ForkNodeParam fp = null;
        Node forkNode = sn.nodes != null && f < sn.nodes.size() ? sn.nodes.get(f) : null;
        if (forkNode != null && sn.nodeparam != null) {
            NodeParam param = sn.nodeparam.get(forkNode);
            if (param instanceof ForkNodeParam) {
                fp = (ForkNodeParam) param;
            }
        }

        // The per-link count, when the refresh has built it.
        if (fp != null && fp.fanOutLink != null && fp.fanOutLink.getNumRows() > 0) {
            Matrix fol = fp.fanOutLink;
            Matrix fop = (fp.fanOutProb != null
                    && fp.fanOutProb.getNumRows() == fol.getNumRows()
                    && fp.fanOutProb.getNumCols() == fol.getNumCols()) ? fp.fanOutProb : null;
            int lo = (r >= 0 && r < fol.getNumCols()) ? r : 0;
            int hi = (r >= 0 && r < fol.getNumCols()) ? r : fol.getNumCols() - 1;
            double best = 0.0;
            for (int c = lo; c <= hi; c++) {
                double acc = 0.0;
                for (int d = 0; d < fol.getNumRows(); d++) {
                    double link = fol.get(d, c);
                    if (link <= 0) {
                        continue;
                    }
                    acc += link * (fop == null ? 1.0 : fop.get(d, c));
                }
                if (acc > best) {
                    best = acc;
                }
            }
            if (best > 0) {
                return (int) Math.round(best);
            }
        }

        // Fallback: out-degree times the node-wide scalar.
        int w = 1;
        if (fp != null && !Double.isNaN(fp.fanOut)) {
            w = Math.max(1, (int) Math.round(fp.fanOut));
        }
        int deg = 0;
        for (int b = 0; b < sn.connmatrix.getNumCols(); b++) {
            if (sn.connmatrix.get(f, b) > 0) {
                deg++;
            }
        }
        return deg * w;
    }
}
