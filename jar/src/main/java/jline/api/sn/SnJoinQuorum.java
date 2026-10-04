package jline.api.sn;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.JoinStrategy;
import jline.lang.nodeparam.JoinNodeParam;
import jline.lang.nodes.Node;

/**
 * Number of sibling tasks a Join node waits for in a given class.
 *
 * <p>Port of matlab/src/api/fj/sn_join_quorum.m.</p>
 */
public final class SnJoinQuorum {
    private SnJoinQuorum() {}

    /**
     * A standard join, an absent declaration, a non-positive quorum and a quorum that is not
     * smaller than the sibling count all return nbranches, i.e. the ordinary AND-join: those are
     * the four ways a join fires only when every sibling has arrived.
     *
     * <p>The count is the one the simulation engines apply (SolverLDES fixes it at FORK time and
     * discards the stragglers when they reach the join), so an analytical solver reading it here
     * charges the same synchronisation event.</p>
     *
     * @param sn        the network structure
     * @param joinNode  the Join node
     * @param jobClass  the class the siblings are matched in
     * @param nbranches the number of siblings the fork emits
     * @return the number of siblings the join waits for
     */
    public static int snJoinQuorum(NetworkStruct sn, Node joinNode, JobClass jobClass, int nbranches) {
        if (sn == null || sn.nodeparam == null || joinNode == null) {
            return nbranches;
        }
        NodeParam param = sn.nodeparam.get(joinNode);
        if (!(param instanceof JoinNodeParam)) {
            return nbranches;
        }
        JoinNodeParam joinParam = (JoinNodeParam) param;
        if (joinParam.joinStrategy == null || joinParam.joinRequired == null) {
            return nbranches;
        }
        JoinStrategy strategy = joinParam.joinStrategy.get(jobClass);
        if (strategy == null || strategy == JoinStrategy.STD) {
            return nbranches;
        }
        Double required = joinParam.joinRequired.get(jobClass);
        if (required == null) {
            return nbranches;
        }
        int q = (int) Math.round(required);
        if (q > 0 && q < nbranches) {
            return q;
        }
        return nbranches;
    }
}
