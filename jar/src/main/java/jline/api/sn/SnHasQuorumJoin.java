package jline.api.sn;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.JoinStrategy;
import jline.lang.nodeparam.JoinNodeParam;

import java.util.Map;

/**
 * Checks if the network has a quorum (k-of-n) join.
 *
 * <p>Port of matlab/src/api/sn/sn_has_quorum_join.m.</p>
 */
public final class SnHasQuorumJoin {
    private SnHasQuorumJoin() {}

    /**
     * True if some Join node declares a non-standard strategy with a positive required count in
     * some class, i.e. it fires before every sibling has arrived. The sibling count is not
     * re-derived here, so a declaration with k &gt;= n reads as a quorum; use
     * {@link SnJoinQuorum} where the branch count is known and the distinction matters, as the
     * fork-join fixed point does.
     *
     * @param sn the network structure
     * @return true if some join declares a positive quorum
     */
    public static boolean snHasQuorumJoin(NetworkStruct sn) {
        if (sn == null || sn.nodeparam == null) {
            return false;
        }
        for (NodeParam param : sn.nodeparam.values()) {
            if (!(param instanceof JoinNodeParam)) {
                continue;
            }
            JoinNodeParam joinParam = (JoinNodeParam) param;
            if (joinParam.joinStrategy == null || joinParam.joinRequired == null) {
                continue;
            }
            for (Map.Entry<JobClass, JoinStrategy> entry : joinParam.joinStrategy.entrySet()) {
                if (entry.getValue() == null || entry.getValue() == JoinStrategy.STD) {
                    continue;
                }
                Double required = joinParam.joinRequired.get(entry.getKey());
                if (required != null && required > 0) {
                    return true;
                }
            }
        }
        return false;
    }
}
