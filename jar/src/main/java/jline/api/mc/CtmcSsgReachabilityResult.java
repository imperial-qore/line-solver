package jline.api.mc;

import jline.lang.NetworkStruct;
import jline.lang.nodes.StatefulNode;
import jline.util.matrix.Matrix;

import java.util.Map;

/**
 * Result data class for CTMC SSG reachability analysis
 */
public final class CtmcSsgReachabilityResult {
    public final Matrix stateSpace;
    public final Matrix stateSpaceAggr;
    public final Matrix stateSpaceHashed;
    public final Map<StatefulNode, Matrix> nodeStateSpace;
    public final NetworkStruct sn;

    public CtmcSsgReachabilityResult(Matrix stateSpace, Matrix stateSpaceAggr, Matrix stateSpaceHashed,
                                     Map<StatefulNode, Matrix> nodeStateSpace, NetworkStruct sn) {
        this.stateSpace = stateSpace;
        this.stateSpaceAggr = stateSpaceAggr;
        this.stateSpaceHashed = stateSpaceHashed;
        this.nodeStateSpace = nodeStateSpace;
        this.sn = sn;
    }

    public Matrix getStateSpace() { return stateSpace; }
    public Matrix getStateSpaceAggr() { return stateSpaceAggr; }
    public Matrix getStateSpaceHashed() { return stateSpaceHashed; }
    public Map<StatefulNode, Matrix> getNodeStateSpace() { return nodeStateSpace; }
    public NetworkStruct getSn() { return sn; }
}
