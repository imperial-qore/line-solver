package jline.api.sn;

import java.util.HashMap;
import java.util.Map;

import jline.lang.NetworkStruct;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.util.matrix.Matrix;

public final class SnGetStateAggr {
    private SnGetStateAggr() {}

    /**
     * Aggregates the state of the network.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return a map of stateful nodes to their aggregated state matrices
     */
    public static Map<StatefulNode, Matrix> snGetStateAggr(NetworkStruct sn) {
        Map<StatefulNode, Matrix> initialState = sn.state;
        Map<StatefulNode, Matrix> initialStateAggr = new HashMap<StatefulNode, Matrix>();
        for (int isf = 0; isf < initialState.size(); isf++) {
            int ind = (int) sn.statefulToNode.get(isf);
            StatefulNode node = sn.stateful.get(isf);
            State.StateMarginalStatistics aggrState =
                    ToMarginal.toMarginalAggr(sn, ind, initialState.get(node), null, null, null, null, null);
            initialStateAggr.put(node, aggrState.nir);
        }
        return initialStateAggr;
    }
}
