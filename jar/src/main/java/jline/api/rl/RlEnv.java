/**
 * @file RL Environment for Queueing Network Routing
 *
 * @since LINE 3.0
 */
package jline.api.rl;

import java.util.List;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.FromMarginal;
import jline.solvers.ssa.SolverSSA;
import jline.io.Ret.SampleResult;
import jline.util.matrix.Matrix;

/**
 * RL environment for queueing network routing decisions.
 */
public class RlEnv {
    public final Network model;
    public final int[] idxOfQueueInNodes;
    public final int[] idxOfSourceInNodes;
    public final int stateSize;
    public final double gamma;
    public final int actionSize;

    public RlEnv(Network model, int[] idxOfQueueInNodes, int[] idxOfSourceInNodes,
                 int stateSize, double gamma) {
        this.model = model;
        this.idxOfQueueInNodes = idxOfQueueInNodes;
        this.idxOfSourceInNodes = idxOfSourceInNodes;
        this.stateSize = stateSize;
        this.gamma = gamma;
        this.actionSize = idxOfQueueInNodes.length;
    }

    public Network getModel() { return model; }
    public int[] getIdxOfQueueInNodes() { return idxOfQueueInNodes; }
    public int[] getIdxOfSourceInNodes() { return idxOfSourceInNodes; }
    public int getStateSize() { return stateSize; }
    public double getGamma() { return gamma; }
    public int getActionSize() { return actionSize; }

    public boolean isInStateSpace(List<Node> nodes) {
        for (int i : idxOfQueueInNodes) {
            Node node = nodes.get(i);
            if (node instanceof StatefulNode) {
                StatefulNode sn = (StatefulNode) node;
                Matrix stateMatrix = sn.getState();
                if (stateMatrix.elementSum() > stateSize) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isInActionSpace(List<Node> nodes) {
        for (int i : idxOfQueueInNodes) {
            Node node = nodes.get(i);
            if (node instanceof StatefulNode) {
                StatefulNode sn = (StatefulNode) node;
                Matrix stateMatrix = sn.getState();
                if (stateMatrix.elementSum() > stateSize - 1) {
                    return false;
                }
            }
        }
        return true;
    }

    public SampleEvent sample() {
        SolverSSA solver = new SolverSSA(model, "verbose", false);
        SampleResult sampleResult = solver.sampleSysAggr(1);
        double t = (sampleResult.t != null && sampleResult.t.length() > 0)
                ? sampleResult.t.get(0) : 0.0;

        int depNode = -1;
        Matrix eventMatrix = sampleResult.event;
        if (eventMatrix != null && eventMatrix.getNumRows() > 0) {
            for (int row = 0; row < eventMatrix.getNumRows(); row++) {
                int nodeIdx = (int) eventMatrix.get(row, 1);
                if (containsInt(idxOfSourceInNodes, nodeIdx) || containsInt(idxOfQueueInNodes, nodeIdx)) {
                    depNode = nodeIdx;
                    break;
                }
            }
            if (depNode == -1 && eventMatrix.getNumRows() > 0) {
                depNode = (int) eventMatrix.get(0, 1);
            }
        }

        return new SampleEvent(t, depNode);
    }

    private static boolean containsInt(int[] array, int value) {
        for (int v : array) {
            if (v == value) return true;
        }
        return false;
    }

    public void update(int[] newState) {
        NetworkStruct sn = model.getStruct(false);
        for (int i = 0; i < idxOfQueueInNodes.length; i++) {
            int nodeIdx = idxOfQueueInNodes[i];
            Node node = model.getNodes().get(nodeIdx);
            if (node instanceof StatefulNode) {
                Matrix marginal = new Matrix(1, 1);
                marginal.set(0, 0, (double) newState[i]);
                Matrix newNodeState = FromMarginal.fromMarginal(sn, nodeIdx, marginal);
                ((StatefulNode) node).setState(newNodeState);
            }
        }
    }

    public void reset() {
        model.reset();
        model.initDefault();
    }
}
