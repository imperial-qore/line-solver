/**
 * @file General RL Environment for Queueing Network Control
 *
 * Implements a general-purpose Reinforcement Learning environment for queueing
 * networks where actions are dispatch/routing decisions at specific nodes.
 * Unlike RlEnv which handles simple source-to-queue routing, this environment
 * supports arbitrary action nodes with configurable action spaces derived from
 * the network topology.
 *
 * Port of: matlab/src/api/rl/rl_env_general.m
 *
 * @since LINE 3.0
 */
package jline.api.rl;

import java.util.ArrayList;
import java.util.HashMap;

import jline.io.Ret.SampleResult;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.constant.EventType;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.State;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;

/**
 * General RL environment for queueing network control decisions.
 */
public class RlEnvGeneral {
    public final Network model;
    public final int[] idxOfQueueInNodes;
    public final int[] idxOfActionNodes;
    public final int stateSize;
    public final double gamma;

    /** Number of queues in the network. */
    public final int nqueues;

    /**
     * Action space map: for each action node index, stores the array of
     * downstream node indices reachable from that node.
     */
    public final HashMap<Integer, int[]> actionSpace;

    public RlEnvGeneral(Network model, int[] idxOfQueueInNodes, int[] idxOfActionNodes,
                        int stateSize, double gamma) {
        this.model = model;
        this.idxOfQueueInNodes = idxOfQueueInNodes;
        this.idxOfActionNodes = idxOfActionNodes;
        this.stateSize = stateSize;
        this.gamma = gamma;
        this.nqueues = idxOfQueueInNodes.length;

        this.actionSpace = new HashMap<Integer, int[]>();
        Matrix connMatrix = model.getConnectionMatrix();
        for (int i : idxOfActionNodes) {
            ArrayList<Integer> reachable = new ArrayList<Integer>();
            int nCols = connMatrix.getNumCols();
            for (int j = 0; j < nCols; j++) {
                if (connMatrix.get(i, j) == 1.0) {
                    reachable.add(j);
                }
            }
            int[] arr = new int[reachable.size()];
            for (int k = 0; k < reachable.size(); k++) {
                arr[k] = reachable.get(k);
            }
            actionSpace.put(i, arr);
        }
    }

    public Network getModel() { return model; }
    public int[] getIdxOfQueueInNodes() { return idxOfQueueInNodes; }
    public int[] getIdxOfActionNodes() { return idxOfActionNodes; }
    public int getStateSize() { return stateSize; }
    public double getGamma() { return gamma; }
    public int getNqueues() { return nqueues; }
    public HashMap<Integer, int[]> getActionSpace() { return actionSpace; }

    /**
     * Checks if the given state vector is within the defined state space.
     */
    public boolean isInStateSpace(int[] state) {
        if (state.length != idxOfQueueInNodes.length) {
            throw new IllegalArgumentException(
                    "State size mismatch: state size=" + state.length
                            + ", required size=" + idxOfQueueInNodes.length);
        }
        for (int i = 0; i < idxOfQueueInNodes.length; i++) {
            if (state[i] > stateSize) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if actions can be taken from the given state.
     */
    public boolean isInActionSpace(int[] state) {
        if (state.length != idxOfQueueInNodes.length) {
            throw new IllegalArgumentException(
                    "State size mismatch: state size=" + state.length
                            + ", required size=" + idxOfQueueInNodes.length);
        }
        for (int i = 0; i < idxOfQueueInNodes.length; i++) {
            if (state[i] > stateSize - 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * Samples the next event from the environment using the SSA solver.
     */
    public GeneralSampleEvent sample() {
        SolverSSA solver = new SolverSSA(model, "verbose", false);
        SampleResult sampleResult = solver.sampleSysAggr(1);
        double dt;
        if (sampleResult.t != null && sampleResult.t.length() > 0) {
            dt = sampleResult.t.get(0);
        } else {
            dt = 0.0;
        }

        int depNode = -1;
        int arvNode = -1;

        Matrix eventMatrix = sampleResult.event;
        if (eventMatrix != null && eventMatrix.getNumRows() > 0) {
            // Parse events from the event matrix
            for (int row = 0; row < eventMatrix.getNumRows(); row++) {
                int nodeIdx = (int) eventMatrix.get(row, 1);
                // Heuristic: first event row is typically the departure,
                // second is the arrival (matches MATLAB event{1}, event{2} ordering)
                if (row == 0) {
                    depNode = nodeIdx;
                } else if (row == 1) {
                    arvNode = nodeIdx;
                }
            }
        }

        return new GeneralSampleEvent(dt, depNode, arvNode, sampleResult);
    }

    /**
     * Updates the model state after an event using the SSA sample result.
     */
    public void update(SampleResult sampleResult) {
        NetworkStruct sn = model.getStruct(false);
        Matrix eventMatrix = sampleResult.event;
        if (eventMatrix == null || eventMatrix.getNumRows() == 0) return;

        // Process each event in the sample
        int jobClass;
        if (eventMatrix.getNumRows() > 0) {
            jobClass = (int) eventMatrix.get(0, 2);
        } else {
            jobClass = 0;
        }

        for (int row = 0; row < eventMatrix.getNumRows(); row++) {
            int nodeIdx = (int) eventMatrix.get(row, 1);
            Node node = model.getNodes().get(nodeIdx);
            if (node instanceof StatefulNode) {
                StatefulNode sNode = (StatefulNode) node;
                Matrix currentState = sNode.getState();
                EventType eventType = (row == 0) ? EventType.DEP : EventType.ARV;
                jline.io.Ret.EventResult result = State.afterEvent(sn, nodeIdx, currentState, eventType, jobClass, true);
                if (result != null && result.outspace != null && result.outspace.getNumRows() > 0) {
                    sNode.setState(result.outspace.getRow(0));
                }
            }
        }
    }

    /**
     * Resets the environment to its initial state.
     */
    public void reset() {
        model.reset();
        model.initDefault();
    }
}
