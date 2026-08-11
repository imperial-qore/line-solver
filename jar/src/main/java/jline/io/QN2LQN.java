/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.layered.*;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Converts a Queueing Network model to a Layered Queueing Network model.
 *
 * <p>This conversion creates an LQN representation of a closed queueing network,
 * mapping stations to processors/tasks/entries/activities and routing to
 * activity precedences (OR-forks, AND-forks, AND-joins). The resulting LQN can
 * be solved by SolverLQNS or SolverLN.</p>
 *
 * <p>Port of MATLAB: matlab/src/io/QN2LQN.m</p>
 */
public class QN2LQN {

    /** Node types that participate in routing precedences */
    private static final Set<NodeType> ROUTING_NODE_TYPES = new HashSet<NodeType>(Arrays.asList(
            NodeType.Queue, NodeType.Delay, NodeType.ClassSwitch,
            NodeType.Router, NodeType.Logger, NodeType.Fork, NodeType.Join
    ));

    /**
     * Convert a Network model to a LayeredNetwork.
     *
     * @param model the queueing network model
     * @return a LayeredNetwork representation
     */
    public static LayeredNetwork convert(Network model) {
        LayeredNetwork lqn = new LayeredNetwork(model.getName());
        NetworkStruct sn = model.getStruct(false);

        int nNodes = sn.nnodes;
        int nClasses = sn.nclasses;
        int nChains = sn.nchains;

        // Pseudo host with INF servers, INF scheduling
        Processor PH = new Processor(lqn, model.getName(), Integer.MAX_VALUE, SchedStrategy.INF);

        // Reference tasks and entries per chain
        Task[] RT = new Task[nChains];
        Entry[] RE = new Entry[nChains];
        for (int c = 0; c < nChains; c++) {
            Matrix inchainMat = sn.inchain.get(c);
            int totalJobs = 0;
            for (int idx = 0; idx < inchainMat.getNumCols(); idx++) {
                totalJobs += (int) sn.njobs.get(0, (int) inchainMat.get(0, idx));
            }
            RT[c] = new Task(lqn, "RefTask_" + (c + 1), totalJobs, SchedStrategy.REF);
            RT[c].on(PH);
            RE[c] = new Entry(lqn, "Chain_" + (c + 1));
            RE[c].on(RT[c]);
        }

        // Create Hosts, Tasks, Entries, and Activities for Queue/Delay nodes
        Processor[] P = new Processor[nNodes];
        Task[] T = new Task[nNodes];
        Entry[][] E = new Entry[nNodes][nClasses];
        Activity[][] A = new Activity[nNodes][nClasses];

        for (int i = 0; i < nNodes; i++) {
            NodeType nodeType = sn.nodetype.get(i);
            if (nodeType == NodeType.Queue || nodeType == NodeType.Delay) {
                int ist = (int) sn.nodeToStation.get(i);
                int nservers = (int) sn.nservers.get(ist);
                SchedStrategy schedStrat = sn.sched.get(sn.stations.get(ist));
                P[i] = new Processor(lqn, sn.nodenames.get(i), nservers, schedStrat);
                T[i] = new Task(lqn, "T_" + sn.nodenames.get(i), Integer.MAX_VALUE, SchedStrategy.INF);
                T[i].on(P[i]);

                for (int r = 0; r < nClasses; r++) {
                    int c = findChain(sn, r);
                    Matrix visits = sn.visits.get(c);
                    if (visits.get(i, r) > 0) {
                        E[i][r] = new Entry(lqn, "E" + (i + 1) + "_" + (r + 1));
                        E[i][r].on(T[i]);
                        Distribution serviceDist = ((jline.lang.nodes.Queue) model.getNodes().get(i)).getService(model.getClasses().get(r));
                        A[i][r] = new Activity(lqn, "Q" + (i + 1) + "_" + (r + 1), serviceDist);
                        A[i][r].on(T[i]).boundTo(E[i][r]).repliesTo(E[i][r]);
                    }
                }
            }
            // ClassSwitch, Router, Logger, Fork, Join, Source, Sink: no host/task/entry needed
        }

        // Create pseudo-activities on reference tasks
        Activity[][][] PA = new Activity[nChains][nNodes][nClasses];
        int[][] boundToRE = new int[nChains][2]; // [node, class] that is boundTo RE
        boolean[] boundToRESet = new boolean[nChains];

        for (int i = 0; i < nNodes; i++) {
            NodeType nodeType = sn.nodetype.get(i);

            if (nodeType == NodeType.ClassSwitch || nodeType == NodeType.Router
                    || nodeType == NodeType.Logger) {
                // Passthrough routing nodes: create Immediate pseudo-activities
                for (int r = 0; r < nClasses; r++) {
                    int c = findChain(sn, r);
                    if (hasIncomingRouting(sn, i, r)) {
                        PA[c][i][r] = new Activity(lqn,
                                "CS_" + (c + 1) + "_" + (i + 1) + "_" + (r + 1),
                                Immediate.getInstance());
                        PA[c][i][r].on(RT[c]);
                    }
                }
            } else if (nodeType == NodeType.Fork || nodeType == NodeType.Join) {
                // Fork/Join nodes: create Immediate pseudo-activities
                for (int r = 0; r < nClasses; r++) {
                    int c = findChain(sn, r);
                    if (hasIncomingRouting(sn, i, r)) {
                        PA[c][i][r] = new Activity(lqn,
                                "FJ_" + (c + 1) + "_" + (i + 1) + "_" + (r + 1),
                                Immediate.getInstance());
                        PA[c][i][r].on(RT[c]);
                    }
                }
            } else if (nodeType == NodeType.Queue || nodeType == NodeType.Delay) {
                for (int r = 0; r < nClasses; r++) {
                    int c = findChain(sn, r);
                    Matrix visits = sn.visits.get(c);
                    if (visits.get(i, r) > 0) {
                        Matrix inchainMat = sn.inchain.get(c);
                        int firstClassInChain = (int) inchainMat.get(0, 0);
                        int refstat = (int) sn.refstat.get(firstClassInChain);
                        if (i == refstat && r == firstClassInChain) {
                            PA[c][i][r] = new Activity(lqn,
                                    "A" + (i + 1) + "_" + (r + 1),
                                    Immediate.getInstance());
                            PA[c][i][r].on(RT[c]).boundTo(RE[c]).synchCall(E[i][r]);
                            boundToRE[c][0] = i;
                            boundToRE[c][1] = r;
                            boundToRESet[c] = true;
                        } else {
                            PA[c][i][r] = new Activity(lqn,
                                    "A" + (i + 1) + "_" + (r + 1),
                                    Immediate.getInstance());
                            PA[c][i][r].on(RT[c]).synchCall(E[i][r]);
                        }
                    }
                }
            }
        }

        // Build OR-fork / AND-fork precedences from routing matrix
        int[][] usedInORFork = new int[nNodes][nClasses];
        for (int c = 0; c < nChains; c++) {
            Matrix inchainMat = sn.inchain.get(c);
            for (int i = 0; i < nNodes; i++) {
                NodeType nodeTypeI = sn.nodetype.get(i);
                if (!ROUTING_NODE_TYPES.contains(nodeTypeI)) {
                    continue;
                }
                for (int ri = 0; ri < inchainMat.getNumCols(); ri++) {
                    int r = (int) inchainMat.get(0, ri);
                    List<String> orforkPrec = new ArrayList<String>();
                    List<Double> orforkProb = new ArrayList<Double>();

                    for (int j = 0; j < nNodes; j++) {
                        NodeType nodeTypeJ = sn.nodetype.get(j);
                        if (!ROUTING_NODE_TYPES.contains(nodeTypeJ)) {
                            continue;
                        }
                        // Skip Join destinations — handled by AND-Join precedences below
                        if (nodeTypeJ == NodeType.Join) {
                            continue;
                        }
                        for (int si = 0; si < inchainMat.getNumCols(); si++) {
                            int s = (int) inchainMat.get(0, si);
                            double pr = sn.rtnodes.get(i * nClasses + r, j * nClasses + s);
                            if (pr > 0 && hasIncomingRouting(sn, i, r)) {
                                if (boundToRESet[c]) {
                                    if (boundToRE[c][0] == j && boundToRE[c][1] == s) {
                                        // Destination is the bound-to-RE activity: create End activity
                                        if (PA[c][i][r] != null) {
                                            Activity endAct = new Activity(lqn,
                                                    "End_" + (c + 1) + "_" + (i + 1) + "_" + (r + 1),
                                                    Immediate.getInstance());
                                            endAct.on(RT[c]);
                                            orforkPrec.add(endAct.getName());
                                            orforkProb.add(pr);
                                        }
                                    } else {
                                        if (PA[c][j][s] != null) {
                                            orforkPrec.add(PA[c][j][s].getName());
                                            orforkProb.add(pr);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (!orforkPrec.isEmpty() && PA[c][i][r] != null) {
                        if (nodeTypeI == NodeType.Fork) {
                            // Fork node: AND-Fork (all branches taken simultaneously)
                            List<String> preActs = new ArrayList<String>();
                            preActs.add(PA[c][i][r].getName());
                            RT[c].addPrecedence(
                                    new ActivityPrecedence(preActs, orforkPrec,
                                            ActivityPrecedenceType.PRE_SEQ,
                                            ActivityPrecedenceType.POST_AND));
                        } else {
                            // OR-Fork (probabilistic choice)
                            Matrix probMatrix = new Matrix(1, orforkPrec.size(), orforkPrec.size());
                            for (int k = 0; k < orforkProb.size(); k++) {
                                probMatrix.set(0, k, orforkProb.get(k));
                            }
                            List<String> preActs = new ArrayList<String>();
                            preActs.add(PA[c][i][r].getName());
                            RT[c].addPrecedence(
                                    new ActivityPrecedence(preActs, orforkPrec,
                                            ActivityPrecedenceType.PRE_SEQ,
                                            ActivityPrecedenceType.POST_OR,
                                            null, probMatrix));
                        }
                        usedInORFork[i][r]++;
                    }
                }
            }
        }

        // AND-Join precedences for Join nodes
        for (int c = 0; c < nChains; c++) {
            Matrix inchainMat = sn.inchain.get(c);
            for (int j = 0; j < nNodes; j++) {
                if (sn.nodetype.get(j) != NodeType.Join) {
                    continue;
                }
                for (int si = 0; si < inchainMat.getNumCols(); si++) {
                    int s = (int) inchainMat.get(0, si);
                    if (PA[c][j][s] == null) {
                        continue;
                    }
                    List<String> joinPre = new ArrayList<String>();
                    for (int i = 0; i < nNodes; i++) {
                        NodeType nodeTypeI = sn.nodetype.get(i);
                        if (!ROUTING_NODE_TYPES.contains(nodeTypeI)) {
                            continue;
                        }
                        for (int ri = 0; ri < inchainMat.getNumCols(); ri++) {
                            int r = (int) inchainMat.get(0, ri);
                            double pr = sn.rtnodes.get(i * nClasses + r, j * nClasses + s);
                            if (pr > 0 && PA[c][i][r] != null) {
                                joinPre.add(PA[c][i][r].getName());
                            }
                        }
                    }
                    if (joinPre.size() > 1) {
                        // Multiple predecessors: AND-Join
                        List<String> postActs = new ArrayList<String>();
                        postActs.add(PA[c][j][s].getName());
                        RT[c].addPrecedence(
                                new ActivityPrecedence(joinPre, postActs,
                                        ActivityPrecedenceType.PRE_AND,
                                        ActivityPrecedenceType.POST_SEQ));
                    } else if (joinPre.size() == 1) {
                        // Single predecessor: Serial precedence
                        List<String> preActs = new ArrayList<String>();
                        preActs.add(joinPre.get(0));
                        List<String> postActs = new ArrayList<String>();
                        postActs.add(PA[c][j][s].getName());
                        RT[c].addPrecedence(
                                new ActivityPrecedence(preActs, postActs));
                    }
                }
            }
        }

        return lqn;
    }

    /**
     * Find the chain index for a given class.
     *
     * @param sn     NetworkStruct
     * @param classR class index (0-based)
     * @return chain index (0-based)
     */
    private static int findChain(NetworkStruct sn, int classR) {
        for (int c = 0; c < sn.nchains; c++) {
            if (sn.chains.get(c, classR) > 0) {
                return c;
            }
        }
        return 0;
    }

    /**
     * Check if any routing leads to node i, class r (i.e., any column entry is positive).
     * MATLAB: any(sn.rtnodes(:, (i-1)*nclasses + r) > 0)
     *
     * @param sn NetworkStruct
     * @param i  node index (0-based)
     * @param r  class index (0-based)
     * @return true if any routing probability to (i,r) is positive
     */
    private static boolean hasIncomingRouting(NetworkStruct sn, int i, int r) {
        int col = i * sn.nclasses + r;
        int nRows = sn.rtnodes.getNumRows();
        for (int row = 0; row < nRows; row++) {
            if (sn.rtnodes.get(row, col) > 0) {
                return true;
            }
        }
        return false;
    }
}
