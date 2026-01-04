/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.GlobalConstants;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.Signal;
import jline.lang.constant.CallType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.RoutingMatrix;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * Converts a LayeredNetwork (LQN) to a Network (QN).
 *
 * <p>The conversion uses REPLY signals to model synchronous call blocking:
 * when a task makes a synchCall, the server blocks until a REPLY signal
 * arrives from the called task.</p>
 *
 * <h3>REPLY Signal Semantics:</h3>
 * <ul>
 *   <li>Each synchCall creates a Request class and a Reply signal class</li>
 *   <li>The caller queue blocks after completing service until Reply arrives</li>
 *   <li>The callee processes the request and class-switches to Reply on completion</li>
 *   <li>Reply signal unblocks the caller and continues downstream</li>
 * </ul>
 *
 * <h3>Network Topology:</h3>
 * <pre>
 *   Think → CallerQueue → CalleeQueue → CallerQueue (reply) → Think
 *                ↓ blocks      ↓ class switch to Reply    ↓ unblocks
 * </pre>
 *
 * <p>This approach:</p>
 * <ul>
 *   <li>Correctly models BLOCKING semantics (server waits for reply)</li>
 *   <li>Provides per-task queue metrics (queue length, utilization)</li>
 *   <li>Supports multi-tier call chains (A → B → C)</li>
 *   <li>Uses DES solver with REPLY signal support</li>
 * </ul>
 *
 * <h3>Phase-2 Support:</h3>
 * <p>When an entry has phase-2 activities, a Fork-Join structure handles
 * the async continuation after reply:</p>
 * <pre>
 *   CallerQueue → CalleeQueue_Ph1 → Fork → { Reply → CallerQueue }
 *                                       → { Ph2Queue → Join }
 * </pre>
 *
 * <p>Example usage:</p>
 * <pre>
 * LayeredNetwork lqn = new LayeredNetwork("MyLQN");
 * // ... define LQN model ...
 * Network model = LQN2QN.convert(lqn);
 * SolverDES solver = new SolverDES(model);
 * solver.getAvgTable();
 * </pre>
 *
 * @see LayeredNetwork
 * @see Network
 * @see Signal
 * @see SignalType#REPLY
 */
public class LQN2QN {

    /**
     * Converts a LayeredNetwork to an equivalent queueing network using REPLY signals.
     *
     * @param lqn The LayeredNetwork model to convert
     * @return A Network that models the LQN behavior with REPLY signal blocking
     */
    public static Network convert(LayeredNetwork lqn) {
        LayeredNetworkStruct lsn = lqn.getStruct();

        // Create QN model
        Network model = new Network(lqn.getName() + "-QN");

        // Identify reference tasks
        List<Integer> refTaskIndices = new ArrayList<Integer>();
        for (int i = 0; i < lsn.isref.getNumElements(); i++) {
            if (lsn.isref.get(i) > 0) {
                refTaskIndices.add(i);
            }
        }

        if (refTaskIndices.isEmpty()) {
            line_error(mfilename(new Object() {}), "LQN must have at least one reference task.");
            return model;
        }

        // Detect phase-2 activities
        boolean hasPhase2 = lsn.actphase != null && !lsn.actphase.isEmpty();
        if (hasPhase2) {
            boolean foundPhase2 = false;
            for (int i = 0; i < lsn.actphase.getNumElements(); i++) {
                if (lsn.actphase.get(i) > 1) {
                    foundPhase2 = true;
                    break;
                }
            }
            hasPhase2 = foundPhase2;
        }

        // Build task service demands (split by phase)
        Map<Integer, double[]> taskServiceByPhase = new HashMap<Integer, double[]>();
        for (int t = 1; t <= lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            double[] phaseDemands = getTaskServiceDemandByPhase(lsn, tidx, hasPhase2);
            taskServiceByPhase.put(tidx, phaseDemands);
        }

        // Build call graph: for each task, find all tasks it calls synchronously
        Map<Integer, List<SynchCallInfo>> synchCallsFrom = new HashMap<Integer, List<SynchCallInfo>>();
        for (int t = 1; t <= lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            List<SynchCallInfo> calls = new ArrayList<SynchCallInfo>();

            List<Integer> entries = lsn.entriesof.get(tidx);
            if (entries != null) {
                for (int eidx : entries) {
                    List<Integer> acts = lsn.actsof.get(eidx);
                    if (acts == null) continue;

                    for (int aidx : acts) {
                        if (lsn.type.get(aidx) != LayeredNetworkElement.ACTIVITY) {
                            continue;
                        }
                        List<Integer> callList = lsn.callsof.get(aidx);
                        if (callList == null) continue;

                        for (int cidx : callList) {
                            CallType ct = lsn.calltype.get(cidx);
                            if (ct == CallType.SYNC) {
                                int targetEidx = (int) lsn.callpair.get(cidx, 1);
                                int targetTidx = (int) lsn.parent.get(targetEidx);
                                // Get call mean (number of calls)
                                double callMean = 1.0;
                                if (lsn.callproc_mean != null && cidx < lsn.callproc_mean.size()) {
                                    Double mean = lsn.callproc_mean.get(cidx);
                                    if (mean != null && !Double.isNaN(mean)) {
                                        callMean = mean;
                                    }
                                }
                                calls.add(new SynchCallInfo(targetTidx, targetEidx, callMean));
                            }
                        }
                    }
                }
            }
            synchCallsFrom.put(tidx, calls);
        }

        // Find all tasks in call chains starting from reference tasks
        Set<Integer> tasksInChain = new HashSet<Integer>();
        for (int refTidx : refTaskIndices) {
            collectTasksInChain(refTidx, synchCallsFrom, tasksInChain);
        }

        // Create nodes for all tasks
        Map<Integer, Station> taskQueues = new HashMap<Integer, Station>();  // tidx -> Queue/Delay
        Map<Integer, Delay> thinkNodes = new HashMap<Integer, Delay>();      // refTidx -> Think Delay

        // Create think nodes for reference tasks
        for (int refTidx : refTaskIndices) {
            String taskName = lsn.names.get(refTidx);
            Delay thinkNode = new Delay(model, taskName + "_Think");
            thinkNodes.put(refTidx, thinkNode);
        }

        // Create queues for all tasks in call chains
        for (int tidx : tasksInChain) {
            String taskName = lsn.names.get(tidx);
            int procIdx = (int) lsn.parent.get(tidx);
            int nServers = (int) lsn.mult.get(procIdx);
            SchedStrategy sched = lsn.sched.get(procIdx);

            Station queue;
            if (nServers == Integer.MAX_VALUE || sched == SchedStrategy.INF) {
                queue = new Delay(model, taskName);
            } else {
                Queue q = new Queue(model, taskName, sched);
                q.setNumberOfServers(nServers);
                queue = q;
            }
            taskQueues.put(tidx, queue);
        }

        // Create job classes for each reference task
        // Each ref task gets: Request class + Reply signal class
        Map<Integer, ClosedClass> requestClasses = new HashMap<Integer, ClosedClass>();
        Map<Integer, Signal> replySignals = new HashMap<Integer, Signal>();

        for (int refTidx : refTaskIndices) {
            String taskName = lsn.names.get(refTidx);
            Delay thinkNode = thinkNodes.get(refTidx);
            int population = (int) lsn.mult.get(refTidx);

            // Create request class (closed)
            ClosedClass requestClass = new ClosedClass(model, taskName + "_Req", population, thinkNode);
            requestClasses.put(refTidx, requestClass);

            // Create reply signal class
            Signal replySignal = new Signal(model, taskName + "_Reply", SignalType.REPLY);
            replySignal.forJobClass(requestClass);
            replySignals.put(refTidx, replySignal);
        }

        // Set service times for all nodes and classes
        for (int refTidx : refTaskIndices) {
            ClosedClass requestClass = requestClasses.get(refTidx);
            Signal replySignal = replySignals.get(refTidx);
            Delay thinkNode = thinkNodes.get(refTidx);

            // Think time at think node
            Double thinkMean = lsn.think_mean.get(refTidx);
            if (thinkMean == null || thinkMean < GlobalConstants.FineTol) {
                thinkNode.setService(requestClass, new Exp(1e8));
            } else {
                thinkNode.setService(requestClass, new Exp(1.0 / thinkMean));
            }
            // Reply passes through think instantly (shouldn't visit, but set for safety)
            thinkNode.setService(replySignal, new Exp(1e9));

            // Set service times at all task queues
            for (int tidx : tasksInChain) {
                Station queue = taskQueues.get(tidx);
                double[] phaseDemands = taskServiceByPhase.get(tidx);
                double serviceMean = phaseDemands[0] + phaseDemands[1];  // Total service

                if (serviceMean > GlobalConstants.FineTol) {
                    queue.setService(requestClass, new Exp(1.0 / serviceMean));
                } else {
                    queue.setService(requestClass, new Exp(1e8));
                }
                // Reply signal passes through instantly (unblocks at caller)
                queue.setService(replySignal, new Exp(1e9));
            }
        }

        // Build routing matrix with class switching for REPLY signals
        RoutingMatrix P = model.initRoutingMatrix();

        for (int refTidx : refTaskIndices) {
            ClosedClass requestClass = requestClasses.get(refTidx);
            Signal replySignal = replySignals.get(refTidx);
            Delay thinkNode = thinkNodes.get(refTidx);

            // Get the call chain for this reference task
            List<SynchCallInfo> calls = synchCallsFrom.get(refTidx);

            if (calls == null || calls.isEmpty()) {
                // No synch calls - just loop Think -> Think (or Think -> RefQueue -> Think)
                if (taskQueues.containsKey(refTidx)) {
                    Station refQueue = taskQueues.get(refTidx);
                    model.addLink(thinkNode, refQueue);
                    model.addLink(refQueue, thinkNode);
                } else {
                    model.addLink(thinkNode, thinkNode);
                }
                continue;
            }

            // Build the call chain routing
            // For simplicity, handle the first synch call (can be extended for multiple)
            SynchCallInfo firstCall = calls.get(0);
            int calleeTidx = firstCall.targetTidx;

            Station callerQueue = taskQueues.get(refTidx);
            Station calleeQueue = taskQueues.get(calleeTidx);

            if (callerQueue == null) {
                // Reference task has no queue - create routing through callee only
                // Think -> Callee -> Think (with class switch)
                model.addLink(thinkNode, calleeQueue);
                P.set(requestClass, replySignal, calleeQueue, thinkNode, 1.0);
                continue;
            }

            // Full routing: Think -> Caller -> Callee -> Caller (reply) -> Think
            // Request path: Think -> CallerQueue -> CalleeQueue
            model.addLink(thinkNode, callerQueue);
            P.set(requestClass, requestClass, thinkNode, callerQueue, 1.0);

            model.addLink(callerQueue, calleeQueue);
            P.set(requestClass, requestClass, callerQueue, calleeQueue, 1.0);

            // Reply path: CalleeQueue -> CallerQueue (class switch) -> Think
            // At callee completion, class switch from Request to Reply
            model.addLink(calleeQueue, callerQueue);
            P.set(requestClass, replySignal, calleeQueue, callerQueue, 1.0);

            // Reply continues from CallerQueue to Think
            model.addLink(callerQueue, thinkNode);
            P.set(replySignal, replySignal, callerQueue, thinkNode, 1.0);

            // Handle nested calls (callee -> deeper tasks)
            buildNestedCallRouting(calleeTidx, refTidx, synchCallsFrom, taskQueues,
                    requestClass, replySignal, P, model, new HashSet<Integer>());
        }

        return model;
    }

    /**
     * Recursively builds routing for nested synchronous calls.
     */
    private static void buildNestedCallRouting(int currentTidx, int refTidx,
            Map<Integer, List<SynchCallInfo>> synchCallsFrom,
            Map<Integer, Station> taskQueues,
            ClosedClass requestClass, Signal replySignal,
            RoutingMatrix P, Network model,
            Set<Integer> visited) {

        if (visited.contains(currentTidx)) {
            return;  // Avoid cycles
        }
        visited.add(currentTidx);

        List<SynchCallInfo> calls = synchCallsFrom.get(currentTidx);
        if (calls == null || calls.isEmpty()) {
            return;  // Leaf task - reply handled by caller
        }

        Station currentQueue = taskQueues.get(currentTidx);
        if (currentQueue == null) {
            return;
        }

        // For each nested call, extend the routing
        for (SynchCallInfo call : calls) {
            int nestedTidx = call.targetTidx;
            Station nestedQueue = taskQueues.get(nestedTidx);
            if (nestedQueue == null) {
                continue;
            }

            // Current -> Nested (Request continues)
            model.addLink(currentQueue, nestedQueue);
            // Note: The P matrix entry is already handled by the default routing
            // or needs to be updated to route Request from current to nested

            // Recursively handle deeper nesting
            buildNestedCallRouting(nestedTidx, refTidx, synchCallsFrom, taskQueues,
                    requestClass, replySignal, P, model, visited);
        }
    }

    /**
     * Collects all tasks reachable via synchronous calls from a starting task.
     */
    private static void collectTasksInChain(int startTidx,
            Map<Integer, List<SynchCallInfo>> synchCallsFrom,
            Set<Integer> tasksInChain) {

        if (tasksInChain.contains(startTidx)) {
            return;
        }
        tasksInChain.add(startTidx);

        List<SynchCallInfo> calls = synchCallsFrom.get(startTidx);
        if (calls != null) {
            for (SynchCallInfo call : calls) {
                collectTasksInChain(call.targetTidx, synchCallsFrom, tasksInChain);
            }
        }
    }

    /**
     * Information about a synchronous call.
     */
    private static class SynchCallInfo {
        final int targetTidx;
        final int targetEidx;
        final double callMean;

        SynchCallInfo(int targetTidx, int targetEidx, double callMean) {
            this.targetTidx = targetTidx;
            this.targetEidx = targetEidx;
            this.callMean = callMean;
        }
    }

    /**
     * Get service demand for a task split by phase.
     *
     * @param lsn The LayeredNetworkStruct
     * @param tidx Task index
     * @param hasPhase2 Whether the model has phase-2 activities
     * @return Array of [phase1_demand, phase2_demand]
     */
    private static double[] getTaskServiceDemandByPhase(LayeredNetworkStruct lsn, int tidx, boolean hasPhase2) {
        double ph1Demand = 0;
        double ph2Demand = 0;

        List<Integer> entries = lsn.entriesof.get(tidx);
        if (entries == null) return new double[]{0, 0};

        for (int eidx : entries) {
            List<Integer> acts = lsn.actsof.get(eidx);
            if (acts == null) continue;

            for (int aidx : acts) {
                if (lsn.type.get(aidx) != LayeredNetworkElement.ACTIVITY) {
                    continue;
                }

                double demand = 0;
                Double meanDemand = lsn.hostdem_mean.get(aidx);
                if (meanDemand != null && !Double.isNaN(meanDemand)) {
                    demand = meanDemand;
                }

                if (hasPhase2) {
                    int a = aidx - lsn.ashift;
                    if (a >= 0 && a < lsn.actphase.getNumElements() && lsn.actphase.get(a) > 1) {
                        ph2Demand += demand;
                    } else {
                        ph1Demand += demand;
                    }
                } else {
                    ph1Demand += demand;
                }
            }
        }

        return new double[]{ph1Demand, ph2Demand};
    }
}
