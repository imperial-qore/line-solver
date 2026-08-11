/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.GlobalConstants;
import jline.lang.ClosedClass;
import jline.lang.ClosedSignal;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Region;
import jline.lang.RoutingMatrix;
import jline.lang.Signal;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.CallType;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Router;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * Converts a LayeredNetwork (LQN) to a Network (QN) using REPLY signals.
 *
 * <p>Java port of {@code matlab/src/io/LQN2QN.m}; the two must stay in step.</p>
 *
 * <h3>Construction</h3>
 * <ul>
 *   <li>One station per host processor (scheduling and multiplicity taken from
 *       the processor). Tasks sharing a processor share the station, as in the
 *       LQN semantics where the processor is the contended resource.</li>
 *   <li>One Delay per reference task, holding its think time.</li>
 *   <li>One closed class per step of the expanded activity graph. A step is an
 *       activity, or one call stage of an activity that issues synchronous
 *       calls. Steps of the reference task chain carry population 0 except the
 *       think class, which carries the reference task multiplicity.</li>
 *   <li>A synchronous call site blocks its caller: the step class has a REPLY
 *       signal bound to it, the token proceeds to the callee, and the callee's
 *       replying activity class-switches to that signal, which returns to the
 *       caller station and unblocks it.</li>
 *   <li>Call multiplicity mean m is unrolled into floor(m) mandatory call
 *       stages plus, if m is not integer, one further stage entered with
 *       probability m-floor(m).</li>
 *   <li>OR-branch and loop probabilities are read from the activity graph
 *       weights. Each call site receives its own copy of the callee subgraph.</li>
 *   <li>AND precedences become Fork and Join nodes, with one Router per branch
 *       since a Fork cannot switch class per output link. A branch tail that
 *       issues a synchronous call is given a merge step, so that it reaches the
 *       Join in an ordinary class rather than as a REPLY signal, which carries
 *       no forked-task identity.</li>
 *   <li>A CacheTask becomes a Cache node. The activity bound to an ItemEntry is
 *       the read step and sits on that node; its two CacheAccess successors
 *       become the hit and the miss class, and since the class switch is
 *       performed by the Cache node itself, the routes leaving it are written
 *       in the successor class.</li>
 *   <li>An asynchronous call is lowered to a non-blocking visit: the caller
 *       does not hold its server for the duration of the call, but it is
 *       serialised behind it, since a closed network has no means of creating
 *       the second token that a truly concurrent send would require.</li>
 *   <li>Entry forwarding splits the reply exits of the forwarding entry: with
 *       the forwarding probability the request is handed to the target entry,
 *       which replies to the original caller, so the forwarder is released
 *       while the caller stays blocked.</li>
 *   <li>Phase-2 activities, the successors of a replying activity, run after
 *       the reply: the replying step's exit routes back to the caller, and
 *       each of its service completions spawns the continuation at the host
 *       station (sn.classspawn). The spawned token walks the phase-2 subgraph
 *       holding only the task's own thread and is destroyed at the chain end,
 *       through a NEGATIVE signal that always misses on a closed chain, or
 *       through the Sink on an open one. A boundary that ends on a call
 *       site is normalised through a merge step at the host station; a
 *       phase 2 that opens with an AND-fork spawns into an immediate head
 *       that feeds the Fork; at an AND-join branch tail the spawned token
 *       inherits the fork identity of the trigger and stands in for it at
 *       the Join; at a cache read the reply is emitted by an immediate
 *       trigger step per hit/miss outcome, whose completion spawns the
 *       matching branch continuation.</li>
 *
 *   <li>An AND-join quorum k of n is applied to the Join node in the class
 *       that entered the Fork; k equal to the branch count is the default
 *       wait-for-all and is left alone. An activity think time becomes an
 *       extra step on a shared ActivityThink delay, in series with the host
 *       demand, so the task keeps its thread for it while its processor is
 *       released.</li>
 * </ul>
 *
 * <p>Not yet represented: delayed-hit retrieval on the cache miss path, the
 * thread pool of a task with an internal AND-fork, task and processor
 * replication with its fan-out, and the setup and delay-off times of a
 * function task. Each is reported through line_warning.</p>
 *
 * @see LayeredNetwork
 * @see Network
 * @see SignalType#REPLY
 */
public class LQN2QN {

    private static final int MAXCALLSTAGES = 20;

    /**
     * Converts a LayeredNetwork to an equivalent queueing network using REPLY signals.
     *
     * @param lqn the LayeredNetwork model to convert
     * @return a Network that models the LQN behaviour with REPLY signal blocking
     */
    public static Network convert(LayeredNetwork lqn) {
        return new LQN2QN(lqn).build();
    }

    /** A routing edge of the step graph. */
    private static class Flow {
        final int from;
        final int to;
        final double prob;
        final boolean fromIsSignal;
        final boolean inTargetClass;

        Flow(int from, int to, double prob, boolean fromIsSignal, boolean inTargetClass) {
            this.from = from;
            this.to = to;
            this.prob = prob;
            this.fromIsSignal = fromIsSignal;
            this.inTargetClass = inTargetClass;
        }
    }

    /** A class switch into the reply signal of a blocking call site. */
    private static class Reply {
        final int calleeExitStep;
        final int callerStep;
        final boolean calleeExitIsSignal;
        final double prob;

        Reply(int calleeExitStep, int callerStep, boolean calleeExitIsSignal, double prob) {
            this.calleeExitStep = calleeExitStep;
            this.callerStep = callerStep;
            this.calleeExitIsSignal = calleeExitIsSignal;
            this.prob = prob;
        }
    }

    /** An exit of a step, with the probability of leaving through it. */
    private static class Port {
        final int step;
        final boolean isSignal;
        double prob;

        Port(int step, boolean isSignal, double prob) {
            this.step = step;
            this.isSignal = isSignal;
            this.prob = prob;
        }
    }

    /** Result of expanding one entry in one call context. */
    private static class EntryResult {
        Integer firstStep;
        final List<Port> replyExits = new ArrayList<Port>();
        final List<Port> terminals = new ArrayList<Port>();
    }

    /** One unrolled call stage of an activity. */
    private static class CallStage {
        final int targetEidx;
        final double prob;
        final boolean isasync;

        CallStage(int targetEidx, double prob, boolean isasync) {
            this.targetEidx = targetEidx;
            this.prob = prob;
            this.isasync = isasync;
        }
    }

    /** Cache read/hit/miss wiring, applied once the classes exist. */
    private static class CacheWiring {
        final Cache node;
        final int readStep;
        final int hitStep;
        final int missStep;
        final DiscreteDistribution itemproc;
        final int nitems;

        CacheWiring(Cache node, int readStep, int hitStep, int missStep,
                    DiscreteDistribution itemproc, int nitems) {
            this.node = node;
            this.readStep = readStep;
            this.hitStep = hitStep;
            this.missStep = missStep;
            this.itemproc = itemproc;
            this.nitems = nitems;
        }
    }

    private final LayeredNetworkStruct lsn;
    private final Network model;

    private final Map<Integer, Station> hostStation = new HashMap<Integer, Station>();
    private final Map<Integer, Boolean> hostIsDelay = new HashMap<Integer, Boolean>();
    private final Map<Integer, Delay> thinkNode = new HashMap<Integer, Delay>();

    // Step arrays, all indexed by step id.
    private final List<Integer> stepAidx = new ArrayList<Integer>();
    private final List<Integer> stepHost = new ArrayList<Integer>();
    private final List<Distribution> stepSvc = new ArrayList<Distribution>();
    private final List<String> stepName = new ArrayList<String>();
    private final List<Boolean> stepBlocks = new ArrayList<Boolean>();
    private final List<Boolean> stepIsThink = new ArrayList<Boolean>();
    private final List<Integer> stepRefTask = new ArrayList<Integer>();
    private final List<Node> stepNode = new ArrayList<Node>();
    private final List<Integer> stepClassOwner = new ArrayList<Integer>();

    private final List<Flow> flow = new ArrayList<Flow>();
    private final List<Reply> reply = new ArrayList<Reply>();
    // spawnPairs rows [triggerStep, targetStep]: each service completion of
    // the trigger class spawns a job of the target class at the same station
    // (sn.classspawn), which is how a phase-2 continuation survives the reply.
    private final List<int[]> spawnPairs = new ArrayList<int[]>();
    // ph2Exits rows [step, isSignal, prob, refTidx]: end of a phase-2 chain.
    // The spawned token is destroyed there: on a closed chain it switches into
    // a NEGATIVE signal aimed at a station nothing visits, whose miss
    // annihilates it; on an open chain it leaves through the Sink.
    private final List<double[]> ph2Exits = new ArrayList<double[]>();
    private final List<Integer> entryStack = new ArrayList<Integer>();
    // Tasks holding a thread during the current expansion; tracks entryStack
    // except across forwarding, which hands the request on and releases the
    // forwarder's thread.
    private final List<Integer> threadStack = new ArrayList<Integer>();
    // Thread-pool tasks holding a thread while a job is at each step.
    private final List<int[]> stepTasks = new ArrayList<int[]>();
    // Tasks whose multiplicity is enforced as a thread pool by a finite
    // capacity region; their calls do not hold the caller's server.
    private boolean[] fcrTask = new boolean[0];

    private final Map<Integer, Cache> cacheNodeOf = new HashMap<Integer, Cache>();
    private final List<CacheWiring> cacheWiring = new ArrayList<CacheWiring>();

    // joinQuorum rows [joinStep, joinAidx]: the quorum is applied once the
    // class that entered the fork exists.
    private final List<int[]> joinQuorum = new ArrayList<int[]>();

    // Shared INF station carrying the activity think times: the task keeps its
    // thread across a think time but its host processor is released.
    private Delay actThinkNode;

    private List<JobClass> stepClass;
    private List<Signal> stepSignal;

    private LQN2QN(LayeredNetwork lqn) {
        this.lsn = lqn.getStruct();
        this.model = new Network(lqn.getName() + "-QN");
    }

    private Network build() {
        List<Integer> refTaskIndices = new ArrayList<Integer>();
        for (int t = 1; t <= lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            if (lsn.isref.get(tidx) > 0) {
                refTaskIndices.add(tidx);
            }
        }
        // open-arrival entries: see _kb/06-solver-catalog.md ("LQN2QN: LQN activity graph -> QN step graph")
        List<Integer> openEntries = new ArrayList<Integer>();
        if (lsn.arrival != null) {
            for (int e = lsn.eshift + 1; e <= lsn.eshift + lsn.nentries; e++) {
                Distribution arv = lsn.arrival.get(e);
                if (arv != null && !Double.isInfinite(arv.getMean())
                        && arv.getMean() > GlobalConstants.FineTol) {
                    openEntries.add(e);
                }
            }
        }

        if (refTaskIndices.isEmpty() && openEntries.isEmpty()) {
            line_error(mfilename(new Object() {}),
                    "LQN must have at least one reference task or open arrival.");
            return model;
        }

        warnUnsupported();
        computeThreadPoolTasks();

        // Stations: one per host processor.
        for (int h = 1; h <= lsn.nhosts; h++) {
            double nservers = lsn.mult.get(h);
            SchedStrategy sched = lsn.sched.get(h);
            if (Double.isInfinite(nservers) || nservers >= Integer.MAX_VALUE || sched == SchedStrategy.INF) {
                hostStation.put(h, new Delay(model, lsn.names.get(h)));
                hostIsDelay.put(h, Boolean.TRUE);
            } else {
                Queue q = new Queue(model, lsn.names.get(h), sched);
                q.setNumberOfServers((int) nservers);
                hostStation.put(h, q);
                hostIsDelay.put(h, Boolean.FALSE);
            }
        }

        // Think delays, one per reference task.
        for (int refTidx : refTaskIndices) {
            thinkNode.put(refTidx, new Delay(model, lsn.names.get(refTidx) + "_Think"));
        }

        // Pass 1: expand the activity graph into a step graph.
        for (int refTidx : refTaskIndices) {
            int thinkStep = addStep(0, 0, null, lsn.names.get(refTidx) + "_Think", false, true, refTidx);
            List<Integer> entries = lsn.entriesof.get(refTidx);
            if (entries == null) {
                continue;
            }
            for (int eidx : entries) {
                EntryResult er = expandEntry(eidx, refTidx);
                if (er.firstStep == null) {
                    continue;
                }
                addRoute(new Port(thinkStep, false, 1.0), er.firstStep, 1.0);
                // reference task cycle closure: see _kb/06-solver-catalog.md ("LQN2QN")
                List<Port> exits = new ArrayList<Port>(er.replyExits);
                exits.addAll(er.terminals);
                for (Port p : exits) {
                    addRoute(p, thinkStep, p.prob);
                }
            }
        }

        // Pass 1b: open arrival chains. Wired to Source and Sink in pass 4
        // once the classes exist. Reference task index 0 marks an open step.
        List<Object[]> openWiring = new ArrayList<Object[]>();   // {eidx, firstStep, List<Port> exits}
        Source srcNode = null;
        Sink snkNode = null;
        if (!openEntries.isEmpty()) {
            srcNode = new Source(model, "Source");
            snkNode = new Sink(model, "Sink");
        }
        for (int eidx : openEntries) {
            EntryResult er = expandEntry(eidx, 0);
            if (er.firstStep == null) {
                line_warning(mfilename(new Object() {}), "Open arrival entry "
                        + lsn.names.get(eidx) + " has no bound activity; ignored.");
                continue;
            }
            List<Port> exits = new ArrayList<Port>(er.replyExits);
            exits.addAll(er.terminals);
            openWiring.add(new Object[]{eidx, er.firstStep, exits});
        }

        // Pass 2: create classes and reply signals.
        int nsteps = stepAidx.size();
        stepClass = new ArrayList<JobClass>();
        stepSignal = new ArrayList<Signal>();
        for (int i = 0; i < nsteps; i++) {
            stepClass.add(null);
            stepSignal.add(null);
        }
        for (int i = 0; i < nsteps; i++) {
            if (stepClassOwner.get(i) != i) {
                // Fork, Join and Router steps carry the job through unchanged.
                continue;
            }
            int refTidx = stepRefTask.get(i);
            if (refTidx == 0) {
                // A step of an open arrival chain travels in an open class.
                stepClass.set(i, new OpenClass(model, stepName.get(i), 0));
            } else {
                int population = stepIsThink.get(i) ? (int) lsn.mult.get(refTidx) : 0;
                stepClass.set(i, new ClosedClass(model, stepName.get(i), population, thinkNode.get(refTidx)));
            }
        }
        for (int i = 0; i < nsteps; i++) {
            stepClass.set(i, stepClass.get(stepClassOwner.get(i)));
        }
        // AND-join quorum, in the class the siblings are matched in.
        for (int[] jq : joinQuorum) {
            applyJoinQuorum((Join) stepNode.get(jq[0]), stepClass.get(jq[0]), jq[1]);
        }
        for (int i = 0; i < nsteps; i++) {
            if (stepBlocks.get(i)) {
                Signal sig = new Signal(model, stepName.get(i) + "_Reply", SignalType.REPLY);
                sig.forJobClass(stepClass.get(i));
                stepSignal.set(i, sig);
            }
        }

        // A Signal installs RAND routing at every node; clear it so that link()
        // only honours the routes set below.
        for (int i = 0; i < nsteps; i++) {
            if (stepSignal.get(i) == null) {
                continue;
            }
            for (Node node : model.getNodes()) {
                if (!(node instanceof Sink)) {
                    node.setRouting(stepSignal.get(i), RoutingStrategy.DISABLED);
                }
            }
        }

        // Spawn bindings for phase-2 continuations: each completion of the
        // trigger class injects a job of the target class at the same station.
        for (int[] sp : spawnPairs) {
            int tgtIdx = model.getClasses().indexOf(stepClass.get(sp[1]));
            if (tgtIdx >= 0) {
                stepClass.get(sp[0]).setSpawnClassIndex(tgtIdx + 1);
            }
        }

        // phase-2 token destructor: see _kb/06-solver-catalog.md ("LQN2QN")
        Queue ph2DumpNode = null;
        Map<Integer, ClosedSignal> ph2DestructorOf = new HashMap<Integer, ClosedSignal>();
        for (double[] pe : ph2Exits) {
            int rft = (int) pe[3];
            if (rft <= 0 || ph2DestructorOf.containsKey(rft)) {
                continue;
            }
            if (ph2DumpNode == null) {
                ph2DumpNode = new Queue(model, "Ph2Sink", SchedStrategy.FCFS);
            }
            ClosedSignal sig = new ClosedSignal(model, "Ph2End_" + lsn.names.get(rft),
                    SignalType.NEGATIVE, thinkNode.get(rft), 0);
            for (Node node : model.getNodes()) {
                if (!(node instanceof Sink)) {
                    node.setRouting(sig, RoutingStrategy.DISABLED);
                }
            }
            ph2DumpNode.setService(sig, new Immediate());
            ph2DestructorOf.put(rft, sig);
        }

        // Pass 3: service times.
        for (int i = 0; i < nsteps; i++) {
            if (stepNode.get(i) != null) {
                // Router-hosted merge step service pairing: see _kb/06-solver-catalog.md ("LQN2QN")
                if (stepClassOwner.get(i) == i && stepNode.get(i) == actThinkNode
                        && actThinkNode != null) {
                    actThinkNode.setService(stepClass.get(i), stepSvc.get(i));
                    continue;
                }
                if (stepClassOwner.get(i) == i && stepNode.get(i) instanceof Router) {
                    if (stepRefTask.get(i) == 0) {
                        // Open chain: no think delay exists, declare the pair at
                        // the caller's host station, which the class never visits.
                        hostStation.get(stepHost.get(i)).setService(stepClass.get(i), new Immediate());
                    } else {
                        thinkNode.get(stepRefTask.get(i)).setService(stepClass.get(i), new Immediate());
                    }
                }
                continue;
            }
            if (stepIsThink.get(i)) {
                int refTidx = stepRefTask.get(i);
                Distribution thinkDist = lsn.think == null ? null : lsn.think.get(refTidx);
                Delay tnode = thinkNode.get(refTidx);
                if (isNonTrivial(thinkDist)) {
                    tnode.setService(stepClass.get(i), thinkDist);
                } else {
                    tnode.setService(stepClass.get(i), new Immediate());
                }
            } else {
                Station station = hostStation.get(stepHost.get(i));
                if (stepSvc.get(i) == null) {
                    station.setService(stepClass.get(i), new Immediate());
                } else {
                    station.setService(stepClass.get(i), stepSvc.get(i));
                }
            }
        }
        // reply signal service declaration: see _kb/06-solver-catalog.md ("LQN2QN")
        for (int i = 0; i < nsteps; i++) {
            if (stepSignal.get(i) == null) {
                continue;
            }
            for (int h = 1; h <= lsn.nhosts; h++) {
                hostStation.get(h).setService(stepSignal.get(i), new Immediate());
            }
            for (Delay tn : thinkNode.values()) {
                tn.setService(stepSignal.get(i), new Immediate());
            }
        }

        // Cache read/hit/miss wiring, now that the classes exist.
        for (CacheWiring cw : cacheWiring) {
            cw.node.setReadItemEntry(stepClass.get(cw.readStep), cw.itemproc, cw.nitems);
            cw.node.setHitClass(stepClass.get(cw.readStep), stepClass.get(cw.hitStep));
            cw.node.setMissClass(stepClass.get(cw.readStep), stepClass.get(cw.missStep));
        }

        // Pass 4: routing.
        RoutingMatrix P = model.initRoutingMatrix();
        for (Flow f : flow) {
            if (f.inTargetClass) {
                P.set(stepClass.get(f.to), stepClass.get(f.to), stationOf(f.from), stationOf(f.to), f.prob);
            } else if (f.fromIsSignal) {
                P.set(stepSignal.get(f.from), stepClass.get(f.to), stationOf(f.from), stationOf(f.to), f.prob);
            } else {
                P.set(stepClass.get(f.from), stepClass.get(f.to), stationOf(f.from), stationOf(f.to), f.prob);
            }
        }
        for (Reply r : reply) {
            if (r.calleeExitIsSignal) {
                // A nested call returns through its own reply signal, which
                // class-switches into the reply signal of the outer call site.
                P.set(stepSignal.get(r.calleeExitStep), stepSignal.get(r.callerStep),
                        stationOf(r.calleeExitStep), stationOf(r.callerStep), r.prob);
            } else {
                P.set(stepClass.get(r.calleeExitStep), stepSignal.get(r.callerStep),
                        stationOf(r.calleeExitStep), stationOf(r.callerStep), r.prob);
            }
        }
        // Phase-2 chain ends: destroy the spawned token.
        for (double[] pe : ph2Exits) {
            int i = (int) pe[0];
            boolean viaSig = pe[1] != 0;
            double p = pe[2];
            int rft = (int) pe[3];
            if (rft == 0) {
                // Open chain: the spawned token leaves through the Sink.
                JobClass ecls = stepClass.get(i);
                P.set(ecls, ecls, stationOf(i), snkNode, p);
            } else if (viaSig) {
                P.set(stepSignal.get(i), ph2DestructorOf.get(rft), stationOf(i), ph2DumpNode, p);
            } else {
                P.set(stepClass.get(i), ph2DestructorOf.get(rft), stationOf(i), ph2DumpNode, p);
            }
        }

        // Open arrival wiring: Source into the first step, exits into the Sink.
        for (Object[] ow : openWiring) {
            int eidx = (Integer) ow[0];
            int firstStep = (Integer) ow[1];
            @SuppressWarnings("unchecked")
            List<Port> exits = (List<Port>) ow[2];
            JobClass firstCls = stepClass.get(firstStep);
            srcNode.setArrival(firstCls, lsn.arrival.get(eidx));
            P.set(firstCls, firstCls, srcNode, stationOf(firstStep), 1.0);
            for (Port ep : exits) {
                // Open chains carry no signals, so every exit is an ordinary class.
                JobClass ecls = stepClass.get(ep.step);
                P.set(ecls, ecls, stationOf(ep.step), snkNode, ep.prob);
            }
        }

        model.link(P);
        buildThreadPoolRegion();

        return model;
    }

    /**
     * Thread pools: one finite capacity region, one linear constraint per
     * task, capping the jobs across the task's step classes (its own steps
     * and those of its nested callees) at the task multiplicity. Membership
     * is by station, so the region spans every station such a class visits.
     * Class coefficients are shared across rows where tasks nest, which keeps
     * intra-region class switches admissible: the sums are left unchanged.
     */
    private void buildThreadPoolRegion() {
        List<Integer> fcrList = new ArrayList<Integer>();
        for (int t = 0; t < fcrTask.length; t++) {
            if (fcrTask[t]) {
                fcrList.add(t);
            }
        }
        if (fcrList.isEmpty()) {
            return;
        }
        List<JobClass> classes = model.getClasses();
        int K = classes.size();
        Matrix A = new Matrix(fcrList.size(), K);
        Matrix b = new Matrix(fcrList.size(), 1);
        List<Node> regionNodes = new ArrayList<Node>();
        boolean anyCoeff = false;
        for (int ti = 0; ti < fcrList.size(); ti++) {
            int tidx = fcrList.get(ti);
            for (int i = 0; i < stepTasks.size(); i++) {
                boolean holds = false;
                for (int ht : stepTasks.get(i)) {
                    if (ht == tidx) {
                        holds = true;
                        break;
                    }
                }
                if (!holds) {
                    continue;
                }
                int cIdx = classes.indexOf(stepClass.get(i));
                if (cIdx >= 0) {
                    A.set(ti, cIdx, 1.0);
                    anyCoeff = true;
                }
                Node nd = stationOf(i);
                if (nd instanceof Station && !regionNodes.contains(nd)) {
                    regionNodes.add(nd);
                }
            }
            b.set(ti, 0, lsn.mult.get(tidx));
        }
        if (!anyCoeff || regionNodes.isEmpty()) {
            return;
        }
        Region fcr = model.addRegion(regionNodes);
        fcr.setLinearConstraints(A, b);
    }

    // ------------------------------------------------------------- step graph

    private Node stationOf(int i) {
        if (stepNode.get(i) != null) {
            return stepNode.get(i);
        }
        if (stepIsThink.get(i)) {
            return thinkNode.get(stepRefTask.get(i));
        }
        return (Node) hostStation.get(stepHost.get(i));
    }

    private int addStep(int aidx, int hidx, Distribution svc, String name,
                        boolean blocks, boolean isthink, int refTidx) {
        stepAidx.add(aidx);
        stepHost.add(hidx);
        stepSvc.add(svc);
        stepName.add(name);
        stepBlocks.add(blocks);
        stepIsThink.add(isthink);
        stepRefTask.add(refTidx);
        stepNode.add(null);
        int id = stepAidx.size() - 1;
        stepClassOwner.add(id);
        // thread-pool holding: see _kb/06-solver-catalog.md ("LQN2QN")
        java.util.TreeSet<Integer> held = new java.util.TreeSet<Integer>();
        for (int ti : threadStack) {
            if (ti >= 0 && ti < fcrTask.length && fcrTask[ti]) {
                held.add(ti);
            }
        }
        int[] heldArr = new int[held.size()];
        int hi = 0;
        for (int ti : held) {
            heldArr[hi++] = ti;
        }
        stepTasks.add(heldArr);
        return id;
    }

    /** A step on a Fork, Join or Router node: no station, no service, no class of its own. */
    private int addAuxStep(Node nodeObj, int ownerStep, String name, int refTidx) {
        int id = addStep(0, 0, null, name, false, false, refTidx);
        stepNode.set(id, nodeObj);
        stepClassOwner.set(id, stepClassOwner.get(ownerStep));
        return id;
    }

    private void addRoute(Port fromPort, int toStep, double prob) {
        flow.add(new Flow(fromPort.step, toStep, prob, fromPort.isSignal, false));
    }

    /**
     * Leaving a Cache node: the switch into the hit or the miss class is made
     * by the node, so the route is declared in the target class.
     */
    private void addCacheRoute(int fromStep, int toStep) {
        flow.add(new Flow(fromStep, toStep, 1.0, false, true));
    }

    // ------------------------------------------------------------- expansion

    /** Expands the activity subgraph bound to an entry, in the current call context. */
    private EntryResult expandEntry(int eidx, int refTidx) {
        EntryResult res = new EntryResult();

        if (entryStack.contains(Integer.valueOf(eidx))) {
            line_warning(mfilename(new Object() {}),
                    "Recursive call cycle at entry " + lsn.names.get(eidx) + " truncated.");
            return res;
        }
        entryStack.add(Integer.valueOf(eidx));
        threadStack.add(Integer.valueOf((int) lsn.parent.get(eidx)));
        try {
            List<Integer> localActs = lsn.actsof.get(eidx);
            if (localActs == null || localActs.isEmpty()) {
                return res;
            }
            // The bound activity is the activity successor of the entry.
            Integer bound = null;
            for (int a : localActs) {
                if (a < lsn.type.getNumElements() && (int) lsn.type.get(a) == LayeredNetworkElement.ACTIVITY
                        && lsn.graph.get(eidx, a) != 0) {
                    bound = a;
                    break;
                }
            }
            if (bound == null) {
                return res;
            }

            EntryWalker walker = new EntryWalker(eidx, refTidx, localActs);
            walker.run(bound, res);

            // forwarding rules: see _kb/06-solver-catalog.md ("LQN2QN")
            List<double[]> fwd = forwardingOf(eidx);
            if (!fwd.isEmpty() && !res.replyExits.isEmpty()) {
                List<Port> ownPorts = new ArrayList<Port>(res.replyExits);
                List<Port> fwdExits = new ArrayList<Port>();
                double pforw = 0.0;
                // forwarder thread release: see _kb/06-solver-catalog.md ("LQN2QN")
                Integer fwdThread = threadStack.remove(threadStack.size() - 1);
                for (double[] f : fwd) {
                    EntryResult fr = expandEntry((int) f[0], refTidx);
                    if (fr.firstStep == null) {
                        continue;
                    }
                    double p = f[1];
                    for (Port op : ownPorts) {
                        addRoute(op, fr.firstStep, op.prob * p);
                    }
                    pforw += p;
                    fwdExits.addAll(fr.replyExits);
                    fwdExits.addAll(fr.terminals);
                }
                threadStack.add(fwdThread);
                // What is left of each of this entry's own ports still replies.
                double residual = Math.max(0.0, 1.0 - pforw);
                for (Port op : res.replyExits) {
                    op.prob = op.prob * residual;
                }
                res.replyExits.addAll(fwdExits);
            }
            return res;
        } finally {
            entryStack.remove(entryStack.size() - 1);
            threadStack.remove(threadStack.size() - 1);
        }
    }

    /** Rows [targetEidx, probability] of the forwarding calls of an entry. */
    private List<double[]> forwardingOf(int eidx) {
        List<double[]> fwd = new ArrayList<double[]>();
        if (lsn.calltype == null || lsn.callpair == null) {
            return fwd;
        }
        for (int cidx = 1; cidx <= lsn.ncalls; cidx++) {
            CallType ct = lsn.calltype.get(cidx);
            if (ct != CallType.FWD || (int) lsn.callpair.get(cidx, 1) != eidx) {
                continue;
            }
            double p = callMean(cidx);
            p = Math.min(Math.max(p, 0.0), 1.0);
            if (p > GlobalConstants.FineTol) {
                fwd.add(new double[]{lsn.callpair.get(cidx, 2), p});
            }
        }
        return fwd;
    }

    /**
     * Walks the intra-task activity graph of one entry, creating steps and
     * expanding every synchronous call site.
     */
    private class EntryWalker {
        private final int eidx;
        private final int refTidx;
        private final List<Integer> localActs;
        private final Map<Integer, Port> visitedExit = new HashMap<Integer, Port>();
        private final Map<Integer, Integer> visitedEntry = new HashMap<Integer, Integer>();
        private final Map<Integer, Integer> joinOf = new HashMap<Integer, Integer>();
        private final List<Integer> forkOwnerStack = new ArrayList<Integer>();
        private boolean sawReply = false;
        private EntryResult res;

        EntryWalker(int eidx, int refTidx, List<Integer> localActs) {
            this.eidx = eidx;
            this.refTidx = refTidx;
            this.localActs = localActs;
        }

        void run(int a0, EntryResult res) {
            this.res = res;
            res.firstStep = walk(a0);

            // chain-end reply-vs-terminal rule: see _kb/06-solver-catalog.md ("LQN2QN")
            if (sawReply) {
                res.replyExits.addAll(res.terminals);
                res.terminals.clear();
            }
        }

        /** Returns the entry step of the activity; its exit port is recorded in visitedExit. */
        private int walk(int aidx) {
            if (visitedEntry.containsKey(aidx)) {
                return visitedEntry.get(aidx);
            }

            int tidx = (int) lsn.parent.get(aidx);

            // The activity bound to an ItemEntry of a CacheTask is the read
            // step: it sits on the Cache node rather than on the processor.
            Cache cacheNode = null;
            if (lsn.iscache != null && tidx < lsn.iscache.getNumElements()
                    && lsn.iscache.get(tidx) != 0 && lsn.graph.get(eidx, aidx) != 0) {
                cacheNode = getCacheNode(tidx);
            }

            int entryStep = makeActivitySteps(aidx, tidx, refTidx, cacheNode);
            Port exitPort = visitedExit.get(aidx);
            visitedEntry.put(aidx, entryStep);

            // reply deferred to phase-2 handling: see _kb/06-solver-catalog.md ("LQN2QN")
            if (repliesHere(aidx)) {
                sawReply = true;
            }

            // Local successors within the same entry.
            List<Integer> succ = new ArrayList<Integer>();
            for (int s : localActs) {
                if (lsn.graph.get(aidx, s) != 0) {
                    succ.add(s);
                }
            }
            if (succ.isEmpty()) {
                res.terminals.add(new Port(exitPort.step, exitPort.isSignal, 1.0));
                return entryStep;
            }

            // phase-2 spawn mechanics: see _kb/06-solver-catalog.md ("LQN2QN")
            if (repliesHere(aidx)) {
                if (cacheNode != null && succ.size() >= 2) {
                    // phase-2 opening at a cache read: see _kb/06-solver-catalog.md ("LQN2QN")
                    int trigH = addStep(aidx, (int) lsn.parent.get(tidx), null,
                            lsn.names.get(aidx) + "_ph2h", false, false, refTidx);
                    int trigM = addStep(aidx, (int) lsn.parent.get(tidx), null,
                            lsn.names.get(aidx) + "_ph2m", false, false, refTidx);
                    addCacheRoute(entryStep, trigH);
                    addCacheRoute(entryStep, trigM);
                    cacheWiring.add(new CacheWiring(cacheNode, entryStep, trigH, trigM,
                            lsn.itemproc == null ? null : lsn.itemproc.get(eidx),
                            (int) lsn.nitems.get(eidx)));
                    res.replyExits.add(new Port(trigH, false, 1.0));
                    res.replyExits.add(new Port(trigM, false, 1.0));
                    List<Integer> savedStack = new ArrayList<Integer>(threadStack);
                    threadStack.clear();
                    threadStack.add(tidx);
                    int nT0 = res.terminals.size();
                    for (int hm = 0; hm < 2; hm++) {
                        int sEntry = walk(succ.get(hm));
                        if (stepNode.get(sEntry) != null) {
                            int hmHead = addStep(aidx, (int) lsn.parent.get(tidx), null,
                                    lsn.names.get(aidx) + "_ph2b" + (hm + 1), false, false, refTidx);
                            addRoute(new Port(hmHead, false, 1.0), sEntry, 1.0);
                            sEntry = hmHead;
                        }
                        spawnPairs.add(new int[]{hm == 0 ? trigH : trigM, sEntry});
                    }
                    while (res.terminals.size() > nT0) {
                        Port t = res.terminals.remove(nT0);
                        ph2Exits.add(new double[]{t.step, t.isSignal ? 1.0 : 0.0, t.prob, refTidx});
                    }
                    threadStack.clear();
                    threadStack.addAll(savedStack);
                    return entryStep;
                }
                // phase-2 lift restricted to branch tail: see _kb/06-solver-catalog.md ("LQN2QN")
                boolean okCtx = cacheNode == null
                        && (forkOwnerStack.isEmpty() || isAndJoinPre(aidx));
                // call-site phase-2 boundary normalization: see _kb/06-solver-catalog.md ("LQN2QN")
                if (okCtx && (exitPort.isSignal
                        || stepNode.get(exitPort.step) instanceof Router)) {
                    int trig = addStep(aidx, (int) lsn.parent.get(tidx), null,
                            lsn.names.get(aidx) + "_ph2t", false, false, refTidx);
                    addRoute(exitPort, trig, 1.0);
                    exitPort = new Port(trig, false, 1.0);
                }
                boolean ph2Spawn = okCtx && !exitPort.isSignal
                        && stepNode.get(exitPort.step) == null;
                if (!ph2Spawn) {
                    line_warning(mfilename(new Object() {}), "Phase-2 activities of "
                            + lsn.names.get(aidx) + " run before the reply: the boundary "
                            + "is not a station departure, a degenerate cache read, or "
                            + "mid-branch inside an AND-fork.");
                } else {
                    res.replyExits.add(new Port(exitPort.step, exitPort.isSignal, 1.0));
                    List<Integer> savedStack = new ArrayList<Integer>(threadStack);
                    threadStack.clear();
                    threadStack.add(tidx);
                    int nT0 = res.terminals.size();
                    List<Integer> posSucc = new ArrayList<Integer>();
                    for (int s : succ) {
                        if (lsn.graph.get(aidx, s) > 0) {
                            posSucc.add(s);
                        }
                    }
                    Integer target = null;
                    if (isAndFork(succ)) {
                        // phase-2 AND-fork spawn: see _kb/06-solver-catalog.md ("LQN2QN")
                        int head = addStep(aidx, (int) lsn.parent.get(tidx), null,
                                lsn.names.get(aidx) + "_ph2", false, false, refTidx);
                        wireAndFork(new Port(head, false, 1.0), succ, aidx);
                        target = head;
                    } else if (isAndJoinPre(aidx)) {
                        // phase-2 AND-join branch-tail spawn: see _kb/06-solver-catalog.md ("LQN2QN")
                        int head = addStep(aidx, (int) lsn.parent.get(tidx), null,
                                lsn.names.get(aidx) + "_ph2", false, false, refTidx);
                        wireAndJoin(new Port(head, false, 1.0), succ.get(0));
                        target = head;
                    } else if (posSucc.size() == 1) {
                        int sEntry = walk(posSucc.get(0));
                        if (stepNode.get(sEntry) == null) {
                            target = sEntry;
                        }
                    }
                    if (target == null) {
                        // phase-2 branching spawn head: see _kb/06-solver-catalog.md ("LQN2QN")
                        int head = addStep(aidx, (int) lsn.parent.get(tidx), null,
                                lsn.names.get(aidx) + "_ph2", false, false, refTidx);
                        for (int s2 : posSucc) {
                            int sEntry = walk(s2);
                            addRoute(new Port(head, false, 1.0), sEntry, lsn.graph.get(aidx, s2));
                        }
                        target = head;
                    }
                    spawnPairs.add(new int[]{exitPort.step, target});
                    while (res.terminals.size() > nT0) {
                        Port t = res.terminals.remove(nT0);
                        ph2Exits.add(new double[]{t.step, t.isSignal ? 1.0 : 0.0, t.prob, refTidx});
                    }
                    threadStack.clear();
                    threadStack.addAll(savedStack);
                    return entryStep;
                }
            }

            if (cacheNode != null) {
                // cache-read hit/miss precedence: see _kb/06-solver-catalog.md ("LQN2QN")
                if (succ.size() < 2) {
                    line_warning(mfilename(new Object() {}), "Cache read " + lsn.names.get(aidx)
                            + " has no hit/miss pair; treated as an ordinary activity.");
                } else {
                    int hEntry = walk(succ.get(0));
                    int mEntry = walk(succ.get(1));
                    addCacheRoute(entryStep, hEntry);
                    addCacheRoute(entryStep, mEntry);
                    cacheWiring.add(new CacheWiring(cacheNode, entryStep, hEntry, mEntry,
                            lsn.itemproc == null ? null : lsn.itemproc.get(eidx),
                            (int) lsn.nitems.get(eidx)));
                    return entryStep;
                }
            }

            if (isAndFork(succ)) {
                wireAndFork(exitPort, succ, aidx);
                return entryStep;
            }

            if (isAndJoinPre(aidx)) {
                wireAndJoin(exitPort, succ.get(0));
                return entryStep;
            }

            for (int s : succ) {
                double p = lsn.graph.get(aidx, s);
                if (p <= 0) {
                    continue;
                }
                int sEntry = walk(s);
                addRoute(exitPort, sEntry, p);
            }
            return entryStep;
        }

        private void wireAndFork(Port from, List<Integer> succ, int aidx) {
            // AND-fork wiring: see _kb/06-solver-catalog.md ("LQN2QN")
            Fork forkNode = new Fork(model, "Fork_" + lsn.names.get(aidx));
            int forkStep = addAuxStep(forkNode, from.step, "Fork_" + lsn.names.get(aidx), refTidx);
            addRoute(from, forkStep, 1.0);
            forkOwnerStack.add(forkStep);
            // replying-branch-first walk order: see _kb/06-solver-catalog.md ("LQN2QN")
            List<Integer> ordered = new ArrayList<Integer>();
            for (int s : succ) {
                if (branchReplies(s)) {
                    ordered.add(s);
                }
            }
            for (int s : succ) {
                if (!branchReplies(s)) {
                    ordered.add(s);
                }
            }
            for (int b = 0; b < ordered.size(); b++) {
                String rname = "Fork_" + lsn.names.get(aidx) + "_" + (b + 1);
                Router routerNode = new Router(model, rname);
                int routerStep = addAuxStep(routerNode, forkStep, rname, refTidx);
                addRoute(new Port(forkStep, false, 1.0), routerStep, 1.0);
                int sEntry = walk(ordered.get(b));
                addRoute(new Port(routerStep, false, 1.0), sEntry, 1.0);
            }
            forkOwnerStack.remove(forkOwnerStack.size() - 1);
        }

        private void wireAndJoin(Port fromPort, int joinAidx) {
            // AND-join wiring: see _kb/06-solver-catalog.md ("LQN2QN")
            if (joinOf.containsKey(joinAidx)) {
                addRoute(fromPort, joinOf.get(joinAidx), 1.0);
                return;
            }
            if (forkOwnerStack.isEmpty()) {
                line_warning(mfilename(new Object() {}), "AND-join at " + lsn.names.get(joinAidx)
                        + " has no enclosing AND-fork; branches are serialised.");
                int sEntry = walk(joinAidx);
                addRoute(fromPort, sEntry, 1.0);
                return;
            }
            int forkOwner = forkOwnerStack.get(forkOwnerStack.size() - 1);
            Join joinNode = new Join(model, "Join_" + lsn.names.get(joinAidx), stepNode.get(forkOwner));
            int joinStep = addAuxStep(joinNode, forkOwner, "Join_" + lsn.names.get(joinAidx), refTidx);
            addRoute(fromPort, joinStep, 1.0);
            int sEntry = walk(joinAidx);
            addRoute(new Port(joinStep, false, 1.0), sEntry, 1.0);
            joinOf.put(joinAidx, joinStep);
            joinQuorum.add(new int[]{joinStep, joinAidx});
        }

        private boolean branchReplies(int a0) {
            // branch-replies search: see _kb/06-solver-catalog.md ("LQN2QN")
            List<Integer> stack = new ArrayList<Integer>();
            Set<Integer> seen = new HashSet<Integer>();
            stack.add(a0);
            while (!stack.isEmpty()) {
                int a = stack.remove(stack.size() - 1);
                if (!seen.add(a)) {
                    continue;
                }
                if (repliesHere(a)) {
                    return true;
                }
                if (isAndJoinPre(a)) {
                    continue;
                }
                for (int s : localActs) {
                    if (lsn.graph.get(a, s) != 0) {
                        stack.add(s);
                    }
                }
            }
            return false;
        }

        private boolean repliesHere(int aidx) {
            if (lsn.replygraph == null || lsn.replygraph.isEmpty()) {
                return false;
            }
            int a = aidx - lsn.ashift;
            int e = eidx - lsn.eshift;
            if (a < 1 || a >= lsn.replygraph.getNumRows() || e < 1 || e >= lsn.replygraph.getNumCols()) {
                return false;
            }
            return lsn.replygraph.get(a, e) != 0;
        }

        /**
         * One step for the host demand, plus one step per unrolled synchronous
         * call stage. The exit port of the activity is left in visitedExit.
         */
        private int makeActivitySteps(int aidx, int tidx, int refTidx, Cache cacheNode) {
            int hidx = (int) lsn.parent.get(tidx);

            if (cacheNode != null) {
                // cache-read step: see _kb/06-solver-catalog.md ("LQN2QN")
                int readStep = addStep(aidx, hidx, null, lsn.names.get(aidx), false, false, refTidx);
                stepNode.set(readStep, cacheNode);
                List<Integer> calls = lsn.callsof.get(aidx);
                if (calls != null && !calls.isEmpty()) {
                    line_warning(mfilename(new Object() {}),
                            "Calls issued by cache read activity " + lsn.names.get(aidx) + " are ignored.");
                }
                if (isNonTrivial(lsn.hostdem == null ? null : lsn.hostdem.get(aidx))) {
                    line_warning(mfilename(new Object() {}),
                            "Host demand of cache read activity " + lsn.names.get(aidx) + " is ignored.");
                }
                visitedExit.put(aidx, new Port(readStep, false, 1.0));
                return readStep;
            }

            Distribution svc = null;
            Distribution d = lsn.hostdem == null ? null : lsn.hostdem.get(aidx);
            if (isNonTrivial(d)) {
                svc = d;
            }

            List<CallStage> callStages = synchCallStages(aidx);
            int entryStep = addStep(aidx, hidx, svc, lsn.names.get(aidx), false, false, refTidx);
            // exit port semantics: see _kb/06-solver-catalog.md ("LQN2QN")
            Port cur = new Port(entryStep, false, 1.0);

            // activity think time: see _kb/06-solver-catalog.md ("LQN2QN")
            Distribution think = actThinkOf(aidx);
            if (think != null) {
                int thinkStep = addStep(aidx, hidx, think, lsn.names.get(aidx) + "_think",
                        false, false, refTidx);
                stepNode.set(thinkStep, actThinkStation());
                addRoute(cur, thinkStep, 1.0);
                cur = new Port(thinkStep, false, 1.0);
            }

            // call-blocking eligibility rules: see _kb/06-solver-catalog.md ("LQN2QN")
            boolean hostBlocks = !Boolean.TRUE.equals(hostIsDelay.get(hidx))
                    && !(tidx < fcrTask.length && fcrTask[tidx]) && refTidx != 0
                    && Double.isFinite(lsn.mult.get(tidx)) && lsn.mult.get(tidx) < Integer.MAX_VALUE
                    && lsn.sched.get(tidx) != SchedStrategy.INF;

            for (int k = 0; k < callStages.size(); k++) {
                CallStage stage = callStages.get(k);
                // async call server release: see _kb/06-solver-catalog.md ("LQN2QN")
                boolean blocks = hostBlocks && !stage.isasync;
                // async call thread release: see _kb/06-solver-catalog.md ("LQN2QN")
                Integer asyncThread = null;
                if (stage.isasync) {
                    asyncThread = threadStack.remove(threadStack.size() - 1);
                }
                EntryResult callee = expandEntry(stage.targetEidx, refTidx);
                if (stage.isasync) {
                    threadStack.add(asyncThread);
                }
                if (callee.firstStep == null) {
                    continue;   // callee not expandable: drop the call, never block
                }
                // dead-end callee return: see _kb/06-solver-catalog.md ("LQN2QN")
                List<Port> calleeReplies = new ArrayList<Port>(callee.replyExits);
                calleeReplies.addAll(callee.terminals);

                // merge-step necessity rules: see _kb/06-solver-catalog.md ("LQN2QN")
                boolean needsMerge = (stage.prob < 1.0) || (k < callStages.size() - 1)
                        || !blocks || isAndJoinPre(aidx);
                int nxt = -1;
                if (needsMerge) {
                    nxt = addStep(aidx, hidx, null, lsn.names.get(aidx) + "_c" + (k + 1) + "_ret",
                            false, false, refTidx);
                    if (!blocks) {
                        // non-blocking merge on a Router: see _kb/06-solver-catalog.md ("LQN2QN")
                        stepNode.set(nxt, new Router(model,
                                lsn.names.get(aidx) + "_c" + (k + 1) + "_ret"));
                    }
                }

                if (blocks) {
                    // first mandatory call bound to service class: see _kb/06-solver-catalog.md ("LQN2QN")
                    int blk;
                    if (k == 0 && stage.prob >= 1.0 && cur.step == entryStep && !cur.isSignal) {
                        blk = entryStep;
                        stepBlocks.set(blk, Boolean.TRUE);
                    } else {
                        blk = addStep(aidx, hidx, null, lsn.names.get(aidx) + "_c" + (k + 1),
                                true, false, refTidx);
                        addRoute(cur, blk, stage.prob);
                        if (stage.prob < 1.0) {
                            addRoute(cur, nxt, 1.0 - stage.prob);
                        }
                    }
                    addRoute(new Port(blk, false, 1.0), callee.firstStep, 1.0);
                    for (Port r : calleeReplies) {
                        reply.add(new Reply(r.step, blk, r.isSignal, r.prob));
                    }
                    if (needsMerge) {
                        addRoute(new Port(blk, true, 1.0), nxt, 1.0);
                        cur = new Port(nxt, false, 1.0);
                    } else {
                        cur = new Port(blk, true, 1.0);
                    }
                } else {
                    addRoute(cur, callee.firstStep, stage.prob);
                    if (stage.prob < 1.0) {
                        addRoute(cur, nxt, 1.0 - stage.prob);
                    }
                    for (Port r : calleeReplies) {
                        addRoute(new Port(r.step, r.isSignal, r.prob), nxt, r.prob);
                    }
                    cur = new Port(nxt, false, 1.0);
                }
            }
            visitedExit.put(aidx, cur);
            return entryStep;
        }
    }

    // ------------------------------------------------------------- utilities

    /**
     * Think time of an activity, or null when it has none. It is a delay in
     * series with the activity's host demand, held at the activity's own task
     * (the thread is kept) but with the host processor released, mirroring lqns.
     */
    private Distribution actThinkOf(int aidx) {
        if (lsn.actthink == null) {
            return null;
        }
        return isNonTrivial(lsn.actthink.get(aidx)) ? lsn.actthink.get(aidx) : null;
    }

    /** Single INF station shared by every activity think time. */
    private Delay actThinkStation() {
        if (actThinkNode == null) {
            actThinkNode = new Delay(model, "ActivityThink");
        }
        return actThinkNode;
    }

    /** One Cache node per CacheTask, named so as not to collide with the processor station. */
    private Cache getCacheNode(int tidx) {
        Cache cnode = cacheNodeOf.get(tidx);
        if (cnode != null) {
            return cnode;
        }
        int[] levelCaps = lsn.itemcap == null ? null : lsn.itemcap.get(tidx);
        if (levelCaps == null || levelCaps.length == 0) {
            levelCaps = new int[]{1};
        }
        Matrix capMatrix = new Matrix(1, levelCaps.length);
        for (int i = 0; i < levelCaps.length; i++) {
            capMatrix.set(0, i, levelCaps[i]);
        }
        cnode = new Cache(model, lsn.names.get(tidx) + "_Cache", (int) lsn.nitems.get(tidx),
                capMatrix, replacementOf(tidx));
        cacheNodeOf.put(tidx, cnode);
        return cnode;
    }

    private ReplacementStrategy replacementOf(int tidx) {
        int id = lsn.replacestrat == null ? 0 : (int) lsn.replacestrat.get(tidx);
        switch (id) {
            case 1: return ReplacementStrategy.FIFO;
            case 2: return ReplacementStrategy.SFIFO;
            case 3: return ReplacementStrategy.LRU;
            case 4: return ReplacementStrategy.HLRU;
            case 5: return ReplacementStrategy.CLIMB;
            case 6: return ReplacementStrategy.QLRU;
            default: return ReplacementStrategy.RR;
        }
    }

    /** An AND-fork is a precedence whose post activities are all marked POST_AND. */
    private boolean isAndFork(List<Integer> succ) {
        if (succ.size() < 2 || lsn.actposttype == null || lsn.actposttype.isEmpty()) {
            return false;
        }
        for (int s : succ) {
            if (s >= lsn.actposttype.getNumElements()
                    || (int) lsn.actposttype.get(s) != ActivityPrecedenceType.ID_POST_AND) {
                return false;
            }
        }
        return true;
    }

    /**
     * Applies the AND-join quorum of a join target to its Join node. A join whose
     * quorum equals its branch count already waits for all branches, which is the
     * default JoinStrategy.STD, so only a genuine quorum k &lt; n is set.
     */
    private void applyJoinQuorum(Join joinNode, JobClass joinClass, int joinAidx) {
        if (joinNode == null || joinClass == null || lsn.actquorum == null
                || joinAidx < 1 || joinAidx >= lsn.actquorum.getNumCols()) {
            return;
        }
        int quorum = (int) lsn.actquorum.get(0, joinAidx);
        int nbranches = countAndJoinBranches(joinAidx);
        if (quorum < 1 || nbranches < 1 || quorum >= nbranches) {
            return;
        }
        joinNode.setStrategy(joinClass, JoinStrategy.Quorum);
        joinNode.setRequired(joinClass, quorum);
    }

    /** Number of branch tails feeding an AND-join, i.e. its PRE_AND predecessors. */
    private int countAndJoinBranches(int joinAidx) {
        if (lsn.graph == null) {
            return 0;
        }
        int count = 0;
        for (int pred = 1; pred < lsn.graph.getNumRows(); pred++) {
            if (pred != joinAidx && lsn.graph.get(pred, joinAidx) != 0 && isAndJoinPre(pred)) {
                count++;
            }
        }
        return count;
    }

    /** An activity marked PRE_AND is one branch tail of an AND-join. */
    private boolean isAndJoinPre(int aidx) {
        return lsn.actpretype != null && !lsn.actpretype.isEmpty()
                && aidx < lsn.actpretype.getNumElements()
                && (int) lsn.actpretype.get(aidx) == ActivityPrecedenceType.ID_PRE_AND;
    }

    private double callMean(int cidx) {
        if (lsn.callproc_mean != null) {
            Double m = lsn.callproc_mean.get(cidx);
            if (m != null && !Double.isNaN(m)) {
                return m;
            }
        }
        if (lsn.callproc != null) {
            DiscreteDistribution d = lsn.callproc.get(cidx);
            if (d != null) {
                return d.getMean();
            }
        }
        return 1.0;
    }

    /** Unrolls the synchronous and asynchronous calls of an activity into stages. */
    private List<CallStage> synchCallStages(int aidx) {
        List<CallStage> stages = new ArrayList<CallStage>();
        List<Integer> calls = lsn.callsof.get(aidx);
        if (calls == null || calls.isEmpty()) {
            return stages;
        }
        for (int cidx : calls) {
            CallType ct = lsn.calltype.get(cidx);
            if (ct != CallType.SYNC && ct != CallType.ASYNC) {
                continue;
            }
            boolean isasync = (ct == CallType.ASYNC);
            int targetEidx = (int) lsn.callpair.get(cidx, 2);
            double m = callMean(cidx);
            int nfull = (int) Math.floor(m + GlobalConstants.FineTol);
            double frac = m - nfull;
            if (nfull > MAXCALLSTAGES) {
                line_warning(mfilename(new Object() {}), "Call multiplicity " + m + " on "
                        + lsn.callnames.get(cidx) + " truncated to " + MAXCALLSTAGES + " stages.");
                nfull = MAXCALLSTAGES;
                frac = 0;
            }
            for (int k = 0; k < nfull; k++) {
                stages.add(new CallStage(targetEidx, 1.0, isasync));
            }
            if (frac > GlobalConstants.FineTol) {
                stages.add(new CallStage(targetEidx, frac, isasync));
            }
        }
        return stages;
    }

    private static boolean isNonTrivial(Distribution d) {
        return d != null && !(d instanceof Immediate) && d.getMean() > GlobalConstants.FineTol;
    }

    private void warnUnsupported() {
        if (lsn.calltype != null) {
            for (int cidx = 1; cidx <= lsn.ncalls; cidx++) {
                if (lsn.calltype.get(cidx) == CallType.ASYNC) {
                    line_warning(mfilename(new Object() {}), "Asynchronous calls are represented by "
                            + "LQN2QN as non-blocking visits: the caller releases its server but "
                            + "remains serialised behind the callee.");
                    break;
                }
            }
        }
        if (lsn.hasretrieval != null && !lsn.hasretrieval.isEmpty()) {
            for (int i = 0; i < lsn.hasretrieval.getNumElements(); i++) {
                if (lsn.hasretrieval.get(i) != 0) {
                    line_warning(mfilename(new Object() {}), "Delayed-hit retrieval on the cache "
                            + "miss path is not represented by LQN2QN.");
                    break;
                }
            }
        }
        if (lsn.repl != null && !lsn.repl.isEmpty()) {
            for (int idx = 1; idx <= lsn.nhosts + lsn.ntasks; idx++) {
                if (idx < lsn.repl.getNumCols() && lsn.repl.get(0, idx) > 1) {
                    line_warning(mfilename(new Object() {}), "Replication of " + lsn.names.get(idx)
                            + " is not represented by LQN2QN: the replicas are collapsed into a "
                            + "single station and their fan-out is ignored.");
                    break;
                }
            }
        }
        if (lsn.isfunction != null && !lsn.isfunction.isEmpty()) {
            for (int idx = 0; idx < lsn.isfunction.getNumElements(); idx++) {
                if (lsn.isfunction.get(idx) != 0) {
                    line_warning(mfilename(new Object() {}), "Setup and delay-off times of function "
                            + "tasks are not represented by LQN2QN: the task is converted as an "
                            + "ordinary always-on station.");
                    break;
                }
            }
        }
    }

    /**
     * Marks the tasks whose multiplicity is a thread pool. Such a task holds
     * one thread per request from entry to reply, also across its nested
     * synchronous calls. Its calls do not hold the caller's server, and a
     * finite capacity region caps the jobs across its step classes at the
     * task multiplicity instead. Excluded, with a warning: a task with an
     * internal AND-fork (its forked siblings would double-count the thread
     * that spawned them), and a task called from inside an AND-fork branch
     * (the fork-join transformation retags branch flows into auxiliary
     * classes, which would silently bypass the task's admission row).
     */
    private void computeThreadPoolTasks() {
        fcrTask = new boolean[lsn.tshift + lsn.ntasks + 1];
        for (int t = 1; t <= lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            if (lsn.isref.get(tidx) > 0) {
                continue;
            }
            if (lsn.iscache != null && !lsn.iscache.isEmpty()
                    && tidx < lsn.iscache.getNumElements() && lsn.iscache.get(tidx) != 0) {
                continue;
            }
            double m = lsn.mult.get(tidx);
            if (Double.isInfinite(m) || m >= Integer.MAX_VALUE || lsn.sched.get(tidx) == SchedStrategy.INF) {
                continue;
            }
            if (taskHasAndFork(tidx)) {
                line_warning(mfilename(new Object() {}), "Multiplicity of task "
                        + lsn.names.get(tidx) + " is not enforced: an AND-fork inside a task "
                        + "cannot be capped by a finite capacity region, whose job count would "
                        + "double-count the forked siblings.");
                continue;
            }
            fcrTask[tidx] = true;
        }

        // Activities inside AND-fork branch bodies, walked from each POST_AND
        // head up to and including the PRE_AND tails.
        List<Integer> branchActs = new ArrayList<Integer>();
        if (lsn.actposttype != null && !lsn.actposttype.isEmpty()) {
            for (int a = lsn.ashift + 1; a <= lsn.ashift + lsn.nacts; a++) {
                if (a >= lsn.actposttype.getNumElements()
                        || (int) lsn.actposttype.get(a) != ActivityPrecedenceType.ID_POST_AND) {
                    continue;
                }
                List<Integer> frontier = new ArrayList<Integer>();
                frontier.add(a);
                while (!frontier.isEmpty()) {
                    int cur = frontier.remove(0);
                    if (branchActs.contains(cur)) {
                        continue;
                    }
                    branchActs.add(cur);
                    if (isAndJoinPre(cur)) {
                        continue;   // branch tail: do not traverse past the join
                    }
                    for (int s = lsn.ashift + 1; s <= lsn.ashift + lsn.nacts; s++) {
                        if (lsn.graph.get(cur, s) != 0
                                && (int) lsn.parent.get(s) == (int) lsn.parent.get(cur)) {
                            frontier.add(s);
                        }
                    }
                }
            }
        }
        if (branchActs.isEmpty()) {
            return;
        }
        List<Integer> front = new ArrayList<Integer>();
        for (int a : branchActs) {
            List<Integer> calls = lsn.callsof.get(a);
            if (calls != null) {
                for (int c : calls) {
                    front.add((int) lsn.parent.get((int) lsn.callpair.get(c, 2)));
                }
            }
        }
        boolean[] shadow = new boolean[fcrTask.length];
        while (!front.isEmpty()) {
            int t = front.remove(0);
            if (t < 0 || t >= shadow.length || shadow[t]) {
                continue;
            }
            shadow[t] = true;
            for (int a = lsn.ashift + 1; a <= lsn.ashift + lsn.nacts; a++) {
                if ((int) lsn.parent.get(a) != t) {
                    continue;
                }
                List<Integer> calls = lsn.callsof.get(a);
                if (calls != null) {
                    for (int c : calls) {
                        front.add((int) lsn.parent.get((int) lsn.callpair.get(c, 2)));
                    }
                }
            }
        }
        for (int t = 0; t < fcrTask.length; t++) {
            if (shadow[t] && fcrTask[t]) {
                fcrTask[t] = false;
                line_warning(mfilename(new Object() {}), "Multiplicity of task "
                        + lsn.names.get(t) + " is not enforced: it is called from inside an "
                        + "AND-fork branch, whose flows the fork-join transformation retags "
                        + "outside the admission constraint.");
            }
        }
    }

    /** True if any activity of the task is the head of an AND-fork branch (POST_AND). */
    private boolean taskHasAndFork(int tidx) {
        if (lsn.actposttype == null || lsn.actposttype.isEmpty()) {
            return false;
        }
        for (int a = lsn.ashift + 1; a <= lsn.ashift + lsn.nacts; a++) {
            if ((int) lsn.parent.get(a) == tidx && a < lsn.actposttype.getNumElements()
                    && (int) lsn.actposttype.get(a) == ActivityPrecedenceType.ID_POST_AND) {
                return true;
            }
        }
        return false;
    }
}
