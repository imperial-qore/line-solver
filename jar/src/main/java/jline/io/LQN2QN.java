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
 * <li>Replication is represented in one of two ways, selected by the
 *       {@code replication} argument. Under materialisation each replica of a
 *       processor is a station of its own and each replica of a task carries
 *       its own copy of the expanded step graph, its own reply signals and its
 *       own admission row; a call from replica i of the caller reaches the
 *       fan-out block {(i*f+k) mod r} of the callee replicas and splits its
 *       call mean uniformly over them, which is deterministic pairing at f=1
 *       and a uniform broadcast at f=r. Under pooling the replicas of a
 *       processor collapse into one station of r times the servers, a
 *       replicated thread pool into one admission row of r times the bound,
 *       and a replicated reference task into one class of r times the
 *       population; that is exact at an infinite-server host and optimistic
 *       elsewhere, since pooled servers share one queue while the replicas
 *       hold r separate ones.</li>
 * </ul>
 *
 * <p>A CacheTask with delayed-hit retrieval gets a retrieval system, which is an
 * ordinary queueing network: one PS fetch station per cache replica, entered and
 * left by the read class, with the Cache node coalescing concurrent misses of the
 * same item. The fetch is what the miss branch does, so the miss activity's host
 * demand moves onto that station; calls issued by the miss activity stay outside
 * the retrieval system and are warned.</p>
 *
 * <p>A SetupTask carries its setup and delay-off times onto its host station as
 * the Queue setup/delay-off pair, per step class: the server shuts down after the
 * delay-off idle period and pays the setup on the next arrival. An infinite-server
 * processor never shuts down, so the pair is dropped there with a warning, as is a
 * setup with no delay-off time.</p>
 *
 * <p>Not yet represented: retrieval on a cache read with phase-2 successors and
 * the thread pool of a task with an internal AND-fork. Each is reported through
 * line_warning.</p>
 *
 * @see LayeredNetwork
 * @see Network
 * @see SignalType#REPLY
 */
public class LQN2QN {

    private static final int MAXCALLSTAGES = 20;
    // Above this many replica subgraph instantiations "auto" pools instead of
    // materialising: the routing matrix is dense in (classes x nodes), so the
    // conversion cost grows quadratically in the instantiation count.
    private static final int MAXREPLINSTANCES = 128;

    /**
     * Converts a LayeredNetwork to an equivalent queueing network using REPLY signals.
     *
     * @param lqn the LayeredNetwork model to convert
     * @return a Network that models the LQN behaviour with REPLY signal blocking
     */
    public static Network convert(LayeredNetwork lqn) {
        return convert(lqn, "auto");
    }

    /**
     * Converts a LayeredNetwork to an equivalent queueing network using REPLY signals.
     *
     * @param lqn         the LayeredNetwork model to convert
     * @param replication how task and processor replication is represented:
     *                    "auto" materialises the replicas while the expansion
     *                    stays within the instantiation budget and pools them
     *                    otherwise, "materialize" always materialises, "pool"
     *                    always pools
     * @return a Network that models the LQN behaviour with REPLY signal blocking
     */
    public static Network convert(LayeredNetwork lqn, String replication) {
        return new LQN2QN(lqn, replication).build();
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
        // Retrieval system of a delayed-hit CacheTask, null when it has none.
        final Queue fetch;
        final Distribution fetchSvc;

        CacheWiring(Cache node, int readStep, int hitStep, int missStep,
                    DiscreteDistribution itemproc, int nitems,
                    Queue fetch, Distribution fetchSvc) {
            this.node = node;
            this.readStep = readStep;
            this.hitStep = hitStep;
            this.missStep = missStep;
            this.itemproc = itemproc;
            this.nitems = nitems;
            this.fetch = fetch;
            this.fetchSvc = fetchSvc;
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
    private final Map<Integer, Queue> fetchNodeOf = new HashMap<Integer, Queue>();
    private final List<CacheWiring> cacheWiring = new ArrayList<CacheWiring>();

    // joinQuorum rows [joinStep, joinAidx]: the quorum is applied once the
    // class that entered the fork exists.
    private final List<int[]> joinQuorum = new ArrayList<int[]>();

    // Shared INF station carrying the activity think times: the task keeps its
    // thread across a think time but its host processor is released.
    private Delay actThinkNode;

    // Replication: an element is r identical copies of itself, either
    // materialised one station and one step-graph copy per replica, or pooled.
    private final String replication;
    private int[] replRaw = new int[0];
    private boolean materialize = false;
    // Stride of the (element, replica) composite keys used for maps and stacks.
    private int rkey = 1;
    // Class and node names must be unique, so a repeated name is disambiguated.
    private final Map<String, Integer> usedNames = new HashMap<String, Integer>();

    private List<JobClass> stepClass;
    private List<Signal> stepSignal;

    private LQN2QN(LayeredNetwork lqn, String replication) {
        this.lsn = lqn.getStruct();
        this.model = new Network(lqn.getName() + "-QN");
        if (!"auto".equals(replication) && !"materialize".equals(replication)
                && !"pool".equals(replication)) {
            line_error(mfilename(new Object() {}),
                    "replication must be \"auto\", \"materialize\" or \"pool\".");
        }
        this.replication = replication;
    }

    /** Composite (element, 0-based replica) key for maps and stacks. */
    private int ekey(int idx, int rep) {
        return idx * rkey + rep;
    }

    /** Element index of a composite key. */
    private int keyIdx(int k) {
        return k / rkey;
    }

    /** Replica of a composite key. */
    private int keyRep(int k) {
        return k % rkey;
    }

    /** Replicas materialised for a host or task index. */
    private int nrep(int idx) {
        return materialize ? replRaw[idx] : 1;
    }

    /** Capacity multiplier carried by a pooled element, 1 when materialised. */
    private int poolFactor(int idx) {
        return materialize ? 1 : replRaw[idx];
    }

    /**
     * Callee replicas reached by one caller replica. An unset fan-out is the
     * smallest value consistent with repl(a)*fanout = repl(b)*fanin.
     */
    private int fanOutOf(int aTidx, int bTidx, int rb) {
        if (rb <= 1) {
            return 1;
        }
        int f = 0;
        if (lsn.fanout != null && aTidx < lsn.fanout.getNumRows() && bTidx < lsn.fanout.getNumCols()) {
            f = (int) lsn.fanout.get(aTidx, bTidx);
        }
        if (f <= 0) {
            int ra = replRaw[aTidx];
            f = rb > ra ? Math.max(1, rb / ra) : 1;
        }
        return Math.min(Math.max(1, f), rb);
    }

    /** Replicas of the callee reached by replica aRep of the caller. */
    private int[] targetReplicas(int aTidx, int aRep, int bTidx) {
        int rb = nrep(bTidx);
        if (rb <= 1) {
            return new int[]{0};
        }
        int f = fanOutOf(aTidx, bTidx, rb);
        int[] reps = new int[f];
        for (int k = 0; k < f; k++) {
            reps[k] = (aRep * f + k) % rb;
        }
        return reps;
    }

    /** Station key of the processor replica running replica trep of a task. */
    private int hostKey(int tidx, int trep) {
        int hidx = (int) lsn.parent.get(tidx);
        return ekey(hidx, trep % nrep(hidx));
    }

    /**
     * Replica 1 keeps the plain name. The suffix stays inside [A-Za-z0-9_]
     * because a class name is a JSON object key, and MATLAB jsondecode mangles
     * any key that is not a valid identifier.
     */
    private String suffixed(String name, int rep) {
        return rep == 0 ? name : name + "_r" + (rep + 1);
    }

    /**
     * A subgraph copied per call site or per replica repeats its names, which
     * Network rejects, so a repeat is disambiguated by occurrence. The suffix is
     * identifier-safe for the same reason as in suffixed().
     */
    private String uniqueName(String base) {
        Integer prev = usedNames.get(base);
        int n = prev == null ? 1 : prev + 1;
        usedNames.put(base, n);
        return n == 1 ? base : base + "_d" + n;
    }

    /** Copies of the task step graphs a materialised expansion would create. */
    private int replInstantiations(List<Integer> refTaskIndices, List<Integer> openEntries) {
        Set<Integer> seeds = new HashSet<Integer>(refTaskIndices);
        for (int e : openEntries) {
            seeds.add((int) lsn.parent.get(e));
        }
        Map<Integer, Integer> memo = new HashMap<Integer, Integer>();
        int total = 0;
        for (int t : seeds) {
            total += replRaw[t] * taskCost(t, new HashSet<Integer>(), memo);
        }
        return total;
    }

    private int taskCost(int tidx, Set<Integer> stack, Map<Integer, Integer> memo) {
        if (stack.contains(tidx)) {
            return 1;   // recursive cycle: truncated anyway
        }
        if (memo.containsKey(tidx)) {
            return memo.get(tidx);
        }
        Set<Integer> deeper = new HashSet<Integer>(stack);
        deeper.add(tidx);
        int total = 1;
        List<Integer> entries = lsn.entriesof.get(tidx);
        if (entries != null) {
            for (int eidx : entries) {
                List<Integer> targets = new ArrayList<Integer>();
                List<Integer> acts = lsn.actsof.get(eidx);
                if (acts != null) {
                    for (int a : acts) {
                        List<Integer> calls = lsn.callsof.get(a);
                        if (calls == null) {
                            continue;
                        }
                        for (int cidx : calls) {
                            CallType ct = lsn.calltype.get(cidx);
                            if (ct == CallType.SYNC || ct == CallType.ASYNC) {
                                targets.add((int) lsn.callpair.get(cidx, 1));
                            }
                        }
                    }
                }
                // A forwarded entry is expanded per replica just as a call is.
                for (double[] f : forwardingOf(eidx)) {
                    targets.add((int) f[0]);
                }
                for (int te : targets) {
                    int b = (int) lsn.parent.get(te);
                    total += fanOutOf(tidx, b, replRaw[b]) * taskCost(b, deeper, memo);
                }
            }
        }
        memo.put(tidx, total);
        return total;
    }

    private Network build() {
        List<Integer> refTaskIndices = new ArrayList<Integer>();
        for (int t = 0; t < lsn.ntasks; t++) {
            int tidx = lsn.tshift + t;
            if (lsn.isref.get(tidx) > 0) {
                refTaskIndices.add(tidx);
            }
        }
        // open-arrival entries: see _kb/06-solver-catalog.md ("LQN2QN: LQN activity graph -> QN step graph")
        List<Integer> openEntries = new ArrayList<Integer>();
        if (lsn.arrival != null) {
            for (int e = lsn.eshift; e < lsn.eshift + lsn.nentries; e++) {
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
        resolveReplication(refTaskIndices, openEntries);
        computeThreadPoolTasks();

        // Stations: one per host processor replica.
        for (int h = 0; h < lsn.nhosts; h++) {
            double nservers = lsn.mult.get(h);
            SchedStrategy sched = lsn.sched.get(h);
            for (int m = 0; m < nrep(h); m++) {
                if (Double.isInfinite(nservers) || nservers >= Integer.MAX_VALUE || sched == SchedStrategy.INF) {
                    hostStation.put(ekey(h, m), new Delay(model, suffixed(lsn.names.get(h), m)));
                    hostIsDelay.put(ekey(h, m), Boolean.TRUE);
                } else {
                    Queue q = new Queue(model, suffixed(lsn.names.get(h), m), sched);
                    // A pooled processor carries the servers of all its replicas.
                    q.setNumberOfServers((int) nservers * poolFactor(h));
                    hostStation.put(ekey(h, m), q);
                    hostIsDelay.put(ekey(h, m), Boolean.FALSE);
                }
            }
        }

        // Think delays, one per reference task replica.
        for (int refTidx : refTaskIndices) {
            for (int m = 0; m < nrep(refTidx); m++) {
                thinkNode.put(ekey(refTidx, m),
                        new Delay(model, suffixed(lsn.names.get(refTidx) + "_Think", m)));
            }
        }

        // Pass 1: expand the activity graph into a step graph. One closed chain
        // per reference task replica: the replicas are separate populations
        // that meet only where they share a station.
        for (int refTidx : refTaskIndices) {
            for (int rep = 0; rep < nrep(refTidx); rep++) {
                int refKey = ekey(refTidx, rep);
                int thinkStep = addStep(0, 0, null,
                        uniqueName(suffixed(lsn.names.get(refTidx) + "_Think", rep)),
                        false, true, refKey);
                List<Integer> entries = lsn.entriesof.get(refTidx);
                if (entries == null) {
                    continue;
                }
                for (int eidx : entries) {
                    EntryResult er = expandEntry(eidx, refKey, rep);
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
            // Each replica of the entry's task receives its own arrival stream.
            for (int rep = 0; rep < nrep((int) lsn.parent.get(eidx)); rep++) {
                EntryResult er = expandEntry(eidx, 0, rep);
                if (er.firstStep == null) {
                    line_warning(mfilename(new Object() {}), "Open arrival entry "
                            + lsn.names.get(eidx) + " has no bound activity; ignored.");
                    continue;
                }
                List<Port> exits = new ArrayList<Port>(er.replyExits);
                exits.addAll(er.terminals);
                openWiring.add(new Object[]{eidx, er.firstStep, exits});
            }
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
            int refKey = stepRefTask.get(i);
            if (refKey == 0) {
                // A step of an open arrival chain travels in an open class.
                stepClass.set(i, new OpenClass(model, stepName.get(i), 0));
            } else {
                // A pooled reference task holds the population of all its replicas.
                int refTidx = keyIdx(refKey);
                int population = stepIsThink.get(i)
                        ? (int) lsn.mult.get(refTidx) * poolFactor(refTidx) : 0;
                stepClass.set(i, new ClosedClass(model, stepName.get(i), population, thinkNode.get(refKey)));
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

        // phase-2 method name destructor: see _kb/06-solver-catalog.md ("LQN2QN")
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
            ClosedSignal sig = new ClosedSignal(model,
                    suffixed("Ph2End_" + lsn.names.get(keyIdx(rft)), keyRep(rft)),
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
                    if (stepRefTask.get(i) == 0) {   // composite key 0 is an open chain
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
                int refTidx = keyIdx(stepRefTask.get(i));
                Distribution thinkDist = lsn.think == null ? null : lsn.think.get(refTidx);
                Delay tnode = thinkNode.get(stepRefTask.get(i));
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
        // A SetupTask is a server that shuts down when idle and pays a setup on
        // the next arrival, which is the Queue setup/delay-off pair at its host.
        Set<Integer> warnedSetup = new HashSet<Integer>();
        for (int i = 0; i < nsteps; i++) {
            if (stepNode.get(i) != null || stepIsThink.get(i) || stepAidx.get(i) == 0) {
                continue;
            }
            int tidx = (int) lsn.parent.get(stepAidx.get(i));
            Distribution[] times = functionTimesOf(tidx);
            if (times == null) {
                continue;
            }
            if (times == NO_DELAYOFF) {
                if (warnedSetup.add(tidx)) {
                    line_warning(mfilename(new Object() {}), "Setup of setup task "
                            + lsn.names.get(tidx) + " is not represented: it has no delay-off "
                            + "time, so its server never shuts down and never sets up again.");
                }
                continue;
            }
            // Delay subclasses Queue, so the infinite server is tested by hostIsDelay.
            if (Boolean.TRUE.equals(hostIsDelay.get(stepHost.get(i)))) {
                if (warnedSetup.add(tidx)) {
                    line_warning(mfilename(new Object() {}), "Setup of setup task "
                            + lsn.names.get(tidx) + " is not represented: its processor is an "
                            + "infinite server, which never shuts down.");
                }
                continue;
            }
            ((Queue) hostStation.get(stepHost.get(i))).setDelayOff(stepClass.get(i),
                    times[0], times[1]);
        }

        // reply signal service declaration: see _kb/06-solver-catalog.md ("LQN2QN")
        for (int i = 0; i < nsteps; i++) {
            if (stepSignal.get(i) == null) {
                continue;
            }
            for (Station hstation : hostStation.values()) {
                hstation.setService(stepSignal.get(i), new Immediate());
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
            if (cw.fetch != null) {
                // Service and routing of the retrieval system are read off the read class.
                cw.fetch.setService(stepClass.get(cw.readStep),
                        cw.fetchSvc == null ? new Immediate() : cw.fetchSvc);
                cw.node.setRetrievalSystem(stepClass.get(cw.readStep),
                        stepClass.get(cw.missStep), new Queue[]{cw.fetch});
            }
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
        // Retrieval systems: the read class circulates cache -> fetch -> cache.
        for (CacheWiring cw : cacheWiring) {
            if (cw.fetch != null) {
                JobClass rcls = stepClass.get(cw.readStep);
                P.set(rcls, rcls, cw.node, cw.fetch, 1.0);
                P.set(rcls, rcls, cw.fetch, cw.node, 1.0);
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
        // One admission row per thread-pool task replica; a pooled task keeps a
        // single row whose bound covers all its replicas.
        java.util.TreeSet<Integer> present = new java.util.TreeSet<Integer>();
        for (int[] held : stepTasks) {
            for (int hk : held) {
                present.add(hk);
            }
        }
        List<Integer> fcrList = new ArrayList<Integer>(present);
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
            int tkey = fcrList.get(ti);
            for (int i = 0; i < stepTasks.size(); i++) {
                boolean holds = false;
                for (int ht : stepTasks.get(i)) {
                    if (ht == tkey) {
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
                for (CacheWiring cw : cacheWiring) {
                    // The caller holds its thread for the whole fetch.
                    if (cw.readStep == i && cw.fetch != null && !regionNodes.contains(cw.fetch)) {
                        regionNodes.add(cw.fetch);
                    }
                }
            }
            b.set(ti, 0, lsn.mult.get(keyIdx(tkey)) * poolFactor(keyIdx(tkey)));
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
                        boolean blocks, boolean isthink, int refKey) {
        stepAidx.add(aidx);
        stepHost.add(hidx);
        stepSvc.add(svc);
        stepName.add(name);
        stepBlocks.add(blocks);
        stepIsThink.add(isthink);
        stepRefTask.add(refKey);
        stepNode.add(null);
        int id = stepAidx.size() - 1;
        stepClassOwner.add(id);
        // Thread holders are (task, replica) keys: a replica has its own pool.
        java.util.TreeSet<Integer> held = new java.util.TreeSet<Integer>();
        for (int tk : threadStack) {
            int ti = keyIdx(tk);
            if (ti >= 0 && ti < fcrTask.length && fcrTask[ti]) {
                held.add(tk);
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
    private int addAuxStep(Node nodeObj, int ownerStep, String name, int refKey) {
        int id = addStep(0, 0, null, name, false, false, refKey);
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

    /**
     * Expands the activity subgraph bound to an entry of replica trep of its
     * task, in the current call context.
     */
    private EntryResult expandEntry(int eidx, int refKey, int trep) {
        EntryResult res = new EntryResult();

        if (entryStack.contains(Integer.valueOf(ekey(eidx, trep)))) {
            line_warning(mfilename(new Object() {}),
                    "Recursive call cycle at entry " + lsn.names.get(eidx) + " truncated.");
            return res;
        }
        entryStack.add(Integer.valueOf(ekey(eidx, trep)));
        threadStack.add(Integer.valueOf(ekey((int) lsn.parent.get(eidx), trep)));
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

            EntryWalker walker = new EntryWalker(eidx, refKey, localActs, trep);
            walker.run(bound, res);

            // forwarding rules: see _kb/06-solver-catalog.md ("LQN2QN")
            List<double[]> fwd = forwardingOf(eidx);
            if (!fwd.isEmpty() && !res.replyExits.isEmpty()) {
                List<Port> ownPorts = new ArrayList<Port>(res.replyExits);
                List<Port> fwdExits = new ArrayList<Port>();
                double pforw = 0.0;
                // forwarder thread release: see _kb/06-solver-catalog.md ("LQN2QN")
                Integer fwdThread = threadStack.remove(threadStack.size() - 1);
                int fwdTidx = (int) lsn.parent.get(eidx);
                for (double[] f : fwd) {
                    // A forwarding call spreads over the reached callee replicas
                    // exactly as a synchronous one does.
                    int[] reps = targetReplicas(fwdTidx, trep, (int) lsn.parent.get((int) f[0]));
                    double p = f[1];
                    int reached = 0;
                    for (int m : reps) {
                        EntryResult fr = expandEntry((int) f[0], refKey, m);
                        if (fr.firstStep == null) {
                            continue;
                        }
                        reached++;
                        for (Port op : ownPorts) {
                            addRoute(op, fr.firstStep, op.prob * p / reps.length);
                        }
                        fwdExits.addAll(fr.replyExits);
                        fwdExits.addAll(fr.terminals);
                    }
                    pforw += p * reached / reps.length;
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
        for (int cidx = 0; cidx < lsn.ncalls; cidx++) {
            CallType ct = lsn.calltype.get(cidx);
            if (ct != CallType.FWD || (int) lsn.callpair.get(cidx, 0) != eidx) {
                continue;
            }
            double p = callMean(cidx);
            p = Math.min(Math.max(p, 0.0), 1.0);
            if (p > GlobalConstants.FineTol) {
                fwd.add(new double[]{lsn.callpair.get(cidx, 1), p});
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
        private final int refKey;
        private final int trep;
        private final List<Integer> localActs;
        private final Map<Integer, Port> visitedExit = new HashMap<Integer, Port>();
        private final Map<Integer, Integer> visitedEntry = new HashMap<Integer, Integer>();
        private final Map<Integer, Integer> joinOf = new HashMap<Integer, Integer>();
        private final List<Integer> forkOwnerStack = new ArrayList<Integer>();
        private boolean sawReply = false;
        private EntryResult res;

        EntryWalker(int eidx, int refKey, List<Integer> localActs, int trep) {
            this.eidx = eidx;
            this.refKey = refKey;
            this.trep = trep;
            this.localActs = localActs;
        }

        /** A class name carries the replica of the task that owns the step. */
        private String sname(String name) {
            return uniqueName(suffixed(name, trep));
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
                cacheNode = getCacheNode(tidx, trep);
            }

            int entryStep = makeActivitySteps(aidx, tidx, refKey, trep, cacheNode);
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
                    int trigH = addStep(aidx, hostKey(tidx, trep), null,
                            sname(lsn.names.get(aidx) + "_ph2h"), false, false, refKey);
                    int trigM = addStep(aidx, hostKey(tidx, trep), null,
                            sname(lsn.names.get(aidx) + "_ph2m"), false, false, refKey);
                    addCacheRoute(entryStep, trigH);
                    addCacheRoute(entryStep, trigM);
                    if (hasRetrieval(tidx)) {
                        line_warning(mfilename(new Object() {}), "Delayed-hit retrieval of "
                                + lsn.names.get(tidx) + " is not represented on a cache read "
                                + "with phase-2 successors.");
                    }
                    cacheWiring.add(new CacheWiring(cacheNode, entryStep, trigH, trigM,
                            lsn.itemproc == null ? null : lsn.itemproc.get(eidx),
                            (int) lsn.nitems.get(eidx), null, null));
                    res.replyExits.add(new Port(trigH, false, 1.0));
                    res.replyExits.add(new Port(trigM, false, 1.0));
                    List<Integer> savedStack = new ArrayList<Integer>(threadStack);
                    threadStack.clear();
                    threadStack.add(ekey(tidx, trep));
                    int nT0 = res.terminals.size();
                    for (int hm = 0; hm < 2; hm++) {
                        int sEntry = walk(succ.get(hm));
                        if (stepNode.get(sEntry) != null) {
                            int hmHead = addStep(aidx, hostKey(tidx, trep), null,
                                    sname(lsn.names.get(aidx) + "_ph2b" + (hm + 1)), false, false, refKey);
                            addRoute(new Port(hmHead, false, 1.0), sEntry, 1.0);
                            sEntry = hmHead;
                        }
                        spawnPairs.add(new int[]{hm == 0 ? trigH : trigM, sEntry});
                    }
                    while (res.terminals.size() > nT0) {
                        Port t = res.terminals.remove(nT0);
                        ph2Exits.add(new double[]{t.step, t.isSignal ? 1.0 : 0.0, t.prob, refKey});
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
                    int trig = addStep(aidx, hostKey(tidx, trep), null,
                            sname(lsn.names.get(aidx) + "_ph2t"), false, false, refKey);
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
                    threadStack.add(ekey(tidx, trep));
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
                        int head = addStep(aidx, hostKey(tidx, trep), null,
                                sname(lsn.names.get(aidx) + "_ph2"), false, false, refKey);
                        wireAndFork(new Port(head, false, 1.0), succ, aidx);
                        target = head;
                    } else if (isAndJoinPre(aidx)) {
                        // phase-2 AND-join branch-tail spawn: see _kb/06-solver-catalog.md ("LQN2QN")
                        int head = addStep(aidx, hostKey(tidx, trep), null,
                                sname(lsn.names.get(aidx) + "_ph2"), false, false, refKey);
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
                        int head = addStep(aidx, hostKey(tidx, trep), null,
                                sname(lsn.names.get(aidx) + "_ph2"), false, false, refKey);
                        for (int s2 : posSucc) {
                            int sEntry = walk(s2);
                            addRoute(new Port(head, false, 1.0), sEntry, lsn.graph.get(aidx, s2));
                        }
                        target = head;
                    }
                    spawnPairs.add(new int[]{exitPort.step, target});
                    while (res.terminals.size() > nT0) {
                        Port t = res.terminals.remove(nT0);
                        ph2Exits.add(new double[]{t.step, t.isSignal ? 1.0 : 0.0, t.prob, refKey});
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
                    Queue fetch = null;
                    Distribution fetchSvc = null;
                    if (hasRetrieval(tidx)) {
                        // The fetch is what the miss branch does, so its demand moves to
                        // the fetch station where concurrent misses coalesce.
                        fetch = getFetchNode(tidx, trep);
                        fetchSvc = stepSvc.get(mEntry);
                        stepSvc.set(mEntry, null);
                        List<Integer> missCalls = lsn.callsof.get(succ.get(1));
                        if (missCalls != null && !missCalls.isEmpty()) {
                            line_warning(mfilename(new Object() {}), "Calls of miss activity "
                                    + lsn.names.get(succ.get(1)) + " stay outside the retrieval "
                                    + "system, so they are not coalesced across concurrent misses.");
                        }
                    }
                    cacheWiring.add(new CacheWiring(cacheNode, entryStep, hEntry, mEntry,
                            lsn.itemproc == null ? null : lsn.itemproc.get(eidx),
                            (int) lsn.nitems.get(eidx), fetch, fetchSvc));
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
            String forkName = sname("Fork_" + lsn.names.get(aidx));
            Fork forkNode = new Fork(model, forkName);
            int forkStep = addAuxStep(forkNode, from.step, forkName, refKey);
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
                String rname = sname("Fork_" + lsn.names.get(aidx) + "_" + (b + 1));
                Router routerNode = new Router(model, rname);
                int routerStep = addAuxStep(routerNode, forkStep, rname, refKey);
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
            String joinName = sname("Join_" + lsn.names.get(joinAidx));
            Join joinNode = new Join(model, joinName, stepNode.get(forkOwner));
            int joinStep = addAuxStep(joinNode, forkOwner, joinName, refKey);
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
            if (a < 0 || a >= lsn.replygraph.getNumRows() || e < 0 || e >= lsn.replygraph.getNumCols()) {
                return false;
            }
            return lsn.replygraph.get(a, e) != 0;
        }

        /**
         * One step for the host demand, plus one step per unrolled synchronous
         * call stage. The exit port of the activity is left in visitedExit.
         */
        private int makeActivitySteps(int aidx, int tidx, int refKey, int trep, Cache cacheNode) {
            int hidx = hostKey(tidx, trep);

            if (cacheNode != null) {
                // cache-read step: see _kb/06-solver-catalog.md ("LQN2QN")
                int readStep = addStep(aidx, hidx, null, sname(lsn.names.get(aidx)), false, false, refKey);
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
            int entryStep = addStep(aidx, hidx, svc, sname(lsn.names.get(aidx)), false, false, refKey);
            // exit port semantics: see _kb/06-solver-catalog.md ("LQN2QN")
            Port cur = new Port(entryStep, false, 1.0);

            // activity think time: see _kb/06-solver-catalog.md ("LQN2QN")
            Distribution think = actThinkOf(aidx);
            if (think != null) {
                int thinkStep = addStep(aidx, hidx, think, sname(lsn.names.get(aidx) + "_think"),
                        false, false, refKey);
                stepNode.set(thinkStep, actThinkStation());
                addRoute(cur, thinkStep, 1.0);
                cur = new Port(thinkStep, false, 1.0);
            }

            // call-blocking eligibility rules: see _kb/06-solver-catalog.md ("LQN2QN")
            boolean hostBlocks = !Boolean.TRUE.equals(hostIsDelay.get(hidx))
                    && !(tidx < fcrTask.length && fcrTask[tidx]) && refKey != 0
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
                // The stage is routed to the callee replicas this caller replica
                // reaches, which share the call mean uniformly.
                List<Integer> calleeFirsts = new ArrayList<Integer>();
                List<Port> calleeReplies = new ArrayList<Port>();
                for (int m : targetReplicas(tidx, trep, (int) lsn.parent.get(stage.targetEidx))) {
                    EntryResult callee = expandEntry(stage.targetEidx, refKey, m);
                    if (callee.firstStep == null) {
                        continue;
                    }
                    calleeFirsts.add(callee.firstStep);
                    // dead-end callee return: see _kb/06-solver-catalog.md ("LQN2QN")
                    calleeReplies.addAll(callee.replyExits);
                    calleeReplies.addAll(callee.terminals);
                }
                if (stage.isasync) {
                    threadStack.add(asyncThread);
                }
                if (calleeFirsts.isEmpty()) {
                    continue;   // callee not expandable: drop the call, never block
                }
                double share = 1.0 / calleeFirsts.size();

                // merge-step necessity rules: see _kb/06-solver-catalog.md ("LQN2QN")
                boolean needsMerge = (stage.prob < 1.0) || (k < callStages.size() - 1)
                        || !blocks || isAndJoinPre(aidx);
                int nxt = -1;
                if (needsMerge) {
                    String retName = sname(lsn.names.get(aidx) + "_c" + (k + 1) + "_ret");
                    nxt = addStep(aidx, hidx, null, retName, false, false, refKey);
                    if (!blocks) {
                        // non-blocking merge on a Router: see _kb/06-solver-catalog.md ("LQN2QN")
                        stepNode.set(nxt, new Router(model, retName));
                    }
                }

                if (blocks) {
                    // first mandatory call bound to service class: see _kb/06-solver-catalog.md ("LQN2QN")
                    int blk;
                    if (k == 0 && stage.prob >= 1.0 && cur.step == entryStep && !cur.isSignal) {
                        blk = entryStep;
                        stepBlocks.set(blk, Boolean.TRUE);
                    } else {
                        blk = addStep(aidx, hidx, null, sname(lsn.names.get(aidx) + "_c" + (k + 1)),
                                true, false, refKey);
                        addRoute(cur, blk, stage.prob);
                        if (stage.prob < 1.0) {
                            addRoute(cur, nxt, 1.0 - stage.prob);
                        }
                    }
                    // Every reached replica replies into the same signal, so the
                    // call site blocks once however many replicas it has.
                    for (int cf : calleeFirsts) {
                        addRoute(new Port(blk, false, 1.0), cf, share);
                    }
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
                    for (int cf : calleeFirsts) {
                        addRoute(cur, cf, stage.prob * share);
                    }
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

    /**
     * Setup and delay-off of a SetupTask, or null when it is always on. Both are
     * needed: without a delay-off the server never shuts down, so it pays the setup
     * once at most and the pair carries no information. A setup without a delay-off
     * returns the NO_DELAYOFF marker so the caller can warn.
     */
    private static final Distribution[] NO_DELAYOFF = new Distribution[0];

    private Distribution[] functionTimesOf(int tidx) {
        if (lsn.hassetup == null || lsn.hassetup.isEmpty()
                || tidx >= lsn.hassetup.getNumElements() || lsn.hassetup.get(tidx) == 0) {
            return null;
        }
        Distribution setup = lsn.setuptime == null ? null : lsn.setuptime.get(tidx);
        if (!isNonTrivial(setup)) {
            return null;
        }
        Distribution delayoff = lsn.delayofftime == null ? null : lsn.delayofftime.get(tidx);
        if (delayoff == null) {
            return NO_DELAYOFF;
        }
        return new Distribution[]{setup, delayoff};
    }

    /** True when the CacheTask coalesces concurrent misses of the same item. */
    private boolean hasRetrieval(int tidx) {
        return lsn.hasretrieval != null && !lsn.hasretrieval.isEmpty()
                && tidx < lsn.hasretrieval.getNumElements() && lsn.hasretrieval.get(tidx) != 0;
    }

    /**
     * The retrieval system of a CacheTask is an ordinary queueing network: one PS
     * fetch station per cache replica, as in the LN cache sublayer.
     */
    private Queue getFetchNode(int tidx, int trep) {
        Queue q = fetchNodeOf.get(ekey(tidx, trep));
        if (q != null) {
            return q;
        }
        q = new Queue(model, suffixed(lsn.names.get(tidx) + "_Cache_Fetch", trep), SchedStrategy.PS);
        fetchNodeOf.put(ekey(tidx, trep), q);
        return q;
    }

    /** One Cache node per CacheTask replica, named so as not to collide with the processor station. */
    private Cache getCacheNode(int tidx, int trep) {
        Cache cnode = cacheNodeOf.get(ekey(tidx, trep));
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
        cnode = new Cache(model, suffixed(lsn.names.get(tidx) + "_Cache", trep),
                (int) lsn.nitems.get(tidx), capMatrix, replacementOf(tidx));
        cacheNodeOf.put(ekey(tidx, trep), cnode);
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
        for (int pred = 0; pred < lsn.graph.getNumRows(); pred++) {
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
            Distribution d = lsn.callproc.get(cidx);
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
            int targetEidx = (int) lsn.callpair.get(cidx, 1);
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

    /**
     * Reads the replication factors and decides whether the replicas are
     * materialised or pooled, warning in the latter case since pooling is an
     * approximation at a queueing host.
     */
    private void resolveReplication(List<Integer> refTaskIndices, List<Integer> openEntries) {
        replRaw = new int[lsn.nhosts + lsn.ntasks + 1];
        java.util.Arrays.fill(replRaw, 1);
        if (lsn.repl != null && !lsn.repl.isEmpty()) {
            for (int idx = 0; idx < lsn.nhosts + lsn.ntasks; idx++) {
                if (idx < lsn.repl.getNumCols()) {
                    replRaw[idx] = Math.max(1, (int) Math.round(lsn.repl.get(0, idx)));
                }
            }
        }
        int firstReplicated = -1;
        rkey = 1;
        for (int idx = 0; idx < lsn.nhosts + lsn.ntasks; idx++) {
            if (replRaw[idx] > 1 && firstReplicated < 0) {
                firstReplicated = idx;
            }
            rkey = Math.max(rkey, replRaw[idx]);
        }
        if (firstReplicated < 0) {
            materialize = false;
            return;
        }
        materialize = "materialize".equals(replication)
                || ("auto".equals(replication)
                    && replInstantiations(refTaskIndices, openEntries) <= MAXREPLINSTANCES);
        if (!materialize) {
            line_warning(mfilename(new Object() {}), "Replication of "
                    + lsn.names.get(firstReplicated) + " is pooled: its replicas become one "
                    + "station of r times the servers, one admission row of r times the bound "
                    + "and one reference class of r times the population, which is exact at an "
                    + "infinite-server host and optimistic elsewhere. Pass \"materialize\" for "
                    + "one station and one step-graph copy per replica.");
        }
    }

    private void warnUnsupported() {
        if (lsn.calltype != null) {
            for (int cidx = 0; cidx < lsn.ncalls; cidx++) {
                if (lsn.calltype.get(cidx) == CallType.ASYNC) {
                    line_warning(mfilename(new Object() {}), "Asynchronous calls are represented by "
                            + "LQN2QN as non-blocking visits: the caller releases its server but "
                            + "remains serialised behind the callee.");
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
        fcrTask = new boolean[lsn.tshift + lsn.ntasks];
        for (int t = 0; t < lsn.ntasks; t++) {
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
            for (int a = lsn.ashift; a < lsn.ashift + lsn.nacts; a++) {
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
                    for (int s = lsn.ashift; s < lsn.ashift + lsn.nacts; s++) {
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
                    front.add((int) lsn.parent.get((int) lsn.callpair.get(c, 1)));
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
            for (int a = lsn.ashift; a < lsn.ashift + lsn.nacts; a++) {
                if ((int) lsn.parent.get(a) != t) {
                    continue;
                }
                List<Integer> calls = lsn.callsof.get(a);
                if (calls != null) {
                    for (int c : calls) {
                        front.add((int) lsn.parent.get((int) lsn.callpair.get(c, 1)));
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
        for (int a = lsn.ashift; a < lsn.ashift + lsn.nacts; a++) {
            if ((int) lsn.parent.get(a) == tidx && a < lsn.actposttype.getNumElements()
                    && (int) lsn.actposttype.get(a) == ActivityPrecedenceType.ID_POST_AND) {
                return true;
            }
        }
        return false;
    }
}
