/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SSJ-based discrete event simulation engine for LayeredNetwork (LQN) models.
 *
 * Mechanically translated from Solver_ssj_ln.kt.
 */
package jline.solvers.ldes.handlers;

import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.CallType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.processes.Distribution;
import jline.lang.processes.DiscreteDistribution;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LNLDESResult;
import jline.util.matrix.Matrix;
import jline.util.Pair;
import umontreal.ssj.simevents.Event;
import umontreal.ssj.simevents.Simulator;
import umontreal.ssj.stat.Tally;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;

public final class Solver_ssj_ln {

    private Solver_ssj_ln() {}

    /**
     * Steady-state LayeredNetwork simulation using SSJ library.
     */
    public static LNLDESResult solver_ssj_ln(LayeredNetworkStruct lsn, SolverOptions options) {
        LNSSJSimulator simulator = new LNSSJSimulator(lsn, options);
        double maxEvents = (double) options.samples;
        simulator.simulate(maxEvents);
        return simulator.getLNLDESResult();
    }

    /**
     * Core discrete event simulation engine for LayeredNetwork models using SSJ library.
     */
    static class LNSSJSimulator {

        // see _kb/09-ldes-and-cache.md (Solver_ssj engine notes: SSJ engine infrastructure)
        private final Simulator ssjSim = new Simulator();

        /**
         * Base class of every event scheduled by this simulator, binding the event
         * to the enclosing run's {@link #ssjSim} rather than to SSJ's JVM-wide
         * default simulator. All event classes below must extend this, not Event.
         */
        private abstract class SimEvent extends Event {
            SimEvent() {
                super(LNSSJSimulator.this.ssjSim);
            }
        }

        private final LayeredNetworkStruct lsn;
        private final SolverOptions options;

        // LQN dimensions
        private final int nhosts;
        private final int ntasks;
        private final int nentries;
        private final int nacts;
        private final int ncalls;

        // Index shifts
        private final int hshift;
        private final int tshift;
        private final int eshift;
        private final int ashift;
        private final int cshift;

        // Random
        private final long seed;
        private final Random random;

        // Host state
        private final int[] hostBusy;
        private final int[] hostServers;
        private SchedStrategy[] hostSchedStrategy;
        private PriorityQueue<HostRequest>[] hostQueues;

        // Task state
        private final int[] taskBusy;
        private final int[] taskMultiplicity;
        private SchedStrategy[] taskSchedStrategy;
        private PriorityQueue<LNRequest>[] taskQueues;

        // PS state
        private List<PSJob>[] hostPSJobs;
        private PSUpdateEvent[] hostPSEvent;

        private final double psTimeSlice = 0.001;

        // Request tracking
        private long nextRequestId = 0L;
        private final Map<Long, LNRequest> activeRequests = new HashMap<Long, LNRequest>();

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Integer, Integer> joinRequiredCount = new HashMap<Integer, Integer>();
        private final Map<Long, JoinProgress> joinProgress = new HashMap<Long, JoinProgress>();
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Long, Integer> pendingWork = new HashMap<Long, Integer>();

        // Sync call tracking
        private final Map<Long, LinkedList<Integer>> pendingSyncCallQueue = new HashMap<Long, LinkedList<Integer>>();
        private final Map<Long, Integer> syncCallActivity = new HashMap<Long, Integer>();
        private final Map<Long, Boolean> syncCallInFlight = new HashMap<Long, Boolean>();

        // Statistics per host
        private final double[] hostBusyTime;
        private final int[] hostCompletions;
        private final double[] lastHostUpdateTime;

        // Statistics per task
        private final double[] taskUtilTime;
        private final double[] taskCompletions;
        private final double[] lastTaskUpdateTime;

        // Statistics per entry
        private final double[] entryThroughput;
        private final Tally[] entryResponseTimeTally;
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final double[] entryResidenceTime;

        // Statistics per activity
        private final double[] activityThroughput;
        private final Tally[] activityServiceTimeTally;
        private final double[] activityBusyTime;
        private final double[] lastActivityUpdateTime;
        private final int[] activityBusy;

        private final Map<Pair<Long, Integer>, Double> activityStartTimes = new HashMap<Pair<Long, Integer>, Double>();

        private final int[] taskQueueLength;
        private final double[] taskQueueLengthTime;
        private final double[] lastTaskQueueUpdateTime;

        private final int[] hostQueueLength;
        private final double[] hostQueueLengthTime;
        private final double[] lastHostQueueUpdateTime;

        private long totalEventCount = 0L;
        private long maxEvents = 0L;
        private long warmupEventThreshold = 0L;
        private boolean warmupDone = false;
        private double warmupEndTime = 0.0;

        private final Map<Integer, Integer> entryBoundActivity = new HashMap<Integer, Integer>();
        private final Map<Integer, Integer> entryReplyActivity = new HashMap<Integer, Integer>();
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Integer, List<Pair<Integer, Double>>> entryForwards = new HashMap<Integer, List<Pair<Integer, Double>>>();
        private final Map<Integer, Integer> activityToEntry = new HashMap<Integer, Integer>();
        private final Map<Integer, List<Integer>> activitySuccessors = new HashMap<Integer, List<Integer>>();
        private final Map<Integer, List<Pair<Integer, Double>>> orForkProbs = new HashMap<Integer, List<Pair<Integer, Double>>>();

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private LNCacheState[] cacheStateByTask;
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Integer, CacheAccess> cacheAccessByDriver = new HashMap<Integer, CacheAccess>();
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Integer, Integer> cacheBranchReply = new HashMap<Integer, Integer>();
        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private final Map<Long, FetchKey> fetcherItem = new HashMap<Long, FetchKey>();
        // Sentinel returned by accessCache when a read is parked as a delayed hit.
        private static final int LN_CACHE_HELD = Integer.MIN_VALUE;

        // Comparators
        private final Comparator<LNRequest> fcfsComparator = new Comparator<LNRequest>() {
            @Override
            public int compare(LNRequest r1, LNRequest r2) {
                return Double.compare(r1.taskArrivalTime, r2.taskArrivalTime);
            }
        };

        private final Comparator<LNRequest> lcfsComparator = new Comparator<LNRequest>() {
            @Override
            public int compare(LNRequest r1, LNRequest r2) {
                return Double.compare(r2.taskArrivalTime, r1.taskArrivalTime);
            }
        };

        private final Comparator<HostRequest> hostFcfsComparator = new Comparator<HostRequest>() {
            @Override
            public int compare(HostRequest h1, HostRequest h2) {
                return Double.compare(h1.queueArrivalTime, h2.queueArrivalTime);
            }
        };

        private final Comparator<HostRequest> hostLcfsComparator = new Comparator<HostRequest>() {
            @Override
            public int compare(HostRequest h1, HostRequest h2) {
                return Double.compare(h2.queueArrivalTime, h1.queueArrivalTime);
            }
        };

        // Debug flag
        private final boolean DEBUG_FORK_JOIN = false;

        public static class PSJob {
            public final HostRequest hostRequest;
            public double remainingService;
            public double lastUpdateTime;

            public PSJob(HostRequest hostRequest, double remainingService, double lastUpdateTime) {
                this.hostRequest = hostRequest;
                this.remainingService = remainingService;
                this.lastUpdateTime = lastUpdateTime;
            }
        }

        public static class LNRequest {
            public final long requestId;
            public final int entryIdx;
            public final CallType callType;
            public final double arrivalTime;
            public double taskArrivalTime;
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            public double serviceStartTime;
            public final long callerRequestId;
            public final int callerTaskIdx;
            public int currentActivityIdx;
            public long parentForkId;
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            public boolean replied;
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            public LNRequest rootRequest;

            public LNRequest(long requestId, int entryIdx, CallType callType, double arrivalTime,
                             double taskArrivalTime, long callerRequestId, int callerTaskIdx,
                             int currentActivityIdx, long parentForkId) {
                this.requestId = requestId;
                this.entryIdx = entryIdx;
                this.callType = callType;
                this.arrivalTime = arrivalTime;
                this.taskArrivalTime = taskArrivalTime;
                this.serviceStartTime = taskArrivalTime;
                this.callerRequestId = callerRequestId;
                this.callerTaskIdx = callerTaskIdx;
                this.currentActivityIdx = currentActivityIdx;
                this.parentForkId = parentForkId;
            }

            public LNRequest getRoot() {
                return rootRequest != null ? rootRequest : this;
            }

            public LNRequest copyWith(long requestId, int currentActivityIdx, long parentForkId) {
                LNRequest r = new LNRequest(requestId, this.entryIdx, this.callType, this.arrivalTime,
                        this.taskArrivalTime, this.callerRequestId, this.callerTaskIdx,
                        currentActivityIdx, parentForkId);
                r.rootRequest = this.getRoot();
                return r;
            }
        }

        public static class HostRequest {
            public final LNRequest request;
            public final int activityIdx;
            public final double serviceTime;
            public final double queueArrivalTime;

            public HostRequest(LNRequest request, int activityIdx, double serviceTime, double queueArrivalTime) {
                this.request = request;
                this.activityIdx = activityIdx;
                this.serviceTime = serviceTime;
                this.queueArrivalTime = queueArrivalTime;
            }
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        public static class JoinProgress {
            public final Set<Integer> arrivedInputs = new HashSet<Integer>();
            public final List<Double> completionTimes = new ArrayList<Double>();
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private static final class LNCacheState {
            final int numItems;
            final int[] levelCapacities;
            final ReplacementStrategy replacementStrategy;
            final LinkedList<Integer>[] levels;
            long totalHits;
            long totalMisses;
            long totalReads;
            // --- delayed-hit retrieval (coalescing) ---
            final boolean hasRetrieval;
            final boolean[] inFlight;                 // [item] -> a fetch for this item is in progress
            final List<HeldLN>[] heldDelayedHits;     // [item] -> reads parked while item is in flight
            long totalDelayedHits;

            @SuppressWarnings("unchecked")
            LNCacheState(int numItems, int[] levelCapacities, ReplacementStrategy rs, boolean hasRetrieval) {
                this.numItems = numItems;
                this.levelCapacities = levelCapacities;
                this.replacementStrategy = rs;
                this.levels = new LinkedList[levelCapacities.length];
                for (int l = 0; l < levelCapacities.length; l++) {
                    this.levels[l] = new LinkedList<Integer>();
                }
                this.hasRetrieval = hasRetrieval;
                if (hasRetrieval) {
                    this.inFlight = new boolean[numItems];
                    this.heldDelayedHits = (List<HeldLN>[]) new List<?>[numItems];
                } else {
                    this.inFlight = null;
                    this.heldDelayedHits = null;
                }
            }
        }

        // A read parked as a delayed hit while its item is being fetched. On release
        // it reads the now-cached item via the hit branch of its own cache access.
        private static final class HeldLN {
            final LNRequest request;
            final CacheAccess access;
            HeldLN(LNRequest request, CacheAccess access) {
                this.request = request;
                this.access = access;
            }
        }

        // Identifies the in-flight fetch triggered by a fetcher request: the cache,
        // the item being fetched, and the cache access that triggered it.
        private static final class FetchKey {
            final LNCacheState cache;
            final int item;
            FetchKey(LNCacheState cache, int item) {
                this.cache = cache;
                this.item = item;
            }
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private static final class CacheAccess {
            final LNCacheState cache;
            final double[] cdf;
            final int hitAidx;
            final int missAidx;
            final int itemEntryIdx;

            CacheAccess(LNCacheState cache, double[] cdf, int hitAidx, int missAidx, int itemEntryIdx) {
                this.cache = cache;
                this.cdf = cdf;
                this.hitAidx = hitAidx;
                this.missAidx = missAidx;
                this.itemEntryIdx = itemEntryIdx;
            }
        }

        @SuppressWarnings("unchecked")
        public LNSSJSimulator(LayeredNetworkStruct lsn, SolverOptions options) {
            this.lsn = lsn;
            this.options = options;
            this.nhosts = lsn.nhosts;
            this.ntasks = lsn.ntasks;
            this.nentries = lsn.nentries;
            this.nacts = lsn.nacts;
            this.ncalls = lsn.ncalls;
            this.hshift = lsn.hshift;
            this.tshift = lsn.tshift;
            this.eshift = lsn.eshift;
            this.ashift = lsn.ashift;
            this.cshift = lsn.cshift;
            this.seed = (long) options.seed;
            this.random = new Random(this.seed);

            this.hostBusy = new int[nhosts];
            this.hostServers = new int[nhosts];

            this.taskBusy = new int[ntasks];
            this.taskMultiplicity = new int[ntasks];

            this.hostBusyTime = new double[nhosts];
            this.hostCompletions = new int[nhosts];
            this.lastHostUpdateTime = new double[nhosts];

            this.taskUtilTime = new double[ntasks];
            this.taskCompletions = new double[ntasks];
            this.lastTaskUpdateTime = new double[ntasks];

            this.entryThroughput = new double[nentries];
            this.entryResidenceTime = new double[nentries];
            this.entryResponseTimeTally = new Tally[nentries];
            for (int i = 0; i < nentries; i++) {
                this.entryResponseTimeTally[i] = new Tally();
            }

            this.activityThroughput = new double[nacts];
            this.activityServiceTimeTally = new Tally[nacts];
            for (int i = 0; i < nacts; i++) {
                this.activityServiceTimeTally[i] = new Tally();
            }
            this.activityBusyTime = new double[nacts];
            this.lastActivityUpdateTime = new double[nacts];
            this.activityBusy = new int[nacts];

            this.taskQueueLength = new int[ntasks];
            this.taskQueueLengthTime = new double[ntasks];
            this.lastTaskQueueUpdateTime = new double[ntasks];

            this.hostQueueLength = new int[nhosts];
            this.hostQueueLengthTime = new double[nhosts];
            this.lastHostQueueUpdateTime = new double[nhosts];

            initialize();
        }

        @SuppressWarnings("unchecked")
        private void initialize() {
            hostSchedStrategy = new SchedStrategy[nhosts];
            for (int i = 0; i < nhosts; i++) hostSchedStrategy[i] = SchedStrategy.FCFS;
            hostPSJobs = new List[nhosts];
            for (int i = 0; i < nhosts; i++) hostPSJobs[i] = new ArrayList<PSJob>();
            hostPSEvent = new PSUpdateEvent[nhosts];

            for (int h = 0; h < nhosts; h++) {
                int hidx = h + 1;
                hostServers[h] = Math.max(1, (int) lsn.mult.get(0, hidx));
                SchedStrategy s = lsn.sched.get(hidx);
                hostSchedStrategy[h] = (s != null) ? s : SchedStrategy.FCFS;
            }

            hostQueues = new PriorityQueue[nhosts];
            for (int h = 0; h < nhosts; h++) {
                hostQueues[h] = new PriorityQueue<HostRequest>(11, getHostComparator(hostSchedStrategy[h]));
            }

            taskSchedStrategy = new SchedStrategy[ntasks];
            for (int i = 0; i < ntasks; i++) taskSchedStrategy[i] = SchedStrategy.FCFS;

            for (int t = 0; t < ntasks; t++) {
                int tidx = tshift + t + 1;
                taskMultiplicity[t] = Math.max(1, (int) lsn.mult.get(0, tidx));
                SchedStrategy sched = lsn.sched.get(tidx);
                if (sched == null) sched = SchedStrategy.FCFS;
                taskSchedStrategy[t] = (sched == SchedStrategy.REF) ? SchedStrategy.FCFS : sched;
            }

            taskQueues = new PriorityQueue[ntasks];
            for (int t = 0; t < ntasks; t++) {
                taskQueues[t] = new PriorityQueue<LNRequest>(11, getTaskComparator(taskSchedStrategy[t]));
            }

            buildEntryActivityMappings();
            buildActivityGraph();
            buildCacheState();
        }

        /**
         * Detect CacheTasks and their ItemEntries from the struct and set up the live
         * cache content and the read-access descriptors.
         *
         * <p>A CacheTask carries {@code lsn.iscache==1} on its task index together with
         * the number of items, per-level capacity vector and replacement strategy. Each
         * ItemEntry on that task carries a discrete popularity distribution
         * ({@code lsn.itemproc}). The read is issued as a synchronous (or asynchronous)
         * call to the ItemEntry, which enters the CacheTask at the entry's bound
         * activity -- the cache-access driver. That driver's POST_CACHE precedence names
         * two continuation activities: the first (lower activity index) is the hit
         * branch, the second the miss branch. Both reply to the ItemEntry.</p>
         */
        private void buildCacheState() {
            if (lsn.iscache == null) return;

            cacheStateByTask = new LNCacheState[ntasks];
            for (int t = 0; t < ntasks; t++) {
                int tidx = tshift + t + 1;
                if (tidx >= lsn.iscache.getNumCols() || lsn.iscache.get(0, tidx) != 1) continue;
                int numItems = (int) lsn.nitems.get(0, tidx);
                if (numItems <= 0) continue;
                int[] cap = (lsn.itemcap != null) ? lsn.itemcap.get(tidx) : null;
                if (cap == null || cap.length == 0) cap = new int[]{1};
                ReplacementStrategy rs = replacementStrategyOf(tidx);
                boolean retr = (lsn.hasretrieval != null && tidx < lsn.hasretrieval.getNumCols()
                        && lsn.hasretrieval.get(0, tidx) == 1);
                cacheStateByTask[t] = new LNCacheState(numItems, cap, rs, retr);
            }

            if (lsn.itemproc == null) return;
            for (int e = 1; e <= nentries; e++) {
                int eidx = eshift + e;
                DiscreteDistribution popularity = lsn.itemproc.get(eidx);
                if (popularity == null) continue;   // not an ItemEntry (or no popularity set)
                int tidx = (int) lsn.parent.get(0, eidx);
                int t = tidx - tshift - 1;
                if (t < 0 || t >= ntasks || cacheStateByTask[t] == null) continue;
                LNCacheState cache = cacheStateByTask[t];

                Integer driver = entryBoundActivity.get(eidx);
                if (driver == null) continue;
                List<Integer> succ = activitySuccessors.get(driver);
                if (succ == null || succ.size() < 2) continue;
                List<Integer> ordered = new ArrayList<Integer>(succ);
                Collections.sort(ordered);   // ascending index: [hit, miss]
                int hitAidx = ordered.get(0);
                int missAidx = ordered.get(1);

                int card = (int) lsn.nitems.get(0, eidx);
                if (card <= 0) card = cache.numItems;
                double[] cdf = buildPopularityCdf(popularity, card);
                if (cdf == null) continue;

                cacheAccessByDriver.put(driver, new CacheAccess(cache, cdf, hitAidx, missAidx, eidx));
                cacheBranchReply.put(hitAidx, eidx);
                cacheBranchReply.put(missAidx, eidx);
            }
        }

        private ReplacementStrategy replacementStrategyOf(int tidx) {
            ReplacementStrategy[] all = ReplacementStrategy.values();
            int ord = (lsn.replacestrat != null && tidx < lsn.replacestrat.getNumCols())
                    ? (int) lsn.replacestrat.get(0, tidx) : ReplacementStrategy.LRU.ordinal();
            if (ord < 0 || ord >= all.length) ord = ReplacementStrategy.LRU.ordinal();
            return all[ord];
        }

        /**
         * Build the cumulative popularity distribution over items 1..card from a discrete
         * distribution, normalized to sum 1. Mirrors the flat-network access sampler
         * (Solver_ssj.createAccessSamplers): evaluate the PMF over item indices 1..card.
         */
        private double[] buildPopularityCdf(DiscreteDistribution popularity, int card) {
            if (popularity == null || popularity.isDisabled()) return null;
            List<Double> itemIndices = new ArrayList<Double>(card);
            for (int i = 1; i <= card; i++) itemIndices.add((double) i);
            Matrix pmf = popularity.evalPMF(itemIndices);
            if (pmf == null || pmf.length() < card) return null;
            double[] cdf = new double[card];
            double sum = 0.0;
            for (int i = 0; i < card; i++) {
                double v = pmf.get(i);
                if (Double.isNaN(v) || v < 0) return null;
                sum += v;
                cdf[i] = sum;
            }
            if (sum <= 0) return null;
            for (int i = 0; i < card; i++) cdf[i] /= sum;
            return cdf;
        }

        private int sampleItem(double[] cdf) {
            double u = random.nextDouble();
            for (int i = 0; i < cdf.length; i++) {
                if (u <= cdf[i]) return i;
            }
            return cdf.length - 1;
        }

        /**
         * Perform a cache read: draw an item from the popularity distribution, test it
         * against the live cache content, update the content per the replacement policy,
         * and return the continuation activity (hit branch on a hit, miss branch on a
         * miss). Mirrors the flat-network {@code processCacheAccess}/{@code
         * handleCacheHit}/{@code handleCacheMiss}, restricted to the layered cache case
         * (no per-item access-cost matrix, no retrieval/delayed-hit system).
         */
        private int accessCache(CacheAccess ca, LNRequest request) {
            LNCacheState cs = ca.cache;
            cs.totalReads++;
            int item = sampleItem(ca.cdf);

            int hitLevel = -1;
            for (int level = 0; level < cs.levels.length; level++) {
                if (cs.levels[level].contains(Integer.valueOf(item))) {
                    hitLevel = level;
                    break;
                }
            }

            if (hitLevel >= 0) {
                cs.totalHits++;
                cacheHit(cs, item, hitLevel);
                return ca.hitAidx;
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            if (cs.hasRetrieval) {
                if (cs.inFlight[item]) {
                    if (cs.heldDelayedHits[item] == null) {
                        cs.heldDelayedHits[item] = new ArrayList<HeldLN>();
                    }
                    cs.heldDelayedHits[item].add(new HeldLN(request, ca));
                    return LN_CACHE_HELD;
                }
                cs.totalMisses++;
                cs.inFlight[item] = true;
                fetcherItem.put(request.requestId, new FetchKey(cs, item));
                return ca.missAidx;
            }

            cs.totalMisses++;
            cacheMiss(cs, item);
            return ca.missAidx;
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private void releaseDelayedHits(LNCacheState cs, int item) {
            cs.inFlight[item] = false;
            cacheMiss(cs, item);
            List<HeldLN> held = cs.heldDelayedHits[item];
            if (held == null || held.isEmpty()) return;
            for (HeldLN h : held) {
                cs.totalDelayedHits++;
                startActivityExecution(h.request, h.access.hitAidx);
            }
            held.clear();
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private void cacheHit(LNCacheState cs, int itemIdx, int hitLevel) {
            ReplacementStrategy rs = cs.replacementStrategy;
            int inew = Math.min(hitLevel + 1, cs.levels.length - 1);   // default: promote one list
            if (inew <= hitLevel) {
                if (rs == ReplacementStrategy.LRU || rs == ReplacementStrategy.HLRU
                        || rs == ReplacementStrategy.QLRU) {
                    cs.levels[hitLevel].remove(Integer.valueOf(itemIdx));
                    cs.levels[hitLevel].addFirst(itemIdx);
                }
                return;
            }
            LinkedList<Integer> lo = cs.levels[hitLevel];
            LinkedList<Integer> hi = cs.levels[inew];
            int kpos = lo.indexOf(Integer.valueOf(itemIdx));
            lo.remove(Integer.valueOf(itemIdx));
            if (hi.size() < cs.levelCapacities[inew]) {
                hi.addFirst(itemIdx);
                return;
            }
            int demoted;
            if (rs == ReplacementStrategy.RR) {
                int r = clampIndex((int) (random.nextDouble() * hi.size()), hi.size());
                demoted = hi.remove(r);
                hi.add(r, Integer.valueOf(itemIdx));
            } else {
                demoted = hi.removeLast();
                hi.addFirst(itemIdx);
            }
            if (rs == ReplacementStrategy.LRU || rs == ReplacementStrategy.SFIFO
                    || rs == ReplacementStrategy.HLRU || rs == ReplacementStrategy.QLRU) {
                lo.addFirst(Integer.valueOf(demoted));
            } else {
                int at = (kpos < 0 || kpos > lo.size()) ? lo.size() : kpos;
                lo.add(at, Integer.valueOf(demoted));
            }
        }

        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
        private void cacheMiss(LNCacheState cs, int itemIdx) {
            int target = 0;   // default: insert into the entry list
            ReplacementStrategy rs = cs.replacementStrategy;
            LinkedList<Integer> list = cs.levels[target];
            int cap = cs.levelCapacities[target];
            if (rs == ReplacementStrategy.RR) {
                if (list.size() >= cap && !list.isEmpty()) {
                    int r = clampIndex((int) (random.nextDouble() * list.size()), list.size());
                    list.set(r, Integer.valueOf(itemIdx));
                } else {
                    list.addFirst(itemIdx);
                }
            } else {
                list.addFirst(itemIdx);
                if (list.size() > cap) {
                    list.removeLast();
                }
            }
        }

        private int clampIndex(int r, int size) {
            if (r < 0) return 0;
            if (r > size - 1) return size - 1;
            return r;
        }

        private Comparator<HostRequest> getHostComparator(SchedStrategy strategy) {
            if (strategy == SchedStrategy.LCFS) return hostLcfsComparator;
            return hostFcfsComparator;
        }

        private Comparator<LNRequest> getTaskComparator(SchedStrategy strategy) {
            if (strategy == SchedStrategy.LCFS) return lcfsComparator;
            return fcfsComparator;
        }

        private boolean isPSScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.PS || strategy == SchedStrategy.DPS || strategy == SchedStrategy.GPS;
        }

        private boolean isINFScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.INF;
        }

        private void buildEntryActivityMappings() {
            for (int e = 1; e <= nentries; e++) {
                int eidx = eshift + e;
                int tidx = (int) lsn.parent.get(0, eidx);
                List<Integer> activities = lsn.actsof.get(tidx);
                if (activities == null) continue;

                for (Integer aidxObj : activities) {
                    int aidx = aidxObj;
                    int rawAidx = aidx - ashift;
                    if (lsn.replygraph != null && rawAidx >= 0 && rawAidx < lsn.replygraph.getNumRows()
                            && e >= 0 && e < lsn.replygraph.getNumCols()) {
                        if (lsn.replygraph.get(rawAidx, e) > 0) {
                            entryReplyActivity.put(eidx, aidx);
                        }
                    }
                }

                if (!activities.isEmpty()) {
                    int firstActivity = findBoundActivityForEntry(eidx, activities);
                    if (firstActivity > 0) {
                        entryBoundActivity.put(eidx, firstActivity);
                    }
                }
            }
            buildForwardingMap();
        }

        /**
         * Collect forwarding calls (CallType.FWD) from the struct into
         * entryForwards: source entry (callpair col 1) -> (target entry (col 2),
         * probability (callproc_mean)). Forwarding is entry-sourced, so it is not
         * carried in callsof and must be gathered separately.
         */
        private void buildForwardingMap() {
            if (lsn.calltype == null) return;
            for (Map.Entry<Integer, CallType> en : lsn.calltype.entrySet()) {
                if (en.getValue() != CallType.FWD) continue;
                int cidx = en.getKey();
                int srcEntry = (int) lsn.callpair.get(cidx, 1);
                int dstEntry = (int) lsn.callpair.get(cidx, 2);
                if (srcEntry <= 0 || dstEntry <= 0) continue;
                Double probObj = lsn.callproc_mean.get(cidx);
                double prob = (probObj != null) ? probObj : 1.0;
                List<Pair<Integer, Double>> lst = entryForwards.get(srcEntry);
                if (lst == null) {
                    lst = new ArrayList<Pair<Integer, Double>>();
                    entryForwards.put(srcEntry, lst);
                }
                lst.add(new Pair<Integer, Double>(dstEntry, prob));
            }
        }

        private int findBoundActivityForEntry(int eidx, List<Integer> activities) {
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            for (Integer aidxObj : activities) {
                int aidx = aidxObj;
                if (lsn.graph.get(eidx, aidx) > 0) {
                    return aidx;
                }
            }
            // Fallback (single-entry / legacy structs without the entry->activity edge):
            // first predecessor-free sequential activity in the task.
            for (Integer aidxObj : activities) {
                int aidx = aidxObj;
                int preType = (int) lsn.actpretype.get(0, aidx);
                if (preType == ActivityPrecedenceType.ID_PRE_SEQ || preType == 0) {
                    boolean hasPredecessor = false;
                    for (Integer otherAidxObj : activities) {
                        int otherAidx = otherAidxObj;
                        if (otherAidx != aidx && lsn.graph.get(otherAidx, aidx) > 0) {
                            hasPredecessor = true;
                            break;
                        }
                    }
                    if (!hasPredecessor) {
                        return aidx;
                    }
                }
            }
            return activities.isEmpty() ? -1 : activities.get(0);
        }

        private int findParentEntry(int aidx, int tidx, List<Integer> taskEntries) {
            if (taskEntries.size() == 1) return taskEntries.get(0);

            for (Integer eidxObj : taskEntries) {
                int eidx = eidxObj;
                Integer boundAct = entryBoundActivity.get(eidx);
                if (boundAct != null && boundAct == aidx) return eidx;
            }

            for (Integer eidxObj : taskEntries) {
                int eidx = eidxObj;
                Integer boundAct = entryBoundActivity.get(eidx);
                if (boundAct == null) continue;
                if (isActivityReachableFrom(boundAct, aidx, tidx)) return eidx;
            }

            return taskEntries.isEmpty() ? -1 : taskEntries.get(0);
        }

        private boolean isActivityReachableFrom(int startAidx, int targetAidx, int tidx) {
            if (startAidx == targetAidx) return true;
            Set<Integer> visited = new HashSet<Integer>();
            LinkedList<Integer> queue = new LinkedList<Integer>();
            queue.add(startAidx);
            while (!queue.isEmpty()) {
                int current = queue.poll();
                if (current == targetAidx) return true;
                if (visited.contains(current)) continue;
                visited.add(current);
                List<Integer> successors = activitySuccessors.get(current);
                if (successors == null) continue;
                for (Integer succ : successors) {
                    int succParent = (int) lsn.parent.get(0, succ);
                    if (succParent == tidx && !visited.contains(succ)) {
                        queue.add(succ);
                    }
                }
            }
            return false;
        }

        private void buildActivityGraph() {
            if (lsn.graph == null) return;

            for (int a = 1; a <= nacts; a++) {
                int aidx = ashift + a;
                List<Integer> successors = new ArrayList<Integer>();

                for (int a2 = 1; a2 <= nacts; a2++) {
                    int aidx2 = ashift + a2;
                    if (lsn.graph.get(aidx, aidx2) > 0) {
                        successors.add(aidx2);
                    }
                }

                if (!successors.isEmpty()) {
                    activitySuccessors.put(aidx, successors);
                }

                int postType = (int) lsn.actposttype.get(0, aidx);
                // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
                if ((postType == ActivityPrecedenceType.ID_POST_OR
                        || postType == ActivityPrecedenceType.ID_POST_LOOP) && successors.size() > 1) {
                    List<Pair<Integer, Double>> probs = new ArrayList<Pair<Integer, Double>>();
                    double totalProb = 0.0;
                    for (Integer succIdx : successors) {
                        double prob = lsn.graph.get(aidx, succIdx);
                        probs.add(new Pair<Integer, Double>(succIdx, prob));
                        totalProb += prob;
                    }
                    if (totalProb > 0 && totalProb != 1.0) {
                        for (int i = 0; i < probs.size(); i++) {
                            probs.set(i, new Pair<Integer, Double>(probs.get(i).getLeft(), probs.get(i).getRight() / totalProb));
                        }
                    }
                    orForkProbs.put(aidx, probs);
                }
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            for (int a = 1; a <= nacts; a++) {
                int aidx = ashift + a;
                if ((int) lsn.actpretype.get(0, aidx) == ActivityPrecedenceType.ID_PRE_AND) {
                    List<Integer> succ = activitySuccessors.get(aidx);
                    if (succ != null && !succ.isEmpty()) {
                        int joinIdx = succ.get(0);
                        Integer c = joinRequiredCount.get(joinIdx);
                        joinRequiredCount.put(joinIdx, (c == null ? 0 : c) + 1);
                    }
                }
            }
        }

        // Per-invocation key for a join activity: rootRequestId x C + joinIdx. joinIdx
        // (an absolute activity index) is far below C, so the encoding is collision-free.
        private static long joinKey(long rootRequestId, int joinIdx) {
            return rootRequestId * 1000003L + joinIdx;
        }

        public void simulate(double maxEventsParam) {
            maxEvents = (long) maxEventsParam;
            warmupEventThreshold = (long) (maxEvents * 0.2);
            ssjSim.init();
            startRefTasks();
            startOpenArrivals();
            ssjSim.start();
        }

        /**
         * Seed the first open (Poisson) arrival for every entry that carries an
         * open-arrival-rate (stored as an exponential inter-arrival in lsn.arrival).
         * Open arrivals drive entries with no closed reference-task caller.
         */
        private void startOpenArrivals() {
            if (lsn.arrival == null) return;
            for (int e = 1; e <= nentries; e++) {
                int eidx = eshift + e;
                @SuppressWarnings("deprecation")
                Distribution arr = lsn.arrival.get(eidx);
                if (arr != null) {
                    new OpenArrival(eidx).schedule(sampleDistribution(arr));
                }
            }
        }

        private void startRefTasks() {
            for (int t = 0; t < ntasks; t++) {
                int tidx = tshift + t + 1;
                SchedStrategy sched = lsn.sched.get(tidx);
                if (sched == SchedStrategy.REF) {
                    List<Integer> entries = lsn.entriesof.get(tidx);
                    if (entries != null && !entries.isEmpty()) {
                        int entryIdx = entries.get(0);
                        @SuppressWarnings("deprecation")
                        Distribution thinkDist = lsn.think.get(tidx);
                        int mult = taskMultiplicity[t];
                        for (int i = 0; i < mult; i++) {
                            double thinkTime = sampleDistribution(thinkDist);
                            new RefTaskThinkComplete(tidx, entryIdx).schedule(thinkTime);
                        }
                    }
                }
            }
        }

        private double sampleDistribution(Distribution dist) {
            if (dist == null) return 0.0;
            if (dist.isImmediate()) return 0.0;
            double[] samples = dist.sample(1, random);
            return (samples.length > 0) ? Math.max(0.0, samples[0]) : 0.0;
        }

        private void updateHostStats(int hostIdx) {
            int h = hostIdx - 1;
            if (h < 0 || h >= nhosts) return;
            double now = ssjSim.time();
            double elapsed = now - lastHostUpdateTime[h];
            if (elapsed > 0 && warmupDone) {
                SchedStrategy strategy = hostSchedStrategy[h];
                int busyServers = isINFScheduling(strategy) ? hostBusy[h] : Math.min(hostBusy[h], hostServers[h]);
                hostBusyTime[h] += (double) busyServers * elapsed;
            }
            lastHostUpdateTime[h] = now;
        }

        private void updateTaskStats(int taskIdx) {
            int t = taskIdx - tshift - 1;
            if (t < 0 || t >= ntasks) return;
            double now = ssjSim.time();
            double elapsed = now - lastTaskUpdateTime[t];
            if (elapsed > 0 && warmupDone) {
                taskUtilTime[t] += (double) taskBusy[t] * elapsed;
            }
            lastTaskUpdateTime[t] = now;
        }

        private void updateActivityStats(int activityIdx) {
            int a = activityIdx - ashift - 1;
            if (a < 0 || a >= nacts) return;
            double now = ssjSim.time();
            double elapsed = now - lastActivityUpdateTime[a];
            if (elapsed > 0 && warmupDone) {
                activityBusyTime[a] += (double) activityBusy[a] * elapsed;
            }
            lastActivityUpdateTime[a] = now;
        }

        private void updateTaskQueueLengthStats(int t) {
            if (t < 0 || t >= ntasks) return;
            double now = ssjSim.time();
            double elapsed = now - lastTaskQueueUpdateTime[t];
            if (elapsed > 0 && warmupDone) {
                taskQueueLengthTime[t] += (double) taskQueueLength[t] * elapsed;
            }
            lastTaskQueueUpdateTime[t] = now;
        }

        private void updateHostQueueLengthStats(int h) {
            if (h < 0 || h >= nhosts) return;
            double now = ssjSim.time();
            double elapsed = now - lastHostQueueUpdateTime[h];
            if (elapsed > 0 && warmupDone) {
                hostQueueLengthTime[h] += (double) hostQueueLength[h] * elapsed;
            }
            lastHostQueueUpdateTime[h] = now;
        }

        private void checkEventCountStop() {
            if (!warmupDone && totalEventCount >= warmupEventThreshold) {
                warmupDone = true;
                warmupEndTime = ssjSim.time();
                resetStatistics();
            }
            if (totalEventCount >= maxEvents) {
                ssjSim.stop();
            }
        }

        private void resetStatistics() {
            for (int h = 0; h < nhosts; h++) {
                hostBusyTime[h] = 0.0;
                hostCompletions[h] = 0;
                lastHostUpdateTime[h] = ssjSim.time();
                hostQueueLengthTime[h] = 0.0;
                lastHostQueueUpdateTime[h] = ssjSim.time();
            }
            for (int t = 0; t < ntasks; t++) {
                taskUtilTime[t] = 0.0;
                taskCompletions[t] = 0.0;
                lastTaskUpdateTime[t] = ssjSim.time();
                taskQueueLengthTime[t] = 0.0;
                lastTaskQueueUpdateTime[t] = ssjSim.time();
            }
            for (int e = 0; e < nentries; e++) {
                entryThroughput[e] = 0.0;
                entryResidenceTime[e] = 0.0;
                entryResponseTimeTally[e].init();
            }
            for (int a = 0; a < nacts; a++) {
                activityThroughput[a] = 0.0;
                activityServiceTimeTally[a].init();
                activityBusyTime[a] = 0.0;
                lastActivityUpdateTime[a] = ssjSim.time();
            }
            activityStartTimes.clear();
            // Cache read/hit/miss/delayed-hit counters (in-flight state is left intact).
            if (cacheStateByTask != null) {
                for (LNCacheState cs : cacheStateByTask) {
                    if (cs == null) continue;
                    cs.totalReads = 0;
                    cs.totalHits = 0;
                    cs.totalMisses = 0;
                    cs.totalDelayedHits = 0;
                }
            }
        }

        private LNRequest createRequest(int entryIdx, CallType callType, long callerRequestId,
                                        int callerTaskIdx, double callMean) {
            long requestId = nextRequestId++;
            LNRequest request = new LNRequest(
                    requestId, entryIdx, callType, ssjSim.time(), ssjSim.time(),
                    callerRequestId, callerTaskIdx, -1, -1
            );
            activeRequests.put(requestId, request);
            return request;
        }

        private LNRequest createRequest(int entryIdx, CallType callType, long callerRequestId, int callerTaskIdx) {
            return createRequest(entryIdx, callType, callerRequestId, callerTaskIdx, 1.0);
        }

        private void arriveAtTask(int tidx, LNRequest request) {
            int t = tidx - tshift - 1;
            if (t < 0 || t >= ntasks) return;

            updateTaskStats(tidx);
            updateTaskQueueLengthStats(t);
            taskQueueLength[t]++;
            request.taskArrivalTime = ssjSim.time();

            SchedStrategy strategy = taskSchedStrategy[t];
            int mult = taskMultiplicity[t];

            if (isINFScheduling(strategy)) {
                taskBusy[t]++;
                startRequestProcessing(tidx, request);
                return;
            }

            if (taskBusy[t] < mult) {
                taskBusy[t]++;
                startRequestProcessing(tidx, request);
            } else {
                taskQueues[t].add(request);
            }
        }

        private void startRequestProcessing(int tidx, LNRequest request) {
            // Service starts now that a task thread is held; entry response time is
            // measured from here (excludes time queued for a busy task thread).
            request.serviceStartTime = ssjSim.time();
            int entryIdx = request.entryIdx;
            Integer boundActivityIdx = entryBoundActivity.get(entryIdx);
            if (boundActivityIdx != null && boundActivityIdx > 0) {
                request.currentActivityIdx = boundActivityIdx;
                startActivityExecution(request, boundActivityIdx);
            } else {
                completeRequest(request);
            }
        }

        private void startActivityExecution(LNRequest request, int aidx) {
            request.currentActivityIdx = aidx;

            int a = aidx - ashift - 1;
            if (a >= 0 && a < nacts) {
                updateActivityStats(aidx);
                activityBusy[a]++;
                activityStartTimes.put(new Pair<Long, Integer>(request.requestId, aidx), ssjSim.time());
            }

            int tidx = (int) lsn.parent.get(0, aidx);
            int hidx = (int) lsn.parent.get(0, tidx);

            @SuppressWarnings("deprecation")
            Distribution hostDemandDist = lsn.hostdem.get(aidx);
            double hostDemand = sampleDistribution(hostDemandDist);

            if (hostDemand > 0) {
                queueForHost(hidx, request, aidx, hostDemand);
            } else {
                proceedAfterHost(request, aidx);
            }
        }

        /**
         * Apply the activity think-time (a pure off-processor delay, distinct from
         * the host demand) before the activity's calls/completion, then proceed.
         * LQN activity think-time contributes to residence but not processor busy time.
         */
        private void proceedAfterHost(LNRequest request, int aidx) {
            @SuppressWarnings("deprecation")
            Distribution thinkDist = (lsn.actthink != null) ? lsn.actthink.get(aidx) : null;
            double think = sampleDistribution(thinkDist);
            if (think > 0) {
                new ActivityThinkComplete(request, aidx).schedule(think);
            } else {
                executeActivityCalls(request, aidx);
            }
        }

        private void queueForHost(int hidx, LNRequest request, int aidx, double serviceTime) {
            int h = hidx - 1;
            if (h < 0 || h >= nhosts) {
                executeActivityCalls(request, aidx);
                return;
            }

            updateHostStats(hidx);

            HostRequest hostRequest = new HostRequest(request, aidx, serviceTime, ssjSim.time());
            SchedStrategy strategy = hostSchedStrategy[h];

            if (isINFScheduling(strategy)) {
                hostBusy[h]++;
                new ActivityHostComplete(hostRequest, hidx).schedule(serviceTime);
            } else if (isPSScheduling(strategy)) {
                startPSJob(h, hidx, hostRequest);
            } else {
                if (hostBusy[h] < hostServers[h]) {
                    hostBusy[h]++;
                    new ActivityHostComplete(hostRequest, hidx).schedule(serviceTime);
                } else {
                    hostQueues[h].add(hostRequest);
                }
            }
        }

        private void startPSJob(int h, int hidx, HostRequest hostRequest) {
            PSJob psJob = new PSJob(hostRequest, hostRequest.serviceTime, ssjSim.time());
            hostPSJobs[h].add(psJob);
            hostBusy[h]++;
            schedulePSUpdate(h, hidx);
        }

        private void schedulePSUpdate(int h, int hidx) {
            if (hostPSEvent[h] != null) {
                hostPSEvent[h].cancel();
            }

            List<PSJob> jobs = hostPSJobs[h];
            if (jobs.isEmpty()) {
                hostPSEvent[h] = null;
                return;
            }

            int numJobs = jobs.size();
            int numServers = hostServers[h];
            double sharePerJob = Math.min(1.0, (double) numServers / (double) numJobs);

            double minCompletionTime = Double.MAX_VALUE;
            double now = ssjSim.time();

            for (PSJob job : jobs) {
                double elapsed = now - job.lastUpdateTime;
                double served = elapsed * sharePerJob;
                job.remainingService = Math.max(0.0, job.remainingService - served);
                job.lastUpdateTime = now;

                double timeToComplete = (sharePerJob > 0) ? job.remainingService / sharePerJob : Double.MAX_VALUE;
                if (timeToComplete < minCompletionTime) {
                    minCompletionTime = timeToComplete;
                }
            }

            if (minCompletionTime < Double.MAX_VALUE && minCompletionTime > 0) {
                PSUpdateEvent event = new PSUpdateEvent(h, hidx);
                event.schedule(minCompletionTime);
                hostPSEvent[h] = event;
            } else if (minCompletionTime == 0.0) {
                processPSCompletions(h, hidx);
            }
        }

        private void processPSCompletions(int h, int hidx) {
            List<PSJob> jobs = hostPSJobs[h];
            int numServers = hostServers[h];
            int numJobs = jobs.size();
            double sharePerJob = (numJobs > 0) ? Math.min(1.0, (double) numServers / (double) numJobs) : 1.0;
            double now = ssjSim.time();

            List<PSJob> completedJobs = new ArrayList<PSJob>();
            for (PSJob job : jobs) {
                double elapsed = now - job.lastUpdateTime;
                double served = elapsed * sharePerJob;
                job.remainingService = Math.max(0.0, job.remainingService - served);
                job.lastUpdateTime = now;

                if (job.remainingService <= 1e-9) {
                    completedJobs.add(job);
                }
            }

            for (PSJob job : completedJobs) {
                jobs.remove(job);
                updateHostStats(hidx);
                hostBusy[h]--;
                if (warmupDone) {
                    hostCompletions[h]++;
                }
                executeActivityCalls(job.hostRequest.request, job.hostRequest.activityIdx);
            }

            if (!jobs.isEmpty()) {
                schedulePSUpdate(h, hidx);
            } else {
                hostPSEvent[h] = null;
            }
        }

        private void completeActivityHostDemand(HostRequest hostRequest, int hidx) {
            int h = hidx - 1;
            if (h < 0 || h >= nhosts) return;

            updateHostStats(hidx);
            hostBusy[h]--;
            if (warmupDone) {
                hostCompletions[h]++;
            }

            if (!hostQueues[h].isEmpty()) {
                HostRequest nextHostRequest = hostQueues[h].poll();
                hostBusy[h]++;
                new ActivityHostComplete(nextHostRequest, hidx).schedule(nextHostRequest.serviceTime);
            }

            proceedAfterHost(hostRequest.request, hostRequest.activityIdx);
        }

        private void executeActivityCalls(LNRequest request, int aidx) {
            List<Integer> calls = lsn.callsof.get(aidx);
            if (calls == null || calls.isEmpty()) {
                completeActivity(request, aidx);
                return;
            }

            LinkedList<Integer> syncCallQueue = new LinkedList<Integer>();

            for (Integer cidx : calls) {
                CallType callType = lsn.calltype.get(cidx);
                if (callType == null) continue;
                int targetEntryIdx = (int) lsn.callpair.get(cidx, 2);
                if (targetEntryIdx <= 0) continue;

                Double meanCallsObj = lsn.callproc_mean.get(cidx);
                double meanCalls = (meanCallsObj != null) ? meanCallsObj : 1.0;
                int numCalls = sampleCallCount(meanCalls);

                for (int c = 0; c < numCalls; c++) {
                    if (callType == CallType.SYNC) {
                        syncCallQueue.add(targetEntryIdx);
                    } else if (callType == CallType.ASYNC) {
                        executeAsyncCall(request, targetEntryIdx, meanCalls);
                    }
                }
            }

            if (!syncCallQueue.isEmpty()) {
                pendingSyncCallQueue.put(request.requestId, syncCallQueue);
                syncCallActivity.put(request.requestId, aidx);
                syncCallInFlight.put(request.requestId, false);
                executeNextSyncCall(request, aidx);
            } else {
                completeActivity(request, aidx);
            }
        }

        private void executeNextSyncCall(LNRequest request, int aidx) {
            LinkedList<Integer> queue = pendingSyncCallQueue.get(request.requestId);
            if (queue == null) return;

            if (queue.isEmpty()) {
                pendingSyncCallQueue.remove(request.requestId);
                syncCallActivity.remove(request.requestId);
                syncCallInFlight.remove(request.requestId);
                completeActivity(request, aidx);
                return;
            }

            int targetEntryIdx = queue.poll();
            syncCallInFlight.put(request.requestId, true);

            executeSyncCall(request, aidx, targetEntryIdx, 1.0);
        }

        private int sampleCallCount(double meanCalls) {
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            if (meanCalls <= 0) return 0;
            int n = (int) Math.floor(meanCalls);
            double frac = meanCalls - n;
            if (frac > 0 && random.nextDouble() < frac) {
                n++;
            }
            return n;
        }

        private void executeSyncCall(LNRequest callerRequest, int callerActivityIdx, int targetEntryIdx, double callMean) {
            int targetTidx = (int) lsn.parent.get(0, targetEntryIdx);
            LNRequest callRequest = createRequest(targetEntryIdx, CallType.SYNC,
                    callerRequest.requestId, (int) lsn.parent.get(0, callerRequest.entryIdx), callMean);
            arriveAtTask(targetTidx, callRequest);
        }

        private void executeAsyncCall(LNRequest callerRequest, int targetEntryIdx, double callMean) {
            int targetTidx = (int) lsn.parent.get(0, targetEntryIdx);
            LNRequest callRequest = createRequest(targetEntryIdx, CallType.ASYNC, -1L, -1, callMean);
            arriveAtTask(targetTidx, callRequest);
        }

        private void handleSyncReply(LNRequest completedRequest) {
            long callerRequestId = completedRequest.callerRequestId;
            if (callerRequestId < 0) return;

            LNRequest callerRequest = activeRequests.get(callerRequestId);
            if (callerRequest == null) return;

            syncCallInFlight.put(callerRequestId, false);

            Integer callerActivityIdx = syncCallActivity.get(callerRequestId);
            if (callerActivityIdx == null) return;

            executeNextSyncCall(callerRequest, callerActivityIdx);
        }

        private void completeActivity(LNRequest request, int aidx) {
            int a = aidx - ashift - 1;
            if (a >= 0 && a < nacts) {
                updateActivityStats(aidx);
                activityBusy[a]--;

                if (warmupDone) {
                    activityThroughput[a] += 1.0;
                    Pair<Long, Integer> startTimeKey = new Pair<Long, Integer>(request.requestId, aidx);
                    Double startTime = activityStartTimes.remove(startTimeKey);
                    if (startTime != null) {
                        double responseTime = ssjSim.time() - startTime;
                        activityServiceTimeTally[a].add(responseTime);
                    }
                } else {
                    activityStartTimes.remove(new Pair<Long, Integer>(request.requestId, aidx));
                }
            }

            totalEventCount++;
            checkEventCountStop();

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            CacheAccess cacheAccess = cacheAccessByDriver.get(aidx);
            if (cacheAccess != null) {
                int nextAidx = accessCache(cacheAccess, request);
                // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
                if (nextAidx != LN_CACHE_HELD) {
                    startActivityExecution(request, nextAidx);
                }
                return;
            }

            List<Integer> successors = activitySuccessors.get(aidx);

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            LNRequest root = request.getRoot();
            if ((request.parentForkId >= 0 || pendingWork.containsKey(root.requestId))
                    && (int) lsn.actpretype.get(0, aidx) == ActivityPrecedenceType.ID_PRE_AND
                    && successors != null && !successors.isEmpty()) {
                int joinIdx = successors.get(0);
                long key = joinKey(root.requestId, joinIdx);
                JoinProgress jp = joinProgress.get(key);
                if (jp == null) {
                    jp = new JoinProgress();
                    joinProgress.put(key, jp);
                }
                jp.arrivedInputs.add(aidx);
                jp.completionTimes.add(ssjSim.time());
                // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
                if (request.rootRequest != null) {
                    activeRequests.remove(request.requestId);
                }
                Integer req = joinRequiredCount.get(joinIdx);
                int required = (req == null) ? 1 : req;
                if (jp.arrivedInputs.size() >= required) {
                    joinProgress.remove(key);
                    // The join's inputs collapse into one continuing line.
                    Integer cur = pendingWork.get(root.requestId);
                    if (cur != null) {
                        pendingWork.put(root.requestId, cur - (required - 1));
                    }
                    startActivityExecution(root, joinIdx);
                }
                return;
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            if (!request.replied && request.parentForkId < 0) {
                boolean atReplyActivity = (entryReplyActivity.containsKey(request.entryIdx)
                        && entryReplyActivity.get(request.entryIdx) == aidx)
                        // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
                        || (cacheBranchReply.containsKey(aidx)
                            && cacheBranchReply.get(aidx) == request.entryIdx);
                boolean phaseBoundary = false;
                int aPhase = (aidx - ashift >= 1 && aidx - ashift <= nacts)
                        ? (int) lsn.actphase.get(0, aidx - ashift) : 1;
                if (aPhase == 1 && successors != null) {
                    for (Integer s : successors) {
                        int sk = s - ashift;
                        if (sk >= 1 && sk <= nacts && (int) lsn.actphase.get(0, sk) == 2) {
                            phaseBoundary = true;
                            break;
                        }
                    }
                }
                if (atReplyActivity || phaseBoundary) {
                    sendEntryReply(request);
                }
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            if (request.replied && !fetcherItem.isEmpty()) {
                FetchKey fk = fetcherItem.remove(request.requestId);
                if (fk != null) {
                    releaseDelayedHits(fk.cache, fk.item);
                }
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            int nsucc = (successors == null) ? 0 : successors.size();
            if (nsucc == 0) {
                checkAndCompleteRequest(request, aidx);
            } else if (nsucc == 1) {
                moveToSuccessorOrJoin(request, successors.get(0));
            } else {
                List<Pair<Integer, Double>> probs = orForkProbs.get(aidx);
                boolean isProbabilistic = (probs != null && !probs.isEmpty());
                if (!isProbabilistic) {
                    // OR-fork / loop edges carry branch probabilities (< 1, sum ~1);
                    // AND-fork edges are all unit weight.
                    for (Integer succ : successors) {
                        double w = lsn.graph.get(aidx, succ);
                        if (w != 1.0 && w > 0.0) {
                            isProbabilistic = true;
                            break;
                        }
                    }
                    if (isProbabilistic) {
                        List<Pair<Integer, Double>> newProbs = new ArrayList<Pair<Integer, Double>>();
                        double totalProb = 0.0;
                        for (Integer succIdx : successors) {
                            double prob = lsn.graph.get(aidx, succIdx);
                            newProbs.add(new Pair<Integer, Double>(succIdx, prob));
                            totalProb += prob;
                        }
                        if (totalProb > 0 && totalProb != 1.0) {
                            for (int i = 0; i < newProbs.size(); i++) {
                                newProbs.set(i, new Pair<Integer, Double>(newProbs.get(i).getLeft(), newProbs.get(i).getRight() / totalProb));
                            }
                        }
                        orForkProbs.put(aidx, newProbs);
                    }
                }
                if (isProbabilistic) {
                    int chosenSuccessor = chooseOrForkBranch(aidx, successors);
                    startActivityExecution(request, chosenSuccessor);
                } else {
                    handleAndFork(request, aidx, successors);
                }
            }
        }

        private void moveToSuccessorOrJoin(LNRequest request, int successorIdx) {
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            startActivityExecution(request, successorIdx);
        }

        private void handleAndFork(LNRequest request, int forkActivityIdx, List<Integer> successors) {
            long parentForkId = nextRequestId++;

            if (DEBUG_FORK_JOIN) {
                String forkName = lsn.hashnames.get(forkActivityIdx);
                StringBuilder succNames = new StringBuilder();
                for (Integer s : successors) {
                    if (succNames.length() > 0) succNames.append(", ");
                    succNames.append(lsn.hashnames.get(s));
                }
                System.out.println("[FORK] time=" + ssjSim.time() + " " + forkName + " -> [" + succNames + "]");
            }

            // One line becomes successors.size() parallel lines.
            long rootId = request.getRoot().requestId;
            Integer cur = pendingWork.get(rootId);
            pendingWork.put(rootId, (cur == null ? 1 : cur) + successors.size() - 1);

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            for (Integer succIdxObj : successors) {
                int succIdx = succIdxObj;
                LNRequest branchRequest = request.copyWith(nextRequestId++, succIdx, parentForkId);
                activeRequests.put(branchRequest.requestId, branchRequest);
                startActivityExecution(branchRequest, succIdx);
            }
        }

        private int chooseOrForkBranch(int aidx, List<Integer> successors) {
            List<Pair<Integer, Double>> probs = orForkProbs.get(aidx);
            if (probs != null && !probs.isEmpty()) {
                double u = random.nextDouble();
                double cumProb = 0.0;
                for (Pair<Integer, Double> pr : probs) {
                    cumProb += pr.getRight();
                    if (u <= cumProb) {
                        return pr.getLeft();
                    }
                }
                return probs.get(probs.size() - 1).getLeft();
            }
            return successors.get(random.nextInt(successors.size()));
        }

        private void checkAndCompleteRequest(LNRequest request, int aidx) {
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            finishLine(request);
        }

        private void finishLine(LNRequest request) {
            LNRequest root = request.getRoot();
            Integer w = pendingWork.get(root.requestId);
            if (w == null) {
                // No fork occurred: a plain sequential invocation, one line only.
                completeRequest(request);
                return;
            }
            int remaining = w - 1;
            if (remaining <= 0) {
                pendingWork.remove(root.requestId);
                completeRequest(root);
            } else {
                pendingWork.put(root.requestId, remaining);
                // Retire this finished branch; the dormant main line stays until the
                // last branch drains and triggers completion on the root.
                if (request.rootRequest != null) {
                    activeRequests.remove(request.requestId);
                }
            }
        }

        /**
         * Unblock the synchronous caller and record entry throughput and service
         * time. Sent at the end of phase 1 (reply activity / phase-1->phase-2
         * boundary); idempotent so a two-phase entry can call it again at full
         * completion without double counting. Entry service time is measured from
         * service start (thread acquisition) to the reply, matching lqsim's entry
         * service-time / phase-utilization convention.
         */
        private void sendEntryReply(LNRequest request) {
            if (request.replied) return;
            request.replied = true;
            int e = request.entryIdx - eshift - 1;
            double responseTime = ssjSim.time() - request.serviceStartTime;
            if (e >= 0 && e < nentries && warmupDone) {
                entryThroughput[e] += 1.0;
                entryResponseTimeTally[e].add(responseTime);
            }
            if (request.callType != CallType.SYNC || request.callerRequestId < 0) {
                return;
            }
            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            List<Pair<Integer, Double>> fwds = entryForwards.get(request.entryIdx);
            if (fwds != null && !fwds.isEmpty()) {
                double u = random.nextDouble();
                double cum = 0.0;
                for (Pair<Integer, Double> fw : fwds) {
                    cum += fw.getRight();
                    if (u < cum) {
                        int targetEidx = fw.getLeft();
                        int targetTidx = (int) lsn.parent.get(0, targetEidx);
                        LNRequest fwdReq = createRequest(targetEidx, CallType.SYNC,
                                request.callerRequestId, request.callerTaskIdx);
                        arriveAtTask(targetTidx, fwdReq);
                        return;   // forwarded; do not reply from this entry
                    }
                }
                // u >= sum of forwarding probabilities: reply normally
            }
            handleSyncReply(request);
        }

        private void completeRequest(LNRequest request) {
            int entryIdx = request.entryIdx;
            int eRes = entryIdx - eshift - 1;

            // Full residence (both phases) drives entry queue length / utilization.
            if (eRes >= 0 && eRes < nentries && warmupDone) {
                entryResidenceTime[eRes] += ssjSim.time() - request.serviceStartTime;
            }

            // Reply now if it was not already sent at a phase-1 boundary (single-phase
            // entries, or a fallback). The task thread is freed below, after phase 2.
            sendEntryReply(request);

            int tidx = (int) lsn.parent.get(0, entryIdx);
            int t = tidx - tshift - 1;
            if (t >= 0 && t < ntasks) {
                updateTaskStats(tidx);
                updateTaskQueueLengthStats(t);
                taskQueueLength[t]--;
                taskBusy[t]--;
                if (warmupDone) {
                    taskCompletions[t] += 1.0;
                }

                if (!taskQueues[t].isEmpty()) {
                    LNRequest nextRequest = taskQueues[t].poll();
                    taskBusy[t]++;
                    startRequestProcessing(tidx, nextRequest);
                }
            }

            if (request.callerRequestId < 0 && request.callerTaskIdx < 0) {
                SchedStrategy sched = lsn.sched.get(tidx);
                if (sched == SchedStrategy.REF) {
                    scheduleNextRefTaskRequest(tidx, entryIdx);
                }
            }

            activeRequests.remove(request.requestId);
        }

        public LNLDESResult getLNLDESResult() {
            LNLDESResult result = new LNLDESResult();
            result.lsn = lsn;

            double simTime = ssjSim.time() - warmupEndTime;
            if (simTime <= 0) {
                return result;
            }

            result.QLN = new Matrix(1, lsn.nidx + 1);
            result.ULN = new Matrix(1, lsn.nidx + 1);
            result.RLN = new Matrix(1, lsn.nidx + 1);
            result.WLN = new Matrix(1, lsn.nidx + 1);
            result.TLN = new Matrix(1, lsn.nidx + 1);
            result.ALN = new Matrix(1, lsn.nidx + 1);

            for (int h = 0; h < nhosts; h++) {
                int hidx = h + 1;
                double util = hostBusyTime[h] / simTime;
                double tput = (double) hostCompletions[h] / simTime;
                result.ULN.set(0, hidx, util);
                result.TLN.set(0, hidx, tput);
            }

            for (int t = 0; t < ntasks; t++) {
                int tidx = tshift + t + 1;
                result.ULN.set(0, tidx, 0.0);
                double tput = taskCompletions[t] / simTime;
                result.TLN.set(0, tidx, tput);
                updateTaskQueueLengthStats(t);
                double avgQLen = taskQueueLengthTime[t] / simTime;
                result.QLN.set(0, tidx, avgQLen);
            }

            for (int e = 0; e < nentries; e++) {
                int eidx = eshift + e + 1;
                double tput = entryThroughput[e] / simTime;
                double respT = (entryResponseTimeTally[e].numberObs() > 0) ? entryResponseTimeTally[e].average() : 0.0;
                result.TLN.set(0, eidx, tput);
                result.RLN.set(0, eidx, respT);
                result.ALN.set(0, eidx, tput);
                // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
                double qLen = entryResidenceTime[e] / simTime;
                result.QLN.set(0, eidx, qLen);
                result.ULN.set(0, eidx, 0.0);
            }

            double[] activityProcUtil = new double[nacts];

            for (int a = 0; a < nacts; a++) {
                int aidx = ashift + a + 1;
                double tput = activityThroughput[a] / simTime;
                result.TLN.set(0, aidx, tput);

                Double meanHostDemObj = lsn.hostdem_mean.get(aidx);
                double meanHostDem = (meanHostDemObj != null) ? meanHostDemObj : 0.0;
                double procUtil = tput * meanHostDem;
                activityProcUtil[a] = procUtil;

                int tidx = (int) lsn.parent.get(0, aidx);
                if (tidx > 0 && tidx <= lsn.nidx) {
                    double currentTaskUtil = result.ULN.get(0, tidx);
                    result.ULN.set(0, tidx, currentTaskUtil + procUtil);

                    List<Integer> taskEntries = lsn.entriesof.get(tidx);
                    if (taskEntries != null && !taskEntries.isEmpty()) {
                        int parentEntryIdx = findParentEntry(aidx, tidx, taskEntries);
                        if (parentEntryIdx > 0) {
                            double currentEntryUtil = result.ULN.get(0, parentEntryIdx);
                            result.ULN.set(0, parentEntryIdx, currentEntryUtil + procUtil);
                        }
                    }
                }

                double respT = (activityServiceTimeTally[a].numberObs() > 0) ? activityServiceTimeTally[a].average() : 0.0;
                result.RLN.set(0, aidx, respT);
                double qLen = tput * respT;
                result.QLN.set(0, aidx, qLen);
            }

            for (int a = 0; a < nacts; a++) {
                int aidx = ashift + a + 1;
                double procUtil = activityProcUtil[a];
                int tidx = (int) lsn.parent.get(0, aidx);
                double taskUtil = (tidx > 0 && tidx <= lsn.nidx) ? result.ULN.get(0, tidx) : 0.0;

                double activityUtil;
                if (taskUtil > 0.0 && taskUtil < 1.0) {
                    activityUtil = Math.min(1.0, procUtil / taskUtil);
                } else if (procUtil > 0.0) {
                    activityUtil = 1.0;
                } else {
                    activityUtil = 0.0;
                }
                result.ULN.set(0, aidx, activityUtil);
            }

            // see _kb/09-ldes-and-cache.md (Layered (LN) engine: Solver_ssj_ln.java notes)
            if (cacheStateByTask != null) {
                for (CacheAccess ca : cacheAccessByDriver.values()) {
                    LNCacheState cs = ca.cache;
                    long reads = cs.totalReads;
                    if (reads <= 0) continue;
                    int tidx = (int) lsn.parent.get(0, ca.itemEntryIdx);
                    result.cacheTaskIdx.add(tidx);
                    result.cacheItemEntryIdx.add(ca.itemEntryIdx);
                    result.cacheHitProb.add((double) cs.totalHits / reads);
                    result.cacheMissProb.add((double) cs.totalMisses / reads);
                    result.cacheDelayedProb.add((double) cs.totalDelayedHits / reads);
                    result.cacheReadRate.add(reads / simTime);
                }
            }

            return result;
        }

        // ==================== Event Classes ====================

        private class RefTaskThinkComplete extends SimEvent {
            private final int tidx;
            private final int entryIdx;

            public RefTaskThinkComplete(int tidx, int entryIdx) {
                this.tidx = tidx;
                this.entryIdx = entryIdx;
            }

            @Override
            public void actions() {
                LNRequest request = createRequest(entryIdx, CallType.ASYNC, -1L, -1);
                arriveAtTask(tidx, request);
            }
        }

        /**
         * An open (Poisson) arrival to an entry: dispatches one open-class request
         * (no caller) to the entry's task, then schedules the next arrival from the
         * entry's exponential inter-arrival distribution.
         */
        private class OpenArrival extends SimEvent {
            private final int entryIdx;

            public OpenArrival(int entryIdx) {
                this.entryIdx = entryIdx;
            }

            @Override
            public void actions() {
                int tidx = (int) lsn.parent.get(0, entryIdx);
                LNRequest request = createRequest(entryIdx, CallType.ASYNC, -1L, -1);
                arriveAtTask(tidx, request);
                @SuppressWarnings("deprecation")
                Distribution arr = lsn.arrival.get(entryIdx);
                if (arr != null) {
                    new OpenArrival(entryIdx).schedule(sampleDistribution(arr));
                }
            }
        }

        private void scheduleNextRefTaskRequest(int tidx, int entryIdx) {
            @SuppressWarnings("deprecation")
            Distribution thinkDist = lsn.think.get(tidx);
            double thinkTime = sampleDistribution(thinkDist);
            new RefTaskThinkComplete(tidx, entryIdx).schedule(thinkTime);
        }

        private class ActivityHostComplete extends SimEvent {
            private final HostRequest hostRequest;
            private final int hidx;

            public ActivityHostComplete(HostRequest hostRequest, int hidx) {
                this.hostRequest = hostRequest;
                this.hidx = hidx;
            }

            @Override
            public void actions() {
                completeActivityHostDemand(hostRequest, hidx);
            }
        }

        private class ActivityThinkComplete extends SimEvent {
            private final LNRequest request;
            private final int aidx;

            public ActivityThinkComplete(LNRequest request, int aidx) {
                this.request = request;
                this.aidx = aidx;
            }

            @Override
            public void actions() {
                executeActivityCalls(request, aidx);
            }
        }

        private class PSUpdateEvent extends SimEvent {
            private final int h;
            private final int hidx;

            public PSUpdateEvent(int h, int hidx) {
                this.h = h;
                this.hidx = hidx;
            }

            @Override
            public void actions() {
                processPSCompletions(h, hidx);
            }
        }
    }
}
