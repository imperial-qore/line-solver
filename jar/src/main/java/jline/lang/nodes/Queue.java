/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.ServiceBinding;
import jline.lang.constant.ServerType;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.sections.*;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.util.Utils.isInf;

/**
 * A queueing station that processes jobs according to various scheduling strategies.
 * 
 * <p>The Queue node is the fundamental service station in queueing networks. It models
 * a system where jobs wait in a buffer and are processed by one or more servers according
 * to a scheduling strategy. Common examples include CPU schedulers, disk I/O queues,
 * network routers, and service desks.</p>
 * 
 * <p>Key features:
 * <ul>
 *   <li>Multiple scheduling strategies: FCFS, PS, LCFS, priority-based, etc.</li>
 *   <li>Single or multi-server configurations</li>
 *   <li>Load-dependent and class-dependent service rates</li>
 *   <li>Support for preemptive and non-preemptive policies</li>
 *   <li>Polling server capabilities with various polling types</li>
 *   <li>Switchover times between job classes</li>
 * </ul>
 * </p>
 * 
 * @see SchedStrategy
 * @see ServiceStation
 * @see Server
 * @since 1.0
 */
public class Queue extends ServiceStation implements Serializable {

    /**
     * Per-class setup time distributions (server activation time).
     * Maps job classes to their setup time distributions.
     */
    protected Map<JobClass, Distribution> setupTimes;
    
    /**
     * Per-class delay-off time distributions (idle time before power-off).
     * Maps job classes to their delay-off time distributions.
     */
    protected Map<JobClass, Distribution> delayOffTimes;

    /**
     * LPS limit: maximum number of jobs that can execute in PS mode for LPS scheduling.
     * Only used when schedStrategy is LPS.
     */
    private double lpsLimit = Double.NaN;

    /**
     * List of server types for heterogeneous multiserver queues.
     * When non-empty, this queue uses heterogeneous server configuration.
     */
    private List<ServerType> serverTypes;

    /**
     * Scheduling policy for heterogeneous servers.
     * Determines how jobs are assigned to server types when multiple compatible types exist.
     */
    private HeteroSchedPolicy heteroSchedPolicy;

    /**
     * Service distributions per server type and job class for heterogeneous queues.
     * Maps ServerType -> JobClass -> Distribution.
     */
    private Map<ServerType, Map<JobClass, Distribution>> heteroServiceDistributions;

    /**
     * Number of servers a job seizes for the whole of its service, per class.
     * An absent entry means one server, the homogeneous default.
     */
    private Map<JobClass, Integer> serverParallelism = new HashMap<JobClass, Integer>();

    /**
     * Set of class indices with immediate feedback enabled.
     * When a job self-loops at this station, it stays in service instead of rejoining queue.
     * null means no classes have immediate feedback, "all" placeholder means all classes.
     */
    private java.util.Set<Integer> immediateFeedbackClasses;

    /**
     * Flag indicating if immediate feedback is enabled for all classes.
     */
    private boolean immediateFeedbackAll = false;

    /**
     * Pass-and-swap (PAS) class compatibility/swap graph (nclasses x nclasses).
     * Entry (r,s) nonzero iff a completing class-r job may take the place of a
     * class-s job (order-independent swap). Undirected; self-loops allowed.
     */
    private Matrix swapGraph = null;

    /**
     * Pass-and-swap (PAS) total service rate function mu(c): maps the ordered
     * state vector c (a 1xN row Matrix of 0-based class indices, c(0) oldest)
     * to the scalar total service rate. The rate of position i is the increment
     * mu(c[0..i]) - mu(c[0..i-1]).
     */
    private SerializableFunction<Matrix, Double> svcRateFun = null;

    /**
     * Time to failure of the server. The clock runs whenever the server is up,
     * busy or idle. null means the server never fails.
     */
    private Distribution breakdownFailure = null;

    /**
     * Repair time of a down server. null means the server never fails.
     */
    private Distribution breakdownRepair = null;

    /**
     * Per-class service distribution used while the server is down. A map with a
     * single null key entry is not used: a single distribution applied to every
     * class is stored under every class by setBreakdown. Empty means the server
     * does not serve at all while down, which is the usual breakdown model.
     */
    private Map<JobClass, Distribution> breakdownDownService = new HashMap<JobClass, Distribution>();

    /**
     * Down-service distribution applied to every class when no per-class map was
     * given. Kept separately because the class list may not be final when
     * setBreakdown is called.
     */
    private Distribution breakdownDownServiceAll = null;

    /**
     * Creates a new queueing station with the specified scheduling strategy.
     * Configures the appropriate server type based on the scheduling strategy.
     * The queue is automatically added to the network model.
     * 
     * @param model the network model to add this queue to
     * @param name the name for this queueing station
     * @param schedStrategy the scheduling strategy to use (e.g., FCFS, PS, LCFS)
     * @throws RuntimeException if the scheduling strategy is not supported
     */
    public Queue(Network model, String name, SchedStrategy schedStrategy) {
        super(name);

        this.setModel(model);
        this.model.addNode(this);
        this.schedStrategy = schedStrategy;
        this.input = new Buffer(model.getClasses());
        this.output = new Dispatcher(model.getClasses());

        this.schedStrategy = schedStrategy;
        this.numberOfServers = 1;
        this.setupTimes = new HashMap<JobClass, Distribution>();
        this.delayOffTimes = new HashMap<JobClass, Distribution>();
        this.serverTypes = new ArrayList<ServerType>();
        this.heteroSchedPolicy = HeteroSchedPolicy.ORDER;
        this.heteroServiceDistributions = new HashMap<ServerType, Map<JobClass, Distribution>>();
        this.immediateFeedbackClasses = new java.util.HashSet<Integer>();
        this.immediateFeedbackAll = false;

        switch (this.schedStrategy) {
            case FCFS:
            case LCFS:
            case SIRO:
            case SEPT:
            case LEPT:
            case SJF:
            case LJF:
            case EDD:
            case SRPT:
            case SRPTPRIO:
            case SETF:
            case FSP:
            case PAS:
            case OI:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new Server(model.getClasses());
                break;
            case INF:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new InfiniteServer(model.getClasses());
                this.numberOfServers = Integer.MAX_VALUE;
                break;
            case HOL:
            case FCFSPRIO:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new Server(model.getClasses());
                break;
            case PS:
            case DPS:
            case GPS:
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
            case LPS:
                this.schedPolicy = SchedStrategyType.PR;
                this.server = new SharedServer(model.getClasses());
                break;
            case LCFSPR:
            case LCFSPRPRIO:
            case FCFSPR:
            case FCFSPRPRIO:
            case EDF:
            case FB:
            case PSJF:
            case LRPT:
                this.schedPolicy = SchedStrategyType.PR;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case LCFSPI:
            case LCFSPIPRIO:
            case FCFSPI:
            case FCFSPIPRIO:
                this.schedPolicy = SchedStrategyType.PNR;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case LCFSPRIO:
                this.schedPolicy = SchedStrategyType.NPPrio;
                this.server = new PreemptiveServer(model.getClasses());
                break;
            case POLLING:
                this.schedPolicy = SchedStrategyType.NP;
                this.server = new PollingServer(model.getClasses());
                break;
            default:
                throw new RuntimeException("Routing Strategy is not supported in JLINE");
        }
    }

    /**
     * Creates a new queueing station with default processor sharing (PS) scheduling.
     * 
     * @param model The network model to add this queue to
     * @param name The name for this queueing station
     */
    public Queue(Network model, String name) {
        this(model, name, SchedStrategy.PS);
    }

    /**
     * Checks whether this queue has a service process configured for the specified job class.
     * 
     * @param jobClass The job class to check
     * @return true if a service process exists for this job class, false otherwise
     */
    public boolean containsJobClass(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Gets the scheduling policy type (preemptive or non-preemptive).
     * 
     * @return The scheduling policy type (PR for preemptive, NP for non-preemptive)
     */
    public SchedStrategyType getSchedPolicy() {
        return this.schedPolicy;
    }

    /**
     * Sets the total service rate function mu(c) of a pass-and-swap (PAS) queue.
     * An order-independent/PAS queue is parameterized by mu(c) as a whole, not
     * by per-class service distributions; this overload of setService accepts
     * the mu(c) function directly.
     *
     * @param muFun maps the ordered state vector c (1xN row Matrix of 0-based
     *              class indices, c(0) oldest) to the scalar total service rate
     */
    public void setService(SerializableFunction<Matrix, Double> muFun) {
        setServiceRateFunction(muFun);
    }

    /**
     * Sets the total service rate function mu(c) of a pass-and-swap (PAS) queue.
     * Also derives a representative per-class rate mu([r]) (single class-r job)
     * so the standard rate/process machinery stays consistent; the authoritative
     * service description remains mu(c).
     *
     * @param muFun the mu(c) service rate function
     */
    public void setServiceRateFunction(SerializableFunction<Matrix, Double> muFun) {
        if (this.schedStrategy != SchedStrategy.PAS && this.schedStrategy != SchedStrategy.OI) {
            throw new RuntimeException("setServiceRateFunction is only applicable to PAS (pass-and-swap) and OI (order-independent) queues.");
        }
        // svcRateFun is copied into QueueNodeParam by refreshLocalVars, so a cached
        // struct would keep the previous mu(c).
        this.svcRateFun = muFun;
        invalidateStruct();
        List<JobClass> classes = this.model.getClasses();
        for (int r = 0; r < classes.size(); r++) {
            Matrix single = new Matrix(1, 1);
            single.set(0, 0, r);
            double rate = muFun.apply(single);
            if (!Double.isInfinite(rate) && rate > 0) {
                this.setService(classes.get(r), new Exp(rate));
            } else {
                this.setService(classes.get(r), new jline.lang.processes.Disabled());
            }
        }
    }

    /**
     * Returns the mu(c) service rate function of a PAS queue, or null.
     *
     * @return the mu(c) function
     */
    public SerializableFunction<Matrix, Double> getServiceRateFunction() {
        return this.svcRateFun;
    }


    /** Result of {@link #checkPermInvariance}. */
    public static class PermCheckResult {
        public boolean ok = true;
        public int[] badc = null;
        public boolean partial = false;
    }

    /**
     * Checks the order-independence (OI) condition on the service rate mu(c):
     * the rate of the job in position j must depend only on the jobs at or
     * ahead of it (positions 1..j) and not on those behind. Since the
     * position-j rate is the prefix increment mu(c[:j]) - mu(c[:j-1]),
     * tail-independence is structural; the substantive requirement is that this
     * increment be independent of the order of the jobs ahead, which (by
     * induction on prefix length) is equivalent to mu(c) being permutation-
     * invariant. Enumerates the reachable multisets (per-class counts &lt;= Nvec,
     * total &lt;= cap) when small, otherwise samples and sets partial=true.
     */
    public PermCheckResult checkPermInvariance(double[] Nvec, double cap) {
        PermCheckResult res = new PermCheckResult();
        final SerializableFunction<Matrix, Double> mu = this.svcRateFun;
        if (mu == null) {
            return res;
        }
        int K = Nvec.length;
        double tol = 1e-9;
        int PERM_ENUM = 5040;      // enumerate all distinct permutations up to this
        int PERM_SAMPLE = 16;      // permutations sampled per multiset above PERM_ENUM
        int LATTICE_BUDGET = 4096;
        int MAXEVAL = 50000;
        boolean hasOpen = false;
        double sumFinite = 0;
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(Nvec[r])) {
                hasOpen = true;
            } else {
                sumFinite += Nvec[r];
            }
        }
        double Lmaxd = (Double.isFinite(cap) && cap >= 0 && cap < 1e18) ? cap : sumFinite;
        int[] ub = new int[K];
        for (int r = 0; r < K; r++) {
            double u = Double.isFinite(Nvec[r]) ? Nvec[r] : Math.min(Lmaxd, 6.0);
            ub[r] = (int) Math.min(u, Lmaxd);
        }
        if (!(Lmaxd >= 2) || !Double.isFinite(Lmaxd)) {
            return res;
        }
        int Lmax = (int) Lmaxd;
        long lat = 1;
        for (int r = 0; r < K; r++) {
            lat *= (ub[r] + 1L);
        }
        boolean exhaustive = !hasOpen && lat <= LATTICE_BUDGET;
        java.util.Random rng = new java.util.Random(0);
        int[] neval = {0};
        if (exhaustive) {
            int[] n = new int[K];
            while (true) {
                int sum = 0, nz = 0;
                for (int r = 0; r < K; r++) {
                    sum += n[r];
                    if (n[r] > 0) nz++;
                }
                if (sum >= 2 && nz >= 2 && sum <= Lmax) {
                    testMultiset(n, K, mu, tol, PERM_ENUM, PERM_SAMPLE, rng, neval, res);
                    if (!res.ok) break;
                    if (neval[0] >= MAXEVAL) {
                        res.partial = true;
                        break;
                    }
                }
                int d = 0;
                while (d < K) {
                    n[d]++;
                    if (n[d] <= ub[d]) break;
                    n[d] = 0;
                    d++;
                }
                if (d >= K) break;
            }
        } else {
            res.partial = true;
            for (int trial = 0; trial < 400; trial++) {
                int len = 2 + rng.nextInt(Math.max(1, Math.min(Lmax, 6) - 1));
                int[] n = new int[K];
                for (int j = 0; j < len; j++) {
                    int r = rng.nextInt(K);
                    if (n[r] < ub[r]) n[r]++;
                }
                int sum = 0, nz = 0;
                for (int r = 0; r < K; r++) {
                    sum += n[r];
                    if (n[r] > 0) nz++;
                }
                if (sum >= 2 && nz >= 2) {
                    testMultiset(n, K, mu, tol, PERM_ENUM, PERM_SAMPLE, rng, neval, res);
                    if (!res.ok || neval[0] >= MAXEVAL) break;
                }
            }
        }
        return res;
    }

    /** Result of {@link #checkRateMonotonicity}. */
    public static class MonoCheckResult {
        public boolean ok = true;
        public int[] badc = null;
        public int badr = -1;
        public boolean partial = false;
    }

    /**
     * Checks the order-independence (OI) condition (1) on the service rate
     * mu(c): the per-job rates must be non-negative, mu(c[:j]) &gt;= mu(c[:j-1])
     * for every microstate and position j.
     *
     * <p>A rate can be permutation-invariant and still fail to parameterize an
     * OI queue. Single-server processor sharing with class-dependent rates,
     * mu(c) = (sum_j mu_{c_j}) / n, is the standard trap: it is flatly
     * invariant under permutations, yet as soon as two classes have different
     * rates its prefix increments go negative -- with mu_hit = 3.0 and
     * mu_miss = 0.7, mu(Hit) = 3.0 while mu(Hit,Miss) = 1.85, so the second job
     * would be served at -1.15.
     *
     * <p>Run this AFTER {@link #checkPermInvariance}: permutation invariance is
     * what makes mu a function of the count vector, and the increments to test
     * are then just mu(n + e_r) - mu(n) over count vectors n and classes r, with
     * no permutation enumeration. Prefixes are non-empty, so mu is never
     * evaluated on an empty microstate, and an increment of exactly zero is
     * accepted -- that is how a class which does not visit this station is
     * expressed. Enumerates the reachable count vectors when small, otherwise
     * samples and sets partial=true.
     */
    public MonoCheckResult checkRateMonotonicity(double[] Nvec, double cap) {
        MonoCheckResult res = new MonoCheckResult();
        final SerializableFunction<Matrix, Double> mu = this.svcRateFun;
        if (mu == null) {
            return res;
        }
        int K = Nvec.length;
        double tol = 1e-9;
        int LATTICE_BUDGET = 4096;
        int MAXEVAL = 50000;
        boolean hasOpen = false;
        double sumFinite = 0;
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(Nvec[r])) {
                hasOpen = true;
            } else {
                sumFinite += Nvec[r];
            }
        }
        double Lmaxd = (Double.isFinite(cap) && cap >= 0 && cap < 1e18) ? cap : sumFinite;
        int[] ub = new int[K];
        for (int r = 0; r < K; r++) {
            double u = Double.isFinite(Nvec[r]) ? Nvec[r] : Math.min(Lmaxd, 6.0);
            ub[r] = (int) Math.min(u, Lmaxd);
        }
        if (!(Lmaxd >= 2) || !Double.isFinite(Lmaxd)) {
            return res;
        }
        int Lmax = (int) Lmaxd;
        long lat = 1;
        for (int r = 0; r < K; r++) {
            lat *= (ub[r] + 1L);
        }
        boolean exhaustive = !hasOpen && lat <= LATTICE_BUDGET;
        java.util.Random rng = new java.util.Random(0);
        int[] neval = {0};
        if (exhaustive) {
            int[] n = new int[K];
            while (true) {
                int sum = 0;
                for (int r = 0; r < K; r++) sum += n[r];
                if (sum >= 1 && sum <= Lmax - 1) {
                    stepPrefix(n, ub, K, mu, tol, neval, res);
                    if (!res.ok) break;
                    if (neval[0] >= MAXEVAL) {
                        res.partial = true;
                        break;
                    }
                }
                int d = 0;
                while (d < K) {
                    n[d]++;
                    if (n[d] <= ub[d]) break;
                    n[d] = 0;
                    d++;
                }
                if (d >= K) break;
            }
        } else {
            res.partial = true;
            for (int trial = 0; trial < 400; trial++) {
                int len = 1 + rng.nextInt(Math.max(1, Math.min(Lmax, 6)));
                int[] n = new int[K];
                for (int j = 0; j < len; j++) {
                    int r = rng.nextInt(K);
                    if (n[r] < ub[r]) n[r]++;
                }
                int sum = 0;
                for (int r = 0; r < K; r++) sum += n[r];
                if (sum >= 1 && sum <= Lmax - 1) {
                    stepPrefix(n, ub, K, mu, tol, neval, res);
                    if (!res.ok || neval[0] >= MAXEVAL) break;
                }
            }
        }
        return res;
    }

    /** No one-job extension of the prefix n may lower the total rate. */
    private static void stepPrefix(int[] n, int[] ub, int K,
            SerializableFunction<Matrix, Double> mu, double tol, int[] neval,
            MonoCheckResult res) {
        double base = mu.apply(rowMatrix(microstate(n, K)));
        neval[0]++;
        for (int r = 0; r < K; r++) {
            if (n[r] >= ub[r]) continue;
            n[r]++;
            int[] cx = microstate(n, K);
            double v = mu.apply(rowMatrix(cx));
            neval[0]++;
            n[r]--;
            if (v - base < -tol * Math.max(1.0, Math.abs(base))) {
                res.ok = false;
                res.badc = cx;
                res.badr = r;
                return;
            }
        }
    }

    /** Canonical microstate of a count vector: classes in index order. */
    private static int[] microstate(int[] n, int K) {
        int len = 0;
        for (int r = 0; r < K; r++) len += n[r];
        int[] c = new int[len];
        int col = 0;
        for (int r = 0; r < K; r++) {
            for (int t = 0; t < n[r]; t++) c[col++] = r;
        }
        return c;
    }

    private static void testMultiset(int[] n, int K, SerializableFunction<Matrix, Double> mu,
            double tol, int PERM_ENUM, int PERM_SAMPLE, java.util.Random rng, int[] neval,
            PermCheckResult res) {
        int len = 0;
        for (int r = 0; r < K; r++) len += n[r];
        int[] c0 = new int[len];
        int col = 0;
        for (int r = 0; r < K; r++) {
            for (int t = 0; t < n[r]; t++) c0[col++] = r;
        }
        double logd = lgammaPos(len + 1);
        for (int r = 0; r < K; r++) {
            if (n[r] > 0) logd -= lgammaPos(n[r] + 1);
        }
        long dcount = Math.round(Math.exp(logd));
        double base = mu.apply(rowMatrix(c0));
        neval[0]++;
        java.util.List<int[]> perms;
        if (dcount <= PERM_ENUM) {
            perms = msPerms(c0);
        } else {
            res.partial = true;
            perms = new java.util.ArrayList<int[]>();
            int[] rev = new int[len];
            for (int i = 0; i < len; i++) rev[i] = c0[len - 1 - i];
            perms.add(rev);
            for (int t = 1; t < PERM_SAMPLE; t++) {
                int[] p = c0.clone();
                for (int i = len - 1; i > 0; i--) {
                    int j = rng.nextInt(i + 1);
                    int tmp = p[i]; p[i] = p[j]; p[j] = tmp;
                }
                perms.add(p);
            }
        }
        for (int[] p : perms) {
            double v = mu.apply(rowMatrix(p));
            neval[0]++;
            if (Math.abs(v - base) > tol * Math.max(1.0, Math.abs(base))) {
                res.ok = false;
                res.badc = c0;
                return;
            }
        }
    }

    private static java.util.List<int[]> msPerms(int[] c0) {
        java.util.List<int[]> out = new java.util.ArrayList<int[]>();
        if (c0.length <= 1) {
            out.add(c0.clone());
            return out;
        }
        java.util.TreeSet<Integer> uniq = new java.util.TreeSet<Integer>();
        for (int x : c0) uniq.add(x);
        for (int u : uniq) {
            int[] rest = new int[c0.length - 1];
            int idx = 0;
            boolean removed = false;
            for (int x : c0) {
                if (!removed && x == u) {
                    removed = true;
                    continue;
                }
                rest[idx++] = x;
            }
            for (int[] sub : msPerms(rest)) {
                int[] perm = new int[c0.length];
                perm[0] = u;
                System.arraycopy(sub, 0, perm, 1, sub.length);
                out.add(perm);
            }
        }
        return out;
    }

    private static Matrix rowMatrix(int[] c) {
        Matrix m = new Matrix(1, c.length);
        for (int i = 0; i < c.length; i++) {
            m.set(0, i, c[i]);
        }
        return m;
    }

    private static double lgammaPos(double x) {
        // Lanczos approximation of ln Gamma(x), x >= 1.
        double[] g = {
            676.5203681218851, -1259.1392167224028, 771.32342877765313,
            -176.61502916214059, 12.507343278686905, -0.13857109526572012,
            9.9843695780195716e-6, 1.5056327351493116e-7
        };
        x -= 1;
        double a = 0.99999999999980993;
        double t = x + 7.5;
        for (int i = 0; i < g.length; i++) {
            a += g[i] / (x + i + 1);
        }
        return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
    }

    /**
     * Sets the class compatibility/swap graph of a pass-and-swap (PAS) queue.
     *
     * @param graph (nclasses x nclasses) adjacency matrix; entry (r,s) nonzero
     *              iff, upon completion of a class-r job, a waiting class-s job
     *              may take its place. Undirected; self-loops allowed.
     */
    public void setSwapGraph(Matrix graph) {
        if (this.schedStrategy == SchedStrategy.OI) {
            throw new RuntimeException("setSwapGraph is not applicable to OI (order-independent) queues; their swap graph is always zero. Use SchedStrategy.PAS to configure a non-trivial swap graph.");
        }
        if (this.schedStrategy != SchedStrategy.PAS) {
            throw new RuntimeException("setSwapGraph is only applicable to PAS (pass-and-swap) queues.");
        }
        int K = this.model.getNumberOfClasses();
        if (graph.getNumRows() != K || graph.getNumCols() != K) {
            throw new RuntimeException("Swap graph must be a " + K + "x" + K + " matrix (nclasses x nclasses).");
        }
        // swapGraph is copied into QueueNodeParam by refreshLocalVars.
        this.swapGraph = graph;
        invalidateStruct();
    }

    /**
     * Returns the (nclasses x nclasses) swap graph of a PAS queue, or null.
     *
     * @return the swap graph matrix
     */
    public Matrix getSwapGraph() {
        return this.swapGraph;
    }

    /**
     * Gets the scheduling strategy used by this queue.
     * 
     * @return The scheduling strategy (e.g., FCFS, PS, LCFS, etc.)
     */
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    /**
     * Gets the scheduling strategy parameter for a specific job class.
     * 
     * <p>For strategies like DPS (Discriminatory Processor Sharing) and GPS
     * (Generalized Processor Sharing), this returns the weight or priority
     * parameter for the job class.</p>
     * 
     * @param jobClass The job class to get the parameter for
     * @return The scheduling parameter value, or 0.0 if not set
     */
    public double getSchedStrategyPar(JobClass jobClass) {
        // For LPS, return the station-wide limit (stored in lpsLimit field)
        // This will be used to populate schedparam Matrix for all classes
        if (this.schedStrategy == SchedStrategy.LPS) {
            return this.lpsLimit;
        }
        return this.schedStrategyPar.getOrDefault(jobClass, 0.0);
    }

    /**
     * Gets the service time distribution for a specific job class.
     * 
     * @param jobClass The job class to get the service distribution for
     * @return The service time distribution for this job class
     */
    public Distribution getService(JobClass jobClass) {
        return this.server.getServiceDistribution(jobClass);
    }

    /**
     * Override, at this queue, the retrieval service rate for a single item of the read
     * class jobinClass in the given cache's retrieval system. The default (when not
     * overridden) is the read class's own service distribution at this queue. item is 0-based.
     *
     * @param cache       the cache whose retrieval system contains this queue
     * @param jobinClass  the read class that creates the retrieval request
     * @param item        the requested item (0-based)
     * @param serviceRate exponential service rate for the item at this queue
     */
    public void setItemServiceRate(Cache cache, JobClass jobinClass, int item, double serviceRate) {
        int jobinClassIdx = jobinClass.getIndex() - 1;
        int retrievalClassIdx = (int) cache.getRetrievalClasses().get(item, jobinClassIdx);
        if (retrievalClassIdx < 0) {
            throw new RuntimeException("No retrieval class defined for the given class/item; call setRetrievalSystem first.");
        }
        JobClass retrievalClass = this.model.getJobClassFromIndex(retrievalClassIdx);
        this.setService(retrievalClass, new Exp(serviceRate));
    }

    /**
     * Prints a summary of this queue's configuration to standard output.
     * 
     * <p>The summary includes the queue name, service processes for each job class,
     * their mean service times and squared coefficients of variation, number of servers,
     * and output routing configuration.</p>
     */
    @Override
    public void printSummary() {
        System.out.format("jline.Queue:\n");
        System.out.format("--Name: %s\n", this.getName());
        System.out.format("--Service Processes:\n");
        for (JobClass jobClass : this.model.getClasses()) {
            System.out.format("----%s: %s (Mean: %g, SCV: %g)\n", jobClass.getName(), this.getServiceProcess(jobClass).toString(), this.getServiceProcess(jobClass).getMean(), this.getServiceProcess(jobClass).getSCV());
        }
        if (isInf(this.getNumberOfServers())) {
            System.out.format("--Number of Servers: Inf\n");
        } else {
            System.out.format("--Number of Servers: %d\n", this.getNumberOfServers());
        }
        this.output.printSummary();
    }

    /**
     * Sets a class-dependent scaling function for service rates.
     * 
     * <p>The function takes a matrix representing the number of jobs of each class
     * at the station and returns the service rate scaling. It may return either a
     * 1x1 matrix, i.e. a chain-independent scaling beta_i(n) shared by every
     * class, or a length-R row vector [beta_{i,1}(n), ..., beta_{i,R}(n)] giving
     * each class its own rate at the same population, which expresses Sauer's
     * chain-dependent service rates mu_{r,i}(n) (Sauer 1983, eq. (40)). This
     * enables modeling of systems where service rates depend on the job mix.</p>
     *
     * @param beta A function that maps job class populations to a service rate scaling
     * @throws RuntimeException if the scheduling strategy doesn't support class dependence
     */
    public void setClassDependence(SerializableFunction<Matrix, Matrix> beta) {
        throw new RuntimeException("Class dependence requires an explicit peak rate: setClassDependence(beta, peakRatePerClass). Pass a 1x1 Matrix (identical peak for every class) or a 1xR per-class vector.");
    }

    /**
     * Sets a class-dependent service-rate scaling and its required peak rate per
     * class. Utilization at the station is reported as U = T*S/peak.
     *
     * @param beta maps the per-class population vector to a service-rate scaling
     * @param peakRatePerClass 1x1 (broadcast) or 1xR peak rate scaling per class
     */
    public void setClassDependence(SerializableFunction<Matrix, Matrix> beta, Matrix peakRatePerClass) {
        if (peakRatePerClass == null || peakRatePerClass.isEmpty()) {
            throw new RuntimeException("Class dependence requires an explicit positive peak rate (peakRatePerClass).");
        }
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                this.setLimitedClassDependence(beta, peakRatePerClass);
                break;
            default:
                throw new RuntimeException("Class-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");

        }
    }

    /**
     * Sets a joint-dependent (non-product-form) service-rate scaling and its
     * required peak rate per class. Unlike setClassDependence, whose handle must
     * express the product-form beta_{i,r}(n_{i,r}) depending on the own-class
     * marginal, eta(ni) may read the joint per-class population vector
     * arbitrarily (e.g. min(ni[0],c)) and is therefore non-product-form: solvers
     * treat it as an approximation with no exactness/uniqueness guarantee.
     * Utilization at the station is reported as U = T*S/peak.
     *
     * @param eta maps the joint per-class population vector to a service-rate scaling
     * @param peakRatePerClass 1x1 (broadcast) or 1xR peak rate scaling per class
     */
    public void setJointDependence(SerializableFunction<Matrix, Matrix> eta, Matrix peakRatePerClass) {
        if (peakRatePerClass == null || peakRatePerClass.isEmpty()) {
            throw new RuntimeException("Joint dependence requires an explicit positive peak rate (peakRatePerClass).");
        }
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                this.setLimitedJointDependence(eta, peakRatePerClass);
                break;
            default:
                throw new RuntimeException("Joint-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");

        }
    }

    /**
     * Sets load-dependent service rate scaling factors.
     *
     * <p>Each element alpha[n] specifies the service rate scaling when there are
     * n jobs at the station. This enables modeling of systems where performance
     * degrades under load.</p>
     *
     * @param alpha A matrix of scaling factors indexed by the number of jobs
     * @throws RuntimeException if the scheduling strategy doesn't support load dependence
     */
    public void setLoadDependence(Matrix alpha) {
        switch (this.schedStrategy) {
            case PS:
            case FCFS:
                this.setLimitedLoadDependence(alpha);
                break;
            default:
                throw new RuntimeException("Load-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) stations.");
        }
    }

    /**
     * Sets the number of servers at this queueing station.
     * 
     * <p>For infinite server (IS) queues, this method has no effect as they
     * always have unlimited servers.</p>
     * 
     * @param numberOfServers The number of parallel servers (must be positive)
     */
    public void setNumberOfServers(int numberOfServers) {
        // Match MATLAB Queue.setNumServers: DPS/GPS do not admit multi-server stations.
        if ((this.schedStrategy == SchedStrategy.DPS || this.schedStrategy == SchedStrategy.GPS)
                && numberOfServers != 1) {
            throw new RuntimeException(String.format(
                    "Cannot use multi-server stations with %s scheduling.",
                    SchedStrategy.toText(this.schedStrategy)));
        }
        if (this.schedStrategy != SchedStrategy.INF) {
            this.numberOfServers = numberOfServers;
            // Station.setNumberOfServers invalidates here and this override did not, so a
            // server count changed after anything had built the struct (a gate call, a
            // getStruct, a solve) left sn.nservers at the OLD value: the model disagreed
            // with itself, the feature recorder reading the node objects and seeing the new
            // count while every sn-based structural predicate saw the old one. Same defect
            // fixed in MATLAB Queue.setNumServers on 2026-09-05.
            invalidateStruct();
        }
    }

    /**
     * Makes the server of this station subject to breakdowns. The server alternates
     * between an UP and a DOWN status: while up it fails after failureDistribution,
     * while down it is restored after repairDistribution.
     *
     * <p>The failure clock runs whenever the server is up, whether or not a job is
     * in service, so a station can fail while idle. Arrivals are unaffected by the
     * server status and keep queueing (subject to the station capacity) while the
     * server is down. A job that is in service when the server fails is not lost:
     * it stays at the station and, service being memoryless in the supported case,
     * resumes when the server is repaired.</p>
     *
     * <p>Only exponential failure and repair distributions are currently expanded
     * into the joint (queue, server status) chain; anything else is rejected here
     * rather than silently approximated.</p>
     *
     * @param failureDistribution     time to failure of an up server
     * @param repairDistribution      repair time of a down server
     * @param downServiceDistribution service distribution used while the server is
     *                                down, applied to every class; null means the
     *                                server does not serve at all while down, which
     *                                is the usual breakdown model
     */
    public void setBreakdown(Distribution failureDistribution, Distribution repairDistribution,
                             Distribution downServiceDistribution) {
        if (failureDistribution == null || repairDistribution == null) {
            throw new IllegalArgumentException("setBreakdown requires a failure and a repair Distribution.");
        }
        if (!(failureDistribution instanceof Exp) || !(repairDistribution instanceof Exp)) {
            throw new RuntimeException("Station '" + this.getName() + "': only exponential failure and repair "
                    + "distributions are supported by setBreakdown. A non-exponential failure or repair process "
                    + "needs its own phase in the joint chain, which is not implemented; use an Environment "
                    + "ensemble for that case.");
        }
        if (failureDistribution.getMean() <= 0 || repairDistribution.getMean() <= 0) {
            throw new IllegalArgumentException(
                    "setBreakdown requires strictly positive failure and repair means.");
        }
        this.breakdownFailure = failureDistribution;
        this.breakdownRepair = repairDistribution;
        this.breakdownDownService.clear();
        this.breakdownDownServiceAll = downServiceDistribution;
        // refreshBreakdown derives sn.hasbreakdown, breakdownMu, repairMu and
        // downServiceRates from these, so a cached struct would answer the unbroken model.
        invalidateStruct();
    }

    /**
     * Makes the server of this station subject to breakdowns, with no service at
     * all while the server is down.
     *
     * @param failureDistribution time to failure of an up server
     * @param repairDistribution  repair time of a down server
     */
    public void setBreakdown(Distribution failureDistribution, Distribution repairDistribution) {
        setBreakdown(failureDistribution, repairDistribution, null);
    }

    /**
     * Sets the degraded service distribution used for one class while the server is
     * down. Overrides the class-independent distribution given to setBreakdown.
     *
     * @param jobClass                the job class
     * @param downServiceDistribution the degraded service distribution
     */
    public void setDownService(JobClass jobClass, Distribution downServiceDistribution) {
        this.breakdownDownService.put(jobClass, downServiceDistribution);
        // read back by refreshBreakdown into sn.downServiceRates
        invalidateStruct();
    }

    /**
     * Time to failure of the server, or null when the server never fails.
     *
     * @return the failure distribution
     */
    public Distribution getBreakdownFailure() {
        return this.breakdownFailure;
    }

    /**
     * Repair time of a down server, or null when the server never fails.
     *
     * @return the repair distribution
     */
    public Distribution getBreakdownRepair() {
        return this.breakdownRepair;
    }

    /**
     * Service distribution used for a class while the server is down, or null when
     * the server does not serve at all while down.
     *
     * @param jobClass the job class
     * @return the degraded service distribution
     */
    public Distribution getDownService(JobClass jobClass) {
        Distribution d = this.breakdownDownService.get(jobClass);
        return (d != null) ? d : this.breakdownDownServiceAll;
    }

    /**
     * Whether this station is subject to server breakdowns.
     *
     * @return true iff both a failure and a repair distribution are configured
     */
    public boolean hasBreakdown() {
        return this.breakdownFailure != null && this.breakdownRepair != null;
    }

    /**
     * Sets the scheduling strategy parameter for a specific job class.
     * 
     * <p>For weighted scheduling strategies (DPS, GPS, etc.), this sets the
     * weight or priority parameter that determines the job class's share of
     * the service capacity.</p>
     * 
     * @param jobClass The job class to set the parameter for
     * @param weight The scheduling parameter value (e.g., weight, priority)
     */
    public void setSchedStrategyPar(JobClass jobClass, double weight) {
        this.schedStrategyPar.put(jobClass, weight);
    }

    /**
     * Sets the maximum number of jobs for LPS scheduling.
     * For LPS: limit is the max number of jobs that can execute in PS mode.
     *
     * @param limit the maximum number of concurrent jobs in PS mode
     * @throws RuntimeException if called on non-LPS queue
     */
    public void setLimit(int limit) {
        if (this.schedStrategy != SchedStrategy.LPS) {
            throw new RuntimeException("setLimit() can only be called on queues with LPS scheduling strategy");
        }
        // Store limit as station-wide parameter for LPS. getSchedStrategyPar returns it
        // for LPS and refreshScheduling folds it into sn.schedparam.
        this.lpsLimit = (double) limit;
        invalidateStruct();
    }

    /**
     * Sets the polling type for this queue (only valid for POLLING scheduling strategy).
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     */
    public void setPollingType(PollingType pollingType) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingType() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingType(pollingType);
        }
        invalidateStructAfterPollingChange();
    }

    /**
     * Sets the polling type for this queue with K value for K-LIMITED (only valid for POLLING scheduling strategy).
     *
     * @param pollingType the polling type (GATED, EXHAUSTIVE, or KLIMITED)
     * @param k the K value for K-LIMITED polling (ignored for other types)
     */
    public void setPollingType(PollingType pollingType, int k) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingType() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingType(pollingType, k);
        }
        invalidateStructAfterPollingChange();
    }

    /**
     * Sets the K value for K-LIMITED polling (only valid for POLLING scheduling strategy with K-LIMITED type).
     *
     * @param k the K value (must be greater than 0)
     */
    public void setPollingK(int k) {
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setPollingK() can only be called on queues with POLLING scheduling strategy");
        }
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setPollingK(k);
        }
        invalidateStructAfterPollingChange();
    }

    /**
     * Sets the switchover time for a job class (only valid for POLLING scheduling strategy).
     *
     * @param jobClass the job class
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchover(JobClass jobClass, Distribution switchoverTime) {
        // Validate input parameters
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (switchoverTime == null) {
            throw new IllegalArgumentException("switchoverTime cannot be null");
        }
        
        // Check scheduling strategy
        if (this.schedStrategy != SchedStrategy.POLLING) {
            throw new RuntimeException("setSwitchover() can only be called on queues with POLLING scheduling strategy");
        }
        
        // Check if job class is valid for this network
        if (!this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }
        
        if (this.server instanceof PollingServer) {
            ((PollingServer) this.server).setSwitchover(jobClass, switchoverTime);
        }
        invalidateStructAfterPollingChange();
    }

    /**
     * Invalidates any cached network struct after a polling parameter change.
     * <p>
     * The polling type, K value and switchover times are all stored on the
     * {@link PollingServer} (the source of truth) and are re-read into the
     * {@link QueueNodeParam} whenever {@code refreshStruct} rebuilds the struct.
     * This method therefore only needs to drop a stale cached struct so it is
     * rebuilt on next access. It deliberately does NOT call {@code getStruct()}:
     * during model construction these setters are frequently invoked before the
     * topology (links and routing) is defined, and forcing a premature build
     * would cache an empty connection matrix, which later breaks chain
     * construction and every solver that relies on it.
     */
    private void invalidateStructAfterPollingChange() {
        if (this.model != null && this.model.getHasStruct()) {
            this.model.resetStruct();
        }
    }

    /**
     * Sets the switchover time from one job class to another (for general scheduling strategies).
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @param switchoverTime the switchover time distribution
     */
    public void setSwitchover(JobClass fromClass, JobClass toClass, Distribution switchoverTime) {
        // Validate input parameters
        if (fromClass == null) {
            throw new IllegalArgumentException("fromClass cannot be null");
        }
        if (toClass == null) {
            throw new IllegalArgumentException("toClass cannot be null");
        }
        if (switchoverTime == null) {
            throw new IllegalArgumentException("switchoverTime cannot be null");
        }
        
        // Check if job classes are valid for this network
        if (!this.model.getClasses().contains(fromClass)) {
            throw new IllegalArgumentException("fromClass is not part of this network");
        }
        if (!this.model.getClasses().contains(toClass)) {
            throw new IllegalArgumentException("toClass is not part of this network");
        }
        
        // For POLLING strategy, delegate to PollingServer but still store in general storage
        if (this.schedStrategy == SchedStrategy.POLLING) {
            if (this.server instanceof PollingServer) {
                // For polling, set the switchover time for the fromClass (ignoring toClass)
                ((PollingServer) this.server).setSwitchover(fromClass, switchoverTime);
            }
        }
        
        // Store in general switchover time storage for all scheduling strategies
        this.setSwitchoverTime(fromClass, toClass, switchoverTime);
        // The POLLING arm above wrote the PollingServer, whose switchover times
        // refreshLocalVars copies into QueueNodeParam; the two-argument setSwitchover
        // already invalidates for exactly that reason.
        invalidateStruct();
    }

    /**
     * Gets the switchover time from one job class to another.
     *
     * @param fromClass the job class to switch from
     * @param toClass the job class to switch to
     * @return the switchover time distribution, or null if not set
     */
    public Distribution getSwitchover(JobClass fromClass, JobClass toClass) {
        return this.getSwitchoverTime(fromClass, toClass);
    }

    /**
     * Gets the switchover time for a job class (POLLING scheduling strategy).
     *
     * @param jobClass the job class
     * @return the switchover time distribution, or null if not set
     */
    public Distribution getSwitchover(JobClass jobClass) {
        if (this.schedStrategy == SchedStrategy.POLLING && this.server instanceof PollingServer) {
            return ((PollingServer) this.server).getSwitchover(jobClass);
        }
        return null;
    }
    
    /**
     * Sets the setup time and delay off time for a job class.
     * This is typically used for function-based tasks that have initialization overhead.
     *
     * @param jobClass the job class
     * @param setupTime the setup time distribution
     * @param delayoffTime the delay off time distribution
     */
    public void setDelayOff(JobClass jobClass, Distribution setupTime, Distribution delayoffTime) {
        // Validate input parameters
        if (jobClass == null) {
            throw new IllegalArgumentException("jobClass cannot be null");
        }
        if (setupTime == null) {
            throw new IllegalArgumentException("setupTime cannot be null");
        }
        if (delayoffTime == null) {
            throw new IllegalArgumentException("delayoffTime cannot be null");
        }
        
        // Check if job class is valid for this network
        if (!this.model.getClasses().contains(jobClass)) {
            throw new IllegalArgumentException("jobClass is not part of this network");
        }
        
        // The values live on the node; refreshing the struct here would sanitize a
        // model that is still being built, e.g. while LineModelIO loads it. Dropping a
        // stale one is not refreshing it: invalidateStruct never calls getStruct, and
        // refreshLocalVars copies these into QueueNodeParam.setupTime/delayoffTime.
        this.setupTimes.put(jobClass, setupTime);
        this.delayOffTimes.put(jobClass, delayoffTime);
        invalidateStruct();
    }
    
    /**
     * Gets the setup time distribution for a job class.
     *
     * @param jobClass the job class
     * @return the setup time distribution, or null if not set
     */
    public Distribution getSetupTime(JobClass jobClass) {
        return this.setupTimes.get(jobClass);
    }
    
    /**
     * Gets the delay-off time distribution for a job class.
     *
     * @param jobClass the job class
     * @return the delay-off time distribution, or null if not set
     */
    public Distribution getDelayOffTime(JobClass jobClass) {
        return this.delayOffTimes.get(jobClass);
    }

    /**
     * Checks if this queue has delay-off times enabled.
     * Delay-off is considered enabled if any job class has both
     * a setup time and a delay-off time distribution configured.
     *
     * @return true if delay-off is enabled, false otherwise
     */
    public boolean isDelayOffEnabled() {
        return !this.setupTimes.isEmpty() && !this.delayOffTimes.isEmpty();
    }

    // ==================== Heterogeneous Server Methods ====================

    /**
     * Adds a server type to this queue for heterogeneous multiserver configuration.
     * <p>
     * When server types are added, the queue becomes a heterogeneous multiserver queue
     * where different server types can have different service rates and serve different
     * subsets of job classes.
     * <p>
     * The total number of servers at this queue becomes the sum of all server type counts.
     *
     * @param serverType the server type to add
     * @throws IllegalArgumentException if serverType is null or already added
     */
    public void addServerType(ServerType serverType) {
        if (serverType == null) {
            throw new IllegalArgumentException("Server type cannot be null");
        }
        if (this.serverTypes.contains(serverType)) {
            throw new IllegalArgumentException("Server type '" + serverType.getName() + "' is already added to this queue");
        }

        // Assign ID and parent
        serverType.setId(this.serverTypes.size());
        serverType.setParentQueue(this);
        this.serverTypes.add(serverType);

        // Initialize service distribution map for this server type
        this.heteroServiceDistributions.put(serverType, new HashMap<JobClass, Distribution>());

        // Update total number of servers
        updateTotalServerCount();
        // This moves numberOfServers as well as the hetero server tables that
        // refreshHeterogeneousServers reads, so any cached struct is stale.
        invalidateStruct();
    }

    /**
     * Updates the total numberOfServers based on all server types.
     */
    private void updateTotalServerCount() {
        if (this.serverTypes.isEmpty()) {
            return;
        }
        int total = 0;
        for (ServerType st : this.serverTypes) {
            total += st.getNumOfServers();
        }
        this.numberOfServers = total;
    }

    /**
     * Gets the list of server types configured for this queue.
     *
     * @return a new list containing the server types
     */
    public List<ServerType> getServerTypes() {
        return new ArrayList<ServerType>(this.serverTypes);
    }

    /**
     * Sets the number of servers that a job of the given class seizes for the whole
     * of its service, JMT's job parallelism ({@code Server.serverNumRequired}).
     * <p>
     * A job waits until n servers are simultaneously free and holds all of them
     * until it completes, so the station serves at most floor(c/n) such jobs at a
     * time. The default is 1.
     *
     * @param jobClass the job class
     * @param n the number of servers required, an integer in [1, c]
     * @throws IllegalArgumentException if n is below 1 or above the server count
     */
    public void setServerParallelism(JobClass jobClass, int n) {
        if (n < 1) {
            throw new IllegalArgumentException("Server parallelism must be a positive integer.");
        }
        if (n > this.numberOfServers) {
            throw new IllegalArgumentException(String.format(
                    "Server parallelism (%d) exceeds the %d servers of station %s, so a job of class %s could never enter service.",
                    n, this.numberOfServers, this.getName(), jobClass.getName()));
        }
        this.serverParallelism.put(jobClass, n);
        // read back by refreshHeterogeneousServers into the struct's parallelism vector
        invalidateStruct();
    }

    /**
     * Gets the number of servers seized by a job of the given class.
     *
     * @param jobClass the job class
     * @return the number of servers required, 1 if unset
     */
    public int getServerParallelism(JobClass jobClass) {
        Integer n = this.serverParallelism.get(jobClass);
        return n == null ? 1 : n;
    }

    /**
     * Checks whether some class seizes more than one server.
     *
     * @return true if job parallelism is declared for at least one class
     */
    public boolean hasServerParallelism() {
        for (Integer n : this.serverParallelism.values()) {
            if (n != null && n > 1) {
                return true;
            }
        }
        return false;
    }

    /**
     * Gets the number of server types configured for this queue.
     *
     * @return the number of server types, or 0 if homogeneous
     */
    public int getNumServerTypes() {
        return this.serverTypes.size();
    }

    /**
     * Checks if this queue is configured as a heterogeneous multiserver queue.
     *
     * @return true if server types are defined, false for homogeneous queue
     */
    public boolean isHeterogeneous() {
        return !this.serverTypes.isEmpty();
    }

    /**
     * Sets the scheduling policy for heterogeneous servers.
     * <p>
     * This policy determines how jobs are assigned to server types when a job's
     * class is compatible with multiple server types.
     *
     * @param policy the heterogeneous scheduling policy
     * @throws IllegalArgumentException if policy is null
     */
    public void setHeteroSchedPolicy(HeteroSchedPolicy policy) {
        if (policy == null) {
            throw new IllegalArgumentException("Heterogeneous scheduling policy cannot be null");
        }
        this.heteroSchedPolicy = policy;
        // read back by refreshHeterogeneousServers into snp.heteroschedpolicy
        invalidateStruct();
    }

    /**
     * Gets the scheduling policy for heterogeneous servers.
     *
     * @return the heterogeneous scheduling policy
     */
    public HeteroSchedPolicy getHeteroSchedPolicy() {
        return this.heteroSchedPolicy;
    }

    /**
     * Sets the service time distribution for a specific job class and server type.
     * <p>
     * This method is used for heterogeneous multiserver queues where different
     * server types may have different service rates for the same job class.
     *
     * @param jobClass the job class
     * @param serverType the server type
     * @param distribution the service time distribution
     * @throws IllegalArgumentException if any parameter is null or serverType is not in this queue
     */
    public void setService(JobClass jobClass, ServerType serverType, Distribution distribution) {
        if (jobClass == null) {
            throw new IllegalArgumentException("Job class cannot be null");
        }
        if (serverType == null) {
            throw new IllegalArgumentException("Server type cannot be null");
        }
        if (distribution == null) {
            throw new IllegalArgumentException("Distribution cannot be null");
        }
        if (!this.serverTypes.contains(serverType)) {
            throw new IllegalArgumentException("Server type '" + serverType.getName() +
                "' is not added to this queue. Call addServerType() first.");
        }

        // Store the heterogeneous service distribution
        Map<JobClass, Distribution> classMap = this.heteroServiceDistributions.get(serverType);
        if (classMap == null) {
            classMap = new HashMap<JobClass, Distribution>();
            this.heteroServiceDistributions.put(serverType, classMap);
        }
        classMap.put(jobClass, distribution);

        // Also ensure the job class is marked as compatible with this server type
        if (!serverType.isCompatible(jobClass)) {
            serverType.addCompatibleClass(jobClass);
        }
        // Unlike the homogeneous ServiceStation.setService, which the in-place
        // refreshRates path maintains, the hetero tables are only read by
        // refreshHeterogeneousServers when the whole struct is rebuilt.
        invalidateStruct();
    }

    /**
     * Gets the service time distribution for a specific job class and server type.
     *
     * @param jobClass the job class
     * @param serverType the server type
     * @return the service time distribution, or null if not set
     */
    public Distribution getService(JobClass jobClass, ServerType serverType) {
        if (serverType == null || !this.heteroServiceDistributions.containsKey(serverType)) {
            return null;
        }
        Map<JobClass, Distribution> classMap = this.heteroServiceDistributions.get(serverType);
        return classMap != null ? classMap.get(jobClass) : null;
    }

    /**
     * Gets all heterogeneous service distributions for this queue.
     *
     * @return the map of server type to job class to distribution
     */
    public Map<ServerType, Map<JobClass, Distribution>> getHeteroServiceDistributions() {
        return this.heteroServiceDistributions;
    }

    /**
     * Gets a server type by its ID.
     *
     * @param id the server type ID
     * @return the server type, or null if not found
     */
    public ServerType getServerType(int id) {
        if (id >= 0 && id < this.serverTypes.size()) {
            return this.serverTypes.get(id);
        }
        return null;
    }

    /**
     * Gets a server type by its name.
     *
     * @param name the server type name
     * @return the server type, or null if not found
     */
    public ServerType getServerType(String name) {
        for (ServerType st : this.serverTypes) {
            if (st.getName().equals(name)) {
                return st;
            }
        }
        return null;
    }

    /**
     * Checks if all job classes have at least one compatible server type.
     * <p>
     * This validation is important because in heterogeneous queues, every job class
     * must be able to be served by at least one server type.
     *
     * @return true if all job classes have at least one compatible server type
     */
    public boolean validateCompatibility() {
        if (!isHeterogeneous()) {
            return true; // Homogeneous queues are always valid
        }

        for (JobClass jobClass : this.model.getClasses()) {
            boolean hasCompatible = false;
            for (ServerType st : this.serverTypes) {
                if (st.isCompatible(jobClass)) {
                    hasCompatible = true;
                    break;
                }
            }
            if (!hasCompatible) {
                return false;
            }
        }
        return true;
    }

    // ==================== Immediate Feedback Methods ====================

    /**
     * Enables or disables immediate feedback for all job classes at this queue.
     * When enabled, jobs that self-loop at this station stay in service instead
     * of rejoining the queue.
     *
     * @param enabled true to enable immediate feedback for all classes, false to disable
     */
    public void setImmediateFeedback(boolean enabled) {
        this.immediateFeedbackAll = enabled;
        if (!enabled) {
            this.immediateFeedbackClasses.clear();
        }
        // refreshStruct folds hasImmediateFeedback(c) into sn.immfeed
        invalidateStruct();
    }

    /**
     * Enables immediate feedback for a specific job class at this queue.
     *
     * @param jobClass the job class to enable immediate feedback for
     */
    public void setImmediateFeedback(JobClass jobClass) {
        if (jobClass != null) {
            this.immediateFeedbackClasses.add(jobClass.getIndex());
            // refreshStruct folds hasImmediateFeedback(c) into sn.immfeed
            invalidateStruct();
        }
    }

    /**
     * Enables immediate feedback for multiple job classes at this queue.
     *
     * @param jobClasses list of job classes to enable immediate feedback for
     */
    public void setImmediateFeedbackForClasses(java.util.List<JobClass> jobClasses) {
        if (jobClasses != null) {
            for (JobClass jc : jobClasses) {
                if (jc != null) {
                    this.immediateFeedbackClasses.add(jc.getIndex());
                }
            }
            // refreshStruct folds hasImmediateFeedback(c) into sn.immfeed
            invalidateStruct();
        }
    }

    /**
     * Checks if immediate feedback is enabled for any class at this queue.
     *
     * @return true if immediate feedback is enabled for at least one class
     */
    public boolean hasImmediateFeedback() {
        return this.immediateFeedbackAll || !this.immediateFeedbackClasses.isEmpty();
    }

    /**
     * Checks if immediate feedback is enabled for a specific class at this queue.
     *
     * @param classId the class index (0-based)
     * @return true if immediate feedback is enabled for the specified class
     */
    public boolean hasImmediateFeedback(int classId) {
        return this.immediateFeedbackAll || this.immediateFeedbackClasses.contains(classId);
    }

    /**
     * Gets the set of class indices with immediate feedback enabled.
     *
     * @return set of class indices, or null if none are set
     */
    public java.util.Set<Integer> getImmediateFeedbackClasses() {
        if (this.immediateFeedbackAll) {
            return null; // Indicates "all"
        }
        return this.immediateFeedbackClasses;
    }

    /**
     * Checks if immediate feedback is enabled for all classes.
     *
     * @return true if immediate feedback is enabled for all classes
     */
    public boolean isImmediateFeedbackAll() {
        return this.immediateFeedbackAll;
    }

}

