package jline.solvers.ctmc;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;
import java.util.Map;
import java.util.Properties;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.Maths;

import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.sparse.FillReducing;
import org.ejml.sparse.csc.factory.DecompositionFactory_DSCC;
import org.ejml.interfaces.decomposition.LUSparseDecomposition_F64;

/**
 * Hardware-aware, profiling-calibrated memory guard for {@link SolverCTMC}.
 *
 * Replaces the historical hard-coded state-space size threshold with a model
 * that (a) probes the memory actually available to this process and (b)
 * calibrates the per-state cost of a sparse steady-state solve by profiling
 * EJML's sparse LU once, caching the fitted power law per machine in the
 * temp directory. Mirrors the MATLAB (ctmc_memory_gate.m / lineGetAvailableMemory.m)
 * and Python-native (memory_guard.py) implementations, sharing the same model.
 *
 * Java 8 compatible: no var, no List.of, no switch expressions.
 */
public final class MemoryGuard {

    private MemoryGuard() {}

    /** Bytes charged per stored factor nonzero (8 value + 8 index, amortized). */
    private static final double BYTES_PER_NZ = 16.0;
    /** Default fraction of available memory the solver may target. */
    public static final double DEFAULT_SAFETY_FRACTION = 0.6;
    /** Lattice sides (n = g*g states) used to fit the fill power law. */
    private static final int CALIB_G1 = 40;
    private static final int CALIB_G2 = 80;
    /** Fallback power law if calibration cannot run. */
    private static final double FALLBACK_ALPHA = BYTES_PER_NZ * 8.0;
    private static final double FALLBACK_BETA = 1.3;
    /** Conservative available memory (bytes) when the probe fails. */
    private static final double FALLBACK_AVAIL_BYTES = 1.0 * 1024 * 1024 * 1024;

    /** Result of a gate decision. */
    public static final class GateResult {
        public final boolean ok;
        public final String message;
        public GateResult(boolean ok, String message) {
            this.ok = ok;
            this.message = message;
        }
    }

    /** Calibrated power-law coefficients. */
    private static final class Calibration {
        double alphaMem;
        double betaMem;
        double alphaT;
        double betaT;
        String sig;
    }

    // in-JVM memo so repeated solves in one session skip disk I/O
    private static Calibration cachedCalibration = null;

    /**
     * Portable memory budget probe (bytes). Because CTMC data is allocated on
     * the JVM heap, the binding constraint is the heap headroom, further capped
     * by free physical memory (queried via the OperatingSystemMXBean, which is
     * implemented on all platforms). Never throws; returns a conservative
     * constant if nothing can be determined.
     */
    public static double getAvailableMemoryBytes() {
        Runtime rt = Runtime.getRuntime();
        double heapHeadroom = (double) rt.maxMemory()
                - (double) (rt.totalMemory() - rt.freeMemory());
        if (heapHeadroom <= 0) {
            heapHeadroom = FALLBACK_AVAIL_BYTES;
        }
        double freePhysical = queryFreePhysicalBytes();
        if (freePhysical > 0) {
            return Math.min(heapHeadroom, freePhysical);
        }
        return heapHeadroom;
    }

    private static double queryFreePhysicalBytes() {
        try {
            java.lang.management.OperatingSystemMXBean bean =
                    ManagementFactory.getOperatingSystemMXBean();
            // see _kb/06-solver-catalog.md for rationale
            String[] candidates = {"getFreeMemorySize", "getFreePhysicalMemorySize"};
            for (int i = 0; i < candidates.length; i++) {
                try {
                    Method m = bean.getClass().getMethod(candidates[i]);
                    m.setAccessible(true);
                    Object v = m.invoke(bean);
                    if (v instanceof Long) {
                        long lv = ((Long) v).longValue();
                        if (lv > 0) {
                            return (double) lv;
                        }
                    }
                } catch (Throwable ignore) {
                    // try next candidate
                }
            }
        } catch (Throwable ignore) {
            // fall through
        }
        return -1;
    }

    private static String machineSignature() {
        String arch = System.getProperty("os.arch", "unknown");
        String os = System.getProperty("os.name", "unknown");
        int cores = Runtime.getRuntime().availableProcessors();
        return arch + "|" + cores + "|" + os;
    }

    private static File cacheFile() {
        String tmp = System.getProperty("java.io.tmpdir", ".");
        return new File(tmp, "line_ctmc_calib_java.properties");
    }

    /**
     * Worst-case log-size of the CTMC state space induced by sn: stars-and-bars
     * job placements per class over the stations that keep no ordered buffer
     * (open classes truncated at the cutoff), times the class-sequence
     * multiplicity of every order-preserving buffer, times the service-phase
     * multiplicity at each station, times one routing pointer per (node,class)
     * doing RROBIN or WRROBIN. Mirrors MATLAB ctmc_state_space_logsize.m and
     * Python state_space_log_size.
     *
     * @param sn      the network structure, after PH conversion
     * @param options solver options carrying the cutoff
     * @return the natural log of the worst-case number of states
     */
    public static double stateSpaceLogSize(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix NK = sn.njobs;

        double cutoffScalar = 0;
        Matrix cutoffMat = options.getCutoffMatrix(M, K);
        for (int ci = 0; ci < cutoffMat.getNumRows(); ci++) {
            for (int cj = 0; cj < cutoffMat.getNumCols(); cj++) {
                double cv = cutoffMat.get(ci, cj);
                if (!Double.isInfinite(cv) && cv > cutoffScalar) {
                    cutoffScalar = cv;
                }
            }
        }
        if (cutoffScalar <= 0) {
            // Same default the CTMC analyzer installs for open/mixed models.
            cutoffScalar = Math.ceil(Math.pow(6000, 1.0 / (M * K)));
        }

        // ORDERED BUFFERS, computed EXACTLY. A station outside the share family
        // keeps the class SEQUENCE of the jobs it holds. Omitting the factor
        // priced gallery_mmap1_multiclass at 726 states against 225840
        // enumerated; BOUNDING it instead of computing it double counts (the
        // sequence term already places the jobs it orders) and refused a working
        // model. The joint sum counts placement and ordering together. A
        // G-network signal never occupies a buffer slot.
        boolean[] isBuffered = new boolean[K];
        int Kb = 0;
        for (int k = 0; k < K; k++) {
            boolean sig = sn.issignal != null && k < sn.issignal.getNumElements()
                    && sn.issignal.get(k) != 0;
            isBuffered[k] = !sig;
            if (isBuffered[k]) {
                Kb++;
            }
        }
        int nOrd = 0;
        if (Kb > 1) {
            for (int i = 0; i < M; i++) {
                SchedStrategy schedI = null;
                if (sn.stations != null && i < sn.stations.size() && sn.sched != null) {
                    schedI = sn.sched.get(sn.stations.get(i));
                }
                if (schedI == SchedStrategy.EXT || isShareScheduling(schedI)) {
                    continue;
                }
                nOrd++;
            }
        }

        double logNstates = 0;
        double[] nkEff = new double[K];
        for (int k = 0; k < K; k++) {
            nkEff[k] = Double.isInfinite(NK.get(k)) ? cutoffScalar : NK.get(k);
        }

        if (nOrd == 0) {
            for (int k = 0; k < K; k++) {
                double mk = Math.max(admittingStations(sn, k, 0, M), 1);
                logNstates += Maths.factln(nkEff[k] + mk - 1) - Maths.factln(mk - 1)
                        - Maths.factln(nkEff[k]);
            }
        } else {
            int mRem = M - nOrd;
            for (int k = 0; k < K; k++) {
                if (!isBuffered[k]) {
                    double mk = Math.max(admittingStations(sn, k, 0, M), 1);
                    logNstates += Maths.factln(nkEff[k] + mk - 1) - Maths.factln(mk - 1)
                            - Maths.factln(nkEff[k]);
                }
            }
            int[] caps = new int[Kb];
            int ci = 0;
            double grid = 1.0;
            for (int k = 0; k < K; k++) {
                if (isBuffered[k]) {
                    caps[ci++] = (int) Math.floor(nkEff[k]);
                    grid *= (Math.floor(nkEff[k]) + 1);
                }
            }
            if (grid <= ORDER_GRID_MAX) {
                java.util.List<int[]> capsList = new java.util.ArrayList<int[]>();
                java.util.List<Double> capTotList = new java.util.ArrayList<Double>();
                for (int i2 = 0; i2 < M; i2++) {
                    SchedStrategy sc = null;
                    if (sn.stations != null && i2 < sn.stations.size() && sn.sched != null) {
                        sc = sn.sched.get(sn.stations.get(i2));
                    }
                    if (sc == SchedStrategy.EXT || isShareScheduling(sc)) {
                        continue;
                    }
                    int[] per = new int[caps.length];
                    int bk = 0;
                    for (int k = 0; k < K && bk < caps.length; k++) {
                        if (!isBuffered[k]) {
                            continue;
                        }
                        int c = caps[bk];
                        if (sn.classcap != null && i2 < sn.classcap.getNumRows()
                                && k < sn.classcap.getNumCols()
                                && Double.isFinite(sn.classcap.get(i2, k))) {
                            c = Math.min(c, (int) Math.floor(sn.classcap.get(i2, k)));
                        }
                        per[bk++] = c;
                    }
                    capsList.add(per);
                    double ct = Double.POSITIVE_INFINITY;
                    if (sn.cap != null && i2 < sn.cap.getNumElements()
                            && Double.isFinite(sn.cap.get(i2)) && sn.cap.get(i2) >= 0) {
                        ct = Math.floor(sn.cap.get(i2));
                    }
                    capTotList.add(ct);
                }
                int[][] capsPer = new int[capsList.size()][];
                double[] capTot = new double[capTotList.size()];
                for (int a2 = 0; a2 < capsList.size(); a2++) {
                    capsPer[a2] = capsList.get(a2);
                    capTot[a2] = capTotList.get(a2);
                }
                logNstates += logOrderedJoint(capsPer, capTot, caps, mRem);
            } else {
                double total = 0;
                for (int c : caps) {
                    total += c;
                }
                double logKb = Math.log(Kb);
                logNstates += nOrd * ((total + 1) * logKb - Math.log(Kb - 1)
                        + Math.log1p(-Math.exp(-(total + 1) * logKb)));
                for (int k = 0; k < K; k++) {
                    if (isBuffered[k]) {
                        double mk = admittingRemaining(sn, k, M);
                        if (mk >= 1) {
                            logNstates += Maths.factln(nkEff[k] + mk - 1)
                                    - Maths.factln(mk - 1) - Maths.factln(nkEff[k]);
                        }
                    }
                }
            }
        }

        if (sn.phasessz != null && !sn.phasessz.isEmpty()) {
            for (int i = 0; i < Math.min(M, sn.phasessz.getNumRows()); i++) {
                SchedStrategy sched = null;
                if (sn.stations != null && i < sn.stations.size() && sn.sched != null) {
                    sched = sn.sched.get(sn.stations.get(i));
                }
                for (int k = 0; k < Math.min(K, sn.phasessz.getNumCols()); k++) {
                    double p = sn.phasessz.get(i, k);
                    if (!Double.isFinite(p) || p <= 1) {
                        continue;
                    }
                    double m;
                    if (sched == SchedStrategy.EXT) {
                        m = 1;
                    } else if (isShareScheduling(sched)) {
                        m = nkEff[k];
                    } else {
                        // nservers may be stored as a row or a column: index linearly.
                        double c = (sn.nservers != null && i < sn.nservers.getNumElements())
                                ? sn.nservers.get(i) : 1;
                        m = Math.min(nkEff[k], c);
                    }
                    if (!Double.isFinite(m)) {
                        m = nkEff[k];
                    }
                    logNstates += Maths.factln(m + p - 1) - Maths.factln(p - 1) - Maths.factln(m);
                }
            }
        }

        if (sn.routing != null && sn.connmatrix != null && !sn.connmatrix.isEmpty()) {
            for (int ind = 0; ind < Math.min(sn.nnodes, sn.connmatrix.getNumRows()); ind++) {
                int nout = 0;
                for (int j = 0; j < sn.connmatrix.getNumCols(); j++) {
                    if (sn.connmatrix.get(ind, j) != 0) {
                        nout++;
                    }
                }
                if (nout <= 1 || sn.nodes == null || ind >= sn.nodes.size()) {
                    continue;
                }
                Map<JobClass, RoutingStrategy> perClass = sn.routing.get(sn.nodes.get(ind));
                if (perClass == null) {
                    continue;
                }
                int nrr = 0;
                for (RoutingStrategy rs : perClass.values()) {
                    if (rs == RoutingStrategy.RROBIN || rs == RoutingStrategy.WRROBIN) {
                        nrr++;
                    }
                }
                if (nrr > 0) {
                    logNstates += nrr * Math.log(nout);
                }
            }
        }
        return logNstates;
    }

    /**
     * How many of stations {@code [from, to)} can hold a job of class k at all.
     *
     * <p>A ZERO per-class capacity means the class is DISABLED there, so it never
     * occupies a slot and the placement term must spread it over the stations that
     * admit it rather than over all M. ld_whittle_bandwidth disables each of its
     * three PS routes for the other two classes; counting all M=4 priced it at
     * C(9,6)^3 = 592704 states, 7852 GB under the quadratic byte model, and the
     * gate refused a model whose true space is 7^3 = 343 and solves at once.
     */
    private static int admittingStations(NetworkStruct sn, int k, int from, int to) {
        int n = 0;
        for (int i = from; i < to; i++) {
            if (sn.classcap != null && i < sn.classcap.getNumRows()
                    && k < sn.classcap.getNumCols() && sn.classcap.get(i, k) == 0) {
                continue;
            }
            n++;
        }
        return n;
    }

    /** As {@link #admittingStations}, restricted to the NON-ordered stations. */
    private static int admittingRemaining(NetworkStruct sn, int k, int M) {
        int n = 0;
        for (int i = 0; i < M; i++) {
            SchedStrategy sc = null;
            if (sn.stations != null && i < sn.stations.size() && sn.sched != null) {
                sc = sn.sched.get(sn.stations.get(i));
            }
            if (!(sc == SchedStrategy.EXT || isShareScheduling(sc))) {
                continue;
            }
            if (sn.classcap != null && i < sn.classcap.getNumRows()
                    && k < sn.classcap.getNumCols() && sn.classcap.get(i, k) == 0) {
                continue;
            }
            n++;
        }
        return n;
    }

    /** Largest (m_1..m_K) box the exact ordered-buffer DP will walk. */
    private static final double ORDER_GRID_MAX = 1.0e6;

    /**
     * Log count of (placement, ordering) configurations over ALL order-preserving
     * stations at once, POPULATION CONSERVED. capsPer[a][k] bounds class k at
     * ordered station a; capTot[a] bounds the buffer TOTAL there, because a finite
     * station capacity is a slot count and not a per-class bound; njobs is the
     * population to share out and mRem share stations take the leftovers.
     *
     * Cutoff truncates an OPEN class's population IN THE NETWORK, exactly as the
     * plain stars-and-bars term treats it, so open classes are conserved too.
     */
    private static double logOrderedJoint(int[][] capsPer, double[] capTot, int[] njobs, int mRem) {
        int kb = njobs.length;
        int[] dims = new int[kb];
        int nstate = 1;
        for (int k = 0; k < kb; k++) {
            dims[k] = njobs[k] + 1;
            nstate *= dims[k];
        }
        double[] L = new double[nstate];
        java.util.Arrays.fill(L, Double.NEGATIVE_INFINITY);
        L[idxOf(njobs, dims)] = 0.0;
        for (int a = 0; a < capsPer.length; a++) {
            double[] ln = new double[nstate];
            java.util.Arrays.fill(ln, Double.NEGATIVE_INFINITY);
            for (int si = 0; si < nstate; si++) {
                if (Double.isInfinite(L[si]) && L[si] < 0) {
                    continue;
                }
                int[] rem = subOf(si, dims);
                int[] av = new int[kb];
                for (int k = 0; k < kb; k++) {
                    av[k] = Math.min(capsPer[a][k], rem[k]);
                }
                int[] m = new int[kb];
                while (true) {
                    int t = 0;
                    for (int k = 0; k < kb; k++) {
                        t += m[k];
                    }
                    if (!(capTot[a] < Double.POSITIVE_INFINITY && t > capTot[a])) {
                        double v = L[si] + Maths.factln(t);
                        for (int k = 0; k < kb; k++) {
                            v -= Maths.factln(m[k]);
                        }
                        int[] nx = new int[kb];
                        for (int k = 0; k < kb; k++) {
                            nx[k] = rem[k] - m[k];
                        }
                        int di = idxOf(nx, dims);
                        if (Double.isInfinite(ln[di]) && ln[di] < 0) {
                            ln[di] = v;
                        } else {
                            double mx = Math.max(ln[di], v);
                            ln[di] = mx + Math.log(Math.exp(ln[di] - mx) + Math.exp(v - mx));
                        }
                    }
                    int pos = 0;
                    while (pos < kb && m[pos] == av[pos]) {
                        m[pos] = 0;
                        pos++;
                    }
                    if (pos == kb) {
                        break;
                    }
                    m[pos]++;
                }
            }
            L = ln;
        }
        double top = Double.NEGATIVE_INFINITY;
        java.util.List<Double> terms = new java.util.ArrayList<Double>();
        for (int si = 0; si < nstate; si++) {
            if (Double.isInfinite(L[si]) && L[si] < 0) {
                continue;
            }
            int[] rem = subOf(si, dims);
            double v = L[si];
            if (mRem >= 1) {
                for (int k = 0; k < kb; k++) {
                    v += Maths.factln(rem[k] + mRem - 1) - Maths.factln(rem[k]) - Maths.factln(mRem - 1);
                }
            } else {
                boolean leftover = false;
                for (int k = 0; k < kb; k++) {
                    if (rem[k] > 0) {
                        leftover = true;
                    }
                }
                if (leftover) {
                    continue;
                }
            }
            terms.add(v);
            if (v > top) {
                top = v;
            }
        }
        if (terms.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }
        double acc = 0;
        for (int i = 0; i < terms.size(); i++) {
            acc += Math.exp(terms.get(i) - top);
        }
        return top + Math.log(acc);
    }

    private static int idxOf(int[] v, int[] dims) {
        int ix = 0;
        int mult = 1;
        for (int k = 0; k < dims.length; k++) {
            ix += v[k] * mult;
            mult *= dims[k];
        }
        return ix;
    }

    private static int[] subOf(int ix, int[] dims) {
        int[] v = new int[dims.length];
        int r = ix;
        for (int k = 0; k < dims.length; k++) {
            v[k] = r % dims[k];
            r /= dims[k];
        }
        return v;
    }

    private static boolean isShareScheduling(SchedStrategy sched) {
        return sched == SchedStrategy.INF || sched == SchedStrategy.PS
                || sched == SchedStrategy.DPS || sched == SchedStrategy.GPS
                || sched == SchedStrategy.PSPRIO || sched == SchedStrategy.DPSPRIO
                || sched == SchedStrategy.GPSPRIO || sched == SchedStrategy.LPS;
    }

    /**
     * Hardware-aware, calibrated pre-gate.
     *
     * @param logNstates natural log of the worst-case state-space size
     * @param force      bypass the hard stop (still warns)
     * @param verbose    print the estimate
     * @param safetyFraction fraction of available memory allowed as budget
     * @return a {@link GateResult}; ok is false only when the predicted
     *         footprint exceeds the budget and force is not set.
     */
    public static GateResult gate(double logNstates, boolean force,
                                  boolean verbose, double safetyFraction) {
        double avail = getAvailableMemoryBytes();
        double budget = safetyFraction * avail;
        Calibration calib = getCalibration(verbose);

        double logPred = Math.log(calib.alphaMem) + calib.betaMem * logNstates;
        double logBudget = Math.log(Math.max(budget, 1.0));
        double predGB = Math.exp(Math.min(logPred, 700.0)) / (1024.0 * 1024 * 1024);
        double budgetGB = budget / (1024.0 * 1024 * 1024);

        if (logPred > logBudget) {
            String msg = String.format(
                    "CTMC predicted peak memory ~%.2f GB exceeds the safe budget "
                    + "~%.2f GB (%.0f%% of %.2f GB available). Reduce the state "
                    + "space (e.g. lower 'cutoff'), use another solver (MVA/NC/FLD), "
                    + "or set force=true to override.",
                    predGB, budgetGB, 100 * safetyFraction,
                    avail / (1024.0 * 1024 * 1024));
            if (!force) {
                return new GateResult(false, msg);
            }
            if (verbose) {
                System.out.println("Warning (forced): " + msg);
            }
            return new GateResult(true, msg);
        }
        if (verbose && logPred > Math.log(Math.max(0.5 * budget, 1.0))) {
            System.out.println(String.format(
                    "CTMC predicted peak memory ~%.2f GB (budget ~%.2f GB).",
                    predGB, budgetGB));
        }
        return new GateResult(true, "");
    }

    private static Calibration getCalibration(boolean verbose) {
        String sig = machineSignature();
        if (cachedCalibration != null && sig.equals(cachedCalibration.sig)) {
            return cachedCalibration;
        }
        Calibration disk = loadCache(sig);
        if (disk != null) {
            cachedCalibration = disk;
            return disk;
        }
        Calibration calib = new Calibration();
        calib.sig = sig;
        try {
            double[] p1 = profilePoint(CALIB_G1);
            double[] p2 = profilePoint(CALIB_G2);
            double[] mem = fitPowerLaw(p1[0], p1[1], p2[0], p2[1]);
            double[] tim = fitPowerLaw(p1[0], Math.max(p1[2], 1e-9),
                    p2[0], Math.max(p2[2], 1e-9));
            if (mem == null) {
                throw new IllegalStateException("degenerate memory fit");
            }
            calib.alphaMem = mem[0];
            calib.betaMem = mem[1];
            calib.alphaT = (tim != null) ? tim[0] : 0.0;
            calib.betaT = (tim != null) ? tim[1] : 1.0;
            storeCache(calib);
            if (verbose) {
                System.out.println(String.format(
                        "CTMC calibration: bytes ~ %.3g*N^%.3f",
                        calib.alphaMem, calib.betaMem));
            }
        } catch (Throwable exc) {
            if (verbose) {
                System.out.println("CTMC calibration failed (" + exc.getMessage()
                        + "); using fallback model");
            }
            calib.alphaMem = FALLBACK_ALPHA;
            calib.betaMem = FALLBACK_BETA;
            calib.alphaT = 0.0;
            calib.betaT = 1.0;
        }
        cachedCalibration = calib;
        return calib;
    }

    /**
     * Factorize the nonsingular block of a g-by-g nearest-neighbour lattice
     * generator (n = g*g states) with EJML's sparse LU. The QBD-like structure
     * is representative of multi-station queueing generators and is cheap and
     * safe to factorize (unlike a random matrix, whose factor is near dense).
     *
     * @return {n, factorBytes, seconds}
     */
    private static double[] profilePoint(int g) {
        DMatrixSparseCSC a = latticeGenerator(g);
        int n = a.getNumRows();
        long t0 = System.nanoTime();
        LUSparseDecomposition_F64<DMatrixSparseCSC> lu =
                DecompositionFactory_DSCC.lu(FillReducing.NONE);
        if (!lu.decompose(a)) {
            throw new IllegalStateException("LU decomposition failed at g=" + g);
        }
        DMatrixSparseCSC lower = lu.getLower(null);
        DMatrixSparseCSC upper = lu.getUpper(null);
        double secs = (System.nanoTime() - t0) / 1.0e9;
        double factnnz = (double) lower.nz_length + (double) upper.nz_length;
        return new double[]{(double) n, BYTES_PER_NZ * factnnz, secs};
    }

    private static DMatrixSparseCSC latticeGenerator(int g) {
        int n = g * g;
        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(n, n, 5 * n);
        double[] rowsum = new double[n];
        for (int r = 0; r < g; r++) {
            for (int c = 0; c < g; c++) {
                int s = r * g + c;
                // right, left, down, up neighbours (unit rates), no wrap-around
                if (c + 1 < g) { addEdge(triplet, rowsum, s, r * g + (c + 1)); }
                if (c - 1 >= 0) { addEdge(triplet, rowsum, s, r * g + (c - 1)); }
                if (r + 1 < g) { addEdge(triplet, rowsum, s, (r + 1) * g + c); }
                if (r - 1 >= 0) { addEdge(triplet, rowsum, s, (r - 1) * g + c); }
            }
        }
        // diagonal = -rowsum (infinitesimal generator)
        for (int s = 0; s < n; s++) {
            triplet.addItem(s, s, -rowsum[s]);
        }
        DMatrixSparseCSC full = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        // drop the last state -> nonsingular block (as ctmc_solve does)
        DMatrixSparseCSC a = new DMatrixSparseCSC(n - 1, n - 1, 0);
        org.ejml.sparse.csc.CommonOps_DSCC.extract(
                full, 0, n - 1, 0, n - 1, a, 0, 0);
        return a;
    }

    private static void addEdge(DMatrixSparseTriplet triplet, double[] rowsum,
                                int from, int to) {
        triplet.addItem(from, to, 1.0);
        rowsum[from] += 1.0;
    }

    private static double[] fitPowerLaw(double x1, double y1, double x2, double y2) {
        if (x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0 || x1 == x2) {
            return null;
        }
        double beta = Math.log(y2 / y1) / Math.log(x2 / x1);
        double alpha = y1 / Math.pow(x1, beta);
        return new double[]{alpha, beta};
    }

    private static Calibration loadCache(String sig) {
        File f = cacheFile();
        if (!f.exists()) {
            return null;
        }
        BufferedReader br = null;
        try {
            Properties props = new Properties();
            br = new BufferedReader(new FileReader(f));
            props.load(br);
            if (!sig.equals(props.getProperty("sig"))) {
                return null;
            }
            Calibration c = new Calibration();
            c.sig = sig;
            c.alphaMem = Double.parseDouble(props.getProperty("alpha_mem"));
            c.betaMem = Double.parseDouble(props.getProperty("beta_mem"));
            c.alphaT = Double.parseDouble(props.getProperty("alpha_t", "0"));
            c.betaT = Double.parseDouble(props.getProperty("beta_t", "1"));
            return c;
        } catch (Throwable ignore) {
            return null;
        } finally {
            if (br != null) {
                try { br.close(); } catch (Throwable ignore) { }
            }
        }
    }

    private static void storeCache(Calibration c) {
        FileWriter fw = null;
        try {
            Properties props = new Properties();
            props.setProperty("sig", c.sig);
            props.setProperty("alpha_mem", Double.toString(c.alphaMem));
            props.setProperty("beta_mem", Double.toString(c.betaMem));
            props.setProperty("alpha_t", Double.toString(c.alphaT));
            props.setProperty("beta_t", Double.toString(c.betaT));
            props.setProperty("timestamp", Long.toString(System.currentTimeMillis()));
            fw = new FileWriter(cacheFile());
            props.store(fw, "LINE CTMC memory-guard calibration");
        } catch (Throwable ignore) {
            // caching is best-effort
        } finally {
            if (fw != null) {
                try { fw.close(); } catch (Throwable ignore) { }
            }
        }
    }
}
