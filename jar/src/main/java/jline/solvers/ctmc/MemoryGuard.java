package jline.solvers.ctmc;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;
import java.util.Properties;

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
