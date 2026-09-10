package jline;

/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Global constants for tolerances and solver configuration.
 * 
 * <p>This class defines system-wide constants used throughout the LINE solver including:
 * <ul>
 *   <li>Numerical tolerances for convergence checks</li>
 *   <li>Special values (infinity, zero threshold)</li>
 *   <li>Version information</li>
 *   <li>Verbosity settings</li>
 * </ul>
 * </p>
 * 
 * <p>The class follows a singleton pattern to ensure consistent configuration
 * across all solver components.</p>
 */
public class GlobalConstants {
    /** Coarse tolerance for approximate comparisons (1e-3) */
    public static final double CoarseTol = 1.0000e-03;
    /** Debug mode flag */
    public static final boolean DummyMode = false;
    /** Fine tolerance for precise comparisons (1e-8) */
    public static final double FineTol = 1.0000e-08;
    /** Value representing immediate service (1/FineTol) */
    public static final double Immediate = 1.0 / FineTol;
    /** Positive infinity constant */
    public static final double Inf = Double.POSITIVE_INFINITY;
    /** Maximum integer value */
    public static final int MaxInt = Integer.MAX_VALUE;
    /** Negative infinity constant */
    public static final double NegInf = Double.NEGATIVE_INFINITY;
    /** LINE solver version */
    public static final String Version = "3.0.7";
    /** Threshold below which values are considered zero (1e-14) */
    public static final double Zero = 1.0000e-14;
    /** Magnitude above which an off-diagonal generator entry counts as an arc (1e-12).
     *  Sign is NOT a criterion: a matrix-exponential generator embeds genuinely negative
     *  off-diagonal entries, and a sign test would read them as absent arcs. */
    public static final double ArcTol = 1.0000e-12;
    /** Global verbosity level for solver output */
    public static VerboseLevel Verbose = initVerbose();

    /**
     * Resolves the initial verbosity from the optional {@code jline.verbose} system
     * property (values: SILENT, STD, DEBUG), falling back to {@link VerboseLevel#STD}
     * when the property is absent or unrecognized. This lets the test harness force
     * SILENT from JVM startup (via surefire {@code -Djline.verbose=SILENT}) so that no
     * gated warning can leak during class static initialization, before any per-test
     * setup runs. Production behavior is unchanged when the property is not set.
     */
    private static VerboseLevel initVerbose() {
        String p = System.getProperty("jline.verbose");
        if (p != null && p.trim().length() > 0) {
            try {
                return VerboseLevel.valueOf(p.trim().toUpperCase());
            } catch (IllegalArgumentException ignored) {
                // unknown value -> fall back to default
            }
        }
        return VerboseLevel.STD;
    }
    private static GlobalConstants single_instance = null;

    /**
     * Returns the singleton instance of GlobalConstants.
     * Creates the instance if it doesn't exist.
     * 
     * @return the singleton GlobalConstants instance
     */
    public static synchronized GlobalConstants getInstance() {
        if (single_instance == null) {
            single_instance = new GlobalConstants();
            // Note: Verbose is already initialized to STD at the field level.
            // Removed redundant setVerbose(STD) call here to avoid overriding
            // any verbose level set by test setup code.
        }
        return single_instance;
    }

    /**
     * Gets the current verbosity level.
     * 
     * @return the current VerboseLevel setting
     */
    public static VerboseLevel getVerbose() {
        return Verbose;
    }
    
    /**
     * Sets the verbosity level for solver output.
     * 
     * @param verbosity the desired verbosity level
     */
    public static void setVerbose(VerboseLevel verbosity) {
        Verbose = verbosity;
    }

    /**
     * Gets the coarse tolerance value.
     * @return coarse tolerance (1e-3)
     */
    public static double getCoarseTol() {
        return CoarseTol;
    }
    
    /**
     * Gets the fine tolerance value.
     * @return fine tolerance (1e-8)
     */
    public static double getFineTol() {
        return FineTol;
    }
    
    /**
     * Gets the zero threshold value.
     * @return zero threshold (1e-14)
     */
    public static double getZero() {
        return Zero;
    }
    
    /**
     * Gets the maximum integer value.
     * @return maximum integer value
     */
    public static int getMaxInt() {
        return MaxInt;
    }
    
    /**
     * Gets the LINE solver version.
     * @return version string
     */
    public static String getVersion() {
        return Version;
    }
}
