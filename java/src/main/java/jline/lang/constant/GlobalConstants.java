package jline.lang.constant;

import jline.solvers.SolverOptions;

/**
 *  Global constants for tolerances and solver configuration
 */
public class GlobalConstants {
    public static final double Zero =  1.0000e-14;
    public static final double CoarseTol = 1.0000e-03;
    public static final double FineTol = 1.0000e-08;
    public static final double Immediate = 100000000;
    public static final String Version = "2.0.31";
    public static final VerboseLevel Verbose = VerboseLevel.STD;
    public static final boolean DummyMode = false;
}
