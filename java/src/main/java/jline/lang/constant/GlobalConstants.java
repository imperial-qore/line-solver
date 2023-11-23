package jline.lang.constant;

/**
 *  Global constants for tolerances and solver configuration
 */
public class GlobalConstants {
    private static GlobalConstants single_instance = null;

    public static final double Zero =  1.0000e-14;
    public static final double CoarseTol = 1.0000e-03;
    public static final double FineTol = 1.0000e-08;
    public static final double Immediate = 100000000;
    public static final String Version = "2.0.31";
    public static VerboseLevel Verbose = VerboseLevel.STD;
    public static final boolean DummyMode = false;

    public void setVerbose(VerboseLevel verbosity) {
        this.Verbose = verbosity;
    }

    public VerboseLevel getVerbose() {
        return this.Verbose;
    }

    // Static method to create a singleton instance of GlobalConstants class
    public static synchronized GlobalConstants getInstance()
    {
        if (single_instance == null) {
            single_instance = new GlobalConstants();
            single_instance.setVerbose(VerboseLevel.STD);
        }

        return single_instance;
    }
}
