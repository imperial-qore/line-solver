package jline.api.qsys;

import jline.util.matrix.Matrix;

/**
 * Result of BMAP/PH/N/N retrial queue analysis.
 */
public class QsysRetrialResult {
    public final double L_orbit;
    public final double N_server;
    public final double L_system;
    public final double Utilization;
    public final double Throughput;
    public final double P_idle;
    public final double P_empty_orbit;
    public final double P_empty_system;
    public final Matrix pi;
    public final int truncLevel;
    /** Relative orbit-truncation error estimate at truncLevel. */
    public final double truncError;
    public final String analyzer;

    public QsysRetrialResult(double L_orbit, double N_server, double L_system, double Utilization,
                             double Throughput, double P_idle, double P_empty_orbit, double P_empty_system,
                             Matrix pi, int truncLevel, String analyzer) {
        this(L_orbit, N_server, L_system, Utilization, Throughput, P_idle, P_empty_orbit, P_empty_system,
                pi, truncLevel, Double.NaN, analyzer);
    }

    public QsysRetrialResult(double L_orbit, double N_server, double L_system, double Utilization,
                             double Throughput, double P_idle, double P_empty_orbit, double P_empty_system,
                             Matrix pi, int truncLevel, double truncError, String analyzer) {
        this.L_orbit = L_orbit;
        this.N_server = N_server;
        this.L_system = L_system;
        this.Utilization = Utilization;
        this.Throughput = Throughput;
        this.P_idle = P_idle;
        this.P_empty_orbit = P_empty_orbit;
        this.P_empty_system = P_empty_system;
        this.pi = pi;
        this.truncLevel = truncLevel;
        this.truncError = truncError;
        this.analyzer = analyzer;
    }

    public double getL_orbit() { return L_orbit; }
    public double getN_server() { return N_server; }
    public double getL_system() { return L_system; }
    public double getUtilization() { return Utilization; }
    public double getThroughput() { return Throughput; }
    public double getP_idle() { return P_idle; }
    public double getP_empty_orbit() { return P_empty_orbit; }
    public double getP_empty_system() { return P_empty_system; }
    public Matrix getPi() { return pi; }
    public int getTruncLevel() { return truncLevel; }
    public double getTruncError() { return truncError; }
    public String getAnalyzer() { return analyzer; }
}
