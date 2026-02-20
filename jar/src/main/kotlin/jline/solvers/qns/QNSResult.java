package jline.solvers.qns;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * QNSResult class stores the results from the QNS solver
 */
public class QNSResult extends SolverResult {
    
    public Matrix QN;  // Queue lengths [stations x classes]
    public Matrix UN;  // Utilizations [stations x classes]
    public Matrix RN;  // Response times [stations x classes]
    public Matrix TN;  // Throughputs [stations x classes]
    public Matrix AN;  // Arrival rates [stations x classes]
    public Matrix WN;  // Residence times [stations x classes]
    public Matrix CN;  // System response times [1 x chains]
    public Matrix XN;  // System throughputs [1 x chains]
    public double runtime;
    public String method;
    public int iter;
    
    /**
     * Default constructor
     */
    public QNSResult() {
        super();
    }
    
    /**
     * Constructor with all parameters
     */
    public QNSResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, 
                    Matrix AN, Matrix WN, Matrix CN, Matrix XN,
                    double runtime, String method, int iter) {
        super();
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.AN = AN;
        this.WN = WN;
        this.CN = CN;
        this.XN = XN;
        this.runtime = runtime;
        this.method = method;
        this.iter = iter;
    }
}