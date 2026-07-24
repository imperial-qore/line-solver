package jline.solvers.wrappers.qns;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * QNSResult class stores the results from the QNS solver.
 * Inherits all metric fields (QN, UN, RN, TN, AN, WN, CN, XN, runtime, method, iter)
 * from SolverResult.
 */
public class QNSResult extends SolverResult {

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
