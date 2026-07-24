package jline.lang;

import jline.util.matrix.Matrix;

import java.util.List;

/**
 * Matrix representation for job class switching probabilities in queueing networks.
 * 
 * <p>ClassSwitchMatrix encapsulates the class switching behavior at nodes where jobs
 * can change their class during service. This is essential for modeling systems where
 * jobs undergo different processing stages or priority changes during their traversal
 * through the network.</p>
 * 
 * <p>Key features:
 * <ul>
 *   <li>Stochastic matrix representation [classes x classes]</li>
 *   <li>Entry (i,j) represents probability of switching from class i to class j</li>
 *   <li>Row sums equal 1.0 for valid probability distributions</li>
 *   <li>Support for both deterministic and probabilistic class switching</li>
 *   <li>Integration with routing matrix computations</li>
 * </ul>
 * </p>
 * 
 * <p>Common applications include modeling priority upgrades/downgrades, 
 * multi-stage service processes, and feedback systems with class-dependent routing.</p>
 * 
 * @see jline.lang.nodes.ClassSwitch
 * @see RoutingMatrix
 * @since 1.0
 */
public class ClassSwitchMatrix {
    private final Matrix m;

    /* --- factory-style constructors that build the inner Matrix --- */
    public ClassSwitchMatrix(double scalar) {
        this.m = Matrix.singleton(scalar);
    }

    public ClassSwitchMatrix(int numRows, int numCols, int arrayLength) {
        this.m = new Matrix(numRows, numCols, arrayLength);
    }

    public ClassSwitchMatrix(int numRows, int numCols) {
        this.m = new Matrix(numRows, numCols);
    }

    public ClassSwitchMatrix(Matrix matrix) {
        this.m = new Matrix(matrix);   // deep-copy to preserve immutability if desired
    }

    public ClassSwitchMatrix(List<Double> array) {
        this.m = new Matrix(array);
    }

    public double get(int r, int c) {
        return m.get(r, c);
    }          // expose whatever is useful

    /**
     * Expose the raw Matrix only if callers genuinely need it
     */
    public Matrix getMatrix() {
        return m;
    }

    public int getNumCols() {
        return m.getNumCols();
    }

    public int getNumRows() {
        return m.getNumRows();
    }

    public int length() {
        return m.length();
    }

    public void set(int jc1, int jc2, double value) {
        m.set(jc1, jc2, value);
    }

    public void set(JobClass jc1, JobClass jc2, double value) {
        m.set(jc1.getIndex()-1, jc2.getIndex()-1, value);
    }

    public void setTo(Matrix m1) {
        m.setTo(m1);
    }

    public double[][] toArray2D() {
        return m.toArray2D();
    }
}
