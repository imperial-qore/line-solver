package jline.opt.variables;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

import java.util.List;

/**
 * Abstract base for line-opt decision variables. Mirrors native-Python
 * {@code line_solver.opt.variables.DecisionVariable}: each variable encodes a
 * tunable model parameter as continuous values in [0, 1] (see
 * {@link #getBounds()}), decodes them to the native domain ({@link #decode}),
 * and applies the decoded value to a per-evaluation model copy
 * ({@link #apply}). Decision variables hold references to the base model's
 * objects; because models are deep-copied per evaluation, objects are
 * re-resolved by name in the target model.
 */
public abstract class DecisionVariable {

    protected final String name;
    protected int dimension = 1;

    protected DecisionVariable(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public int getDimension() {
        return dimension;
    }

    /** Bounds per encoded dimension, each {low, high}. Standard is {0,1}. */
    public abstract double[][] getBounds();

    /** Decode the encoded slice (length {@link #getDimension()}) to a value. */
    public abstract Object decode(double[] x);

    /** Apply a decoded value to a (copied) model. */
    public abstract void apply(Network model, Object value);

    /** Type identifier used by decomposition and the gradient dispatcher. */
    public abstract String getVariableType();

    /**
     * LQN layer name(s) this variable perturbs, or null for flat-network
     * variables. Used by layer freezing (both explicit {@code frozen_layers} and
     * adaptive {@code auto_freeze}): a variable whose layer set intersects the
     * frozen set is held fixed. The model argument is the LayeredNetwork the
     * variable is applied to. LQN variable subclasses override this.
     */
    public List<String> getLayer(Object model) {
        return null;
    }

    /**
     * The variable's current (decoded) value in the given model, or null. Used
     * by layer freezing to hold a variable at the model's existing parameter
     * value. Flat and non-introspectable variables return null. LQN variable
     * subclasses override this.
     */
    public Object currentValue(Object model) {
        return null;
    }

    // ---- shared resolution helpers ----------------------------------------

    protected static JobClass resolveClass(Network model, JobClass jobclass) {
        String target = jobclass.getName();
        for (JobClass c : model.getClasses()) {
            if (c.getName().equals(target)) {
                return c;
            }
        }
        return null;
    }

    protected static Node resolveNode(Network model, Node node) {
        String target = node.getName();
        for (Node c : model.getNodes()) {
            if (c.getName().equals(target)) {
                return c;
            }
        }
        return null;
    }

    /** Node-by-node binary connection (adjacency) matrix. */
    protected static Matrix connectionMatrix(Network model) {
        return model.getConnectionMatrix();
    }

    protected static double[][] unitBounds(int dim) {
        double[][] b = new double[dim][2];
        for (int i = 0; i < dim; i++) {
            b[i][0] = 0.0;
            b[i][1] = 1.0;
        }
        return b;
    }

    protected static int indexOfNode(List<Node> nodes, String nameToFind) {
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i).getName().equals(nameToFind)) {
                return i;
            }
        }
        return -1;
    }
}
