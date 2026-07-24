package jline.opt.variables;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.nodes.Node;

import java.util.List;

/**
 * Optimize the routing probabilities of a job class from a source node to a
 * list of target nodes. Uses stick-breaking encoding (dimension = targets - 1)
 * so probabilities sum to 1. On apply, default routing is reconstructed from
 * the connection matrix for every class (a uniform split over outgoing links,
 * reproducing LINE's link() semantics), then the overridden (class, source) row
 * is replaced with the decoded probabilities. Mirrors native-Python
 * {@code RoutingProbabilities}.
 */
public class RoutingProbabilities extends DecisionVariable {

    private final JobClass jobclass;
    private final Node source;
    private final List<Node> targets;

    public RoutingProbabilities(JobClass jobclass, Node source, List<Node> targets) {
        this(jobclass, source, targets,
                jobclass.getName() + "_routing_from_" + source.getName());
    }

    public RoutingProbabilities(JobClass jobclass, Node source, List<Node> targets, String name) {
        super(name);
        this.jobclass = jobclass;
        this.source = source;
        this.targets = targets;
        this.dimension = Math.max(1, targets.size() - 1);
    }

    public JobClass getJobClass() {
        return jobclass;
    }

    public Node getSource() {
        return source;
    }

    public List<Node> getTargets() {
        return targets;
    }

    public double[][] getBounds() {
        return unitBounds(dimension);
    }

    public Object decode(double[] x) {
        int n = targets.size();
        if (n == 1) {
            return new double[]{1.0};
        }
        double[] probs = new double[n];
        double remaining = 1.0;
        for (int i = 0; i < n - 1; i++) {
            probs[i] = remaining * x[i];
            remaining -= probs[i];
        }
        probs[n - 1] = remaining;
        return probs;
    }

    public void apply(Network model, Object value) {
        // Override only the source node's outgoing routing for this class, via
        // setProbRouting, leaving all other routes intact. This works whether
        // the model was built with link() or addLink() (model.link() is
        // rejected after addLink()).
        double[] probs = (double[]) value;
        JobClass cls = resolveClass(model, jobclass);
        if (cls == null) {
            return;
        }
        Node src = resolveNode(model, source);
        if (src == null) {
            return;
        }
        for (int i = 0; i < targets.size(); i++) {
            Node target = resolveNode(model, targets.get(i));
            if (target != null) {
                src.setProbRouting(cls, target, probs[i]);
            }
        }
    }

    public String getVariableType() {
        return "routing";
    }
}
