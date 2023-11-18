package jline.lang;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Node;

/**
 * Class modelling the output section of a Node
 */
public class OutputStrategy implements Serializable {
    private final JobClass jobClass;
    private RoutingStrategy routingStrategy;
    private double probability;
    private Node destination;


    public static List<RoutingStrategy> legalStrategies = Arrays.asList(RoutingStrategy.DISABLED, RoutingStrategy.PROB, RoutingStrategy.RAND);

    public OutputStrategy (JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        this.jobClass = jobClass;
        this.routingStrategy = routingStrategy;
        this.destination = destination;
        this.probability = probability;

        if (!legalStrategies.contains(routingStrategy)) {
            throw new RuntimeException("Unsupported Routing Strategy!");
        }
    }

    public OutputStrategy (JobClass jobClass, RoutingStrategy routingStrategy) {
        this(jobClass, routingStrategy, null, 1);
    }

    public JobClass getJobClass() {
        return this.jobClass;
    }

    public RoutingStrategy getRoutingStrategy() {
        return this.routingStrategy;
    }

    public void setRoutingStrategy(RoutingStrategy routingStrategy) {
        this.routingStrategy = routingStrategy;
    }

    public double getProbability() {
        return this.probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    public Node getDestination() {
        return this.destination;
    }

    public void setDestination(Node destination) {
        this.destination = destination;
    }
}
