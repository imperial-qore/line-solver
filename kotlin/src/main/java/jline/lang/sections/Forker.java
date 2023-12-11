package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.OutputStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Node;
import jline.solvers.ssa.events.ForkEvent;

import java.util.List;

/**
 * Output section that forks incoming jobs into sibling tasks
 */
public class Forker extends OutputSection {
	public double taskPerLink;
    protected List<JobClass> jobClasses;

    public Forker(List<JobClass> customerClasses) {
        super("Forker");
        this.jobClasses = customerClasses;
        this.taskPerLink = 1.0;
        this.initDispatcherJobClasses(customerClasses);
    }

    private void initDispatcherJobClasses(List<JobClass> customerClasses){
        for(JobClass jobClass : customerClasses){
            this.outputStrategies.add(new OutputStrategy(jobClass, RoutingStrategy.RAND));
        }
    }

    @Override
    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
        for (OutputStrategy outputStrategy : this.outputStrategies) {
            if ((outputStrategy.getJobClass() == jobClass) && (outputStrategy.getDestination() == destination)) {
                outputStrategy.setRoutingStrategy(routingStrategy);
                outputStrategy.setProbability(probability);
                this.probabilityUpdate();
                return;
            }
        }

        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy, destination, probability);
        outputStrategies.add(outputStrategy);
        outputEvents.put(outputStrategy, new ForkEvent(this, destination, jobClass));
        this.probabilityUpdate();
    }
}
