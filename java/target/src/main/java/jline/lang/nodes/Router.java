package jline.lang.nodes;

import java.io.Serializable;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.constant.ServiceStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

/**
 * A node that routes jobs without imposing any delay
 */
public class Router extends Node implements HasSchedStrategy, Serializable {
    protected double cap;
    protected int numberOfServers;
    protected SchedStrategyType schedPolicy;
    protected SchedStrategy schedStrategy;

    public Router(Network model, String name) {
        super(name);

        List<JobClass> jobClasses = model.getClasses();
        this.server = new ServiceTunnel();
        this.input = new Buffer(jobClasses);
        this.output = new Dispatcher(jobClasses);
        this.cap = Double.POSITIVE_INFINITY;

        this.schedPolicy = SchedStrategyType.NP;
        this.schedStrategy = SchedStrategy.FCFS;
        this.numberOfServers = 1;
        this.setModel(model);
        this.model.addNode(this);
    }

    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    public Distribution getServiceProcess(JobClass jobClass) {
        return this.server.getServiceDistribution(jobClass);
    }

    public void setService(JobClass jobClass, Distribution distribution) {
        this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, distribution));
    }

    public void setSchedPolicy(SchedStrategyType schedPolicy) {
        this.schedPolicy = schedPolicy;
    }

    public double minRate() {
        return Double.POSITIVE_INFINITY;
    }

    public double maxRate() {
        return 0.0;
    }
    public double avgRate() {
        return 0.0;
    }
    public int rateCt() {return 0;}
}
