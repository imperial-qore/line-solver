/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.constant.JoinStrategy;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Input section of a join node
 */
public class Joiner extends InputSection implements Serializable {
    public Map<JobClass, JoinStrategy> joinStrategy;
    public Map<JobClass, Double> joinRequired;
    public Map<JobClass, JobClass> joinJobClasses;
    protected List<JobClass> jobClasses;

    public Joiner(List<JobClass> customerClasses) {
        super("Joiner");
        this.joinJobClasses = new HashMap<>();
        this.jobClasses = customerClasses;
        this.joinStrategy = new HashMap<>();
        this.joinRequired = new HashMap<>();

        for (JobClass jobclass : this.jobClasses) {
            this.joinStrategy.put(jobclass, JoinStrategy.STD);
            this.joinRequired.put(jobclass, -1.0);
            this.joinJobClasses.put(jobclass, jobclass);
        }
    }

    public void setRequired(JobClass jobClass, double njobs) {
        this.joinRequired.put(jobClass, njobs);
    }

    public void setStrategy(JobClass jobClass, JoinStrategy joinStrategy) {
        this.joinStrategy.put(jobClass, joinStrategy);
    }

//    @Override
//    public void setOutputStrategy(JobClass jobClass, RoutingStrategy routingStrategy, Node destination, double probability) {
//        for (OutputStrategy outputStrategy : this.outputStrategies) {
//            if ((outputStrategy.getJobClass() == jobClass) && (outputStrategy.getDestination() == destination)) {
//                outputStrategy.setRoutingStrategy(routingStrategy);
//                outputStrategy.setProbability(probability);
//                this.probabilityUpdate();
//                return;
//            }
//        }
//
//        OutputStrategy outputStrategy = new OutputStrategy(jobClass, routingStrategy, destination, probability);
//        outputStrategies.add(outputStrategy);
//        outputEvents.put(outputStrategy, new JoinOutputEvent(this, destination, jobClass));
//        this.probabilityUpdate();
//    }
}
