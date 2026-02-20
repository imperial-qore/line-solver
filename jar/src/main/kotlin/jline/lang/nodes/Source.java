/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.ServiceBinding;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.sections.Dispatcher;
import jline.lang.sections.RandomSource;
import jline.lang.sections.ServiceTunnel;
import jline.lang.workflow.Workflow;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An abstraction of the external world jobs in open classes come from
 */
public class Source extends Station implements Serializable {
    protected List<ServiceBinding> arrivalProcess;
    protected SchedStrategy schedStrategy;

    /**
     * Creates a new source node in the specified network.
     * Initializes arrival processes for all job classes with Disabled distribution.
     * 
     * @param model the network model to add this source to
     * @param name the name for this source node
     */
    public Source(Network model, String name) {
        super(name);
        this.numberOfServers = 1;

        List<JobClass> jobClasses = model.getClasses();
        this.output = new Dispatcher(jobClasses);
        this.server = new ServiceTunnel();
        this.input = new RandomSource(jobClasses);
        this.schedStrategy = SchedStrategy.EXT;
        this.setModel(model);
        this.model.addNode(this);
        this.arrivalProcess = new ArrayList<ServiceBinding>();

        for (JobClass jobClass : jobClasses) {
            this.classCap.put(jobClass, Integer.MAX_VALUE);
            this.setArrival(jobClass, new Disabled());
        }
    }

    /**
     * Checks if this source has an arrival process configured for the specified job class.
     * 
     * @param jobClass the job class to check
     * @return true if an arrival process exists for this job class, false otherwise
     */
    public boolean containsJobClass(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.arrivalProcess) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Gets the arrival distribution for a specific job class.
     * 
     * @param jobClass the job class to query
     * @return the arrival distribution, or Disabled if none configured
     */
    public Distribution getArrivalDistribution(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.arrivalProcess) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return serviceProcess.getDistribution();
            }
        }
        return new Disabled();
    }

    /**
     * Gets the arrival process distribution for a specific job class.
     * This is an alias for getArrivalDistribution.
     * 
     * @param jobClass the job class to query
     * @return the arrival distribution, or Disabled if none configured
     */
    public final Distribution getArrivalProcess(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.arrivalProcess) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return serviceProcess.getDistribution();
            }
        }

        return new Disabled();
    }

    /**
     * Gets the scheduling strategy for this source node.
     * Sources always use EXT (external) scheduling strategy.
     * 
     * @return the scheduling strategy (always SchedStrategy.EXT)
     */
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    @Override
    public void printSummary() {
        System.out.format("jline.Source:\n");
        System.out.format("--Name: %s\n", this.getName());
        System.out.format("--Arrival Processes:\n");
        for (JobClass jobClass : this.model.getClasses()) {
            System.out.format("----%s: %s\n", jobClass.getName(), this.getArrivalDistribution(jobClass).toString());
        }

        this.output.printSummary();
    }

    /**
     * Removes the arrival process for a specific job class.
     * 
     * @param jobClass the job class whose arrival process should be removed
     */
    protected void removeArrivalProcess(JobClass jobClass) {
        Iterator<ServiceBinding> serviceProcessIterator = this.arrivalProcess.iterator();
        while (serviceProcessIterator.hasNext()) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcessIterator.next().getJobClass().getIndex() == jobClass.getIndex()) {
                serviceProcessIterator.remove();
            }
        }
        ((RandomSource) this.input).removeServiceProcess(jobClass);
    }

    /**
     * Sets the arrival distribution for a specific job class.
     * If distribution is null or Disabled, sets class capacity to 0.
     * If a Workflow is provided, it will be converted to a PH distribution.
     *
     * @param jobClass the job class to configure
     * @param distribution the arrival distribution or workflow to set
     */
    public void setArrival(JobClass jobClass, Object distribution) {
        // If Workflow, convert to PH distribution
        Distribution actualDistribution;
        if (distribution instanceof Workflow) {
            actualDistribution = ((Workflow) distribution).toPH();
        } else if (distribution instanceof Distribution) {
            actualDistribution = (Distribution) distribution;
        } else {
            throw new IllegalArgumentException("distribution must be a Distribution or Workflow object");
        }

        ServiceBinding arrivalProcess = new ServiceBinding(jobClass, ServiceStrategy.LI, actualDistribution);
        //this.input.setServiceProcess(arrivalProcess);
        this.setArrivalProcess(arrivalProcess);
        if ((actualDistribution == null) || (actualDistribution instanceof Disabled)) {
            this.classCap.put(jobClass, 0);
        } else {
            this.classCap.put(jobClass, Integer.MAX_VALUE);
        }
    }

    /**
     * Sets the arrival process using a service binding.
     * Removes any existing arrival process for the same job class.
     * 
     * @param arrivalProcess the service binding defining the arrival process
     */
    public void setArrivalProcess(ServiceBinding arrivalProcess) {
        removeArrivalProcess(arrivalProcess.getJobClass());
        this.arrivalProcess.add(arrivalProcess);
    }

}
