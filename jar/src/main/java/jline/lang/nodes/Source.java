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
    // The shared marked process driving the marked classes (null if none). Typed as
    // Distribution rather than MarkedMAP because MMAPt and MPHt are marked too and do
    // not subclass it: they carry a schedule, so they deliberately stay outside the
    // Markovian hierarchy for the same reason MAPt does.
    protected jline.lang.processes.Distribution markedProcess;
    protected List<JobClass> markedClasses; // classes bound to marks 1..K of markedProcess (null if none)
    // Per-class batch-size law: each arrival epoch releases a batch of this many
    // jobs instead of one. Null entry (the default) means single arrivals.
    protected java.util.Map<JobClass, jline.lang.processes.DiscreteDistribution> arrivalBatch
            = new java.util.HashMap<JobClass, jline.lang.processes.DiscreteDistribution>();

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

    @Override
    public void removeJobClass(JobClass jobClass) {
        super.removeJobClass(jobClass);
        this.removeArrivalProcess(jobClass);
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
     * Sets a batch-size law for a class, turning each arrival epoch of that
     * class into the simultaneous release of a batch of jobs.
     *
     * <p>The interarrival distribution set by {@link #setArrival} continues to
     * govern the spacing of the epochs; this governs how many jobs each epoch
     * releases. A {@code Geometric(a)} interarrival time with a
     * {@code Geometric(beta)} batch size is the Geo^X arrival process of a
     * discrete-time queue, whose analytical counterpart is
     * {@code jline.api.dqsys.Dqsys_geoxgeo1}.
     *
     * <p>The batch size must be supported on {1,2,...}: an epoch that releases
     * no job is not an arrival epoch. Batch laws are therefore rejected here if
     * they can return zero, rather than being silently clamped.
     *
     * @param jobClass  the job class to configure
     * @param batchSize the batch-size law, or null to restore single arrivals
     */
    public void setArrivalBatch(JobClass jobClass, jline.lang.processes.DiscreteDistribution batchSize) {
        if (batchSize == null) {
            this.arrivalBatch.remove(jobClass);
            // sn.arrivalbatch is populated only by the full struct build, so a batch set or
            // cleared after anything had built the struct never reached it: the recorder read
            // the node object and emitted BatchArrival while every sn-based predicate saw
            // none. Same defect as Queue.setNumberOfServers, fixed 2026-09-05.
            invalidateStruct();
            return;
        }
        if (batchSize.getMean() < 1.0) {
            throw new IllegalArgumentException("Arrival batch size for class '"
                    + jobClass.getName() + "' has mean " + batchSize.getMean()
                    + "; a batch must carry at least one job, so its support must be {1,2,...}");
        }
        this.arrivalBatch.put(jobClass, batchSize);
        invalidateStruct();
    }

    /**
     * @param jobClass the job class to query
     * @return the batch-size law bound to the class, or null for single arrivals
     */
    public jline.lang.processes.DiscreteDistribution getArrivalBatch(JobClass jobClass) {
        return this.arrivalBatch.get(jobClass);
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

    /**
     * Binds a MarkedMAP with K marks to K open classes: mark k emits jobs of
     * class classes.get(k-1), with all marks driven by one shared modulating
     * chain. Mirrors MATLAB Source.setMarkedArrival.
     *
     * @param mmap    the shared marked arrival process (M3A layout)
     * @param classes the K distinct open classes, ordered by mark index
     */
    public void setMarkedArrival(jline.lang.processes.MarkedMAP mmap, List<JobClass> classes) {
        bindMarkedArrival(mmap, mmap.getNumberOfTypes(), classes);
    }

    /**
     * Binds a time-inhomogeneous marked arrival process with K marks to K open classes.
     *
     * <p>Same contract as the MarkedMAP overload: mark k emits jobs of class
     * classes.get(k-1), all marks share one modulating chain, and only the first marked
     * class carries the arrival stream. The matrices vary with the wall clock, which
     * changes nothing about the binding.
     *
     * @param mmapt   the shared marked schedule
     * @param classes the K distinct open classes, ordered by mark index
     */
    public void setMarkedArrival(jline.lang.processes.MMAPt mmapt, List<JobClass> classes) {
        bindMarkedArrival(mmapt, mmapt.getNumberOfTypes(), classes);
    }

    /**
     * Binds a time-inhomogeneous marked phase-type process with K marks to K open classes.
     *
     * @param mpht    the shared marked schedule
     * @param classes the K distinct open classes, ordered by mark index
     */
    public void setMarkedArrival(jline.lang.processes.MPHt mpht, List<JobClass> classes) {
        bindMarkedArrival(mpht, mpht.getNumberOfTypes(), classes);
    }

    /**
     * Binds a BATCH marked time-inhomogeneous process with K marks to K open classes.
     *
     * <p>Same contract as every other overload, with one addition that matters downstream:
     * an epoch of a BMMAPt releases a BATCH, and the batch is HOMOGENEOUS in its mark, so
     * all of its jobs enter the class bound to that mark. The batch size comes from the
     * process itself, which is why {@code Source.setArrivalBatch} must not also be set on
     * the same class: the two are mutually exclusive by construction, exactly as they are
     * for a BMAP.
     *
     * @param bmmapt  the shared batch marked schedule
     * @param classes the K distinct open classes, ordered by mark index
     */
    public void setMarkedArrival(jline.lang.processes.BMMAPt bmmapt, List<JobClass> classes) {
        bindMarkedArrival(bmmapt, bmmapt.getNumberOfTypes(), classes);
    }

    /** The one binding path behind every setMarkedArrival overload. */
    private void bindMarkedArrival(jline.lang.processes.Distribution mmap, int K,
                                   List<JobClass> classes) {
        if (classes == null || classes.size() != K) {
            throw new IllegalArgumentException(String.format(
                    "The %s has %d types but %d classes were supplied.",
                    mmap.getName(), K, classes == null ? 0 : classes.size()));
        }
        java.util.Set<Integer> seen = new java.util.HashSet<Integer>();
        for (JobClass cls : classes) {
            if (!(cls instanceof jline.lang.OpenClass)) {
                throw new IllegalArgumentException("setMarkedArrival requires open classes.");
            }
            if (!seen.add(cls.getIndex())) {
                throw new IllegalArgumentException("setMarkedArrival requires distinct classes for the marks.");
            }
        }
        // see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        for (JobClass cls : classes) {
            this.setArrival(cls, mmap);
        }
        this.markedProcess = mmap;
        this.markedClasses = new ArrayList<JobClass>(classes);
    }

    /**
     * Drop the marked-arrival binding, leaving whatever per-class arrival
     * processes are currently set.
     * <p>
     * Used where a transformation replaces the shared marked stream by per-class
     * renewal streams (MAPQN2RENV freezes each modulating phase into its own
     * Exp), after which a struct still advertising markidx would claim a marked
     * arrival the model no longer has.
     */
    public void clearMarkedArrival() {
        this.markedProcess = null;
        this.markedClasses = null;
    }

    /**
     * @return the shared marked process bound via setMarkedArrival, or null. It is a
     *         MarkedMAP, an MMAPt or an MPHt.
     */
    public jline.lang.processes.Distribution getMarkedProcess() {
        return this.markedProcess;
    }

    /**
     * @return the classes bound to marks 1..K of the marked process, or null
     */
    public List<JobClass> getMarkedClasses() {
        return this.markedClasses;
    }

}
