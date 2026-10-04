/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Distribution;

import java.io.Serializable;

/**
 * A class for associating job classes, service strategies and distributions
 */
public class ServiceBinding implements Serializable {
    private final JobClass jobClass;
    private final ServiceStrategy serviceStrategy;
    private Distribution distribution;

    public ServiceBinding(JobClass jobClass, ServiceStrategy serviceStrat) {
        this.jobClass = jobClass;
        this.serviceStrategy = serviceStrat;
    }

    public ServiceBinding(JobClass jobClass, ServiceStrategy serviceStrat, Distribution distribution) {
        this(jobClass, serviceStrat);
        this.distribution = distribution;
    }

    public final Distribution getDistribution() {
        return this.distribution;
    }

    public final JobClass getJobClass() {
        return this.jobClass;
    }

    public final ServiceStrategy getServiceStrategy() {
        return this.serviceStrategy;
    }
}
