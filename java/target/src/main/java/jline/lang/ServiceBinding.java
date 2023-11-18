package jline.lang;

import java.io.Serializable;

import jline.lang.distributions.Distribution;
import jline.lang.*;
import jline.lang.constant.ServiceStrategy;

public class ServiceBinding implements Serializable {
    private final JobClass jobClass;
    private final ServiceStrategy serviceStrategy;
    private Distribution distribution;
    
    public ServiceBinding(JobClass jobClass, ServiceStrategy serviceStrategy) {
        this.jobClass = jobClass;
        this.serviceStrategy = serviceStrategy;
    }
    public ServiceBinding(JobClass jobClass, ServiceStrategy serviceStrategy, Distribution distribution) {
        this (jobClass, serviceStrategy);
        this.distribution = distribution;
    }


    public final JobClass getJobClass() {
        return this.jobClass;
    }
    public final Distribution getDistribution() { return this.distribution; }
}
