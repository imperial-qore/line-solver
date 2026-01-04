/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * A section offering a service
 */
public class ServiceSection extends Section implements Serializable {
    protected double numberOfServers;
    protected Map<JobClass, ServiceBinding> serviceProcesses;

    public ServiceSection(String className) {
        super(className);

        this.serviceProcesses = new HashMap<JobClass, ServiceBinding>();
    }

    public boolean containsJobClass(JobClass jobClass) {
        // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
        for (JobClass key : this.serviceProcesses.keySet()) {
            if (key.getIndex() == jobClass.getIndex()) {
                return true;
            }
        }
        return false;
    }

    public Distribution getServiceDistribution(JobClass jobClass) {
        // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
        for (Map.Entry<JobClass, ServiceBinding> entry : this.serviceProcesses.entrySet()) {
            if (entry.getKey().getIndex() == jobClass.getIndex()) {
                return entry.getValue().getDistribution();
            }
        }
        return new Disabled();
    }

    public ServiceBinding getServiceProcess(JobClass jobClass) {
        // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
        for (Map.Entry<JobClass, ServiceBinding> entry : this.serviceProcesses.entrySet()) {
            if (entry.getKey().getIndex() == jobClass.getIndex()) {
                return entry.getValue();
            }
        }
        return null;
    }

    public void removeServiceProcess(JobClass jobClass) {
        // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
        JobClass keyToRemove = null;
        for (JobClass key : this.serviceProcesses.keySet()) {
            if (key.getIndex() == jobClass.getIndex()) {
                keyToRemove = key;
                break;
            }
        }
        if (keyToRemove != null) {
            this.serviceProcesses.remove(keyToRemove);
        }
    }

    public void setServiceProcesses(ServiceBinding serviceProcess) {
        this.serviceProcesses.put(serviceProcess.getJobClass(), serviceProcess);
    }
}