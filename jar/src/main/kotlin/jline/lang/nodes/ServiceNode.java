/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * A station with a service process
 */
public abstract class ServiceNode extends StatefulNode implements Serializable {
    protected List<ServiceBinding> serviceProcesses;
    protected SchedStrategy schedStrategy;
    protected SchedStrategyType schedPolicy;
    protected HashMap<JobClass, Double> schedStrategyPar;


    public ServiceNode(String name) {
        super(name);
        this.serviceProcesses = new ArrayList<ServiceBinding>();
        this.schedStrategyPar = new HashMap<JobClass, Double>();
    }

    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    public Distribution getServiceProcess(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            // Use index comparison instead of object identity to handle Signal resolution
            // (Signal -> OpenSignal/ClosedSignal creates new object with same index)
            if (serviceProcess.getJobClass().getIndex() == jobClass.getIndex()) {
                return serviceProcess.getDistribution();
            }
        }

        return new Disabled();
    }

    private boolean hasJobClass(JobClass jobClass) {
        return this.schedStrategyPar.containsKey(jobClass);
    }

    private void removeServiceProcess(JobClass jobClass) {
        this.schedStrategyPar.remove(jobClass);
        this.serviceProcesses.removeIf(serviceProcess -> serviceProcess.getJobClass().getIndex() == jobClass.getIndex());
        this.server.removeServiceProcess(jobClass);
    }

    public void setService(JobClass jobClass, Distribution distribution) {
        setService(jobClass, distribution, 1);
    }

    public void setService(JobClass jobClass, Distribution distribution, double weight) {
        boolean resetState = false;
        if (this.hasJobClass(jobClass)) {
            resetState = true;
            this.removeServiceProcess(jobClass);
            //clearState();
        }

        if (distribution.isImmediate()) {
            this.serviceProcesses.add(new ServiceBinding(jobClass, ServiceStrategy.LI, Immediate.getInstance()));
            if (this.server != null) {
                this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, Immediate.getInstance()));
            }
        } else {
            this.serviceProcesses.add(new ServiceBinding(jobClass, ServiceStrategy.LI, distribution));
            if (this.server != null) {
                this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, distribution));
            }
        }

        this.schedStrategyPar.put(jobClass, weight);

        if (resetState) {
            this.model.setInitialized(false);
        }

    }
}
