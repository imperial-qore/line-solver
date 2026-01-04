/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.nodes;

import jline.lang.JobClass;
import jline.lang.SelfLoopingClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ServiceStrategy;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.processes.Immediate;
import jline.lang.workflow.Workflow;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static jline.io.InputOutputKt.line_warning;
import static jline.io.InputOutputKt.mfilename;

/**
 * A station with a service process
 */
public abstract class ServiceStation extends Station implements Serializable {
    protected List<ServiceBinding> serviceProcesses;

    public ServiceStation(String name) {
        super(name);
        this.serviceProcesses = new ArrayList<ServiceBinding>();
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
        this.classCap.remove(jobClass);
        this.schedStrategyPar.remove(jobClass);
        this.serviceProcesses.removeIf(serviceProcess -> serviceProcess.getJobClass().getIndex() == jobClass.getIndex());
        this.server.removeServiceProcess(jobClass);
    }

    public void setService(JobClass jobClass, Distribution distribution) {
        setService(jobClass, (Object) distribution, 1);
    }

    public void setService(JobClass jobClass, Distribution distribution, double weight) {
        setService(jobClass, (Object) distribution, weight);
    }

    public void setService(JobClass jobClass, Object distribution, double weight) {
        // If Workflow, convert to PH distribution
        Distribution actualDistribution;
        if (distribution instanceof Workflow) {
            actualDistribution = ((Workflow) distribution).toPH();
        } else if (distribution instanceof Distribution) {
            actualDistribution = (Distribution) distribution;
        } else {
            throw new IllegalArgumentException("distribution must be a Distribution or Workflow object");
        }

        boolean resetState = false;
        if (this.hasJobClass(jobClass)) {
            resetState = true;
            this.removeServiceProcess(jobClass);
            //clearState();
        }

        if (jobClass instanceof SelfLoopingClass && jobClass.getReferenceStation().getNodeIndex() != this.getNodeIndex() && !(actualDistribution instanceof Disabled)) {
            line_warning(mfilename(new Object() {
            }), "For a self-looping class, service cannot be set on stations other than the reference station of that class.");
        }

        if (actualDistribution.isImmediate()) {
            this.serviceProcesses.add(new ServiceBinding(jobClass, ServiceStrategy.LI, Immediate.getInstance()));
            if (this.server != null) {
                this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, Immediate.getInstance()));
            }
        } else {
            this.serviceProcesses.add(new ServiceBinding(jobClass, ServiceStrategy.LI, actualDistribution));
            if (this.server != null) {
                this.server.setServiceProcesses(new ServiceBinding(jobClass, ServiceStrategy.LI, actualDistribution));
            }
        }

        this.classCap.put(jobClass, Integer.MAX_VALUE);
        this.schedStrategyPar.put(jobClass, weight);

        if (resetState) {
            this.model.setInitialized(false);
        }

        // Default to Drop for finite capacity, WaitingQueue otherwise
        if (this.cap < Integer.MAX_VALUE && !Double.isInfinite(this.cap)) {
            this.dropRule.put(jobClass, DropStrategy.Drop);
        } else {
            this.dropRule.put(jobClass, DropStrategy.WaitingQueue);
        }
    }
}
