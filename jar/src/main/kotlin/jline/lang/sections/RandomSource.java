/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.ServiceStrategy;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Input buffer of a Source
 */
public class RandomSource extends InputSection implements Serializable {
    protected List<ServiceBinding> serviceProcesses;

    public RandomSource(List<JobClass> jobClasses) {
        super("RandomSource");
        serviceProcesses = new ArrayList<ServiceBinding>();

        for (JobClass jobClass : jobClasses) {
            serviceProcesses.add(new ServiceBinding(jobClass, ServiceStrategy.LI));
        }
    }

    public void removeServiceProcess(JobClass jobClass) {
        Iterator<ServiceBinding> serviceProcessIterator = this.serviceProcesses.iterator();
        while (serviceProcessIterator.hasNext()) {
            // Use index comparison to handle Signal resolution (Signal -> OpenSignal/ClosedSignal)
            if (serviceProcessIterator.next().getJobClass().getIndex() == jobClass.getIndex()) {
                serviceProcessIterator.remove();
            }
        }
    }

    @Override
    public void setServiceProcess(ServiceBinding serviceProcess) {
        removeServiceProcess(serviceProcess.getJobClass());
        serviceProcesses.add(serviceProcess);
    }
}
