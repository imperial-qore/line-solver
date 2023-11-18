package jline.lang.sections;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.ServiceStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

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
    
    @Override
    public void setServiceProcess(ServiceBinding serviceProcess) {
        removeServiceProcess(serviceProcess.getJobClass());
        serviceProcesses.add(serviceProcess);
    }

    public void removeServiceProcess(JobClass jobClass) {
        Iterator<ServiceBinding> serviceProcessIterator = this.serviceProcesses.iterator();
        while (serviceProcessIterator.hasNext()) {
            if (serviceProcessIterator.next().getJobClass() == jobClass) {
                serviceProcessIterator.remove();
            }
        }
    }
}
