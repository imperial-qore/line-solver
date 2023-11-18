package jline.lang.sections;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import jline.lang.*;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

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

    public ServiceBinding getServiceProcess(JobClass jobClass) {
        return this.serviceProcesses.get(jobClass);
    }

    public Distribution getServiceDistribution(JobClass jobClass) {
        if (!this.serviceProcesses.containsKey(jobClass)) {
            return new Immediate();
        }
        return this.serviceProcesses.get(jobClass).getDistribution();
    }

    public void setServiceProcesses(ServiceBinding serviceProcess) {
        this.serviceProcesses.put(serviceProcess.getJobClass(), serviceProcess);
    }
    
    public boolean containsJobClass(JobClass jobClass) {
    	return this.serviceProcesses.containsKey(jobClass);
    }
    
    public void removeServiceProcess(JobClass jobClass) {
    	this.serviceProcesses.remove(jobClass);
    }
}