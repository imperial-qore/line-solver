package jline.lang.nodes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ServiceStrategy;
import jline.lang.distributions.*;
import jline.lang.sections.*;

/**
 * An abstraction of the external world jobs in open classes come from
 */
public class Source extends Station implements HasSchedStrategy, Serializable {
    protected List<ServiceBinding> serviceProcesses;
    protected SchedStrategy schedStrategy;

    public Source(String name) {
        super(name);
        this.numberOfServers = 1;
    }
    public Source(Network model, String name) {
        this(name);

        List<JobClass> jobClasses = model.getClasses();
        this.output = new Dispatcher(jobClasses);
        this.server = new ServiceTunnel();
        this.input = new RandomSource(jobClasses);
        this.schedStrategy = SchedStrategy.EXT;
        this.setModel(model);
        this.model.addNode(this);
        this.serviceProcesses = new ArrayList<ServiceBinding>();

        for (JobClass jobClass : jobClasses) {
            this.classCap.put(jobClass, Double.POSITIVE_INFINITY);
            this.setArrival(jobClass, new Disabled());
        }
    }

    public void setArrival(JobClass jobClass, Distribution distribution) {
        ServiceBinding arrivalProcess = new ServiceBinding(jobClass, ServiceStrategy.LI, distribution);
        //this.input.setServiceProcess(arrivalProcess);
        this.setServiceProcess(arrivalProcess);
        if ((distribution == null) || (distribution instanceof Disabled)){
            this.classCap.put(jobClass, 0.0);
        } else {
            this.classCap.put(jobClass, Double.POSITIVE_INFINITY);
        }
    }

    protected void removeServiceProcess(JobClass jobClass) {
        Iterator<ServiceBinding> serviceProcessIterator = this.serviceProcesses.iterator();
        while (serviceProcessIterator.hasNext()) {
            if (serviceProcessIterator.next().getJobClass() == jobClass) {
                serviceProcessIterator.remove();
            }
        }
        ((RandomSource)this.input).removeServiceProcess(jobClass);
    }

    public void setServiceProcess(ServiceBinding arrivalProcess) {
        removeServiceProcess(arrivalProcess.getJobClass());
        serviceProcesses.add(arrivalProcess);
    }

    public final Distribution getServiceProcess(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            if (serviceProcess.getJobClass() == jobClass) {
                return serviceProcess.getDistribution();
            }
        }

        return new Disabled();
    }

    public Distribution getArrivalDistribution(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            if (serviceProcess.getJobClass() == jobClass) {
                return serviceProcess.getDistribution();
            }
        }
        return new Disabled();
    }
    
    public boolean containsJobClass(JobClass jobClass) {
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            if (serviceProcess.getJobClass() == jobClass) {
                return true;
            }
        }
        return false;
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
    public SchedStrategy getSchedStrategy() {
        return this.schedStrategy;
    }

    public double minRate() {
        double acc = Double.POSITIVE_INFINITY;
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            double dRate = serviceProcess.getDistribution().getRate();
            if (dRate != 0) {
                acc = Math.min(acc, dRate);
            }
        }
        return acc;
    }
    public double maxRate() {
        double acc = 0;
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            if (serviceProcess.getDistribution().getRate() == Double.POSITIVE_INFINITY) {
                continue;
            }
            acc = Math.max(acc, serviceProcess.getDistribution().getRate());
        }
        return acc;
    }
    public double avgRate() {
        double acc = 0;
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            double accVal = serviceProcess.getDistribution().getRate();
            if ((accVal == Double.POSITIVE_INFINITY) || (accVal == 0)) {
                continue;
            }
            acc += accVal;
        }
        return acc;
    }

    public int rateCt() {
        int acc = 0;
        for (ServiceBinding serviceProcess : this.serviceProcesses) {
            double accVal = serviceProcess.getDistribution().getRate();
            if ((accVal == Double.POSITIVE_INFINITY) || (accVal == 0)) {
                continue;
            }
            acc++;
        }
        return acc;
    }
}
