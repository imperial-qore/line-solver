package jline.lang.sections;

import jline.lang.JobClass;
import jline.lang.ServiceBinding;
import jline.lang.constant.ServiceStrategy;
import jline.lang.distributions.Exp;

import java.io.Serializable;
import java.util.List;

/**
 * A preemptive service section
 */
public class PreemptiveServer extends ServiceSection implements Serializable {
    public PreemptiveServer(List<JobClass> customerClasses) {
        super("PreemptiveServer");
        this.numberOfServers = 1;
    }

    private void initServers(List<JobClass> customerClasses){
        for (JobClass jobClass : customerClasses) {
            this.serviceProcesses.put(jobClass, new ServiceBinding(jobClass, ServiceStrategy.LI, new Exp(0)));
        }
    }
}
