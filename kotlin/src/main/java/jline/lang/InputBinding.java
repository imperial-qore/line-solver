package jline.lang;

import java.io.Serializable;

import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategyType;

public class InputBinding implements Serializable {
    protected JobClass jobClass;
    protected SchedStrategyType schedPolicy;
    protected DropStrategy dropStrategy;
    public InputBinding(JobClass jobClass, SchedStrategyType schedPolicy, DropStrategy dropStrategy) {
        this.jobClass = jobClass;
        this.schedPolicy = schedPolicy;
        this.dropStrategy = dropStrategy;
    }

    public final JobClass getJobClass() {
        return this.jobClass;
    }
}
