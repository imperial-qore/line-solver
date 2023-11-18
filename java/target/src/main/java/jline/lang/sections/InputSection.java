package jline.lang.sections;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SchedStrategyType;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

/**
 * Input section of a station
 */
public class InputSection extends Section  implements Serializable {
    protected SchedStrategyType schedPolicy;
    protected List<InputBinding> inputJobProcesses;
    protected List<SchedStrategy> inputJobClasses;

    public InputSection(String className) {
        super(className);
        this.inputJobProcesses = new ArrayList<InputBinding>();
        this.inputJobClasses = new ArrayList<SchedStrategy>();
    }

    public void setInputJobProcess(InputBinding process) {
        removeInputJobProcess(process.getJobClass());
        inputJobProcesses.add(process);
    }

    protected void removeInputJobProcess(JobClass jobClass) {
        Iterator<InputBinding> inputJobProcessIter = this.inputJobProcesses.iterator();
        while (inputJobProcessIter.hasNext()) {
            if (inputJobProcessIter.next().getJobClass() == jobClass) {
                inputJobProcessIter.remove();
            }
        }
    }

    public void setInputJobClasses(int index, SchedStrategy schedStrategy) {
        int diff = index - this.inputJobClasses.size() + 1;
        while (diff > 0){
            this.inputJobClasses.add(null);
            diff--;
        }
        this.inputJobClasses.set(index, schedStrategy);
    }


    public void setServiceProcess(ServiceBinding serviceProcess) {

    }
}
