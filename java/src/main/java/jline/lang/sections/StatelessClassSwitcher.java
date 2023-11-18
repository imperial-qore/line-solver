package jline.lang.sections;

import java.io.Serializable;
import java.util.List;
import java.util.function.Function;

import jline.lang.*;
import jline.util.Matrix;

/**
 * A class switcher that does not keep a local state
 */
public class StatelessClassSwitcher extends ClassSwitcher implements Serializable {
    public StatelessClassSwitcher(List<JobClass> jobClasses, Matrix csMatrix) {
        super(jobClasses, "StatelessClassSwitcher");
        
        this.csFun = (input) -> csMatrix.get(input.r, input.s);
    }

    public void updateClasses(List<JobClass>  jobClasses) {
        this.jobClasses = jobClasses;
    }
    
    public void updateClassSwitch(Matrix csMatrix) {
    	this.csFun = (input) -> csMatrix.get(input.r, input.s);
    }
    public List<JobClass> getJobClasses() {
        return this.jobClasses;
    }
}
