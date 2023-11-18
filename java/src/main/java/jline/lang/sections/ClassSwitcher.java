package jline.lang.sections;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;

import jline.lang.*;
import jline.util.SerializableFunction;
import jline.util.Matrix;

/**
 * A job class switcher based on a probability table
 */
public class ClassSwitcher extends ServiceSection implements Serializable {
	protected SerializableFunction<CSFunInput, Double> csFun;
    protected List<JobClass> jobClasses;
    
    public ClassSwitcher(List<JobClass> jobClasses, String name) {
        super(name);

        this.jobClasses = jobClasses;
        this.numberOfServers = 1;
        this.serviceProcesses = new HashMap<JobClass, ServiceBinding>();
    }
    
    public double applyCsFun(int r, int s) {
    	return this.csFun.apply(new CSFunInput(r, s, null, null));
    }

    public static class CSFunInput {
        /* Currently this class is designed for the usage of statelessclassswitcher
         * state and statedep are not used
         * This class is the input of function (csFun)
         * Currently the csFun is the csMatrix for class switch node
         * r is the row index and s is the column index
         * The return value is the matrix.get(r,s)
         * Now there are only two places using this class.
         * 	One is creating csMask in refreshChains
         * 	Another is obtaining the csMatrix of particular class swicther in getRoutingMatrix
         * If future modification is done, please modify the aforementioned code as well
         */
        public int r;
        public int s;
        public Matrix state;
        public Matrix statedep;

        public CSFunInput(int r, int s, Matrix state, Matrix statedep) {
            this.r = r;
            this.s = s;
            this.state = state;
            this.statedep = statedep;
        }
    }

}
