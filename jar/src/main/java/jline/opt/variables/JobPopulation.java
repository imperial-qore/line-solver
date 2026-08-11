package jline.opt.variables;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;

/**
 * Optimize the fixed circulating population of a closed class. Encodes an
 * integer job count in [minJobs, maxJobs]. Mirrors native-Python
 * {@code JobPopulation}.
 */
public class JobPopulation extends DecisionVariable {

    private final ClosedClass jobclass;
    private final int minJobs;
    private final int maxJobs;

    public JobPopulation(ClosedClass jobclass, int minJobs, int maxJobs) {
        this(jobclass, minJobs, maxJobs, jobclass.getName() + "_population");
    }

    public JobPopulation(ClosedClass jobclass, int minJobs, int maxJobs, String name) {
        super(name);
        this.jobclass = jobclass;
        this.minJobs = minJobs;
        this.maxJobs = maxJobs;
    }

    public ClosedClass getJobClass() {
        return jobclass;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        double continuous = minJobs + x[0] * (maxJobs - minJobs);
        int v = (int) Math.round(continuous);
        if (v < minJobs) {
            v = minJobs;
        }
        if (v > maxJobs) {
            v = maxJobs;
        }
        return v;
    }

    public void apply(Network model, Object value) {
        int v = ((Number) value).intValue();
        JobClass cls = resolveClass(model, jobclass);
        if (cls instanceof ClosedClass) {
            ((ClosedClass) cls).setPopulation(v);
        }
    }

    public String getVariableType() {
        return "job_population";
    }
}
