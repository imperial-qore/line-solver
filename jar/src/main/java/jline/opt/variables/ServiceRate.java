package jline.opt.variables;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.nodes.Node;
import jline.lang.nodes.ServiceStation;
import jline.lang.processes.Exp;
import jline.opt.results.SensitivityData;

/**
 * Optimize the exponential processing rate of a station for a job class.
 * Encodes a continuous rate in [minRate, maxRate]. This is a continuous
 * (differentiable) variable, so it exposes {@link #paramKey()} and
 * {@link #decodeJacobian} for the analytic-gradient path. Mirrors native-Python
 * {@code ServiceRate}.
 */
public class ServiceRate extends DecisionVariable {

    private final ServiceStation station;
    private final JobClass jobclass;
    private final double minRate;
    private final double maxRate;

    public ServiceRate(ServiceStation station, JobClass jobclass, double minRate, double maxRate) {
        this(station, jobclass, minRate, maxRate,
                station.getName() + "_" + jobclass.getName() + "_rate");
    }

    public ServiceRate(ServiceStation station, JobClass jobclass, double minRate,
                       double maxRate, String name) {
        super(name);
        this.station = station;
        this.jobclass = jobclass;
        this.minRate = minRate;
        this.maxRate = maxRate;
    }

    public ServiceStation getStation() {
        return station;
    }

    public JobClass getJobClass() {
        return jobclass;
    }

    public double[][] getBounds() {
        return unitBounds(1);
    }

    public Object decode(double[] x) {
        return minRate + x[0] * (maxRate - minRate);
    }

    public void apply(Network model, Object value) {
        double rate = ((Number) value).doubleValue();
        JobClass cls = resolveClass(model, jobclass);
        if (cls == null) {
            cls = jobclass;
        }
        for (Node node : model.getNodes()) {
            if (node.getName().equals(station.getName()) && node instanceof ServiceStation) {
                ((ServiceStation) node).setService(cls, new Exp(rate));
                return;
            }
        }
    }

    public String getVariableType() {
        return "service_rate";
    }

    /** Sensitivity parameter key this variable controls. */
    public String paramKey() {
        return SensitivityData.paramKey(station.getName(), jobclass.getName());
    }

    /** d(decoded rate)/d(encoded x): constant slope of the linear map. */
    public double decodeJacobian(double[] x) {
        return maxRate - minRate;
    }
}
