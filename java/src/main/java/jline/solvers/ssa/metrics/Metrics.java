package jline.solvers.ssa.metrics;

import jline.lang.Network;
import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.state.SSAStateMatrix;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Metrics {
    protected List<Metric> stateMetrics;
    protected List<Metric> eventMetrics;
    protected boolean useMSER5;
    protected boolean useR5;
    protected int r5value;
    protected boolean record;

    public Metrics() {
        this.stateMetrics = new ArrayList<Metric>();
        this.eventMetrics = new ArrayList<Metric>();
        this.useMSER5 = false;
        this.useR5 = false;
        this.r5value = 19;
        this.record = true;
    }

    public void MSER5() {
        this.useMSER5 = true;
        for (Metric metric : this.stateMetrics) {
            metric.configureMSER5();
        }

        for (Metric metric : this.eventMetrics) {
            metric.configureMSER5();
        }
    }

    public void R5(int k) {
        this.useR5 = true;
        this.r5value = k;
        for (Metric metric : this.stateMetrics) {
            metric.configureR5(k);
        }

        for (Metric metric : this.eventMetrics) {
            metric.configureR5(k);
        }
    }

    public void setRecord(boolean record) {
        this.record = record;

        for (Metric metric : this.stateMetrics) {
            metric.setRecord(record);
        }

        for (Metric metric : this.eventMetrics) {
            metric.setRecord(record);
        }
    }

    public void addMetric(Metric metric) {
        if (metric.useMatrix) {
            this.stateMetrics.add(metric);
        } else {
            this.eventMetrics.add(metric);
        }

        if (this.useMSER5) {
            metric.configureMSER5();
        } else if (this.useR5) {
            metric.configureR5(this.r5value);
        }

        if (!this.record) {
            metric.setRecord(false);
        }
    }

    public Metric getMetricByName(String name) {
        for (Metric metric : this.stateMetrics) {
            if (metric.getName().equals(name)) {
                return metric;
            }
        }
        for (Metric metric : this.eventMetrics) {
            if (metric.getName().equals(name)) {
                return metric;
            }
        }

        return new DummyMetric(-1,-1,-1);
    }

    public String getMetricStringByName(String name) {
        Object metric;
        metric = this.getMetricByName(name).getMetric();
        String outString;
        if (metric instanceof Double) {
            outString = String.format("%,.5f", metric);
        } else {
            outString = metric.toString();
        }

        return outString;
    }

    public double getMetricValueByName(String name) {
        Object metric;
        metric = this.getMetricByName(name).getMetric();
        String outString;
        if (metric instanceof Double) {
            return (double)metric;
        } else if (metric instanceof Integer){
            return (double)metric;
        } else if (metric instanceof Long){
            return (double)metric;
        }

        return 0.0;
    }

    public List<String> getAllMetricNames() {
        List<String> outList = new ArrayList<String>();

        for (Metric metric : this.stateMetrics) {
            outList.add(metric.getName());
        }
        for (Metric metric : this.eventMetrics) {
            outList.add(metric.getName());
        }
        return outList;
    }

    private List<Metric> allMetrics() {
        return Stream.concat(this.stateMetrics.stream(), this.eventMetrics.stream()).collect(Collectors.toList());
    }

    public static void outputSummary(Network network, Metrics[][] metricList) {
        System.out.println("Metrics Summary");
        System.out.println("------------------------------------");
        int nStateful = network.getNumberOfStatefulNodes();
        int nClasses = network.getNumberOfClasses();
        List<String> metricNames = new ArrayList<String>();
        List<String> metricShortNames = new ArrayList<String>();
        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (Metric metric : metricList[i][j].allMetrics()) {
                    if (metric == null) {
                        continue;
                    }
                    String metricName = metric.getName();
                    if (!metricNames.contains(metricName)) {
                        metricNames.add(metricName);
                        metricShortNames.add(metric.getShortName());
                    }
                }
            }
        }

        System.out.format("%-10s\t%-10s", "Node", "Class");
        for (String metricName : metricShortNames) {
            System.out.format("\t%s    ", metricName);
        }

        System.out.format("\n");

        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                System.out.format("%-10s\t%-10s", network.getNodeByStatefulIndex(i).getName(),
                        network.getClassByIndex(j).getName());
                for (String metricName : metricNames) {
                    System.out.format("\t%s",metricList[i][j].getMetricStringByName(metricName));
                }
                System.out.format("\n");
            }
        }
    }

    public static void saveSummary(String filename, Network network, Metrics[][] metricList, List<String> additionallines) {
        FileWriter writer;
        try {
            writer = new FileWriter(filename);
            int nStateful = network.getNumberOfStatefulNodes();
            int nClasses = network.getNumberOfClasses();
            List<String> metricNames = new ArrayList<String>();
            List<String> metricShortNames = new ArrayList<String>();
            for (int i = 0; i < nStateful; i++) {
                for (int j = 0; j < nClasses; j++) {
                    for (Metric metric : metricList[i][j].allMetrics()) {
                        String metricName = metric.getName();
                        if (!metricNames.contains(metricName)) {
                            metricNames.add(metricName);
                            metricShortNames.add(metric.getShortName());
                        }
                    }
                }
            }

            writer.append("Node,Class,");
            for (String metricName : metricShortNames) {
                writer.append(metricName);
                writer.append(",");
            }

            writer.append("\n");

            for (int i = 0; i < nStateful; i++) {
                for (int j = 0; j < nClasses; j++) {
                    writer.append(network.getNodeByStatefulIndex(i).getName());
                    writer.append(",");
                    writer.append(network.getClassByIndex(j).getName());
                    writer.append(",");
                    for (String metricName : metricNames) {
                        writer.append(String.format("%f", metricList[i][j].getMetricValueByName(metricName)));
                        writer.append(",");
                    }
                    writer.append("\n");
                }
            }

            for (String additionalline : additionallines) {
                writer.append(additionalline);
                writer.append("\n");
            }

            writer.flush();
            writer.close();
        } catch(java.io.IOException e) {
        }
    }

    public void fromStateMatrix(double t, SSAStateMatrix networkState) {
        for (Metric metric : this.stateMetrics) {
            metric.fromStateMatrix(t, networkState);
        }
    }

    public void fromEvent(double t, Event e) {
        for (Metric metric : this.eventMetrics) {
            metric.fromEvent(t, e);
        }
    }


    public void fromEventN(double t, Event e, int n) {
        for (Metric metric : this.eventMetrics) {
            metric.fromEvent(t, e, n);
        }
    }
    public void taper(double t) {
        for (Metric metric : this.stateMetrics) {
            metric.taper(t);
        }
        for(Metric metric : this.eventMetrics) {
            metric.taper(t);
        }
    }

    public static String getlineHeader(Network network, Metrics[][] metricList) {
        int nStateful = network.getNumberOfStatefulNodes();
        int nClasses = network.getNumberOfClasses();

        List<String> metricNames = new ArrayList<String>();
        List<String> metricShortNames = new ArrayList<String>();
        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (Metric metric : metricList[i][j].allMetrics()) {
                    String metricName = metric.getName();
                    if (!metricNames.contains(metricName)) {
                        metricNames.add(metricName);
                        metricShortNames.add(metric.getShortName());
                    }
                }
            }
        }

        String outString = "";

        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (String metricName : metricShortNames) {
                    outString += network.getNodeByStatefulIndex(i).getName() + "_" + network.getClassByIndex(j).getName() + "_" + metricName + ",";
                }
            }
        }

        return outString.substring(0, outString.length()-1); // remove last comma
    }

    public static String getlineSummary(Network network, Metrics[][] metricList) {
        int nStateful = network.getNumberOfStatefulNodes();
        int nClasses = network.getNumberOfClasses();
        List<String> metricNames = new ArrayList<String>();
        List<String> metricShortNames = new ArrayList<String>();
        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (Metric metric : metricList[i][j].allMetrics()) {
                    String metricName = metric.getName();
                    if (!metricNames.contains(metricName)) {
                        metricNames.add(metricName);
                        metricShortNames.add(metric.getShortName());
                    }
                }
            }
        }

        String outString = network.getName() + ",";

        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (String metricName : metricNames) {
                    outString += String.format("%f", metricList[i][j].getMetricValueByName(metricName)) + ",";
                }
            }
        }
        return outString.substring(0, outString.length()-1); // remove last comma
    }

    public static List<Double> getlineValues(Network network, Metrics[][] metricList) {
        int nStateful = network.getNumberOfStatefulNodes();
        int nClasses = network.getNumberOfClasses();
        List<String> metricNames = new ArrayList<String>();
        List<String> metricShortNames = new ArrayList<String>();
        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (Metric metric : metricList[i][j].allMetrics()) {
                    String metricName = metric.getName();
                    if (!metricNames.contains(metricName)) {
                        metricNames.add(metricName);
                        metricShortNames.add(metric.getShortName());
                    }
                }
            }
        }

        List<Double> outList = new ArrayList<Double>();

        for (int i = 0; i < nStateful; i++) {
            for (int j = 0; j < nClasses; j++) {
                for (String metricName : metricNames) {
                    outList.add(metricList[i][j].getMetricValueByName(metricName));
                }
            }
        }
        return outList;
    }
}
