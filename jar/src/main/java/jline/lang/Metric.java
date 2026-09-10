/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.MetricType;
import jline.lang.nodes.Station;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

/**
 * Constants for specifying a Metric
 */
public class Metric {
    public String type;
    public Station station;
    public JobClass jobClass;
    public boolean isDisabled;
    private double alfa;
    private int analyzedSamples;
    private String className;
    private int discardedSamples;
    private double lowerLimit;
    private int maxSamples;
    private double meanValue;
    private MetricType metricType;
    private String nodeType;
    private double precision;
    private String stationName;
    private boolean successful;
    private double upperLimit;

    public Metric() {

    }

    public Metric(NamedNodeMap attributes) {
        for (int j = 0; j < attributes.getLength(); j++) {
            Node istAttr = attributes.item(j);
            String nodeName = istAttr.getNodeName();
            String nodeValue = istAttr.getNodeValue();
            switch (nodeName) {
                case "alfa":
                    this.alfa = Double.parseDouble(nodeValue);
                    break;
                case "analyzedSamples":
                    this.analyzedSamples = Integer.parseInt(nodeValue);
                    break;
                case "class":
                    this.className = nodeValue;
                    break;
                case "discardedSamples":
                    this.discardedSamples = Integer.parseInt(nodeValue);
                    break;
                case "lowerLimit":
                    this.lowerLimit = Double.parseDouble(nodeValue);
                    break;
                case "maxSamples":
                    this.maxSamples = Integer.parseInt(nodeValue);
                    break;
                case "meanValue":
                    this.meanValue = Double.parseDouble(nodeValue);
                    break;
                case "measureType":
                    this.metricType = MetricType.toMetricType(nodeValue);
                    break;
                case "nodeType":
                    this.nodeType = nodeValue;
                    break;
                case "precision":
                    this.precision = Double.parseDouble(nodeValue);
                    break;
                case "station":
                    this.stationName = nodeValue;
                    break;
                case "successful":
                    this.successful = Boolean.parseBoolean(nodeValue);
                    break;
                case "upperLimit":
                    this.upperLimit = Double.parseDouble(nodeValue);
                    break;
            }
        }
    }

    public double getAlfa() {
        return alfa;
    }

    public int getAnalyzedSamples() {
        return analyzedSamples;
    }

    public String getClassName() {
        return className;
    }

    public int getDiscardedSamples() {
        return discardedSamples;
    }

    public double getLowerLimit() {
        return lowerLimit;
    }

    public int getMaxSamples() {
        return maxSamples;
    }

    public double getMeanValue() {
        return meanValue;
    }

    public MetricType getMetricType() {
        return metricType;
    }

    public String getNodeType() {
        return nodeType;
    }

    public double getPrecision() {
        return precision;
    }

    public String getStationName() {
        return stationName;
    }

    public double getUpperLimit() {
        return upperLimit;
    }

    public boolean isSuccessful() {
        return successful;
    }

    public void setClassName(String className) {
        this.className = className;
    }

    public void setStationName(String stationName) {
        this.stationName = stationName;
    }

    /**
     * Overwrites the mean. Needed by the JMVA reader, which must rescale a
     * chain-level measure into the per-class value it stands for.
     */
    public void setMeanValue(double meanValue) {
        this.meanValue = meanValue;
        this.lowerLimit = meanValue;
        this.upperLimit = meanValue;
    }

    /**
     * Overwrites the metric type. The JMVA reader uses it to restate JMVA's
     * "Residence time" as a per-visit response time once it has divided the
     * value by the visit count.
     */
    public void setMetricType(MetricType metricType) {
        this.metricType = metricType;
    }

    /**
     * Sets the analysed-sample count. JMVA is analytical and reports no sample
     * counts, so its reader must stamp the saturated value: the consumer in
     * SolverJMT.getResults() drops any closed-class metric whose
     * analyzedSamples does not exceed the chain population, a JSIM
     * sample-sufficiency test that would otherwise discard every analytical
     * result. MATLAB writes Inf here for the same reason.
     */
    public void setAnalyzedSamples(int analyzedSamples) {
        this.analyzedSamples = analyzedSamples;
    }

    /** Sets the node type ("station" or "region"). */
    public void setNodeType(String nodeType) {
        this.nodeType = nodeType;
    }
}

