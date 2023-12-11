package jline.lang.distributions;

import java.io.Serializable;

import jline.util.Pair;

abstract public class ContinuousDistribution extends Distribution implements Serializable {
    public ContinuousDistribution(String name, int numParam, Pair<Double,Double> support) {
        super(name, numParam, support);
    }
}
