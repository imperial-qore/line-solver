package jline.lang.distributions;

import java.io.Serializable;

import jline.util.Pair;

public abstract class DiscreteDistribution extends Distribution implements Serializable {
    public DiscreteDistribution(String name, int numParam, Pair<Double,Double> support) {
        super(name, numParam, support);
    }
}
