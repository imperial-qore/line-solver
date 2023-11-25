package jline.lang.processes;

import jline.lang.distributions.APH;
import jline.lang.distributions.Distribution;
import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.math3.random.EmpiricalDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class Replayer extends Distribution {
    private String fileName;
    EmpiricalDistribution data;

    public Replayer(Object data) {
        super("Replayer", 1, new Pair<Double, Double>(0.0, Double.POSITIVE_INFINITY));
        this.data = new EmpiricalDistribution();
        if (data instanceof String) {
            this.fileName = (String) data;
            try {
                this.data.load(new File(fileName));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } else if (data instanceof Matrix) {
            this.data.load(((Matrix) data).toArray1D());
        }
    }

    @Override
    public Matrix sample(long n) {
        return null;
    }

    @Override
    public Matrix sample(long n, Random random) {
        return null;
    }

    @Override
    public double getMean() {
        return this.data.getNumericalMean();
    }

    @Override
    public double getRate() {
        return 1 / this.getMean();
    }

    @Override
    public double getSCV() {
        double mean = this.getMean();
        double var = this.getVar();
        return var / mean / mean;
    }

    @Override
    public double getVar() {
        return this.data.getNumericalVariance();
    }

    // TODO: not implemented yet
    @Override
    public double getSkew() {
        return Double.NaN;
    }

    @Override
    public double evalCDF(double t) {
        return this.data.cumulativeProbability(t);
    }

    // TODO: not implemented yet
    @Override
    public double evalLST(double s) {
        return Double.NaN;
    }

    // TODO: in MATLAB this also uses the skewness
    public static APH fitAPH(String fileName) {
        Replayer replayer = new Replayer(fileName);
        return APH.fitMeanAndSCV(replayer.getMean(), replayer.getSCV());
    }
}
