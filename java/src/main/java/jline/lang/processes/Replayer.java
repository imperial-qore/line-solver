package jline.lang.processes;

import jline.lang.distributions.APH;
import jline.lang.distributions.Distribution;
import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.io.*;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;

public class Replayer extends Distribution {
    private String fileName;
    double[] data;
    EmpiricalDistribution ecdf;

    public Replayer(Object data) {
        super("Replayer", 1, new Pair<Double, Double>(0.0, Double.POSITIVE_INFINITY));

        if (data instanceof String) {
            this.fileName = (String) data;
                File file = new File(fileName);
            Scanner scan = null;
            try {
                scan = new Scanner(file);
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
            scan.useDelimiter("\n");

            Vector<Double> content = new Vector<Double>();
            while(scan.hasNext())
            {
                content.add(Double.parseDouble(scan.next()));
            }
            this.data = new double[content.size()];
            for (int i = 0; i < content.size(); i++) {
                this.data[i] = content.get(i);
            }
        } else if (data instanceof Matrix) {
            this.data = ((Matrix) data).toArray1D();
        }
        this.ecdf = null;
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
        Mean mean = new Mean();
        double meanValue = mean.evaluate(data);
        return meanValue;
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
        Variance var = new Variance();
        double varValue = var.evaluate(data);
        return varValue;
    }

    // TODO: not implemented yet
    @Override
    public double getSkew() {
        Skewness skewness = new Skewness();
        double skewnessValue = skewness.evaluate(data);
        return skewnessValue;
    }

    public String getFileName() {
        return fileName;
    }

    @Override
    public double evalCDF(double t) {
        if (this.ecdf == null) {
            this.ecdf = new EmpiricalDistribution();
            this.ecdf.load(this.data);
        }
        return this.ecdf.cumulativeProbability(t);
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
