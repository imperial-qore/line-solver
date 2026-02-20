/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;

/**
 * A distribution that replays empirical data values from trace files.
 * 
 * <p>The Replayer distribution reads real-world trace data from files and
 * reproduces the exact sequence of observations. This is essential for
 * trace-driven simulation, workload characterization, and model validation
 * against real data patterns.</p>
 * 
 * <p>Key capabilities:
 * <ul>
 *   <li>Reading trace data from text files (one value per line)</li>
 *   <li>Sequential replay of empirical observations</li>
 *   <li>Empirical CDF computation from trace data</li>
 *   <li>Statistical moment calculation (mean, variance, skewness)</li>
 *   <li>Support for both file paths and direct data arrays</li>
 * </ul>
 * </p>
 * 
 * <p>Common applications include reproducing real arrival processes,
 * service time patterns, and validating analytical models against
 * measurements from production systems.</p>
 * 
 * @see ReplayerNodeParam
 * @see EmpiricalCDF
 * @since 1.0
 */
public class Replayer extends Distribution {
    double[] data;
    EmpiricalDistribution ecdf;
    private String fileName;

    public Replayer(Object data) {
        super("Replayer", 1, new Pair<Double, Double>(0.0, Inf));

        if (data instanceof String) {
            this.fileName = (String) data;
            File file = new File(fileName);
            Scanner scan;
            try {
                scan = new Scanner(file);
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
            Vector<Double> content = new Vector<Double>();
            while (scan.hasNextLine()) {
                String line = scan.nextLine().trim();
                if (!line.isEmpty()) {
                    try {
                        content.add(Double.parseDouble(line));
                    } catch (NumberFormatException e) {
                        throw new RuntimeException("Invalid number format in trace file '" + fileName + "': " + line, e);
                    }
                }
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
    public double evalCDF(double t) {
        if (this.ecdf == null) {
            this.ecdf = new EmpiricalDistribution();
            this.ecdf.load(this.data);
        }
        return this.ecdf.cumulativeProbability(t);
    }

    @Override
    public double evalLST(double s) {
        // For empirical data, compute LST numerically as the average of e^(-s*x) over all data points
        double sum = 0.0;
        for (double x : data) {
            sum += Math.exp(-s * x);
        }
        return sum / data.length;
    }

    public APH fitAPH() {
        return APH.fit(this.getMean(), this.getSCV(), this.getSkewness());
    }

    public String getFileName() {
        return fileName;
    }

    @Override
    public double getMean() {
        Mean mean = new Mean();
        double meanValue = mean.evaluate(data);
        return meanValue;
    }

    @Override
    public double getRate() {
        return 1.0 / this.getMean();
    }

    @Override
    public double getSCV() {
        double mean = this.getMean();
        double var = this.getVar();
        return var / mean / mean;
    }

    @Override
    public double getSkewness() {
        // This uses the same bias of matlab if the latter calls skewness(data,0)
        Skewness skewness = new Skewness();
        double skewnessValue = skewness.evaluate(data);
        return skewnessValue;
    }

    @Override
    public double getVar() {
        Variance var = new Variance();
        double varValue = var.evaluate(data);
        return varValue;
    }

    @Override
    public double[] sample(int n) {
        return null;
    }

    @Override
    public double[] sample(int n, Random random) {
        return null;
    }

}