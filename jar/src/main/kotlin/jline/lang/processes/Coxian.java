/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mam.Map_meanKt;
import jline.api.mam.Map_normalizeKt;
import jline.api.mam.Map_sampleKt;
import jline.api.mam.Map_scvKt;
import jline.GlobalConstants;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A general Coxian distribution with n phases.
 */
@SuppressWarnings("unchecked")
public class Coxian extends Markovian {

    public Coxian(Matrix mu0, Matrix phi0) {
        this(mu0.toList1D(), phi0.toList1D());
    }

    public Coxian(List<Double> mu0, List<Double> phi0) {
        super("Coxian", 1);

        if (!checkParameter(mu0, phi0))
            line_error(mfilename(new Object() {
            }), "Invalid input parameters to Coxian distribution. The exit probability vector must end with value 1.0.");

        this.setParam(1, "mu", mu0);
        this.setParam(2, "phi", phi0);
        nPhases = mu0.size();

        MatrixCell rep = new MatrixCell();
        if (this.getNumberOfPhases() == 2) {
            double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
            double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
            double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
            Matrix matrix1 = new Matrix(2, 2, 4);
            Matrix matrix2 = new Matrix(2, 2, 4);
            matrix1.set(0, 0, -mu1);
            matrix1.set(0, 1, (1 - phi1) * mu1);
            matrix1.set(1, 1, -mu2);
            matrix2.set(0, 0, phi1 * mu1);
            matrix2.set(1, 0, mu2);
            MatrixCell D = Map_normalizeKt.map_normalize(matrix1, matrix2);
            rep.set(0, D.get(0));
            rep.set(1, D.get(1));
        } else {
            Matrix mu = getMu();
            Matrix phi = getPhi();
            //diag(mu(1:end-1).*(1-phi(1:end-1)),1)
            Matrix expression1 = new Matrix(mu.getNumRows(), mu.getNumRows(), mu.getNumRows() - 1);
            for (int i = 0; i < mu.getNumRows() - 1; i++)
                expression1.set(i, i + 1, mu.get(i, 0) * (1 - phi.get(i, 0)));

            //diag(-mu)
            Matrix tmp = mu.copy();
            tmp.changeSign();
            Matrix expression2 = Matrix.diag(tmp.getNonZeroValues());

            //diag(-mu)+diag(mu(1:end-1).*(1-phi(1:end-1)),1)
            Matrix expression3 = expression1.add(1, expression2);

            //phi.*mu
            Matrix expression4 = new Matrix(mu.getNumRows(), 1);
            for (int i = 0; i < mu.getNumRows(); i++)
                expression4.set(i, 0, mu.get(i, 0) * (phi.get(i, 0)));

            //zeros(length(mu),length(mu)-1)
            Matrix expression5 = new Matrix(mu.getNumRows(), mu.getNumRows() - 1);

            //[phi.*mu,zeros(length(mu),length(mu)-1)]
            Matrix expression6 = new Matrix(0, 0, 0);
            Matrix.concatColumns(expression4, expression5, expression6);

            rep = Map_normalizeKt.map_normalize(expression3, expression6);
        }
        setProcess(rep);
    }

    public static Coxian fitCentral(double mean, double var, double skew) {
        Cox2 cx2 = Cox2.fitCentral(mean, var, skew);
        Coxian cx = new Coxian(cx2.getMu(), cx2.getPhi());
        double scv = var / mean / mean;
        if (Math.abs(1 - Map_scvKt.map_scv(cx.D(0), cx.D(1)) / scv) > 0.01) {
            cx = Coxian.fitMeanAndSCV(mean, scv);
        }
        cx.immediate = mean < GlobalConstants.CoarseTol;
        return cx;
    }

    public static Coxian fitMeanAndSCV(double mean, double var, double skew) {
        Cox2 cx2 = Cox2.fitCentral(mean, var, skew);
        Coxian cx = new Coxian(cx2.getMu(), cx2.getPhi());
        double scv = var / mean / mean;
        if (Math.abs(1 - Map_scvKt.map_scv(cx.D(0), cx.D(1)) / scv) > 0.01) {
            cx = Coxian.fitMeanAndSCV(mean, scv);
        }
        cx.immediate = mean < GlobalConstants.CoarseTol;
        return cx;
    }

    // Fit a Coxian distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
    public static Coxian fitMeanAndSCV(double mean, double SCV) {

        double n;
        List<Double> mu = new LinkedList<>();
        List<Double> phi = new LinkedList<>();
        double lambda;

        if (SCV >= 1 - GlobalConstants.CoarseTol && SCV <= 1 + GlobalConstants.CoarseTol) {
            n = 1;
            mu.add(1.0 / mean);
            phi.add(1.0);
        } else if (SCV > 0.5 + GlobalConstants.CoarseTol && SCV < 1 - GlobalConstants.CoarseTol) {
            phi.add(0.0);
            phi.add(0.0);
            n = 2;
            mu.add(2 / mean / (1 + FastMath.sqrt(1 + (2 * (SCV - 1)))));
            mu.add(2 / mean / (1 - FastMath.sqrt(1 + (2 * (SCV - 1)))));
        } else if (SCV <= 0.5 + GlobalConstants.CoarseTol) {
            n = FastMath.ceil(1.0 / SCV);
            lambda = n / mean;
            for (int i = 0; i < n; i++) {
                mu.add(lambda);
                phi.add(0.0);
            }
        } else { // SCV > 1 + GlobalConstants.CoarseTol
            n = 2;
            // transform hyperexp into coxian
            mu.add(2 / mean);
            mu.add((2 / mean) / (2 * SCV));
            phi.add(1 - (mu.get(1) / mu.get(0)));
            phi.add(1.0);
        }

        phi.set((int) n - 1, 1.0);
        Coxian cx = new Coxian(mu, phi);
        cx.immediate = mean < GlobalConstants.CoarseTol;
        return cx;
    }

    private boolean checkParameter(List<Double> mu, List<Double> phi) {
        if (mu.size() != phi.size()) {
            return false;
        } else {
            return !((Math.abs(phi.get(phi.size() - 1) - 1) > GlobalConstants.CoarseTol) && !Utils.isInf(phi.get(phi.size() - 1)) && !Double.isNaN(phi.get(phi.size() - 1)));
        }
    }

    @Override
    public double evalCDF(double t) {
        return super.evalCDF(t);
    }

    @Override
    public double evalLST(double s) {
        return super.evalLST(s);
    }

    @Override
    public double getMean() {
        if (this.getNumberOfPhases() == 2) {
            double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
            double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
            double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
            return 1 / mu1 + (1 - phi1) / mu2;
        } else {
            MatrixCell rep = this.getProcess();
            return Map_meanKt.map_mean(rep.get(0), rep.get(1));
        }
    }

    public Matrix getMu() {
        List<Double> mu = (List<Double>) this.getParam(1).getValue();
        Matrix res = new Matrix(mu.size(), 1, mu.size());
        if (this.getNumberOfPhases() == 2) {
            res.set(0, 0, mu.get(0));
            res.set(1, 0, mu.get(1));
        } else {
            for (int i = 0; i < mu.size(); i++)
                res.set(i, 0, mu.get(i));
        }
        return res;
    }

    @Override
    public long getNumberOfPhases() {
        return ((List<Double>) this.getParam(1).getValue()).size();
    }

    public Matrix getPhi() {
        List<Double> phi = (List<Double>) this.getParam(2).getValue();
        Matrix res = new Matrix((int) this.getNumberOfPhases(), 1, (int) this.getNumberOfPhases());
        if (this.getNumberOfPhases() == 2) {
            res.set(0, 0, phi.get(0));
            res.set(1, 0, 1);
        } else {
            for (int i = 0; i < phi.size(); i++)
                res.set(i, 0, phi.get(i));
        }
        return res;
    }

    @Override
    public double getRate() {
        return 1 / getMean();
    }

    @Override
    public double getSCV() {
        if (this.getNumberOfPhases() == 2) {
            double mu1 = ((List<Double>) this.getParam(1).getValue()).get(0);
            double mu2 = ((List<Double>) this.getParam(1).getValue()).get(1);
            double phi1 = ((List<Double>) this.getParam(2).getValue()).get(0);
            double mean = 1 / mu1 + (1 - phi1) / mu2;
            double var = ((2 * mu2 * (mu1 - mu1 * phi1)) / (mu1 + mu2 - mu1 * phi1) + (2 * mu1 * mu2 * phi1) / (mu1 + mu2 - mu1 * phi1)) / (mu1 * mu1 * ((mu2 * (mu1 - mu1 * phi1)) / (mu1 + mu2 - mu1 * phi1) + (mu1 * mu2 * phi1) / (mu1 + mu2 - mu1 * phi1))) - (1 / mu1 - (phi1 - 1) / mu2) * (1 / mu1 - (phi1 - 1) / mu2) - (((phi1 - 1) / (mu2 * mu2) + (phi1 - 1) / (mu1 * mu2)) * ((2 * mu2 * (mu1 - mu1 * phi1)) / (mu1 + mu2 - mu1 * phi1) + (2 * mu1 * mu2 * phi1) / (mu1 + mu2 - mu1 * phi1))) / ((mu2 * (mu1 - mu1 * phi1)) / (mu1 + mu2 - mu1 * phi1) + (mu1 * mu2 * phi1) / (mu1 + mu2 - mu1 * phi1));
            return var / FastMath.pow(mean, 2);
        } else {
            MatrixCell rep = this.getProcess();
            Matrix D0 = rep.get(0);
            Matrix D1 = rep.get(1);
            MAP map = new MAP(D0, D1);
            return map.getSCV();
        }
    }

    @Override
    public double getSkewness() {
        return super.getSkewness();
    }

    @Override
    public double getVar() {
        return this.getSCV() * Math.pow(this.getMean(), 2);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, null);
    }

    @Override
    public double[] sample(int n, Random random) {
        return Map_sampleKt.map_sample(D(0), D(1), n, random);
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getMean()
     */
    public double mean() {
        return getMean();
    }
    
    /**
     * Kotlin-style property alias for getRate()
     */
    public double rate() {
        return getRate();
    }
    
    /**
     * Kotlin-style property alias for getSCV()
     */
    public double scv() {
        return getSCV();
    }
    
    /**
     * Kotlin-style property alias for getSkewness()
     */
    public double skewness() {
        return getSkewness();
    }
    
    /**
     * Kotlin-style property alias for getVar()
     */
    public double var() {
        return getVar();
    }
    
    /**
     * Kotlin-style property alias for getMu()
     */
    public Matrix mu() {
        return getMu();
    }
    
    /**
     * Kotlin-style property alias for getPhi()
     */
    public Matrix phi() {
        return getPhi();
    }
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numberOfPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Kotlin-style property alias for getNumberOfPhases()
     */
    public long numPhases() {
        return getNumberOfPhases();
    }
}
