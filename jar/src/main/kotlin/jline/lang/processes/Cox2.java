/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * A Coxian distribution with 2 phases.
 */
@SuppressWarnings("unchecked")
public class Cox2 extends Coxian {
    public Cox2(double mu1, double mu2, double phi1) {
        super(Arrays.asList(mu1, mu2), Arrays.asList(phi1, 1.0));
    }

    public static Cox2 fitCentral(double mean, double var, double skew) {
        // Fit the distribution from first three central moments (mean, variance, skewness)
        double scv = var / mean / mean;
        double e1 = mean;
        double e2 = (1 + scv) * FastMath.pow(e1, 2);
        double e3 = -(2 * FastMath.pow(e1, 3) - 3 * e1 * e2 - skew * FastMath.pow(e2 - FastMath.pow(e1, 2), 3.0 / 2.0));

        //consider the two possible solutions
        double phi = (6 * FastMath.pow(e1, 3) - 6 * e2 * e1 + e3) / (-6 * FastMath.pow(e1, 3) + 3 * e2 * e1);
        double mu11 = (2 * (e3 - 3 * e1 * e2)) / (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3) +
                (3 * e1 * e2 - e3 + sqrt(24 * FastMath.pow(e1, 3) * e3 - 27 * FastMath.pow(e1, 2) * FastMath.pow(e2, 2) - 18 * e1 * e2 * e3 + 18 * FastMath.pow(e2, 3) + FastMath.pow(e3, 2))) /
                        (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3);
        double mu12 = (2 * (e3 - 3 * e1 * e2)) / (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3) -
                (e3 - 3 * e1 * e2 + sqrt(24 * FastMath.pow(e1, 3) * e3 - 27 * FastMath.pow(e1, 2) * FastMath.pow(e2, 2) - 18 * e1 * e2 * e3 + 18 * FastMath.pow(e2, 3) + FastMath.pow(e3, 2))) /
                        (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3);
        double mu21 = -(3 * e1 * e2 - e3 + sqrt(24 * FastMath.pow(e1, 3) * e3 - 27 * FastMath.pow(e1, 2) * FastMath.pow(e2, 2) - 18 * e1 * e2 * e3 + 18 * FastMath.pow(e2, 3) + FastMath.pow(e3, 2))) /
                (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3);
        double mu22 = (e3 - 3 * e1 * e2 + sqrt(24 * FastMath.pow(e1, 3) * e3 - 27 * FastMath.pow(e1, 2) * FastMath.pow(e2, 2) - 18 * e1 * e2 * e3 + 18 * FastMath.pow(e2, 3) + FastMath.pow(e3, 2))) /
                (-3 * FastMath.pow(e2, 2) + 2 * e1 * e3);

        Cox2 cx = null;
        if (phi >= 0 && phi <= 1 && mu11 >= 0 && mu21 >= 0) {
            // if the first solution is feasible
            cx = new Cox2(mu11, mu21, phi);
        } else if (phi >= 0 && phi <= 1 && mu12 >= 0 && mu22 >= 0) {
            // if the second solution is feasible
            cx = new Cox2(mu12, mu22, phi);
        } else {
            //line_warning(mfilename, 'Cox2.fitCentral: Third moment could not be fitted exactly.\n');
            // fit is not feasible
            if (scv >= 0.5) {
                // line_warning(mfilename, 'Infeasible combination of central moments, fitting only mean and squared coefficient of variation.');
                cx = Cox2.fitMeanAndSCV(mean, scv);
            } else {
                // line_warning(mfilename, 'Infeasible combination of central moments, fitting only mean.');
                cx = Cox2.fitMean(mean);
            }
        }
        return cx;
    }

    public static Cox2 fitMean(double mean) {
        double phi = 1.0 - GlobalConstants.CoarseTol;
        double l0 = 1.0 / mean;
        double l1 = 1.0 / mean;
        return new Cox2(l0, l1, phi);
    }

    public static Cox2 fitMeanAndSCV(double mean, double SCV) {
        double phi;
        double l0;
        double l1;
        if (SCV < 1 && SCV >= 0.5) {
            l0 = 2 / mean / (1 + sqrt(1 + 2 * (SCV - 1)));
            l1 = 2 / mean / (1 - sqrt(1 + 2 * (SCV - 1)));
            phi = 0.0;
        } else if (SCV == 1.0) {
            l0 = 1.0 / mean;
            l1 = 1.0 / mean;
            phi = 1.0;
        } else {
            l0 = 2.0 / mean;
            l1 = l0 / 2 / SCV;
            phi = 1.0 - l1 / l0;
        }
        return new Cox2(l0, l1, phi);
    }
}