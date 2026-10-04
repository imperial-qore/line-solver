/**
 * @file Continuous-time Markov chain transient analysis
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import odesolver.LSODA;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Ctmc_transient {
    private Ctmc_transient() {}

    public static Pair<double[], List<double[]>> ctmc_transient(Matrix Q, Matrix pi0, double t0, double t1) {
        return ctmc_transient(Q, pi0, t0, t1, null);
    }

    public static Pair<double[], List<double[]>> ctmc_transient(final Matrix Q,
                                                                final Matrix pi0,
                                                                double t0,
                                                                double t1,
                                                                Double timestep) {
        LSODA lsoda = new LSODA(0.0, 0.0, 1.0e-6, 1.0e-6, 12, 5);
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return pi0.length();
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] ydot)
                    throws MaxCountExceededException, DimensionMismatchException {
                for (int i = 0; i < Q.getNumCols(); i++) {
                    for (int j = 0; j < Q.getNumRows(); j++) {
                        ydot[i] += y[j] * Q.get(j, i);
                    }
                }
            }
        };
        double[] y = new double[pi0.length()];

        if (timestep != null && timestep > 0) {
            // grid built on an index, never by accumulating the step: repeated addition
            // drifts and made the loop bound miss or duplicate the endpoint
            List<Double> timePoints = new ArrayList<Double>();
            int nsteps = (int) Math.floor((t1 - t0) / timestep + 1e-9);
            for (int i = 0; i <= nsteps; i++) {
                timePoints.add(t0 + i * timestep);
            }
            if (timePoints.get(timePoints.size() - 1) < t1) {
                timePoints.add(t1);
            }

            // t0 carries pi0 itself; integrating over a zero-length first interval is
            // what threw inside LSODA
            List<double[]> piResults = new ArrayList<double[]>();
            double[] cur = pi0.toArray1D();
            piResults.add(cur.clone());
            for (int i = 1; i < timePoints.size(); i++) {
                lsoda.integrate(ode, timePoints.get(i - 1), cur, timePoints.get(i), y);
                cur = y.clone();
                piResults.add(cur);
            }

            double[] tArr = new double[timePoints.size()];
            for (int i = 0; i < tArr.length; i++) {
                tArr[i] = timePoints.get(i);
            }
            return new Pair<double[], List<double[]>>(tArr, piResults);
        } else {
            lsoda.integrate(ode, t0, pi0.toArray1D(), t1, y);
            List<double[]> pi = new ArrayList<double[]>();
            for (Double[] doubleArray : lsoda.getYvec()) {
                double[] arr = new double[doubleArray.length];
                for (int k = 0; k < doubleArray.length; k++) {
                    arr[k] = doubleArray[k].doubleValue();
                }
                pi.add(arr);
            }
            double[] tArr = new double[lsoda.getTvec().size()];
            for (int k = 0; k < tArr.length; k++) {
                tArr[k] = lsoda.getTvec().get(k).doubleValue();
            }
            return new Pair<double[], List<double[]>>(tArr, pi);
        }
    }

    public static Pair<double[], List<double[]>> ctmc_transient(Matrix Q, Matrix pi0, double t1) {
        return ctmc_transient(Q, pi0, 0.0, t1);
    }

    public static Pair<double[], List<double[]>> ctmc_transient(Matrix Q, double t1) {
        double[] pi0 = new double[Q.length()];
        Arrays.fill(pi0, 1.0 / Q.length());
        return ctmc_transient(Q, new Matrix(pi0), 0.0, t1);
    }
}
