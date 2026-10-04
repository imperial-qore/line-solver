/**
 * @file NPFQN Non-Exponential Approximation
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class Npfqn_nonexp_approx {
    private Npfqn_nonexp_approx() {}

    /**
     * Approximates non-product-form queueing networks using the specified method.
     */
    public static Ret.npfqnNonexpApprox npfqn_nonexp_approx(String method,
                                                            NetworkStruct sn,
                                                            Matrix ST,
                                                            Matrix V,
                                                            Matrix SCV,
                                                            Matrix Tin,
                                                            Matrix Uin,
                                                            Matrix gamma,
                                                            Matrix nservers) {
        int M = sn.nstations;
        Matrix rho = new Matrix(M, 1);
        rho.zero();
        Matrix scva = Matrix.ones(M, 1);
        Matrix scvs = Matrix.ones(M, 1);
        Matrix eta = Matrix.ones(M, 1);
        Matrix T = Tin.copy();
        Matrix U = Uin.copy();
        Matrix STLocal = ST.copy();
        Matrix gammaLocal = gamma.copy();
        Matrix nserversLocal = nservers.copy();

        if ("default".equals(method) || "none".equals(method) || "hmva".equals(method)) {
            return new Ret.npfqnNonexpApprox(STLocal, gammaLocal, nserversLocal, rho, scva, scvs, eta);
        } else if ("interp".equals(method)) {
            int ist = 0;
            while (ist < M) {
                Matrix nnzClasses = new Matrix(1, ST.getNumCols());
                nnzClasses.zero();
                {
                    int j = 0;
                    while (j < STLocal.getNumCols()) {
                        if (Double.isFinite(STLocal.get(ist, j)) && Double.isFinite(SCV.get(ist, j))) {
                            nnzClasses.set(0, j, 1);
                        }
                        j++;
                    }
                }
                int j = 0;
                while (j < nnzClasses.getNumElements()) {
                    if (nnzClasses.get(j) > 0) {
                        rho.set(ist, rho.get(ist) + U.get(ist, j));
                    }
                    j++;
                }
                if (nnzClasses.elementSum() != 0.0) {
                    SchedStrategy s = sn.sched.get(sn.stations.get(ist));
                    if (s == SchedStrategy.FCFS) {
                        Matrix STinnz = new Matrix(1, (int) nnzClasses.elementSum());
                        Matrix SCVinnz = STinnz.copy();
                        Matrix Tinnz = STinnz.copy();
                        int tempj = 0;
                        int jj = 0;
                        while (jj < nnzClasses.getNumElements()) {
                            if (nnzClasses.get(jj) > 0) {
                                STinnz.set(tempj, STLocal.get(ist, jj));
                                SCVinnz.set(tempj, SCV.get(ist, jj));
                                Tinnz.set(tempj, T.get(ist, jj));
                                tempj++;
                            }
                            jj++;
                        }
                        if (STinnz.elementMax() - STinnz.elementMin() > 0
                                || SCVinnz.elementMax() > 1 + GlobalConstants.FineTol
                                || SCVinnz.elementMin() < 1 - GlobalConstants.FineTol) {
                            scva.set(ist, 1.0);
                            scvs.set(ist, SCVinnz.mult(Tinnz.transpose()).toDouble() / Tinnz.elementSum());
                            gammaLocal.set(ist, (Math.pow(rho.get(ist), nserversLocal.get(ist)) + rho.get(ist)) / 2);
                            if (scvs.get(ist) > 1 - 1e-6 && scvs.get(ist) < 1 + 1e-6 && nserversLocal.get(ist) == 1.0) {
                                eta.set(ist, rho.get(ist));
                            } else {
                                eta.set(ist, FastMath.exp(-2 * (1 - rho.get(ist)) / (scvs.get(ist) + scva.get(ist) * rho.get(ist))));
                            }
                            int order = 8;
                            double ai = FastMath.pow(rho.get(ist), order);
                            double bi = FastMath.pow(rho.get(ist), order);
                            int k = 0;
                            while (k < nnzClasses.getNumElements()) {
                                if (nnzClasses.get(k) > 0 && sn.rates.get(ist, k) > 0) {
                                    STLocal.set(ist, k,
                                            FastMath.max(0.0, 1 - ai) * STLocal.get(ist, k)
                                                    + ai * (bi * eta.get(ist) + FastMath.max(0.0, 1 - bi) * gammaLocal.get(ist))
                                                    * (nserversLocal.get(ist) / Tinnz.elementSum()));
                                }
                                k++;
                            }
                            nserversLocal.set(ist, 1.0);
                        }
                    }
                }
                ist++;
            }
        } else {
            throw new IllegalArgumentException("Unknown approximation method: " + method);
        }

        return new Ret.npfqnNonexpApprox(STLocal, gammaLocal, nserversLocal, rho, scva, scvs, eta);
    }
}
