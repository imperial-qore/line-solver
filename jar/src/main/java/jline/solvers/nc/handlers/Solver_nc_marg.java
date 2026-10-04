package jline.solvers.nc.handlers;

import org.apache.commons.math3.util.FastMath;

import jline.api.mam.Map_pie;
import jline.api.pfqn.ld.Pfqn_ncld;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_nc_marg {
    private Solver_nc_marg() {}

    /**
     * Convenience overload accepting a boxed {@code Double} {@code lG} (which may be null,
     * meaning the normalizing constant is unknown and will be computed). Delegates to
     * {@link #solver_nc_marg(NetworkStruct, SolverOptions, double)}.
     */
    public static SolverNC.SolverNCMargReturn solver_nc_marg(NetworkStruct sn, SolverOptions options, Double lG) {
        return solver_nc_marg(sn, options, lG != null ? lG.doubleValue() : Double.NaN);
    }

    public static SolverNC.SolverNCMargReturn solver_nc_marg(NetworkStruct sn, SolverOptions options, double lG) {
        int M = sn.nstations;
        int K = sn.nclasses;
        java.util.Map state = sn.state;
        Matrix S = sn.nservers;
        Matrix V = new Matrix(sn.nstateful, K);
        for (int i = 0; i < sn.visits.size(); i++) {
            V = V.add(1.0, sn.visits.get(i));
        }
        Matrix rates = sn.rates;
        Matrix ST = rates.copy();
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                ST.set(i, j, 1.0 / ST.get(i, j));
            }
        }
        ST.removeNaN();

        Ret.snGetDemands ret = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = ret.Dchain;
        Matrix STchain = ret.STchain;
        Matrix Nchain = ret.Nchain;

        long startTimeMillis = System.nanoTime();
        M = STchain.getNumRows();
        K = STchain.getNumCols();

        Matrix mu = new Matrix(0, (int) Nchain.elementSum());
        for (int i = 0; i < M; i++) {
            Matrix tmp = new Matrix(1, (int) Nchain.elementSum());
            if (Utils.isInf(S.get(i))) {
                for (int j = 0; j < tmp.length(); j++) {
                    tmp.set(j, (double) (j + 1));
                }
            } else {
                for (int j = 0; j < tmp.length(); j++) {
                    tmp.set(j, FastMath.min((double) (j + 1), S.get(i)));
                }
            }
            mu = Matrix.concatRows(mu, tmp, null);
        }

        if (Double.isNaN(lG)) {
            Matrix Z_tmp = Nchain.copy();
            Z_tmp.fill(0.0);
            lG = (double) Pfqn_ncld.pfqn_ncld(Lchain, Nchain, Z_tmp, mu, options).lG;
        }

        double G = FastMath.exp(lG);
        Matrix lPr = new Matrix(sn.nstations, 1);
        lPr.fill(0.0);

        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            int isf = (int) sn.stationToStateful.get(ist);
            State.StateMarginalStatistics ret1 = ToMarginal.toMarginal(
                    sn, ind, (Matrix) state.get(sn.stateful.get(isf)), null, null, null, null, null);
            Matrix nirvec = ret1.nir;
            Matrix sivec = ret1.sir;
            java.util.List<Matrix> kirvec = ret1.kir;

            // matrix.elementMin negative check disabled per Java source comments
            Matrix nivec_chain = nirvec.mult(sn.chains.transpose());

            Matrix Lchain_tmp = new Matrix(0, Lchain.getNumCols());
            Matrix mu_tmp = new Matrix(0, mu.getNumCols());
            for (int i = 0; i < sn.nstations; i++) {
                if (i != ist) {
                    Matrix Lchain_row_i = Matrix.extractRows(Lchain, i, i + 1, null);
                    Matrix mu_row_i = Matrix.extractRows(mu, i, i + 1, null);
                    Lchain_tmp = Matrix.concatRows(Lchain_tmp, Lchain_row_i, null);
                    mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null);
                }
            }
            Matrix Nchain_tmp = Nchain.copy();
            for (int i = 0; i < Nchain_tmp.length(); i++) {
                Nchain_tmp.set(i, Nchain_tmp.get(i) - nivec_chain.get(i));
            }
            Matrix Zchain_tmp = Nchain.copy();
            Zchain_tmp.fill(0.0);

            double lG_minus_i = (double) Pfqn_ncld.pfqn_ncld(Lchain_tmp, Nchain_tmp, Zchain_tmp, mu_tmp, options).lG;
            double lF_i = 0.0;

            SchedStrategy schedStrat = sn.sched.get(sn.stations.get(ist));
            if (schedStrat == SchedStrategy.FCFS) {
                for (int r = 0; r < K; r++) {
                    MatrixCell PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
                    if (!PHr.isEmpty()) {
                        Matrix kir = new Matrix(1, kirvec.size());
                        for (int i = 0; i < kir.length(); i++) {
                            kir.set(i, kirvec.get(i).get(0, r));
                        }
                        if (kir.length() > 1) {
                            throw new RuntimeException("solver_nc_marg: Cannot return state probability "
                                    + "because the product-form solution requires exponential service times at FCFS nodes.");
                        }

                        if (Math.abs(ST.get(ist, r) - Matrix.extractRows(ST, ist, ist + 1, null).elementMax()) > 1e-6) {
                            throw new RuntimeException("solver_nc_marg: Cannot return state probability "
                                    + "because the product-form solution requires identical service times at FCFS nodes.");
                        }
                    }
                }
                boolean allZero = true;
                for (int i = 0; allZero && i < sivec.getNumRows(); i++) {
                    for (int j = 0; allZero && j < sivec.getNumCols(); j++) {
                        if (Math.abs(sivec.get(i, j)) > 1e-6) {
                            allZero = false;
                        }
                    }
                }
                if (!allZero) {
                    double sumLog = 0.0;
                    for (int r = 0; r < K; r++) {
                        sumLog += nirvec.get(0, r) * Math.log(V.get(ist, r));
                    }
                    double sum_kirvec = 0.0;
                    for (int i = 0; i < kirvec.size(); i++) {
                        sum_kirvec += kirvec.get(i).elementSum();
                    }
                    Matrix mu_row_ist = new Matrix(1, (int) sum_kirvec);
                    Matrix.extract(mu, ist, ist + 1, 0, (int) sum_kirvec, mu_row_ist, 0, 0);
                    for (int i = 0; i < mu_row_ist.length(); i++) {
                        mu_row_ist.set(i, FastMath.log(mu_row_ist.get(i)));
                    }
                    lF_i += (sumLog - mu_row_ist.elementSum());
                } else {
                    lF_i = 0.0;
                }
            } else if (schedStrat == SchedStrategy.PS) {
                for (int r = 0; r < K; r++) {
                    MatrixCell PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
                    if (!PHr.isEmpty()) {
                        Matrix kir = new Matrix(1, kirvec.size());
                        for (int i = 0; i < kir.length(); i++) {
                            kir.set(i, kirvec.get(i).get(0, r));
                        }

                        Matrix PHr_tmp = PHr.get(0);
                        for (int i = 0; i < PHr_tmp.getNumRows(); i++) {
                            for (int j = 0; j < PHr_tmp.getNumCols(); j++) {
                                PHr_tmp.set(i, j, -1.0 * PHr_tmp.get(i, j));
                            }
                        }
                        PHr_tmp = PHr_tmp.inv();
                        Matrix Ar = Map_pie.map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                        Matrix kir_tmp = Ar.copy();
                        for (int i = 0; i < kir_tmp.getNumRows(); i++) {
                            for (int j = 0; j < kir_tmp.getNumCols(); j++) {
                                kir_tmp.set(i, j, kir.get(i, j) * FastMath.log(V.get(ist, r) * kir_tmp.get(i, j)));
                            }
                        }

                        lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum());
                    }
                }

                double sum_kirvec = 0.0;
                for (int i = 0; i < kirvec.size(); i++) {
                    sum_kirvec += kirvec.get(i).elementSum();
                }
                Matrix mu_row_ist = new Matrix(1, (int) sum_kirvec);
                Matrix.extract(mu, ist, ist + 1, 0, (int) sum_kirvec, mu_row_ist, 0, 0);
                for (int i = 0; i < mu_row_ist.getNumCols(); i++) {
                    mu_row_ist.set(i, FastMath.log(mu_row_ist.get(i)));
                }
                lF_i += (Maths.factln(sum_kirvec) - mu_row_ist.elementSum());
            } else if (schedStrat == SchedStrategy.INF) {
                for (int r = 0; r < K; r++) {
                    MatrixCell PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
                    if (!PHr.isEmpty()) {
                        Matrix kir = new Matrix(1, kirvec.size());
                        for (int i = 0; i < kir.length(); i++) {
                            kir.set(i, kirvec.get(i).get(0, r));
                        }
                        Matrix PHr_tmp = PHr.get(0);
                        for (int i = 0; i < PHr_tmp.getNumRows(); i++) {
                            for (int j = 0; j < PHr_tmp.getNumCols(); j++) {
                                PHr_tmp.set(i, j, -1.0 * PHr_tmp.get(i, j));
                            }
                        }
                        PHr_tmp = PHr_tmp.inv();
                        Matrix Ar = Map_pie.map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                        Matrix kir_tmp = Ar.copy();
                        for (int i = 0; i < kir_tmp.getNumRows(); i++) {
                            for (int j = 0; j < kir_tmp.getNumCols(); j++) {
                                kir_tmp.set(i, j, kir.get(i, j) * FastMath.log(V.get(ist, r) * kir_tmp.get(i, j)));
                            }
                        }

                        lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum());
                    }
                }
            }

            lPr.set(ist, lF_i + lG_minus_i - lG);
        }
        long endTimeMillis = System.nanoTime();
        double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;
        lPr.removeNaN();
        return new SolverNC.SolverNCMargReturn(lPr, G, lG, runtime);
    }
}
