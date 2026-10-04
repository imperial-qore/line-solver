/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import jline.api.mam.Map_mean;
import jline.api.mam.Map_pie;
import jline.io.InputOutput;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_ldqbd_transient {
    private Solver_mam_ldqbd_transient() {}

    public static TransientResult solver_mam_ldqbd_transient(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        if (K != 1) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Transient QBD method requires a single-class model.");
            return emptyTransientResult(M, K);
        }

        double N = sn.njobs.get(0, 0);
        if (!Double.isInfinite(N)) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Transient QBD method requires an open model.");
            return emptyTransientResult(M, K);
        }

        int sourceIdx = -1;
        int queueIdx = -1;
        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.EXT) {
                sourceIdx = i;
            } else if (sched == SchedStrategy.FCFS) {
                queueIdx = i;
            }
        }

        if (sourceIdx < 0 || queueIdx < 0) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Transient QBD method requires exactly one Source and one FCFS Queue.");
            return emptyTransientResult(M, K);
        }

        double lambda = sn.rates.get(sourceIdx, 0);
        Station queueStation = sn.stations.get(queueIdx);
        JobClass jobClass = sn.jobclasses.get(0);
        java.util.Map<JobClass, MatrixCell> queueProcMap = sn.proc.get(queueStation);
        if (queueProcMap == null) {
            throw new RuntimeException("No service process for queue station");
        }
        MatrixCell PH_queue = queueProcMap.get(jobClass);
        if (PH_queue == null) {
            throw new RuntimeException("No service process for queue station");
        }
        int nServers = (int) sn.nservers.get(queueIdx, 0);
        int bufCap = (int) sn.cap.get(queueIdx, 0);

        Matrix D0 = PH_queue.get(0);
        if (D0 == null) {
            throw new RuntimeException("No D0 matrix in service process");
        }
        int nPhases;
        boolean isPH;
        double mu;

        if (D0.getNumRows() == 1 && D0.getNumCols() == 1) {
            mu = -D0.get(0, 0);
            nPhases = 1;
            isPH = false;
        } else {
            nPhases = D0.getNumRows();
            isPH = true;
            mu = 0.0;
        }

        if (isPH && nServers > 1) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Transient QBD with PH service supports single-server only.");
            return emptyTransientResult(M, K);
        }

        Matrix D1 = null;
        Matrix alpha = null;
        Matrix t_exit = null;
        if (isPH) {
            D1 = PH_queue.get(1);
            if (D1 == null) {
                throw new RuntimeException("No D1 matrix in PH service process");
            }
            alpha = Map_pie.map_pie(PH_queue);
            t_exit = D0.mult(Matrix.ones(nPhases, 1)).scale(-1.0);
        }

        double T_start = options.timespan[0];
        double T_end = options.timespan[1];
        double T_duration = T_end - T_start;

        int Cap;
        if (Double.isInfinite((double) bufCap) || bufCap > 1000000) {
            double muEff;
            if (isPH) {
                double meanSvc = Map_mean.map_mean(PH_queue);
                muEff = (meanSvc > 0) ? 1.0 / meanSvc : 1.0;
            } else {
                muEff = mu;
            }
            double rho = lambda / (nServers * muEff);
            double tol = options.tol;
            if (rho < 1.0) {
                Cap = Math.min(10000,
                        Math.max(200, (int) Math.ceil(-Math.log(tol) / (-Math.log(rho)))));
            } else {
                Cap = 10000;
            }
        } else {
            Cap = bufCap;
        }

        Matrix Q;
        int dim;

        if (!isPH) {
            dim = Cap + 1;
            Q = new Matrix(dim, dim);
            for (int n = 0; n <= Cap; n++) {
                double dep = Math.min(n, nServers) * mu;
                double arr = (n < Cap) ? lambda : 0.0;
                if (n > 0) {
                    Q.set(n, n - 1, dep);
                }
                if (n < Cap) {
                    Q.set(n, n + 1, arr);
                }
                Q.set(n, n, -(dep + arr));
            }
        } else {
            dim = 1 + Cap * nPhases;
            Q = new Matrix(dim, dim);

            Q.set(0, 0, -lambda);
            for (int j = 0; j < nPhases; j++) {
                Q.set(0, 1 + j, lambda * alpha.get(0, j));
            }

            for (int n = 1; n <= Cap; n++) {
                int rowStart = 1 + (n - 1) * nPhases;

                for (int i = 0; i < nPhases; i++) {
                    for (int j = 0; j < nPhases; j++) {
                        Q.set(rowStart + i, rowStart + j, D0.get(i, j));
                    }
                    if (n < Cap) {
                        Q.set(rowStart + i, rowStart + i, Q.get(rowStart + i, rowStart + i) - lambda);
                    }
                }

                if (n < Cap) {
                    int nextStart = rowStart + nPhases;
                    for (int i = 0; i < nPhases; i++) {
                        Q.set(rowStart + i, nextStart + i, lambda);
                    }
                }

                if (n == 1) {
                    for (int i = 0; i < nPhases; i++) {
                        Q.set(rowStart + i, 0, t_exit.get(i, 0));
                    }
                } else {
                    int prevStart = rowStart - nPhases;
                    for (int i = 0; i < nPhases; i++) {
                        for (int j = 0; j < nPhases; j++) {
                            Q.set(rowStart + i, prevStart + j, D1.get(i, j));
                        }
                    }
                }
            }
        }

        int nTimePoints = Math.min(101, Math.max(11, (int) Math.round(T_duration * 10)));
        double dt = T_duration / (nTimePoints - 1);
        double[] times = new double[nTimePoints];
        for (int i = 0; i < nTimePoints; i++) {
            times[i] = T_start + i * dt;
        }

        Matrix pi_t = new Matrix(1, dim);
        pi_t.set(0, 0, 1.0);

        Matrix eQdt = Q.scale(dt).expm();

        double[] queue_lengths = new double[nTimePoints];
        double[] util_values = new double[nTimePoints];
        double[] tput_values = new double[nTimePoints];

        Matrix pi_current = pi_t;
        for (int t_idx = 0; t_idx < nTimePoints; t_idx++) {
            double q = 0.0;
            double u = 0.0;
            double tput = 0.0;

            for (int n = 0; n <= Cap; n++) {
                double p_n;
                if (!isPH) {
                    p_n = pi_current.get(0, n);
                } else {
                    if (n == 0) {
                        p_n = pi_current.get(0, 0);
                    } else {
                        int idx_s = 1 + (n - 1) * nPhases;
                        double sum = 0.0;
                        for (int j = 0; j < nPhases; j++) {
                            sum += pi_current.get(0, idx_s + j);
                        }
                        p_n = sum;
                    }
                }
                q += n * p_n;
                if (n >= 1) {
                    u += ((double) Math.min(n, nServers) / nServers) * p_n;
                    if (!isPH) {
                        tput += Math.min(n, nServers) * mu * p_n;
                    } else {
                        int idx_s = 1 + (n - 1) * nPhases;
                        for (int j = 0; j < nPhases; j++) {
                            tput += pi_current.get(0, idx_s + j) * t_exit.get(j, 0);
                        }
                    }
                }
            }
            queue_lengths[t_idx] = q;
            util_values[t_idx] = u;
            tput_values[t_idx] = tput;

            if (t_idx < nTimePoints - 1) {
                pi_current = pi_current.mult(eQdt);
            }
        }

        Matrix[][] Qt = new Matrix[M][K];
        Matrix[][] Ut = new Matrix[M][K];
        Matrix[][] Tt = new Matrix[M][K];

        Matrix qResult = new Matrix(nTimePoints, 2);
        Matrix uResult = new Matrix(nTimePoints, 2);
        Matrix tResult = new Matrix(nTimePoints, 2);
        for (int t = 0; t < nTimePoints; t++) {
            qResult.set(t, 0, queue_lengths[t]);
            qResult.set(t, 1, times[t]);
            uResult.set(t, 0, util_values[t]);
            uResult.set(t, 1, times[t]);
            tResult.set(t, 0, tput_values[t]);
            tResult.set(t, 1, times[t]);
        }
        Qt[queueIdx][0] = qResult;
        Ut[queueIdx][0] = uResult;
        Tt[queueIdx][0] = tResult;

        return new TransientResult(Qt, Ut, Tt);
    }

    private static TransientResult emptyTransientResult(int M, int K) {
        return new TransientResult(new Matrix[M][K], new Matrix[M][K], new Matrix[M][K]);
    }
}
