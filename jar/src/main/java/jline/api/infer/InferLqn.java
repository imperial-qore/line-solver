package jline.api.infer;

import jline.VerboseLevel;
import jline.lang.layered.Activity;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Task;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ln.LNOptions;
import jline.solvers.ln.SolverLN;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * LQN parameter identification via an Extended Kalman Filter.
 *
 * <p>Estimates hidden Layered Queueing Network parameters (activity host
 * demands, task think times) from measured performance data (response times,
 * utilizations, throughputs), following Zheng, Yang, Woodside, Litoiu, Iszlai,
 * "Tracking Time-Varying Parameters in Software Systems with Extended Kalman
 * Filters", CASCON 2005 (equations 1-9). The parameter is modelled as a
 * zero-mean random walk a_k = a_{k-1} + w and the measurement as z_k = h(a_k)+v,
 * where h is the LQN performance model. This is the JAR counterpart of the
 * MATLAB {@code infer_lqn} family and the native-Python
 * {@code line_solver.inference.infer_lqn}; results match across codebases.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London. All rights reserved.
 */
public class InferLqn {

    private static final double EPS = 2.220446049250313e-16;

    private InferLqn() {
    }

    /** Finite-difference sensitivity matrix H = dh/da together with h0 = h(a). */
    public static class JacobianResult {
        public final Matrix H;
        public final Matrix h0;

        public JacobianResult(Matrix H, Matrix h0) {
            this.H = H;
            this.h0 = h0;
        }
    }

    /**
     * Apply a parameter vector to a LayeredNetwork in place.
     *
     * @param model     the LQN whose parameters are set
     * @param paramSpec list of parameters to set
     * @param a         (np x 1) parameter vector
     */
    public static void setParams(LayeredNetwork model, List<ParamSpec> paramSpec, Matrix a) {
        if (a.getNumRows() * a.getNumCols() != paramSpec.size()) {
            throw new RuntimeException("Length of parameter vector does not match paramSpec.");
        }
        for (int i = 0; i < paramSpec.size(); i++) {
            ParamSpec s = paramSpec.get(i);
            double val = a.get(i);
            switch (s.type) {
                case HOSTDEM: {
                    Activity act = findActivity(model, s.name);
                    if (act == null) {
                        throw new RuntimeException("Activity '" + s.name + "' not found.");
                    }
                    act.setHostDemand(val);
                    break;
                }
                case THINK: {
                    Task tsk = findTask(model, s.name);
                    if (tsk == null) {
                        throw new RuntimeException("Task '" + s.name + "' not found.");
                    }
                    tsk.setThinkTime(val);
                    break;
                }
                default:
                    throw new RuntimeException("Unknown parameter type.");
            }
        }
        // invalidate and regenerate the cached layered struct so the next solve
        // re-reads the mutated model objects
        model.getStruct(true);
    }

    /** Read the current values of the parameters named in paramSpec. */
    public static Matrix getParams(LayeredNetwork model, List<ParamSpec> paramSpec) {
        Matrix a0 = new Matrix(paramSpec.size(), 1);
        for (int i = 0; i < paramSpec.size(); i++) {
            ParamSpec s = paramSpec.get(i);
            switch (s.type) {
                case HOSTDEM: {
                    Activity act = findActivity(model, s.name);
                    if (act == null) {
                        throw new RuntimeException("Activity '" + s.name + "' not found.");
                    }
                    a0.set(i, 0, act.getHostDemandMean());
                    break;
                }
                case THINK: {
                    Task tsk = findTask(model, s.name);
                    if (tsk == null) {
                        throw new RuntimeException("Task '" + s.name + "' not found.");
                    }
                    a0.set(i, 0, tsk.getThinkTimeMean());
                    break;
                }
                default:
                    throw new RuntimeException("Unknown parameter type.");
            }
        }
        return a0;
    }

    /** Extract the observation vector selected by obsSpec from a solved table. */
    public static Matrix getObs(LayeredNetworkAvgTable table, List<ObsSpec> obsSpec) {
        List<String> names = table.getNodeNames();
        List<Double> qlen = table.getQLen();
        List<Double> util = table.getUtil();
        List<Double> respt = table.getRespT();
        List<Double> tput = table.getTput();
        Matrix z = new Matrix(obsSpec.size(), 1);
        for (int i = 0; i < obsSpec.size(); i++) {
            ObsSpec s = obsSpec.get(i);
            int idx = names.indexOf(s.name);
            if (idx < 0) {
                throw new RuntimeException("Element '" + s.name + "' not found in the LQN.");
            }
            double v;
            switch (s.metric) {
                case RESPT:
                    v = respt.get(idx);
                    break;
                case UTIL:
                    v = util.get(idx);
                    break;
                case TPUT:
                    v = tput.get(idx);
                    break;
                case QLEN:
                    v = qlen.get(idx);
                    break;
                default:
                    throw new RuntimeException("Unknown metric.");
            }
            z.set(i, 0, v);
        }
        return z;
    }

    /**
     * Forward finite-difference sensitivity matrix of h at a.
     *
     * <p>Column i is (h(a + d_i) - h0)/d_i with d_i = fdStep*max(|a_i|, fdFloor).
     * hfun is evaluated np+1 times. This is the approximate sensitivity matrix
     * H_k used in the EKF update.</p>
     */
    public static JacobianResult jacobian(Function<Matrix, Matrix> hfun, Matrix a,
                                          double fdStep, double fdFloor) {
        int np = a.getNumRows();
        Matrix h0 = hfun.apply(a);
        int no = h0.getNumRows();
        Matrix H = new Matrix(no, np);
        for (int i = 0; i < np; i++) {
            double d = fdStep * Math.max(Math.abs(a.get(i, 0)), fdFloor);
            Matrix ap = a.copy();
            ap.set(i, 0, ap.get(i, 0) + d);
            Matrix hi = hfun.apply(ap);
            for (int r = 0; r < no; r++) {
                H.set(r, i, (hi.get(r, 0) - h0.get(r, 0)) / d);
            }
        }
        return new JacobianResult(H, h0);
    }

    /**
     * Extended Kalman Filter tracking a hidden parameter vector across Z.
     *
     * @param hfun maps a parameter vector to a predicted observation z = h(a)
     * @param a0   (np x 1) initial estimate
     * @param P0   (np x np) initial covariance
     * @param Z    (no x nsteps) measurement matrix, one column per step
     * @param Q    (np x np) drift covariance
     * @param R    (no x no) measurement covariance
     * @param opt  options (fdStep, fdFloor, clampPositive, aTrue, verbose)
     */
    public static InferLqnResult ekf(Function<Matrix, Matrix> hfun, Matrix a0, Matrix P0,
                                     Matrix Z, Matrix Q, Matrix R, InferLqnOptions opt) {
        if (opt == null) {
            opt = new InferLqnOptions();
        }
        int np = a0.getNumRows();
        int no = Z.getNumRows();
        int nsteps = Z.getNumCols();

        Matrix ahat = new Matrix(np, nsteps);
        Matrix eHist = new Matrix(no, nsteps);
        Matrix zpredHist = new Matrix(no, nsteps);
        List<Matrix> Phist = new ArrayList<Matrix>();

        Matrix a = a0.copy();
        Matrix P = P0.copy();
        Matrix I = Matrix.eye(np);

        for (int k = 0; k < nsteps; k++) {
            // (1) predict: zero-mean drift; project covariance (eq 5)
            Matrix aPred = a;
            Matrix Ppred = P.add(Q);
            // (2,4) predicted measurement and sensitivity matrix
            JacobianResult jr = jacobian(hfun, aPred, opt.fdStep, opt.fdFloor);
            Matrix H = jr.H;
            Matrix zpred = jr.h0;
            // (3) prediction error
            Matrix zk = new Matrix(no, 1);
            for (int r = 0; r < no; r++) {
                zk.set(r, 0, Z.get(r, k));
            }
            Matrix e = zk.sub(zpred);
            // (6) Kalman gain (suboptimal because h is nonlinear)
            Matrix S = H.mult(Ppred).mult(H.transpose()).add(R);
            Matrix K = Ppred.mult(H.transpose()).mult(S.inv());
            // (4-update) improved estimate
            a = aPred.add(K.mult(e));
            if (opt.clampPositive) {
                for (int r = 0; r < np; r++) {
                    if (a.get(r, 0) < opt.fdFloor) {
                        a.set(r, 0, opt.fdFloor);
                    }
                }
            }
            // (7) covariance update; symmetrize for numerical stability
            P = I.sub(K.mult(H)).mult(Ppred);
            P = P.add(P.transpose()).scale(0.5);

            for (int r = 0; r < np; r++) {
                ahat.set(r, k, a.get(r, 0));
            }
            for (int r = 0; r < no; r++) {
                eHist.set(r, k, e.get(r, 0));
                zpredHist.set(r, k, zpred.get(r, 0));
            }
            Phist.add(P.copy());
            if (opt.verbose) {
                System.out.printf("[InferLqn.ekf] step %d/%d  ||e||=%.4g%n", k + 1, nsteps, norm(e));
            }
        }

        InferLqnResult res = new InferLqnResult();
        res.ahat = ahat;
        res.e = eHist;
        res.zpred = zpredHist;
        res.P = P;
        res.Phist = Phist;
        res.Er = rms(eHist);
        if (opt.aTrue != null) {
            double acc = 0.0;
            int cnt = 0;
            for (int r = 0; r < np; r++) {
                double at = opt.aTrue.get(r, 0);
                for (int k = 0; k < nsteps; k++) {
                    double d = ahat.get(r, k) - at;
                    acc += d * d;
                    cnt++;
                }
            }
            res.Ea = Math.sqrt(acc / cnt);
        } else {
            res.Ea = null;
        }
        return res;
    }

    /**
     * Identify hidden LQN parameters from measured performance data.
     *
     * <p>Estimates the parameters named in paramSpec (activity host demands
     * and/or task think times) of the LayeredNetwork model from the measurement
     * sequence Z (no x nsteps) using an EKF over the observation model defined by
     * obsSpec. A single measurement column with opt.QFac = 0 reduces to one-shot
     * least-squares calibration. On return the model has the final estimate
     * applied.</p>
     */
    public static InferLqnResult inferLqn(final LayeredNetwork model, final List<ParamSpec> paramSpec,
                                          final List<ObsSpec> obsSpec, Matrix Z, InferLqnOptions opt) {
        if (opt == null) {
            opt = new InferLqnOptions();
        }
        int no = obsSpec.size();
        if (Z.getNumRows() != no) {
            throw new RuntimeException("Row count of Z must equal obsSpec size.");
        }

        Matrix a0 = (opt.a0 != null) ? opt.a0.copy() : getParams(model, paramSpec);
        int np = a0.getNumRows();

        double gammaT;
        if (opt.gammaT != null) {
            gammaT = opt.gammaT;
        } else if (opt.T != null && opt.Tstar != null) {
            gammaT = opt.T / opt.Tstar;
        } else {
            gammaT = 1.0;
        }

        Matrix Q = (opt.Q != null) ? opt.Q : buildQ(a0, opt.QFac, opt.cvA);   // eq 9a
        Matrix R = (opt.R != null) ? opt.R : buildR(Z, opt.RFac, gammaT);     // eq 9b
        Matrix P0 = (opt.P0 != null) ? opt.P0 : buildP0(a0);

        final InferLqnOptions.ObservationSolver solver = opt.solver;
        Function<Matrix, Matrix> hfun = new Function<Matrix, Matrix>() {
            @Override
            public Matrix apply(Matrix a) {
                setParams(model, paramSpec, a);
                LayeredNetworkAvgTable table;
                if (solver != null) {
                    table = solver.evaluate(model);
                } else {
                    SolverOptions so = new LNOptions();
                    so.verbose = VerboseLevel.SILENT;
                    table = (LayeredNetworkAvgTable) new SolverLN(model, so).getEnsembleAvg();
                }
                return getObs(table, obsSpec);
            }
        };

        InferLqnResult res = ekf(hfun, a0, P0, Z, Q, R, opt);
        res.a0 = a0;
        res.Q = Q;
        res.R = R;
        res.P0 = P0;

        // apply the final estimate to the model
        Matrix aLast = new Matrix(np, 1);
        int nsteps = res.ahat.getNumCols();
        for (int r = 0; r < np; r++) {
            aLast.set(r, 0, res.ahat.get(r, nsteps - 1));
        }
        setParams(model, paramSpec, aLast);
        return res;
    }

    // ---------------- helpers ----------------

    private static Activity findActivity(LayeredNetwork model, String name) {
        for (Activity a : model.getActivities().values()) {
            if (a.getName().equals(name)) {
                return a;
            }
        }
        return null;
    }

    private static Task findTask(LayeredNetwork model, String name) {
        for (Task t : model.getTasks().values()) {
            if (t.getName().equals(name)) {
                return t;
            }
        }
        return null;
    }

    private static Matrix buildQ(Matrix a0, double QFac, double cvA) {
        int np = a0.getNumRows();
        double[] d = new double[np];
        for (int i = 0; i < np; i++) {
            double q = QFac * Math.abs(a0.get(i, 0)) * cvA;
            d[i] = Math.max(q * q, EPS);
        }
        return Matrix.diag(d);
    }

    private static Matrix buildR(Matrix Z, double RFac, double gammaT) {
        int no = Z.getNumRows();
        int nsteps = Z.getNumCols();
        double[] d = new double[no];
        for (int i = 0; i < no; i++) {
            double zbar = 0.0;
            for (int k = 0; k < nsteps; k++) {
                zbar += Z.get(i, k);
            }
            zbar /= nsteps;
            double r = (RFac * Math.abs(zbar)) / 1.96;
            d[i] = Math.max((r * r) / gammaT, EPS);
        }
        return Matrix.diag(d);
    }

    private static Matrix buildP0(Matrix a0) {
        int np = a0.getNumRows();
        double[] d = new double[np];
        for (int i = 0; i < np; i++) {
            double p = 0.5 * Math.abs(a0.get(i, 0));
            d[i] = Math.max(p * p, EPS);
        }
        return Matrix.diag(d);
    }

    private static double norm(Matrix v) {
        double s = 0.0;
        int n = v.getNumRows() * v.getNumCols();
        for (int i = 0; i < n; i++) {
            double x = v.get(i);
            s += x * x;
        }
        return Math.sqrt(s);
    }

    private static double rms(Matrix M) {
        int rows = M.getNumRows();
        int cols = M.getNumCols();
        double s = 0.0;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double x = M.get(r, c);
                s += x * x;
            }
        }
        return Math.sqrt(s / (rows * cols));
    }
}
