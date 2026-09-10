package jline.opt.sensitivity;

import jline.api.mc.Ctmc_sens;
import jline.api.mc.Ctmc_solve;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Node;
import jline.lang.nodes.ServiceStation;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.opt.results.SensitivityData;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Analytic d(QLen)/d(rate) for an arbitrary Markovian model, via the
 * generator-derivative equation of Trivedi and Bobbio (2017), Eq. (9.81).
 * Mirrors MATLAB {@code opt.sens.ctmcSensitivities}.
 *
 * <p>This is the fallback for the cases the differentiated-MVA primitive cannot
 * reach: open, multiserver, and non-unit-visit networks, for which
 * {@link OpenSensitivity} and {@link ClosedSensitivity} return {@code null}. It
 * is exact wherever {@link SolverCTMC} is exact, and correspondingly it is
 * limited by the state-space size rather than by the product-form assumptions.</p>
 *
 * <p>Only queue-length sensitivities are produced. Utilization, response time,
 * and throughput are reward rates whose definition involves the solved metrics
 * themselves, so they need the reward-derivative term of Eq. (9.83) and are not
 * covered here. Returns {@code null} when the model has no CTMC representation of
 * tractable size or no perturbable exponential rate.</p>
 */
final class CtmcSensitivity {

    private CtmcSensitivity() {
    }

    static SensitivityData compute(Network model) {
        NetworkStruct sn = model.getStruct();
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix rates = sn.rates;

        SolverCTMC solver = new SolverCTMC(model);
        Matrix Q;
        Matrix spaceAggr;
        try {
            Q = solver.getGenerator().infGen;
            spaceAggr = solver.getStateSpaceAggr();
        } catch (RuntimeException e) {
            return null;   // state space unavailable or over the memory gate
        }
        if (Q == null || spaceAggr == null || spaceAggr.isEmpty()) {
            return null;
        }
        int n = Q.length();
        Matrix pi = Ctmc_solve.ctmc_solve(Q);
        SolverOptions options = solver.getOptions();

        // One parameter per finite positive exponential service rate. The rate
        // setter substitutes an Exp of the perturbed rate, which is a perturbation
        // of theta only where the nominal process is itself exponential; for an
        // Erlang or Coxian service the substitution would change the distribution
        // family and the difference quotient would not be dQ/dtheta, so those
        // stations are skipped rather than reported wrong.
        List<int[]> params = new ArrayList<int[]>();   // {stationIdx, classIdx, nodeIdx}
        List<Double> paramValue = new ArrayList<Double>();
        List<String> paramKey = new ArrayList<String>();
        for (int ist = 0; ist < M; ist++) {
            Station station = sn.stations.get(ist);
            Map<JobClass, ProcessType> procRow = sn.procid != null ? sn.procid.get(station) : null;
            for (int k = 0; k < K; k++) {
                double rate = rates.get(ist, k);
                if (Double.isNaN(rate) || Double.isInfinite(rate) || rate <= 0) {
                    continue;
                }
                if (procRow == null || procRow.get(sn.jobclasses.get(k)) != ProcessType.EXP) {
                    continue;
                }
                int nodeIdx = (int) sn.stationToNode.get(ist);
                params.add(new int[]{ist, k, nodeIdx});
                paramValue.add(rate);
                paramKey.add(SensitivityData.paramKey(sn.nodenames.get(nodeIdx), sn.classnames.get(k)));
            }
        }

        SensitivityData sens = new SensitivityData();
        boolean any = false;
        for (int p = 0; p < params.size(); p++) {
            int nodeIdx = params.get(p)[2];
            int classIdx = params.get(p)[1];
            double theta = paramValue.get(p);
            double h = Math.max(Math.abs(theta), 1.0) * 1e-6;

            Matrix Qp = perturbedGenerator(model, options, nodeIdx, classIdx, theta + h);
            Matrix Qm = perturbedGenerator(model, options, nodeIdx, classIdx, theta - h);
            if (Qp == null || Qm == null || Qp.length() != n || Qm.length() != n) {
                continue;   // perturbation changed the state space; skip this parameter
            }
            Matrix dQ = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    dQ.set(i, j, (Qp.get(i, j) - Qm.get(i, j)) / (2.0 * h));
                }
            }
            Matrix dpi = Ctmc_sens.ctmc_sens(Q, dQ, pi);   // 1 x n

            String pkey = paramKey.get(p);
            for (int ist = 0; ist < M; ist++) {
                String nameI = sn.nodenames.get((int) sn.stationToNode.get(ist));
                for (int r = 0; r < K; r++) {
                    int col = ist * K + r;
                    if (col >= spaceAggr.getNumCols()) {
                        continue;
                    }
                    double val = 0.0;
                    boolean finite = true;
                    for (int s = 0; s < n; s++) {
                        double rv = spaceAggr.get(s, col);
                        if (Double.isNaN(rv) || Double.isInfinite(rv)) {
                            finite = false;   // Source stations carry an infinite population
                            break;
                        }
                        val += dpi.get(0, s) * rv;
                    }
                    if (!finite) {
                        continue;
                    }
                    sens.add("QLen", SensitivityData.metricKey(nameI, sn.classnames.get(r)), pkey, val);
                    any = true;
                }
            }
        }
        return any ? sens : null;
    }

    /**
     * Rebuild the generator with the service rate of node {@code nodeIdx} class
     * {@code classIdx} set to {@code value}, on a copy of the model so the
     * caller's model is left untouched. The distribution is replaced by an
     * exponential of the requested rate, so this is only meaningful where the
     * nominal service is itself exponential.
     *
     * <p>The hard refresh is required, not defensive: {@code setService}
     * deliberately leaves the cached struct in place, so a copy that inherited a
     * built struct would report the old rate and the difference quotient would
     * silently come out as zero.</p>
     */
    private static Matrix perturbedGenerator(Network model, SolverOptions options,
                                             int nodeIdx, int classIdx, double value) {
        Network modelCopy = model.copy();
        Node node = modelCopy.getNodes().get(nodeIdx);
        if (!(node instanceof ServiceStation)) {
            return null;
        }
        JobClass cls = modelCopy.getClasses().get(classIdx);
        ((ServiceStation) node).setService(cls, new Exp(value));
        modelCopy.refreshStruct(true);
        SolverCTMC solverCopy = new SolverCTMC(modelCopy, options);
        try {
            return solverCopy.getGenerator().infGen;
        } catch (RuntimeException e) {
            return null;
        }
    }
}
