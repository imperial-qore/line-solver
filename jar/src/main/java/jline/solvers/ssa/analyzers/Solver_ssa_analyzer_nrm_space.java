package jline.solvers.ssa.analyzers;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SSAResult;
import jline.solvers.ssa.handlers.Solver_ssa_nrm_space;
import jline.solvers.ssa.handlers.Solver_ssa_nrm_space.SolverSSAResultNRMSpace;
import jline.util.matrix.ColumnView;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.lang.JobClass;
import jline.lang.nodes.Station;

import java.util.Map;

public final class Solver_ssa_analyzer_nrm_space {
    private Solver_ssa_analyzer_nrm_space() {}

    public static SSAResult solver_ssa_analyzer_nrm_space(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        Map<Station, Map<JobClass, MatrixCell>> PH = sn.proc;
        Map<Station, SchedStrategy> sched = sn.sched;

        Matrix XN = new Matrix(1, K);
        XN.fill(Double.NaN);
        Matrix UN = new Matrix(M, K);
        UN.fill(Double.NaN);
        Matrix QN = new Matrix(M, K);
        QN.fill(Double.NaN);
        Matrix RN = new Matrix(M, K);
        RN.fill(Double.NaN);
        Matrix TN = new Matrix(M, K);
        TN.fill(Double.NaN);
        Matrix CN = new Matrix(1, K);
        CN.fill(Double.NaN);

        SolverSSAResultNRMSpace result = Solver_ssa_nrm_space.solver_ssa_nrm_space(sn, options);
        Matrix pi = result.pi;
        Matrix space = result.outspace;
        Matrix depRates = result.depRates;

        if (pi.getNumRows() > 1) {
            pi = pi.transpose();
        }

        for (int k = 0; k < K; k++) {
            int refnd = (int) sn.stationToNode.get((int) sn.refstat.get(k));
            ColumnView colView = depRates.getColumnView(refnd * K + k);
            XN.set(0, k, pi.multColumnView(colView));
        }

        Matrix rates = sn.rates;
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int k = 0; k < K; k++) {
                ColumnView dcolView = depRates.getColumnView(ind * K + k);
                ColumnView scolView = space.getColumnView(ind * K + k);
                TN.set(ist, k, pi.multColumnView(dcolView));
                QN.set(ist, k, pi.multColumnView(scolView));
            }

            SchedStrategy sp = sched.get(sn.stations.get(ist));
            if (sp == SchedStrategy.INF || sp == SchedStrategy.EXT) {
                for (int k = 0; k < K; k++) {
                    UN.set(ist, k, QN.get(ist, k));
                }
            } else if (sp == SchedStrategy.PS || sp == SchedStrategy.DPS || sp == SchedStrategy.GPS
                    || sp == SchedStrategy.PSPRIO || sp == SchedStrategy.DPSPRIO || sp == SchedStrategy.GPSPRIO
                    || sp == SchedStrategy.FCFS || sp == SchedStrategy.LCFS) {
                for (int k = 0; k < K; k++) {
                    Map<JobClass, MatrixCell> stationPH = PH.get(sn.stations.get(ist));
                    MatrixCell phEntry = stationPH != null ? stationPH.get(sn.jobclasses.get(k)) : null;
                    if (phEntry != null && !phEntry.isEmpty()) {
                        ColumnView dcolView = depRates.getColumnView(ind * K + k);
                        double throughput = pi.multColumnView(dcolView);
                        UN.set(ist, k, throughput / rates.get(ist, k) / S.get(ist));
                    }
                }
            }
        }

        for (int k = 0; k < K; k++) {
            for (int ist = 0; ist < M; ist++) {
                if (TN.get(ist, k) > 0) {
                    RN.set(ist, k, QN.get(ist, k) / TN.get(ist, k));
                } else {
                    RN.set(ist, k, 0.0);
                }
            }
            CN.set(0, k, NK.get(k) / XN.get(0, k));
        }

        QN.apply(Double.NaN, 0.0, "equal");
        UN.apply(Double.NaN, 0.0, "equal");
        RN.apply(Double.NaN, 0.0, "equal");
        XN.apply(Double.NaN, 0.0, "equal");
        TN.apply(Double.NaN, 0.0, "equal");
        CN.apply(Double.NaN, 0.0, "equal");

        return new SSAResult(QN, UN, RN, TN, CN, XN, null, null, sn);
    }
}
