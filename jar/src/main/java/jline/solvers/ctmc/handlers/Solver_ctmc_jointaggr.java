package jline.solvers.ctmc.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.api.mc.Ctmc_solve_reducible;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.MatFileUtils;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public class Solver_ctmc_jointaggr {

    public static SolverCTMC.SolverCtmcJointResult solver_ctmc_jointaggr(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        String fname = "";
        long Tstart = System.nanoTime();

        ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
        Matrix Q = solverCTMCResult.getQ();
        Matrix SS = solverCTMCResult.getStateSpace();
        Matrix SSq = solverCTMCResult.getStateSpaceAggr();
        sn = solverCTMCResult.getSn();

        if (options.keep) {
            try {
                MatFileUtils.ensureWorkspaceDirectoryExists();
                fname = MatFileUtils.genFilename("workspace");
                MatFileUtils.saveCTMCWorkspace(SS, Q, SSq, fname);
            } catch (Exception e) {
                InputOutput.line_warning("solver_ctmc_jointaggr", "Could not save workspace to .mat file: %s", e.getMessage());
                fname = "";
            }
        }

        // Use ctmc_solve_reducible to handle reducible CTMCs (matches MATLAB)
        Pair<Matrix, List<List<Integer>>> sol = Ctmc_solve_reducible.ctmc_solve_reducible(Q);
        Matrix pi = sol.getLeft();
        for (int row = 0; row < pi.getNumRows(); row++) {
            for (int col = 0; col < pi.getNumCols(); col++) {
                if (pi.get(row, col) < GlobalConstants.Zero) {
                    pi.set(row, col, 0.0);
                }
            }
        }

        List<Double> nvec = new ArrayList<Double>();

        // Build nvec using ToMarginal to get nir (marginal job counts)
        for (int i = 0; i < sn.nstations; i++) {
            int nodeIdx = (int) sn.stationToNode.get(i);
            if (sn.isstateful.get(nodeIdx) != 0.0) {
                int isf = (int) sn.stationToStateful.get(i);
                Matrix nodeState = sn.state.get(sn.stateful.get(isf));
                if (nodeState != null) {
                    jline.lang.state.State.StateMarginalStatistics margStats = ToMarginal.toMarginal(sn, nodeIdx, nodeState, null, null, null, null, null);
                    Matrix nir = margStats.nir;

                    for (int col = 0; col < nir.getNumCols(); col++) {
                        nvec.add(nir.get(0, col));
                    }
                }
            }
        }

        // Convert nvec to matrix for comparison
        Matrix nvecMatrix = new Matrix(1, nvec.size());
        for (int i = 0; i < nvec.size(); i++) {
            nvecMatrix.set(0, i, nvec.get(i));
        }

        // Find matching rows in SSq (aggregate state space) using findrows logic
        double Pnir = 0.0;
        for (int row = 0; row < SSq.getNumRows(); row++) {
            Matrix rowMatrix = SSq.getRow(row);
            boolean matches = true;
            if (nvecMatrix.getNumCols() != rowMatrix.getNumCols()) {
                matches = false;
            } else {
                for (int col = 0; col < nvecMatrix.getNumCols(); col++) {
                    if (Math.abs(nvecMatrix.get(0, col) - rowMatrix.get(0, col)) > 1e-10) {
                        matches = false;
                        break;
                    }
                }
            }
            if (matches) {
                Pnir += pi.get(row);
            }
        }

        long Tstop = System.nanoTime();
        double runtime = (Tstop - Tstart) / 1000000000.0;

        Matrix result = new Matrix(1, 1);
        result.set(0, 0, Pnir);
        return new SolverCTMC.SolverCtmcJointResult(result, runtime, fname);
    }
}
