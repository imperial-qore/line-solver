package jline.solvers.ctmc.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mc.Ctmc_solve;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.MatFileUtils;
import jline.util.matrix.Matrix;

public class Solver_ctmc_joint {
    private final SolverCTMC solverCTMC;

    public Solver_ctmc_joint(SolverCTMC solverCTMC) {
        this.solverCTMC = solverCTMC;
    }

    public static SolverCTMC.SolverCtmcJointResult solver_ctmc_joint(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;
        String fname = "";
        long T0 = System.nanoTime();
        ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
        Matrix SS = solverCTMCResult.getStateSpace();
        Matrix Q = solverCTMCResult.getQ();
        if (options.keep) {
            try {
                MatFileUtils.ensureWorkspaceDirectoryExists();
                fname = MatFileUtils.genFilename("workspace");
                MatFileUtils.saveCTMCWorkspace(SS, Q, null, fname);
            } catch (Exception e) {
                InputOutput.line_warning("solver_ctmc_joint", "Could not save workspace to .mat file: %s", e.getMessage());
                fname = "";
            }
        }
        Matrix pi = jline.solvers.ctmc.CtmcStationary.solve(Q, SS, sn, options);
        for (int row = 0; row < pi.getNumRows(); row++) {
            for (int col = 0; col < pi.getNumCols(); col++) {
                if (pi.get(row, col) < 0) {
                    pi.set(row, col, 0);
                }
            }
        }
        List<Matrix> statevecList = new ArrayList<Matrix>();
        // Iterate over nodes (not stations) to match MATLAB solver_ctmc_joint.m
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                int requiredlength = sn.space.get(sn.stateful.get(isf)).getNumCols();
                Matrix stateMatrix = sn.state.get(sn.stateful.get(isf));
                int currentlength = stateMatrix == null ? 0 : stateMatrix.length();
                Matrix state_i = new Matrix(1, requiredlength);
                int numZeros = requiredlength - currentlength;
                for (int col = 0; col < numZeros; col++) {
                    state_i.set(0, col, 0);
                }
                for (int col = 0; col < currentlength; col++) {
                    state_i.set(0, numZeros + col, stateMatrix.get(0, col));
                }
                statevecList.add(state_i);
            }
        }
        int totalCols = 0;
        for (Matrix m : statevecList) totalCols += m.getNumCols();
        Matrix statevec = new Matrix(1, totalCols);
        int currentCol = 0;
        for (Matrix matrix : statevecList) {
            for (int col = 0; col < matrix.getNumCols(); col++) {
                statevec.set(0, currentCol++, matrix.get(0, col));
            }
        }
        List<Integer> PnirIndex = Matrix.findRows(SS, statevec);
        Matrix Pnir = new Matrix(PnirIndex.size(), 1);
        for (int i = 0; i < PnirIndex.size(); i++) {
            int stateIdx = PnirIndex.get(i);
            Pnir.set(i, 0, pi.get(0, stateIdx));
        }

        double runtime = (System.nanoTime() - T0) / 1000000000.0;
        return new SolverCTMC.SolverCtmcJointResult(Pnir, runtime, fname);
    }
}
