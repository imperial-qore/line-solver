package jline.solvers.ctmc.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.MatFileUtils;
import jline.util.matrix.Matrix;

public class Solver_ctmc_marg {
    private final SolverCTMC solverCTMC;

    public Solver_ctmc_marg(SolverCTMC solverCTMC) {
        this.solverCTMC = solverCTMC;
    }

    public static Matrix solver_ctmc_marg(NetworkStruct sn, SolverOptions options) {
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
                MatFileUtils.saveCTMCWorkspace(SS, Q, null, fname);
            } catch (Exception e) {
                InputOutput.line_warning("solver_ctmc_marg", "Could not save workspace to .mat file: %s", e.getMessage());
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
        // see _kb/06-solver-catalog.md for rationale
        double pisum = pi.sumSubMatrix(0, pi.getNumRows(), 0, pi.getNumCols());
        for (int row = 0; row < pi.getNumRows(); row++) {
            for (int col = 0; col < pi.getNumCols(); col++) {
                pi.set(row, col, pi.get(row, col) / pisum);
            }
        }
        Map<Integer, Integer> statesz = new HashMap<Integer, Integer>();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                if (!statesz.containsKey(isf)) {
                    statesz.put(isf, sn.space.get(sn.stateful.get(isf)).getNumCols());
                }
            }
        }
        List<Integer> cstatesz = new ArrayList<Integer>();
        cstatesz.add(0);
        int cumulativeSum = 0;
        for (int i = 0; i < statesz.size(); i++) {
            cumulativeSum += statesz.get(i);
            cstatesz.add(cumulativeSum);
        }

        Matrix Pnir = new Matrix(1, sn.nstations);
        Pnir.zero();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                int ist = (int) sn.nodeToStation.get(ind);
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
                // see _kb/06-solver-catalog.md for rationale
                int colStart = cstatesz.get(isf);
                int colStop = colStart + state_i.getNumCols();
                Matrix SS_slice = new Matrix(SS.getNumRows(), state_i.getNumCols());
                Matrix.extract(SS, 0, SS.getNumRows(), colStart, colStop, SS_slice, 0, 0);

                List<Integer> rowmatched = Matrix.findRows(SS_slice, state_i);
                double sum = 0.0;
                for (Integer row : rowmatched) {
                    // see _kb/06-solver-catalog.md for rationale
                    sum += pi.get(row);
                }
                Pnir.set(0, ist, sum);
            }
        }
        long Tstop = System.nanoTime();
        double runtime = (Tstop - Tstart) / 1000000000.0;
        return Pnir;
    }
}
