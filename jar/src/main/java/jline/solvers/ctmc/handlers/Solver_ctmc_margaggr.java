package jline.solvers.ctmc.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.api.mc.Ctmc_solve_reducible;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.VerboseLevel;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.ResultCTMCMargAggr;
import jline.util.MatFileUtils;
import jline.util.matrix.Matrix;

public final class Solver_ctmc_margaggr {
    private Solver_ctmc_margaggr() {}

    public static ResultCTMCMargAggr solver_ctmc_margaggr(NetworkStruct sn, SolverOptions options) {
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
                // genFilename already ends in .mat, and the MATLAB twin routes this
                // through line_printf, so it stays silent at VerboseLevel.SILENT.
                if (options.verbose != VerboseLevel.SILENT) {
                    System.out.println("\nCTMC generator and state space saved in: " + fname);
                }
            } catch (Exception e) {
                InputOutput.line_warning("solver_ctmc_margaggr",
                        "Could not save workspace to .mat file: %s", e.getMessage());
                fname = "";
            }
        }

        jline.util.Pair<Matrix, java.util.List<java.util.List<Integer>>> ctmcRes =
                Ctmc_solve_reducible.ctmc_solve_reducible(Q);
        Matrix pi = ctmcRes.getFirst();
        for (int row = 0; row < pi.getNumRows(); row++) {
            for (int col = 0; col < pi.getNumCols(); col++) {
                if (pi.get(row, col) < GlobalConstants.Zero) {
                    pi.set(row, col, 0.0);
                }
            }
        }

        List<Integer> statesz = new ArrayList<Integer>();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                statesz.add(sn.space.get(sn.stateful.get(isf)).getNumCols());
            }
        }

        List<Integer> cstatesz = new ArrayList<Integer>();
        cstatesz.add(0);
        int cumulativeSum = 0;
        for (Integer size : statesz) {
            cumulativeSum += size;
            cstatesz.add(cumulativeSum);
        }

        Matrix Pnir = new Matrix(1, sn.nstations);
        Pnir.zero();

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                int isf = (int) sn.nodeToStateful.get(ind);
                int ist = (int) sn.nodeToStation.get(ind);

                Matrix nodeState = sn.state.get(sn.stateful.get(isf));
                if (nodeState != null) {
                    int stateLength = sn.space.get(sn.stateful.get(isf)).getNumCols();

                    State.StateMarginalStatistics margStats =
                            ToMarginal.toMarginal(sn, ind, nodeState, null, null, null, null, null);
                    Matrix nivec = margStats.nir;

                    Pnir.set(0, ist, 0.0);
                    for (int s = 0; s < SS.getNumRows(); s++) {
                        int colStart = cstatesz.get(isf);
                        int colEnd = cstatesz.get(isf) + stateLength;
                        Matrix stateSlice = Matrix.extract(SS, s, s + 1, colStart, colEnd);

                        State.StateMarginalStatistics sivecStats =
                                ToMarginal.toMarginal(sn, ind, stateSlice, null, null, null, null, null);
                        Matrix sivec = sivecStats.nir;

                        boolean matches = true;
                        if (sivec.getNumCols() != nivec.getNumCols()) {
                            matches = false;
                        } else {
                            for (int col = 0; col < sivec.getNumCols(); col++) {
                                if (Math.abs(sivec.get(0, col) - nivec.get(0, col)) > 1e-10) {
                                    matches = false;
                                    break;
                                }
                            }
                        }

                        if (matches) {
                            Pnir.set(0, ist, Pnir.get(0, ist) + pi.get(s));
                        }
                    }
                }
            }
        }

        long Tstop = System.nanoTime();
        double runtime = (Tstop - Tstart) / 1000000000.0;

        return new ResultCTMCMargAggr(Pnir, pi, runtime, fname);
    }
}
