package jline.api.mc;

import jline.lang.NetworkStruct;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * CTMC State Space Generator for Reachability Analysis.
 */
public final class Ctmc_ssg_reachability {
    private Ctmc_ssg_reachability() {}

    public static CtmcSsgReachabilityResult ctmc_ssg_reachability(NetworkStruct sn, SolverOptions options) {
        Matrix cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses);
        State.StateSpaceGeneratorResult spaceResult = State.spaceGenerator(sn, cutoffMatrix, options);

        Matrix stateSpace = spaceResult.SS;
        Matrix stateSpaceHashed = spaceResult.SSh;
        Map<StatefulNode, Matrix> nodeStateSpace = spaceResult.sn.space;
        sn.space = nodeStateSpace;
        sn.spaceHash = spaceResult.ST.spaceHash;

        options.config.hide_immediate = true;

        int nclasses = sn.nclasses;
        int A = sn.sync.size();
        Matrix stateSpaceAggr = new Matrix(stateSpaceHashed.getNumRows(), stateSpaceHashed.getNumCols());
        stateSpaceAggr.zero();

        for (int a = 0; a < A; a++) {
            Map<Integer, Matrix> stateCell = new HashMap<Integer, Matrix>();
            for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
                Matrix state = stateSpaceHashed.getRow(s);

                for (int ind = 0; ind < sn.nnodes; ind++) {
                    if (sn.isstateful.get(ind) > 0) {
                        int isf = (int) sn.nodeToStateful.get(ind);
                        stateCell.put(isf, sn.space.get(sn.stateful.get(isf)).getRow((int) state.get(isf)));

                        if (sn.isstation.get(ind) > 0) {
                            int ist = (int) sn.nodeToStation.get(ind);
                            State.StateMarginalStatistics marginalResult =
                                    ToMarginal.toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null);
                            Matrix nir = marginalResult.nir;

                            int startCol = ist * nclasses;
                            int endCol = (ist + 1) * nclasses;

                            if (stateSpaceAggr.getNumCols() < endCol) {
                                stateSpaceAggr.expandMatrix(stateSpaceAggr.getNumRows(), endCol,
                                        stateSpaceAggr.getNumNonZeros());
                            }

                            for (int c = 0; c < nir.length(); c++) {
                                stateSpaceAggr.set(s, startCol + c, nir.get(c));
                            }
                        }
                    }
                }
            }
        }

        return new CtmcSsgReachabilityResult(stateSpace, stateSpaceAggr, stateSpaceHashed, nodeStateSpace, sn);
    }
}
