package jline.api.mc;

import jline.lang.NetworkStruct;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Ctmc_ssg {
    private Ctmc_ssg() {}

    public static SolverCTMC.CtmcSsgResult ctmc_ssg(NetworkStruct sn, SolverOptions options) {
        // Get properly dimensioned cutoff matrix
        Matrix cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses);
        State.StateSpaceGeneratorResult spaceResult = State.spaceGenerator(sn, cutoffMatrix, options);
        java.util.Map<jline.lang.nodes.StatefulNode, Matrix> nodeStateSpace = spaceResult.sn.space;
        sn.space = nodeStateSpace;
        sn.spaceHash = spaceResult.ST.spaceHash;
        Matrix stateSpace = spaceResult.SS;
        Matrix stateSpaceHashed = spaceResult.SSh;

        int nclasses = sn.nclasses;

        Matrix stateSpaceAggr = new Matrix(stateSpaceHashed.getNumRows(), stateSpaceHashed.getNumCols());

        MatrixCell stateCell = new MatrixCell();
        for (int s = 0; s < stateSpaceHashed.getNumRows(); s++) {
            Matrix state = stateSpaceHashed.getRow(s);
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.isstateful.get(ind) == 1.0) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    int nodeStateIdx = (int) state.get(isf);
                    Matrix nodeState = sn.space.get(sn.stateful.get(isf)).getRow(nodeStateIdx);
                    stateCell.set(isf, nodeState);
                    if (sn.isstation.get(ind) == 1.0) {
                        int ist = (int) sn.nodeToStation.get(ind);
                        State.StateMarginalStatistics marginalResult =
                                ToMarginal.toMarginal(sn, ind, stateCell.get(isf), null, null, null, null, null);
                        Matrix nir = marginalResult.nir;

                        int startCol = ist * nclasses;
                        stateSpaceAggr.expandMatrix(stateSpaceAggr.getNumRows(),
                                (ist + 1) * nclasses,
                                stateSpaceAggr.getNumNonZeros());
                        for (int i = 0; i < nir.getNumCols(); i++) {
                            stateSpaceAggr.set(s, startCol + i, nir.get(i));
                        }
                    }
                }
            }
        }
        return new SolverCTMC.CtmcSsgResult(stateSpace, stateSpaceAggr, stateSpaceHashed, nodeStateSpace, sn);
    }
}
