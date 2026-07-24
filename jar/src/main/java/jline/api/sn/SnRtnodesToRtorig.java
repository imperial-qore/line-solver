package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.api.mc.Dtmc_stochcomp;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class SnRtnodesToRtorig {
    private SnRtnodesToRtorig() {}

    /**
     * Converts routing matrices from nodes to original format, specifically handling class switching nodes.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return a pair containing rtorigcell (K x K cell of csshift x csshift matrices) and rtorig
     */
    public static Pair<List<List<Matrix>>, Matrix> snRtnodesToRtorig(NetworkStruct sn) {
        int K = sn.nclasses;
        Matrix rtnodes = sn.rtnodes;

        // Find CS shift point - the index before the first class switching node
        int csshift = sn.nnodes;
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodenames.get(ind).startsWith("CS_")) {
                csshift = ind;
                break;
            }
        }

        // Build columns to keep for stochastic complement
        List<Integer> colToKeep = new ArrayList<Integer>();
        for (int ind = 0; ind < csshift; ind++) {
            for (int k = 0; k < K; k++) {
                colToKeep.add(ind * K + k);
            }
        }

        // Compute stochastic complement to remove class switching nodes
        Matrix rtorig = Dtmc_stochcomp.dtmc_stochcomp(rtnodes, colToKeep);

        // Initialize the cell array equivalent - K x K matrices, each of size csshift x csshift
        List<List<Matrix>> rtorigcell = new ArrayList<List<Matrix>>(K);
        for (int r = 0; r < K; r++) {
            List<Matrix> row = new ArrayList<Matrix>(K);
            for (int s = 0; s < K; s++) {
                row.add(new Matrix(csshift, csshift));
            }
            rtorigcell.add(row);
        }

        // Replace NaNs with 0 for cache routing probabilities
        for (int i = 0; i < rtorig.getNumRows(); i++) {
            for (int j = 0; j < rtorig.getNumCols(); j++) {
                if (Double.isNaN(rtorig.get(i, j))) {
                    rtorig.set(i, j, 0.0);
                }
            }
        }

        // Populate cell array with routing probabilities between classes
        for (int ind = 0; ind < csshift; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Sink) {
                for (int jnd = 0; jnd < csshift; jnd++) {
                    for (int r = 0; r < K; r++) {
                        for (int s = 0; s < K; s++) {
                            rtorigcell.get(r).get(s).set(ind, jnd, rtorig.get((ind * K) + r, (jnd * K) + s));
                        }
                    }
                }
            }
        }

        return new Pair<List<List<Matrix>>, Matrix>(rtorigcell, rtorig);
    }
}
