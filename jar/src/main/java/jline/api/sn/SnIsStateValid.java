/**
 * Stochastic Network State Validation Utility
 *
 * Validates the consistency and feasibility of queueing network states by checking
 * job conservation laws, capacity constraints, and scheduling discipline requirements.
 * Essential for ensuring model correctness before numerical analysis.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.util.matrix.Matrix;

public final class SnIsStateValid {
    private SnIsStateValid() {}

    /**
     * Checks if the network state is valid.
     *
     * @param sn the NetworkStruct object for the queueing network model
     * @return true if the state is valid, false otherwise
     */
    public static boolean snIsStateValid(NetworkStruct sn) {
        NetworkStruct snTmp = sn.copy();
        Matrix nir = new Matrix(snTmp.nstations, snTmp.nclasses);
        Matrix sir = new Matrix(snTmp.nstations, snTmp.nclasses);

        for (int ist = 0; ist < snTmp.nstations; ist++) {
            // statefulIdx = (int) snTmp.stationToStateful.get(ist); // unused, kept for parity
            if (snTmp.state.get(snTmp.stations.get(ist)).getNumRows() > 1) {
                if (snTmp.stateprior.get(snTmp.stations.get(ist)).elementMax() < 1 - GlobalConstants.FineTol) {
                    line_warning(mfilename(new Object()),
                            "isStateValid will ignore some states of station %d, define a unique initial state to address this problem.\n",
                            ist);
                }
                Matrix initialState = new Matrix(1, snTmp.state.get(snTmp.stations.get(ist)).getNumCols());
                Matrix.extractRows(snTmp.state.get(snTmp.stations.get(ist)), 0, 1, initialState);
                snTmp.state.put(snTmp.stations.get(ist), initialState);
            }

            State.StateMarginalStatistics stats = ToMarginal.toMarginal(snTmp,
                    (int) snTmp.stationToNode.get(0, ist),
                    snTmp.state.get(snTmp.stations.get(ist)),
                    null,
                    null,
                    null,
                    null,
                    null);

            for (int i = 0; i < snTmp.nclasses; i++) {
                nir.set(ist, i, stats.nir.get(0, i));
                sir.set(ist, i, stats.sir.get(0, i));
            }
        }

        return State.isValid(snTmp, nir, sir);
    }
}
