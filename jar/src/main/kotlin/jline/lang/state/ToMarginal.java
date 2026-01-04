package jline.lang.state;

import static jline.GlobalConstants.Inf;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

public class ToMarginal implements Serializable {
    /**
     * Computes marginal statistics for a stateful node in a Stochastic Petri Net.
     * 
     * This method extracts marginal statistics from the global state space, focusing on
     * a specific node and computing its local state information including job populations,
     * server allocations, and phase distributions. It handles different node types
     * (transitions, places, caches) and their specific state representations.
     * 
     * <p>For transition nodes with multi-phase service, it computes the distribution
     * of jobs across different service phases. For station nodes, it aggregates
     * queue lengths and server utilizations. The method is essential for event
     * processing in SPN models where local node state must be extracted from
     * the global state space.
     * 
     * @param sn Network structure containing the complete SPN model definition
     * @param ind Index of the node for which marginal statistics are computed
     * @param state_i Current state vector for this node from the global state space
     * @param phasesz Number of phases for each service mode (for transition nodes)
     * @param phaseshift Cumulative phase counts for indexing (for transition nodes)
     * @param space_buf Buffer space matrix containing queue state definitions
     * @param space_srv Server space matrix containing server allocation states
     * @param space_var Variable space matrix containing phase variable definitions
     * 
     * @return StateMarginalStatistics containing:
     *         <ul>
     *         <li>ni: Total number of jobs at the node</li>
     *         <li>nir: Number of jobs per class</li>
     *         <li>sir: Number of jobs in service per class</li>
     *         <li>kir: Per-phase job distribution for each class</li>
     *         </ul>
     * 
     * @see State.StateMarginalStatistics for the returned data structure
     * @see #toMarginalAggr for aggregated marginal computation variant
     */
    public static State.StateMarginalStatistics toMarginal(NetworkStruct sn, int ind, Matrix state_i, Matrix phasesz, Matrix phaseshift, Matrix space_buf, Matrix space_srv, Matrix space_var) {
        /* ------------------------------------------------------------------
         *  Cached (non-station) stateful node
         * ------------------------------------------------------------------ */
        if (sn.isstation.get(ind, 0) == 0 && sn.isstateful.get(ind, 0) > 0) {

            /* --- Transition node (multi-phase service) -------------------- */
            if (sn.nodetype.get(ind) == NodeType.Transition) {

                /* number of service modes for this transition */
                int R = ((TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).nmodes;
                int numRows = state_i.getNumRows();

                // kir[mode] is a matrix of size numRows × numPhases for that mode
                // MATLAB: kir = zeros(size(state_i,1), R, max(phasesz))
                Matrix sir = new Matrix(numRows, R);
                List<Matrix> kir = new ArrayList<Matrix>();

                for (int r = 0; r < R; r++) {
                    int numPh = (int) phasesz.get(0, r);
                    Matrix kir_r = new Matrix(numRows, numPh);

                    for (int k = 0; k < numPh; k++) {
                        int col = (int) phaseshift.get(0, r) + k;
                        for (int row = 0; row < numRows; row++) {
                            double val = space_srv.get(row, col);
                            kir_r.set(row, k, val);
                            sir.set(row, r, sir.get(row, r) + val);
                        }
                    }
                    kir.add(kir_r);
                }

                Matrix nir = sir.copy();
                Matrix ni = nir.sumRows();

                return new State.StateMarginalStatistics(ni, nir, sir, kir);
            }

            /* --- Cache node ------------------------------------------------ */
            if (sn.nodetype.get(ind) == NodeType.Cache) {
                // For cache nodes, the state vector contains:
                // [job counts for each class, cache item states]
                // The job counts are the first sn.nclasses elements
                int R = sn.nclasses;
                Matrix nir = new Matrix(1, R);
                double niVal = 0.0;
                
                // Extract job counts per class (first R elements of state_i)
                for (int r = 0; r < R && r < state_i.length(); r++) {
                    double v = state_i.get(r);
                    nir.set(0, r, v);
                    niVal += v;
                }
                
                Matrix ni = new Matrix(1, 1);
                ni.set(0, 0, niVal);
                
                Matrix sir = nir.copy();  // for cache, all jobs are "in service"
                List<Matrix> kir = new ArrayList<Matrix>();
                kir.add(sir);
                
                return new State.StateMarginalStatistics(ni, nir, sir, kir);
            }

            /* --- Generic stateful (non-transition) node ------------------- */
            double nvarsSum = sn.nvars.sumRows(ind);
            int nCols = (int) (state_i.length() - nvarsSum);

            // Guard against negative column count
            if (nCols < 0) {
                line_error(mfilename(new Object[]{}), "Warning: nCols calculated as negative: " + nCols + 
                                 " (state_i.length=" + state_i.length() + ", nvarsSum=" + nvarsSum + ")");
                nCols = 0;
            }

            Matrix nir = new Matrix(1, Math.max(1, nCols));
            double niVal = 0.0;
            for (int i = 0; i < nCols && i < state_i.length(); i++) {
                double v = state_i.get(i);
                nir.set(0, i, v);
                niVal += v;
            }
            Matrix ni = new Matrix(1, 1);
            ni.set(0, 0, niVal);

            Matrix sir = nir;                // all jobs are “in service”
            List<Matrix> kir = new ArrayList<Matrix>();
            kir.add(sir);

            return new State.StateMarginalStatistics(ni, nir, sir, kir);
        }


        int R = sn.nclasses;
        int ist = (int) sn.nodeToStation.get(ind);

        if (phasesz == null) {
            phasesz = new Matrix(1, sn.phasessz.getNumCols());
            Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), phasesz, 0, 0);
        }
        if (phaseshift == null) {
            phaseshift = new Matrix(1, sn.phaseshift.getNumCols());
            Matrix.extract(sn.phaseshift, ist, ist + 1, 0, sn.phaseshift.getNumCols(), phaseshift, 0, 0);
        }

        boolean isExponential = (phasesz.elementMax() == 1);

        if (space_var == null) {
            int col = (int) sn.nvars.sumRows(ind);
            
            if (col == 0) {
                // No variables - create empty matrix
                space_var = new Matrix(state_i.getNumRows(), 0);
            } else {
                // Validate bounds to prevent matrix extraction errors
                int varStartCol = state_i.getNumCols() - col;
                if (varStartCol < 0 || col < 0 || varStartCol >= state_i.getNumCols()) {
                    // Handle edge case where variable space dimensions are invalid
                    space_var = new Matrix(state_i.getNumRows(), Math.max(0, col));
                    space_var.zero();
                } else {
                    space_var = new Matrix(state_i.getNumRows(), col);
                    Matrix.extract(state_i, 0, state_i.getNumRows(), varStartCol, state_i.getNumCols(), space_var, 0, 0);
                }
            }
        }
        if (space_srv == null) {
            int sumPhasesz = (int) phasesz.elementSum();
            int sumNvars = (int) sn.nvars.sumRows(ind);
            
            // Validate bounds to prevent matrix extraction errors
            int srvStartCol = state_i.getNumCols() - sumPhasesz - sumNvars;
            int srvEndCol = state_i.getNumCols() - sumNvars;
            
            if (srvStartCol < 0 || srvEndCol > state_i.getNumCols() || srvStartCol >= srvEndCol) {
                // Handle edge case where state_i dimensions don't match expected layout
                // This can happen with degenerate models or edge cases in state space generation
                space_srv = new Matrix(state_i.getNumRows(), Math.max(0, sumPhasesz));
                space_srv.zero();
            } else {
                space_srv = new Matrix(state_i.getNumRows(), sumPhasesz);
                Matrix.extract(state_i, 0, state_i.getNumRows(), srvStartCol, srvEndCol, space_srv, 0, 0);
            }
        }
        if (space_buf == null) {
            int col = state_i.getNumCols() - (int) (phasesz.elementSum() + sn.nvars.sumRows(ind));
            
            // Validate bounds to prevent matrix extraction errors
            if (col < 0 || col > state_i.getNumCols()) {
                // Handle edge case where calculated buffer size is invalid
                space_buf = new Matrix(state_i.getNumRows(), Math.max(0, col));
                space_buf.zero();
            } else {
                space_buf = new Matrix(state_i.getNumRows(), col);
                Matrix.extract(state_i, 0, state_i.getNumRows(), 0, col, space_buf, 0, 0);
            }
        }

        Matrix nir = new Matrix(state_i.getNumRows(), R);
        Matrix sir = new Matrix(state_i.getNumRows(), R);
        List<Matrix> kir = new ArrayList<Matrix>();
        if (isExponential) {
            sir = space_srv;
            kir.add(space_srv);
        } else {
            // Initialize kir
            for (int i = 0; i < phasesz.elementMax(); i++)
                kir.add(new Matrix(state_i.getNumRows(), R));

            for (int r = 0; r < R; r++) {
                for (int k = 0; k < phasesz.get(r); k++) {
                    // kir(:,r,k) = space_srv(:,phaseshift(r)+k);
                    Matrix tmp_kir = kir.get(k);
                    Matrix.extract(space_srv, 0, space_srv.getNumRows(), (int) phaseshift.get(r) + k, (int) phaseshift.get(r) + k + 1, tmp_kir, 0, r);

                    // sir(:,r) = sir(:,r) + kir(:,r,k);
                    for (int i = 0; i < sir.getNumRows(); i++)
                        sir.set(i, r, sir.get(i, r) + tmp_kir.get(i, r));
                }
            }
        }

        switch (sn.sched.get(sn.stations.get(ist))) {
            case INF:
            case PS:
            case DPS:
            case GPS:
            case LPS:
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
                nir = sir.copy();
                break;
            case EXT:
                nir.fill(Inf);
                break;
            case FCFS:
            case HOL:
            case FCFSPRIO:
            case LCFS:
            case LCFSPRIO:
                for (int r = 0; r < R; r++) {
                    Matrix sumval;
                    if (space_buf.getNumRows() == 0 || space_buf.getNumCols() == 0 || (space_buf.getNumRows() == 1 && space_buf.getNumCols() == 1 && space_buf.value() == 0)) {
                        sumval = new Matrix(1, 1);
                        sumval.set(0, 0, 0);
                    } else {
                        // +1 since loop starts from 0 but classes start from 1
                        sumval = space_buf.countEachRow(r + 1);
                    }
                    for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r) + sumval.get(i));
                }
                break;
            case LCFSPI:
            case LCFSPR:
            case LCFSPIPRIO:
            case LCFSPRPRIO:
                if (space_buf.length() > 1) {
                    // space_buf = space_buf(1:2:end);
                    Matrix sub_space_buf = new Matrix(1, (space_buf.getNumCols() * space_buf.getNumRows() + 1) / 2);
                    for (int i = 0; i < sub_space_buf.getNumCols(); i++)
                        sub_space_buf.set(0, i, space_buf.get(i * 2));

                    for (int r = 0; r < R; r++) {
                        Matrix sumval = sub_space_buf.countEachRow(r);
                        for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r) + sumval.get(i));
                    }
                } else {
                    nir = sir.copy();
                }
                break;
            case SIRO:
            case SEPT:
            case LEPT:
            case SRPT:
            case SRPTPRIO:
                for (int r = 0; r < R; r++) {
                    for (int i = 0; i < sir.getNumRows(); i++)
                        nir.set(i, r, sir.get(i, r) + space_buf.get(i, r));
                }
                break;
            case POLLING:
                // For polling servers, treat like FCFS with buffer counting
                for (int r = 0; r < R; r++) {
                    Matrix sumval;
                    if (space_buf.getNumRows() == 0 || space_buf.getNumCols() == 0 || (space_buf.getNumRows() == 1 && space_buf.getNumCols() == 1 && space_buf.value() == 0)) {
                        sumval = new Matrix(1, 1);
                        sumval.set(0, 0, 0);
                    } else {
                        // +1 since loop starts from 0 but classes start from 1
                        sumval = space_buf.countEachRow(r + 1);
                    }
                    for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r) + sumval.get(i));
                }
                break;
            default:
                for (int r = 0; r < R; r++) {
                    for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r));
                }
        }

        if (sn.nodetype.get(ind) != NodeType.Place) {
            for (int r = 0; r < R; r++) {
                if (Double.isNaN(sn.rates.get(ist, r))) {
                    for (int i = 0; i < nir.getNumRows(); i++) {
                        nir.remove(i, r);
                        sir.remove(i, r);
                    }
                    for (int k = 0; k < phasesz.get(r); k++) {
                        for (int j = 0; j < kir.size(); j++) {
                            kir.get(j).remove(r, k);
                        }
                    }
                }
            }
        }

        Matrix ni = nir.sumRows();

        return new State.StateMarginalStatistics(ni, nir, sir, kir);
    }

    /**
     * Computes aggregated marginal statistics for a stateful node using specified aggregation weights.
     * 
     * This method is a specialized variant of {@link #toMarginal} that applies aggregation weights
     * to compute weighted marginal statistics. It is particularly useful in event processing
     * where multiple state configurations need to be combined with specific probabilities or
     * weights, such as during transition enabling calculations or place token aggregation.
     * 
     * <p>The aggregation weights K and Ks determine how different components of the state
     * space are combined to produce the final marginal statistics. This is essential for
     * handling complex state dependencies in SPN models where the effective state at a
     * node depends on weighted contributions from multiple sources.
     * 
     * @param sn Network structure containing the complete SPN model definition
     * @param ind Index of the node for which aggregated marginal statistics are computed
     * @param state_i Current state vector for this node from the global state space
     * @param K Aggregation weight matrix for combining state components
     * @param Ks Cumulative aggregation weights for proper normalization
     * @param space_buf Buffer space matrix containing queue state definitions
     * @param space_srv Server space matrix containing server allocation states
     * @param space_var Variable space matrix containing phase variable definitions
     * 
     * @return StateMarginalStatistics containing weighted/aggregated:
     *         <ul>
     *         <li>ni: Total weighted number of jobs at the node</li>
     *         <li>nir: Weighted number of jobs per class</li>
     *         <li>sir: Weighted number of jobs in service per class</li>
     *         <li>kir: Weighted per-phase job distribution for each class</li>
     *         </ul>
     * 
     * @see State.StateMarginalStatistics for the returned data structure
     * @see #toMarginal for the non-aggregated variant
     * @see State#handleEnableEvent for usage in transition enabling calculations
     */
    public static State.StateMarginalStatistics toMarginalAggr(NetworkStruct sn, int ind, Matrix state_i, Matrix K, Matrix Ks, Matrix space_buf, Matrix space_srv, Matrix space_var) {

        int ist = (int) sn.nodeToStation.get(ind);
        int R = sn.nclasses;

        if (sn.isstation.get(ind, 0) == 0 && sn.isstateful.get(ind, 0) > 0) {
            // For non-station stateful nodes, return the sum of all state elements except the last nvars elements
            double nvarsSum = sn.nvars.sumRows(ind);
            int nCols = (int) (state_i.getNumCols() - nvarsSum);
            
            // Guard against negative column count
            if (nCols < 0) {
                nCols = 0;
            }
            
            Matrix nir = new Matrix(state_i.getNumRows(), Math.max(1, nCols));
            Matrix ni = new Matrix(state_i.getNumRows(), 1);
            
            for (int i = 0; i < state_i.getNumRows(); i++) {
                double niVal = 0.0;
                for (int j = 0; j < nCols && j < state_i.getNumCols(); j++) {
                    double v = state_i.get(i, j);
                    if (j < nir.getNumCols()) {
                        nir.set(i, j, v);
                    }
                    niVal += v;
                }
                ni.set(i, 0, niVal);
            }
            
            return new State.StateMarginalStatistics(ni, nir, null, null);
        }


        if (K == null) {
            K = new Matrix(1, sn.phasessz.getNumCols());
            Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), K, 0, 0);
        }
        if (Ks == null) {
            Ks = new Matrix(1, sn.phaseshift.getNumCols());
            Matrix.extract(sn.phaseshift, ist, ist + 1, 0, sn.phaseshift.getNumCols(), Ks, 0, 0);
        }

        if (space_var == null) {
            int col = (int) sn.nvars.sumRows(ind);
            if (col == 0) {
                // No variables - create empty matrix
                space_var = new Matrix(state_i.getNumRows(), 0);
            } else {
                space_var = new Matrix(state_i.getNumRows(), col);
                Matrix.extract(state_i, 0, state_i.getNumRows(), state_i.getNumCols() - col, state_i.getNumCols(), space_var, 0, 0);
            }
        }

        if (space_srv == null) {
            int sumK = (int) K.elementSum();
            int sumNvars = (int) sn.nvars.sumRows(ind);
            int startCol = state_i.getNumCols() - sumK;
            int endCol = state_i.getNumCols() - sumNvars;
            int numCols = endCol - startCol;
            space_srv = new Matrix(state_i.getNumRows(), numCols);
            Matrix.extract(state_i, 0, state_i.getNumRows(), startCol, endCol, space_srv, 0, 0);
        }

        if (space_buf == null) {
            int col = state_i.getNumCols() - (int) (K.elementSum());
            space_buf = new Matrix(state_i.getNumRows(), col);
            Matrix.extract(state_i, 0, state_i.getNumRows(), 0, col, space_buf, 0, 0);
        }
        Matrix nir = new Matrix(state_i.getNumRows(), R);
        nir.zero();
        for (int r = 0; r < R; r++) {
            for (int k = 0; k < K.get(r); k++) {
                Matrix tmp_kir = new Matrix(state_i.getNumRows(), 1);
                Matrix.extract(space_srv, 0, space_srv.getNumRows(), (int) Ks.get(r) + k, (int) Ks.get(r) + k + 1, tmp_kir, 0, 0);
                for (int i = 0; i < nir.getNumRows(); i++) {
                    nir.set(i, r, nir.get(i, r) + tmp_kir.get(i));
                }
            }
        }

        // MATLAB LINE does not handle LCFSPR, INF, PS, DPS, GPS cases
        // nir: class-r jobs at the station
        switch (sn.sched.get(sn.stations.get(ist))) {
            case EXT:
                nir.fill(Inf);
                break;
            case FCFS:
            case HOL:
            case FCFSPRIO:
            case LCFS:
            case LCFSPRIO:
                for (int r = 0; r < R; r++) {
                    Matrix sumval;
                    if (space_buf.getNumRows() == 0 || space_buf.getNumCols() == 0 || (space_buf.getNumRows() == 1 && space_buf.getNumCols() == 1 && space_buf.value() == 0)) {
                        sumval = new Matrix(1, 1);
                        sumval.set(0, 0, 0);
                    } else {
                        // +1 since loop starts from 0 but classes start from 1
                        sumval = space_buf.countEachRow(r + 1);
                    }
                    for (int i = 0; i < nir.getNumRows(); i++) nir.set(i, r, nir.get(i, r) + sumval.get(i));
                }
                break;
            case SIRO:
            case SEPT:
            case LEPT:
            case SRPT:
            case SRPTPRIO:
                for (int r = 0; r < R; r++) {
                    for (int i = 0; i < nir.getNumRows(); i++)
                        nir.set(i, r, nir.get(i, r) + space_buf.get(i, r));
                }
                break;
            default:
                // do nothing, no-op case
        }

        for (int r = 0; r < R; r++) {
            if (Double.isNaN(sn.rates.get(ist, r)) && sn.nodetype.get(ind) != NodeType.Place) { // if station disabled
                for (int i = 0; i < nir.getNumRows(); i++) {
                    nir.remove(i, r);
                }
            }
        }

        Matrix ni = nir.sumRows();
        return new State.StateMarginalStatistics(ni, nir, null, null);

    }
}