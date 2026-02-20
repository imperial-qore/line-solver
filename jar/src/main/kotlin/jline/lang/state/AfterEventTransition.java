package jline.lang.state;

import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.List;
import java.util.Map;

public class AfterEventTransition implements Serializable {
    static Ret.EventResult afterEventTransition(NetworkStruct sn, int ind, EventType event, int jobClass, boolean isSimulation,
                                                Matrix inspace, Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                                int M, int R,
                                                int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                                double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        switch (event) {
            case ENABLE:
                // no-op this is a global event
                outspace = inspace.copy();
                outrate = new Matrix(outspace.getNumRows(), 1);
                outrate.fill(0);
                outprob = new Matrix(outspace.getNumRows(), 1);
                outprob.ones();
                break;
            case FIRE:
                // no-op this is a global event
                outspace = inspace.copy();
                outrate = new Matrix(outspace.getNumRows(), 1);
                outrate.fill(0);
                outprob = new Matrix(outspace.getNumRows(), 1);
                outprob.ones();
                break;
            case PHASE:
                int mode = jobClass; // in a Transition, jobClass is interpreted as the mode
                outspace = new Matrix(0, 0);
                outrate = new Matrix(0, 0);
                outprob = new Matrix(0, 0);

                // Check bounds first - mode must be valid
                if (mode < 0 || mode >= K.getNumCols()) {
                    // Invalid mode index - return empty to indicate no phase transition possible
                    return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0));
                }

                // For single-phase modes (exponential), no phase transitions are possible
                if (K.get(mode) <= 1) {
                    // Return empty to indicate no phase transition is possible
                    return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0));
                }

                State.StateMarginalStatistics transitionStats = ToMarginal.toMarginal(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
                Matrix niTrans = transitionStats.ni;
                Matrix nirTrans = transitionStats.nir;
                List<Matrix> kirTrans = transitionStats.kir;

                // ni, nir, kir is the count of *enabled* servers while they progress their execution through the phases
                if (nirTrans.get(mode) > 0) {
                    // Check if we're in the last phase and there's no next phase to transition to
                    boolean hasValidTransitions = false;
                    
                    for (int k = 0; k < K.get(mode); k++) {
                        Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                            int colIndex = (int) (Ks.get(mode) + k);
                            if (spaceSrv.get(row, colIndex) > 0) {
                                en.set(row, 0, 1);
                            } else {
                                en.set(row, 0, 0);
                            }
                        }

                        boolean anyEn = false;
                        for (int row = 0; row < en.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                                anyEn = true;
                                break;
                            }
                        }

                        if (anyEn) {
                            for (int kdest = 0; kdest < K.get(mode); kdest++) {
                                if (kdest != k) { // setdiff(1:K(mode),k) equivalent
                                    double rateValue = 0;
                                    Matrix spaceSrvK = new Matrix(0, 0);
                                    Matrix spaceBufK = new Matrix(0, 0);
                                    Matrix spaceVarK = new Matrix(0, 0);

                                    // Extract enabled rows and track original row indices
                                    java.util.List<Integer> enabledRows = new java.util.ArrayList<Integer>();
                                    for (int row = 0; row < en.getNumRows(); row++) {
                                        if (en.get(row, 0) == 1) {
                                            enabledRows.add(row);
                                            Matrix spaceSrvRow = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                            Matrix spaceBufRow = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                            Matrix spaceVarRow = Matrix.extractRows(spaceVar, row, row + 1, null);

                                            if (spaceSrvK.isEmpty()) {
                                                spaceSrvK = spaceSrvRow.copy();
                                                spaceBufK = spaceBufRow.copy();
                                                spaceVarK = spaceVarRow.copy();
                                            } else {
                                                spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvRow, null);
                                                spaceBufK = Matrix.concatRows(spaceBufK, spaceBufRow, null);
                                                spaceVarK = Matrix.concatRows(spaceVarK, spaceVarRow, null);
                                            }
                                        }
                                    }

                                    // Update server state: remove from phase k, add to phase kdest
                                    for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                        spaceSrvK.set(row, (int) (Ks.get(mode) + k), spaceSrvK.get(row, (int) (Ks.get(mode) + k)) - 1);
                                        spaceSrvK.set(row, (int) (Ks.get(mode) + kdest), spaceSrvK.get(row, (int) (Ks.get(mode) + kdest)) + 1);
                                    }

                                    // Calculate rate using firing process - per row using kir
                                    TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                                    MatrixCell firingProcCell = null;
                                    if (transParam.firingproc != null) {
                                        // Find the Mode key that corresponds to this mode index
                                        for (Map.Entry<jline.lang.Mode, MatrixCell> entry : transParam.firingproc.entrySet()) {
                                            if (entry.getKey().getIndex() - 1 == mode) { // Mode indices are 1-based
                                                firingProcCell = entry.getValue();
                                                break;
                                            }
                                        }
                                    }

                                    // Compute per-row rates using kir values from original rows
                                    // MATLAB: rate = firingproc{mode}{1}(k,kdest) * kir(:,mode,k)
                                    Matrix rateMatrix = new Matrix(enabledRows.size(), 1);
                                    for (int i = 0; i < enabledRows.size(); i++) {
                                        int origRow = enabledRows.get(i);
                                        double kirValue = kirTrans.get(mode).get(origRow, k);

                                        if (firingProcCell != null && firingProcCell.get(0) != null) {
                                            Matrix firingProc = firingProcCell.get(0);
                                            double firingProcRate = firingProc.get(k, kdest);
                                            rateValue = firingProcRate * kirValue;
                                        } else {
                                            rateValue = kirValue;
                                        }

                                        // MATLAB: outrate = nir(mode) .* rate
                                        double finalRate = nirTrans.get(origRow, mode) * rateValue;
                                        rateMatrix.set(i, 0, finalRate);
                                    }

                                    // Build output states
                                    Matrix leftBottomTrans = Matrix.concatColumns(spaceBufK, spaceSrvK, null);
                                    Matrix bottomTrans = Matrix.concatColumns(leftBottomTrans, spaceVarK, null);
                                    outspace = Matrix.concatRows(outspace, bottomTrans, null);

                                    outrate = Matrix.concatRows(outrate, rateMatrix, null);

                                    Matrix probMatrix = new Matrix(spaceBufK.getNumRows(), 1);
                                    probMatrix.ones();
                                    outprob = Matrix.concatRows(outprob, probMatrix, null);
                                    hasValidTransitions = true;
                                }
                            }
                        }
                    }
                    
                    // If no valid transitions were found, return empty result (will be interpreted as rate 0)
                    if (!hasValidTransitions && outspace.isEmpty()) {
                        // Return empty matrices to indicate rate of 0
                        return new Ret.EventResult(inspace.copy(), Matrix.zeros(1, 1), Matrix.ones(1, 1));
                    }

                    if (isSimulation && eventCache.isEnabled()) {
                        eventCache.put(key, new Ret.EventResult(outspace, outrate, outprob));
                    }

                    if (isSimulation) {
                        if (outspace.getNumRows() > 1) {
                            Matrix tot_rate = outrate.sumCols();
                            Matrix cum_rate = Matrix.scaleMult(outrate.cumsumViaCol(), 1.0 / tot_rate.value());
                            int firing_ctr = -1;
                            double rand = Maths.rand();
                            for (int row = 0; row < cum_rate.getNumRows(); row++) {
                                if (rand > cum_rate.get(row, 0)) {
                                    firing_ctr = row;
                                }
                            }
                            firing_ctr++;
                            outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                            double outrate_val = outrate.sumCols().value();
                            outrate = new Matrix(1, 1);
                            outrate.set(0, 0, outrate_val);
                            outprob = Matrix.extractRows(outprob, firing_ctr, firing_ctr + 1, null);
                        }
                    }
                }
                break;
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}