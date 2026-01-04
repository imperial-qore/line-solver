package jline.lang.state;

import jline.lang.GlobalSync;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class AfterGlobalEvent implements Serializable {
    /**
     * Processes a global event in a Stochastic Petri Net (SPN) and computes the resulting state space.
     * 
     * This is the main entry point for handling global synchronization events in SPN models,
     * particularly for transition firings that affect multiple places simultaneously. The method
     * orchestrates the complete event processing workflow including both ENABLE and FIRE phases
     * of transition execution.
     * 
     * <p>The method handles two main types of global events:
     * <ul>
     * <li><b>ENABLE events:</b> Check if transitions can be enabled based on token availability
     *     in input places and generate all possible enabling combinations</li>
     * <li><b>FIRE events:</b> Execute transition firings, consume tokens from input places,
     *     produce tokens in output places, and update server/phase states</li>
     * </ul>
     * 
     * <p>For simulation mode, the method supports stochastic selection of outcomes when multiple
     * firing possibilities exist, using proper probability distributions for realistic behavior.
     * 
     * @param sn Network structure containing the complete SPN model definition including
     *           nodes, arcs, job classes, and transition parameters
     * @param ind Node index of the transition being processed (must be a Transition node)
     * @param glspace Global state space represented as a list of state matrices, one for each
     *                stateful node in the network (places and transitions with servers)
     * @param glevent Global synchronization event containing active events (the main transition
     *                firing) and passive events (input/output place interactions)
     * @param isSimulation If true, enables stochastic selection for simulation; if false,
     *                     returns all possible outcomes for analytical computation
     * 
     * @return AfterGlobalEventResult containing:
     *         <ul>
     *         <li>outglspace: Updated global state space after event processing</li>
     *         <li>outrates: Transition rates for each resulting state</li>
     *         <li>outprobs: Probabilities for each resulting state transition</li>
     *         </ul>
     * 
     * @throws IllegalArgumentException if the node is not a valid transition or if the
     *                                  global event structure is malformed
     * @throws IllegalStateException if the transition cannot be enabled or fired from
     *                               the current state
     * 
     * @see GlobalSync for event structure and synchronization semantics
     * @see State#handleEnableEvent for ENABLE phase processing details
     * @see State#handleFireEvent for FIRE phase processing details
     * @see TransitionNodeParam for transition-specific configuration
     */
    public static AfterGlobalEventResult afterGlobalEvent(NetworkStruct sn, int ind, List<Matrix> glspace, GlobalSync glevent, boolean isSimulation) {
        List<Matrix> outglspace = new ArrayList<Matrix>();
        for (Matrix m : glspace) {
            outglspace.add(m.copy());
        }
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(0, 0);

        int R = sn.nclasses;

        if (sn.nodetype.get(ind) == NodeType.Transition && glevent.getActive().size() > 0) {
            int isf = (int) sn.nodeToStateful.get(ind);
            Matrix inspace = glspace.get(isf);
            double V = sn.nvars.getRow(ind).elementSum();

            // Extract state components
            int varStartCol = (int) (inspace.getNumCols() - V);
            if (varStartCol < 0) {
                throw new IllegalStateException("Invalid varStartCol calculation: inspace.getNumCols()=" + 
                    inspace.getNumCols() + ", V=" + V + ", resulting in negative varStartCol=" + varStartCol);
            }
            Matrix spaceVar = Matrix.extractColumns(inspace, varStartCol, inspace.getNumCols(), null);

            if (sn.nodetype.get(ind) == NodeType.Transition) {
                TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                Matrix fK = transParam.firingphases;
                int nmodes = transParam.nmodes;

                // Handle NaN firingphases (non-phase-type distributions like Pareto)
                // Infer phase count from D0 matrix size
                if (fK.hasNaN()) {
                    fK = new Matrix(1, nmodes);
                    for (int m = 0; m < nmodes; m++) {
                        if (transParam.firingproc != null && transParam.firingproc.containsKey(m)
                            && transParam.firingproc.get(m) != null && !transParam.firingproc.get(m).isEmpty()) {
                            Matrix D0 = transParam.firingproc.get(m).get(0);
                            if (D0 != null) {
                                fK.set(0, m, D0.getNumRows());
                            } else {
                                fK.set(0, m, 1);
                            }
                        } else {
                            fK.set(0, m, 1);
                        }
                    }
                }

                Matrix fKs = new Matrix(1, fK.getNumCols() + 1);
                fKs.set(0, 0, 0);
                for (int i = 0; i < fK.getNumCols(); i++) {
                    fKs.set(0, i + 1, fKs.get(0, i) + fK.get(0, i));
                }
                int mode = glevent.getActive().get(0).getMode();
                
                // Validate bounds for spaceBuf extraction
                if (nmodes > inspace.getNumCols()) {
                    throw new IllegalStateException("Cannot extract spaceBuf: nmodes=" + nmodes + 
                        " exceeds inspace.getNumCols()=" + inspace.getNumCols());
                }
                Matrix spaceBuf = Matrix.extractColumns(inspace, 0, nmodes, null);
                
                // Validate bounds for spaceSrv extraction
                int srvEndCol = (int) (nmodes + fK.elementSum());
                if (srvEndCol > inspace.getNumCols()) {
                    throw new IllegalStateException("Cannot extract spaceSrv: end column " + srvEndCol + 
                        " (nmodes=" + nmodes + " + fK.elementSum()=" + fK.elementSum() + 
                        ") exceeds inspace.getNumCols()=" + inspace.getNumCols());
                }
                Matrix spaceSrv = Matrix.extractColumns(inspace, nmodes, srvEndCol, null);

                State.EventHandleResult result = null;
                switch (glevent.getActive().get(0).getEvent()) {
                    case ENABLE:
                        result = State.handleEnableEvent(sn, ind, glevent, glspace, outglspace, inspace, spaceBuf, spaceSrv,
                                spaceVar, fK, fKs, mode, transParam, R);
                        break;

                    case FIRE:
                        result = State.handleFireEvent(sn, ind, glevent, glspace, outglspace, inspace, spaceBuf, spaceSrv,
                                spaceVar, fK, fKs, mode, transParam, R);
                        break;
                }
                if (result != null) {
                    outspace = result.outspace;
                    outrate = result.outrate;
                    outprob = result.outprob;
                }
            }
        }

        // Handle simulation case - select one outcome randomly
        if (isSimulation && outspace.getNumRows() > 1) {
            double totRate = outrate.elementSum();
            Matrix cumRate = outrate.cumsumViaCol();
            cumRate = Matrix.scaleMult(cumRate, 1.0 / totRate);
            double rand = Maths.rand();
            int firingCtr = 0;
            for (int i = 0; i < cumRate.getNumRows(); i++) {
                if (rand > cumRate.get(i, 0)) {
                    firingCtr = i + 1;
                }
            }
            if (firingCtr >= outspace.getNumRows()) {
                firingCtr = outspace.getNumRows() - 1;
            }
            outspace = Matrix.extractRows(outspace, firingCtr, firingCtr + 1, null);
            outrate = new Matrix(1, 1);
            outrate.set(0, 0, totRate);
            outprob = Matrix.extractRows(outprob, firingCtr, firingCtr + 1, null);

            // Update outglspace with selected state
            if (!outspace.isEmpty()) {
                int isf = (int) sn.nodeToStateful.get(ind);
                outglspace.set(isf, outspace);
            }
        }
        
        // Ensure outrate and outprob are not empty
        if (outrate.isEmpty()) {
            outrate = new Matrix(1, 1);
            outrate.set(0, 0, 0.0);
        }
        if (outprob.isEmpty()) {
            outprob = new Matrix(1, 1);
            outprob.set(0, 0, 0.0);
        }

        return new AfterGlobalEventResult(outglspace, outrate, outprob);
    }

    /**
     * Result container for global event processing in Stochastic Petri Net models.
     * 
     * This class encapsulates the complete outcome of processing a global synchronization
     * event, containing the updated state space and the associated rates and probabilities
     * for all possible resulting states after event execution.
     * 
     * @see AfterGlobalEvent#afterGlobalEvent for the method that produces this result
     */
    public static class AfterGlobalEventResult {
        /**
         * Updated global state space after event processing.
         * Each matrix in the list represents the state space for a specific stateful node
         * in the network after the global event has been processed.
         */
        public final List<Matrix> outglspace;
        
        /**
         * Transition rates matrix for each resulting state.
         * Each row corresponds to a possible outcome state and contains the rate
         * at which that transition occurs (may be immediate for logical transitions).
         */
        public final Matrix outrate;
        
        /**
         * Probability matrix for each resulting state transition.
         * Each row corresponds to a possible outcome state and contains the probability
         * of that specific outcome occurring when the event is fired.
         */
        public final Matrix outprob;

        /**
         * Constructs a new result container for global event processing.
         * 
         * @param outglspace Updated global state space after event processing
         * @param outrate Transition rates for each resulting state
         * @param outprob Probabilities for each resulting state transition
         */
        public AfterGlobalEventResult(List<Matrix> outglspace, Matrix outrate, Matrix outprob) {
            this.outglspace = outglspace;
            this.outrate = outrate;
            this.outprob = outprob;
        }
    }
}