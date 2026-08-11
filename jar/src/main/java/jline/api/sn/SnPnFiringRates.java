/**
 * @file Recover per-mode transition firing rates of a Petri net from the Place throughputs
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.util.matrix.Matrix;

/**
 * Recovers per-mode transition firing rates from the Place throughputs.
 *
 * The firing rates of a Petri net are not carried by the network structure, but
 * they are determined by the Place throughputs together with the net structure.
 * Writing x for the vector of per-mode firing rates, two families of equations
 * hold at steady state, for every Place p and class k:
 *
 * <ul>
 * <li>departure: the sum over the modes consuming (p,k) of x, weighted by the
 * input arc multiplicity when tputIsTokens is true and unweighted when it is
 * false, equals TN(p,k)</li>
 * <li>balance: the sum over all modes of x times (produced minus consumed)
 * equals zero</li>
 * </ul>
 *
 * The system is solved in least squares. That is deliberate: an exact solver
 * supplies throughputs that satisfy it exactly and the fit is then the exact
 * answer, whereas a simulator supplies estimates that satisfy it only up to
 * sampling error and the least-squares fit is the right estimator there. A
 * residual test would reject every simulated run.
 *
 * Port of matlab/src/api/sn/sn_pn_firing_rates.m
 */
public final class SnPnFiringRates {
    private SnPnFiringRates() {}

    /**
     * Outcome of the recovery. The rates field is null when the firing rates
     * cannot be recovered, in which case the caller keeps the value it had.
     */
    public static final class Ret {
        /** Firing rate per (transition, mode) pair, as a column vector. */
        public final Matrix rates;
        /** Tokens consumed, indexed [mode][place][class]. */
        public final double[][][] consumed;
        /** Tokens produced, indexed [mode][place][class]. */
        public final double[][][] produced;
        /** Node indices of the Places, in the order used by the arrays above. */
        public final List<Integer> placeNodes;

        public Ret(Matrix rates, double[][][] consumed, double[][][] produced, List<Integer> placeNodes) {
            this.rates = rates;
            this.consumed = consumed;
            this.produced = produced;
            this.placeNodes = placeNodes;
        }
    }

    /**
     * Recovers the per-mode firing rates from the Place throughputs.
     *
     * @param sn network structure
     * @param TN average throughputs at stations
     * @param tputIsTokens true when TN counts tokens, false when it counts firing events
     * @return the recovery outcome, with a null rates field when undetermined
     */
    public static Ret snPnFiringRates(NetworkStruct sn, Matrix TN, boolean tputIsTokens) {
        Ret undetermined = new Ret(null, null, null, new ArrayList<Integer>());

        int R = sn.nclasses;
        if (TN == null || TN.isEmpty()) {
            return undetermined;
        }

        List<Integer> placeNodes = new ArrayList<Integer>();
        List<Integer> transNodes = new ArrayList<Integer>();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt == NodeType.Place) {
                placeNodes.add(ind);
            } else if (nt == NodeType.Transition) {
                transNodes.add(ind);
            }
        }
        if (placeNodes.isEmpty() || transNodes.isEmpty()) {
            return undetermined;
        }

        // see _kb/03-api-layer.md for rationale
        for (int ind = 0; ind < sn.nnodes; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt == NodeType.Source || nt == NodeType.Sink) {
                return undetermined;
            }
        }
        List<Integer> statefulNodes = new ArrayList<Integer>();
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind) > 0) {
                statefulNodes.add(ind);
            }
        }
        for (int pp = 0; pp < placeNodes.size(); pp++) {
            int sfp = statefulNodes.indexOf(placeNodes.get(pp));
            if (sfp < 0) {
                return undetermined;
            }
            for (int sfj = 0; sfj < statefulNodes.size(); sfj++) {
                if (sfj == sfp) {
                    continue;
                }
                if (sn.nodetype.get(statefulNodes.get(sfj)) == NodeType.Transition) {
                    continue;
                }
                for (int r = 0; r < R; r++) {
                    for (int k = 0; k < R; k++) {
                        if (sn.rt.get(sfp * R + r, sfj * R + k) > 0
                                || sn.rt.get(sfj * R + r, sfp * R + k) > 0) {
                            return undetermined;
                        }
                    }
                }
            }
        }

        // Enumerate the (transition, mode) pairs: a mode is what carries a
        // firing rate, and a transition may hold several.
        List<Integer> modeTrans = new ArrayList<Integer>();
        List<Integer> modeIdx = new ArrayList<Integer>();
        List<Boolean> modeTimed = new ArrayList<Boolean>();
        for (int tt = 0; tt < transNodes.size(); tt++) {
            int ind = transNodes.get(tt);
            NodeParam param = sn.nodeparam.get(sn.nodes.get(ind));
            if (!(param instanceof TransitionNodeParam)) {
                return undetermined;
            }
            TransitionNodeParam tparam = (TransitionNodeParam) param;
            for (int m = 0; m < tparam.nmodes; m++) {
                modeTrans.add(ind);
                modeIdx.add(m);
                // see _kb/03-api-layer.md for rationale
                boolean timed = true;
                if (tparam.timing != null && tparam.timing.size() > m) {
                    timed = tparam.timing.get(m) != TimingStrategy.IMMEDIATE;
                }
                modeTimed.add(Boolean.valueOf(timed));
            }
        }
        int nModes = modeTrans.size();
        if (nModes == 0) {
            return undetermined;
        }

        int nPlaces = placeNodes.size();
        double[][][] consumed = new double[nModes][nPlaces][R];
        double[][][] produced = new double[nModes][nPlaces][R];
        for (int mm = 0; mm < nModes; mm++) {
            TransitionNodeParam tparam = (TransitionNodeParam) sn.nodeparam.get(sn.nodes.get(modeTrans.get(mm)));
            Matrix enab = tparam.enabling.get(modeIdx.get(mm));
            Matrix fire = tparam.firing.get(modeIdx.get(mm));
            for (int pp = 0; pp < nPlaces; pp++) {
                int pind = placeNodes.get(pp);
                for (int k = 0; k < R; k++) {
                    consumed[mm][pp][k] = Math.max(0.0, enab.get(pind, k));
                    produced[mm][pp][k] = Math.max(0.0, fire.get(pind, k));
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        int nEq = 2 * nPlaces * R;
        Matrix Afull = new Matrix(nEq, nModes);
        Matrix bfull = new Matrix(nEq, 1);
        int row = 0;
        int nMeasured = 0;
        for (int pp = 0; pp < nPlaces; pp++) {
            int ist = (int) sn.nodeToStation.get(placeNodes.get(pp));
            for (int k = 0; k < R; k++) {
                boolean anyTimed = false;
                for (int mm = 0; mm < nModes; mm++) {
                    double coeff = 0.0;
                    if (modeTimed.get(mm).booleanValue()) {
                        if (tputIsTokens) {
                            coeff = consumed[mm][pp][k];
                        } else {
                            coeff = consumed[mm][pp][k] > 0 ? 1.0 : 0.0;
                        }
                    }
                    if (coeff != 0.0) {
                        anyTimed = true;
                    }
                    Afull.set(row, mm, coeff);
                }
                if (anyTimed) {
                    bfull.set(row, 0, ist >= 0 ? TN.get(ist, k) : 0.0);
                    row++;
                    nMeasured++;
                }

                for (int mm = 0; mm < nModes; mm++) {
                    Afull.set(row, mm, produced[mm][pp][k] - consumed[mm][pp][k]);
                }
                bfull.set(row, 0, 0.0);
                row++;
            }
        }

        // see _kb/03-api-layer.md for rationale
        if (nMeasured == 0) {
            return undetermined;
        }

        Matrix A = Matrix.extract(Afull, 0, row, 0, nModes);
        Matrix b = Matrix.extract(bfull, 0, row, 0, 1);

        Matrix xfit = A.pinv().mult(b);

        // see _kb/03-api-layer.md for rationale
        double xmax = 0.0;
        for (int mm = 0; mm < nModes; mm++) {
            xmax = Math.max(xmax, Math.abs(xfit.get(mm, 0)));
        }
        double negTol = -1e-6 * Math.max(1.0, xmax);
        for (int mm = 0; mm < nModes; mm++) {
            if (xfit.get(mm, 0) < negTol) {
                return undetermined;
            }
        }

        return new Ret(xfit, consumed, produced, placeNodes);
    }
}
