/**
 * @file Whether the fluid Petri route can answer a model, and why not.
 *
 * Port of {@code matlab/src/solvers/FLD/fluid_petri_applicable.m}.
 *
 * @since LINE 3.0
 */
package jline.solvers.fluid.petri;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;

/**
 * Whether the fluid Petri route can answer this model, and why not.
 *
 * <p>A QUEUEING STATION IS THE ONE STRUCTURAL EXCLUSION: a net whose tokens also visit a Queue or a
 * Delay is two formalisms at once, and LINE has no reference semantics for the hand-off. Saying so
 * by name is the point of this class -- the drift would otherwise be built over the Petri
 * coordinates alone and integrate the queueing part of the model as though it were not there.
 */
public final class PetriApplicable {
    private PetriApplicable() {}

    /** The verdict and, when negative, the element that decides it. */
    public static final class Verdict {
        public final boolean ok;
        public final String reason;

        Verdict(boolean ok, String reason) {
            this.ok = ok;
            this.reason = reason;
        }
    }

    /**
     * @param sn      the model
     * @param options the solver options, whose {@code config.hide_immediate} is consulted
     * @return whether the fluid Petri drift can be built, and the refusal otherwise
     */
    public static Verdict applicable(NetworkStruct sn, SolverOptions options) {
        final int I = sn.nodetype.size();
        boolean hasTransition = false;
        for (int i = 0; i < I; i++) {
            if (sn.nodetype.get(i) == NodeType.Transition) {
                hasTransition = true;
                break;
            }
        }
        if (!hasTransition) {
            return new Verdict(false,
                    "the model has no Transition node, so it is not a Petri net");
        }

        for (int ind = 0; ind < I; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt != NodeType.Place && nt != NodeType.Transition && nt != NodeType.Source
                    && nt != NodeType.Sink) {
                return new Verdict(false, "node " + sn.nodenames.get(ind) + " is a " + nt
                        + ". The fluid Petri route solves the marking of a Petri net, and a model "
                        + "that also holds queueing stations is two formalisms at once with no "
                        + "reference semantics for the hand-off; use SolverCTMC, SolverJMT, "
                        + "SolverSSA or SolverLDES");
            }
        }

        // A queueing place declares a service process, which is what turns it
        // into a station with an embedded queue and a depository.
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Place) {
                continue;
            }
            int ist = (int) sn.nodeToStation.get(ind);
            if (ist < 0 || sn.rates == null || ist >= sn.rates.getNumRows()) {
                continue;
            }
            for (int k = 0; k < sn.nclasses; k++) {
                double v = sn.rates.get(ist, k);
                if (!Double.isNaN(v) && v > 0) {
                    return new Verdict(false, "place " + sn.nodenames.get(ind)
                            + " is a QUEUEING place (it declares a service process), whose "
                            + "embedded queue this drift does not carry; use SolverLDES");
                }
            }
        }

        if (options != null && options.config != null && options.config.hide_immediate) {
            return new Verdict(false,
                    "options.config.hide_immediate eliminates the immediate transitions from the "
                    + "event set, but the fluid Petri route needs them: it solves their firing "
                    + "flows as algebraic unknowns");
        }
        return new Verdict(true, "");
    }
}
