package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;

/**
 * Decides whether a model lives on a discrete (slotted) time scale.
 *
 * <p>A model is discrete-time when every enabled interarrival and service law
 * is lattice-valued on a common slot length and at least one of them is
 * intrinsically discrete. The lattice families are Geometric (support
 * {1,2,...} slots), DMAP, DiscreteUniform with integral bounds, and Det whose
 * value is a positive integral number of slots. Immediate is deliberately not
 * one: a zero interval is not a point of {d,2d,...}, the same refusal the LDES
 * slotted engine makes in slotSnap.
 *
 * <p>The test runs on procid, rates and scv rather than on proc, because the
 * struct refresh may already have replaced a lattice law by a continuous
 * surrogate. procid keeps the requested family and (mean, SCV) identify the
 * member of it exactly for every family above.
 *
 * <p>MATLAB twin: sn_is_discrete_time.m
 */
public final class SnIsDiscreteTime {
    private SnIsDiscreteTime() {}

    /** Outcome of the test: the verdict, the slot length and why it failed. */
    public static final class Result {
        public final boolean isDiscreteTime;
        public final double slotLength;
        public final boolean hasLattice;
        public final boolean hasContinuous;
        public final boolean hasDmap;
        public final String reason;

        public Result(boolean isDiscreteTime, double slotLength, boolean hasLattice,
                      boolean hasContinuous, boolean hasDmap, String reason) {
            this.isDiscreteTime = isDiscreteTime;
            this.slotLength = slotLength;
            this.hasLattice = hasLattice;
            this.hasContinuous = hasContinuous;
            this.hasDmap = hasDmap;
            this.reason = reason;
        }
    }

    public static Result snIsDiscreteTime(NetworkStruct sn, SolverOptions options) {
        double tol = 1e-8;
        String timescale = "auto";
        double slotLength = 1.0;

        if (options != null && options.config != null) {
            Object ts = options.config.get("timescale");
            if (ts instanceof String) {
                timescale = ((String) ts).toLowerCase();
            }
            Object sl = options.config.get("slotlength");
            if (sl instanceof Number) {
                slotLength = ((Number) sl).doubleValue();
            }
        }
        if (!"auto".equals(timescale) && !"discrete".equals(timescale) && !"continuous".equals(timescale)) {
            throw new RuntimeException("config.timescale must be 'auto', 'discrete' or 'continuous'.");
        }
        if (slotLength <= 0 || Double.isInfinite(slotLength) || Double.isNaN(slotLength)) {
            throw new RuntimeException("config.slotlength must be a positive finite scalar.");
        }
        if ("continuous".equals(timescale)) {
            return new Result(false, slotLength, false, false, false, "");
        }

        boolean hasLattice = false;
        boolean hasContinuous = false;
        boolean hasDmap = false;
        String reason = "";

        for (int ist = 0; ist < sn.nstations; ist++) {
            Station station = sn.stations.get(ist);
            for (int r = 0; r < sn.nclasses; r++) {
                ProcessType procType = null;
                if (sn.procid.containsKey(station)) {
                    procType = sn.procid.get(station).get(sn.jobclasses.get(r));
                }
                if (procType == null || procType == ProcessType.DISABLED) {
                    continue;
                }
                double rate = sn.rates.get(ist, r);
                if (Double.isNaN(rate) || rate <= 0) {
                    continue;
                }
                double meanSlots = 1.0 / (rate * slotLength);

                if (procType == ProcessType.GEOMETRIC) {
                    hasLattice = true;
                    if (meanSlots < 1 - tol) {
                        hasContinuous = true;
                        reason = "Geometric at station " + ist + " class " + r + " has mean "
                                + meanSlots + " slots, below the one-slot minimum of its support.";
                    }
                } else if (procType == ProcessType.DMAP) {
                    hasLattice = true;
                    hasDmap = true;
                } else if (procType == ProcessType.DUNIFORM) {
                    hasLattice = true;
                    double varSlots = sn.scv.get(ist, r) * meanSlots * meanSlots;
                    double width = Math.sqrt(Math.max(0.0, 12 * varSlots + 1)) - 1;
                    long lo = Math.round(meanSlots - width / 2);
                    if (lo < 1) {
                        hasContinuous = true;
                        reason = "DiscreteUniform at station " + ist + " class " + r
                                + " is not contained in {1,2,...}.";
                    }
                } else if (procType == ProcessType.DET) {
                    if (Math.abs(meanSlots - Math.round(meanSlots)) <= tol * Math.max(1.0, meanSlots)
                            && Math.round(meanSlots) >= 1) {
                        hasLattice = true;
                    } else {
                        // a Det off the lattice is what makes the model continuous
                        hasContinuous = true;
                    }
                } else {
                    hasContinuous = true;
                }
            }
        }

        if ("discrete".equals(timescale)) {
            if (hasLattice && hasContinuous) {
                throw new RuntimeException("config.timescale='discrete' was requested but the model "
                        + "mixes lattice and non-lattice laws. " + reason);
            }
            if (!hasLattice) {
                throw new RuntimeException("config.timescale='discrete' was requested but no "
                        + "interarrival or service law is lattice-valued on a slot of " + slotLength + ".");
            }
            return new Result(true, slotLength, hasLattice, hasContinuous, hasDmap, "");
        }

        boolean isDT = hasLattice && !hasContinuous;

        if (!isDT && hasDmap) {
            // A DMAP has no continuous-time reading: its (D0,D1) are probability
            // matrices, so the CTMC machinery would compute inv(-D0) where the
            // law needs inv(I-D0) and return a wrong number in silence.
            throw new RuntimeException("The model mixes a DMAP with continuous-time laws. A DMAP is "
                    + "only defined on a slotted time scale, so no solver can interpret this model. "
                    + reason);
        }

        return new Result(isDT, slotLength, hasLattice, hasContinuous, hasDmap, reason);
    }
}
