/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.api.sn.SnMapModulation;
import jline.api.sn.SnMapModulation.MapModulation;
import jline.lang.Environment;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.GlobalConstants;
import jline.lang.nodes.ServiceStation;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

import java.util.List;

/**
 * Markov-modulated image of a network with MAP/MMPP/MMAP arrival or service
 * processes as a queueing network in a random environment.
 *
 * <p>Every non-renewal process is a point process modulated by the CTMC with
 * generator Q = D0 + D1, whose conditional intensity in phase k is
 * lambda(k) = sum_j D1(k,j). The transformation freezes each phase into an
 * environment stage in which the process is the Poisson process of that
 * intensity, i.e. an exponential arrival or service time, and lets the
 * environment switch stages at the rates of Q. With P modulated processes the
 * stage set is the Cartesian product of their phase spaces and the environment
 * generator is the Kronecker sum of the individual Q's, so only one process
 * changes phase at a time, as in the original model.
 *
 * <p>The image is exact in structure for an MMPP (diagonal D1): the modulating
 * chain, its stationary distribution and the phase-conditional intensities are
 * all preserved. For a general MAP the phase jumps that occur AT an event epoch
 * (off-diagonal D1) are aggregated into Q and their correlation with the event
 * stream is lost, so the image matches the modulating chain and the conditional
 * intensities but not the full inter-event autocorrelation.
 *
 * <p>Populations are carried across stage switches unchanged (identity reset),
 * as a phase switch moves no job.
 *
 * <p>Mirrors matlab/src/io/map2renv.m.
 */
public class MAPQN2RENV {

    /** Default cap on the number of environment stages. */
    public static final int DEFAULT_MAX_STAGES = 64;

    /** The environment image together with the shape of the transformation. */
    public static final class RenvImage {
        /** The random-environment model. */
        public final Environment env;
        /** Number of stages, i.e. the product of the phase orders. */
        public final int nstages;
        /** Phase order of each modulated process. */
        public final int[] orders;
        /** True when every modulated process was an MMPP (diagonal D1). */
        public final boolean isMMPP;
        /** Longest mean stage sojourn of the environment. */
        public final double maxHoldTime;

        RenvImage(Environment env, int nstages, int[] orders, boolean isMMPP, double maxHoldTime) {
            this.env = env;
            this.nstages = nstages;
            this.orders = orders;
            this.isMMPP = isMMPP;
            this.maxHoldTime = maxHoldTime;
        }
    }

    /**
     * Retained name for the transformation now implemented by
     * {@link #map2renv(Network, SolverOptions)}, which generalizes it from a
     * single MMPP2 service process to any number of MAP, MMPP2 or MMAP arrival
     * and service processes of arbitrary phase order.
     *
     * @param model Network with at least one MAP/MMPP2/MMAP process
     * @return Environment model with exponential rates modulated by the phases
     */
    public static Environment mapqn2renv(Network model) {
        return map2renv(model, null);
    }

    /**
     * @param model   Network with at least one MAP/MMPP2/MMAP process
     * @param options solver options; config.map_env_maxstages caps the stage count
     * @return Environment model with exponential rates modulated by the phases
     */
    public static Environment map2renv(Network model, SolverOptions options) {
        return map2renvImage(model, options).env;
    }

    /**
     * @param model   Network with at least one MAP/MMPP2/MMAP process
     * @param options solver options; config.map_env_maxstages caps the stage count
     * @return the environment image and the shape of the transformation
     */
    public static RenvImage map2renvImage(Network model, SolverOptions options) {
        if (model == null) {
            throw new RuntimeException("map2renv requires a Network model.");
        }
        int maxStages = DEFAULT_MAX_STAGES;
        if (options != null && options.config != null && options.config.map_env_maxstages > 0) {
            maxStages = options.config.map_env_maxstages;
        }

        NetworkStruct sn = model.getStruct(false);
        List<MapModulation> mods = SnMapModulation.snMapModulation(sn);
        if (mods.isEmpty()) {
            throw new RuntimeException(
                    "The model declares no MAP, MMPP2 or MMAP process, so it has no random-environment image.");
        }

        int P = mods.size();
        int[] orders = new int[P];
        int nstages = 1;
        boolean isMMPP = true;
        for (int p = 0; p < P; p++) {
            orders[p] = mods.get(p).order;
            nstages *= orders[p];
            isMMPP = isMMPP && mods.get(p).isMMPP;
        }
        if (nstages > maxStages) {
            StringBuilder ord = new StringBuilder();
            for (int p = 0; p < P; p++) {
                ord.append(p > 0 ? " " : "").append(orders[p]);
            }
            throw new RuntimeException(String.format(
                    "The random-environment image of this model has %d stages (phase orders [%s]), above the "
                            + "options.config.map_env_maxstages cap of %d. Reduce the order of the modulating "
                            + "processes or raise the cap.", nstages, ord.toString(), maxStages));
        }

        // Stage s enumerates the phase tuples in column-major order, phaseOf[s][p]
        // is the phase of process p in stage s.
        int[][] phaseOf = new int[nstages][P];
        for (int s = 0; s < nstages; s++) {
            int rem = s;
            for (int p = 0; p < P; p++) {
                phaseOf[s][p] = rem % orders[p];
                rem = rem / orders[p];
            }
        }

        Environment envModel = new Environment(model.getName() + "_renv", nstages);
        for (int s = 0; s < nstages; s++) {
            String name = stageName(phaseOf[s]);
            envModel.addStage(s, name, "item", buildStage(model, mods, phaseOf[s], name));
        }

        // Kronecker sum of the phase generators: a transition changes the phase of
        // one process only, at the rate that process assigns to it.
        double[] exitRate = new double[nstages];
        for (int s = 0; s < nstages; s++) {
            for (int p = 0; p < P; p++) {
                Matrix Qp = mods.get(p).phaseGenerator();
                int k = phaseOf[s][p];
                int stride = 1;
                for (int q = 0; q < p; q++) {
                    stride *= orders[q];
                }
                for (int l = 0; l < orders[p]; l++) {
                    if (l == k || Qp.get(k, l) <= GlobalConstants.Zero) {
                        continue;
                    }
                    int t = s + (l - k) * stride;
                    envModel.addTransition(s, t, new Exp(Qp.get(k, l)));
                    exitRate[s] += Qp.get(k, l);
                }
            }
        }
        double maxHold = 0;
        for (int s = 0; s < nstages; s++) {
            if (exitRate[s] > 0) {
                maxHold = Math.max(maxHold, 1.0 / exitRate[s]);
            }
        }

        envModel.init();
        return new RenvImage(envModel, nstages, orders, isMMPP, maxHold);
    }

    private static String stageName(int[] phases) {
        StringBuilder sb = new StringBuilder("Phase");
        for (int p = 0; p < phases.length; p++) {
            sb.append('_').append(phases[p] + 1);
        }
        return sb.toString();
    }

    /**
     * Copy of the base model in which every modulated process is the exponential
     * process of its phase-conditional intensity.
     */
    private static Network buildStage(Network model, List<MapModulation> mods, int[] phases, String stageName) {
        Network stageNet = model.copy();
        stageNet.setName(model.getName() + "_" + stageName);
        List<Station> stations = stageNet.getStations();
        List<JobClass> classes = stageNet.getClasses();
        for (int p = 0; p < mods.size(); p++) {
            MapModulation mod = mods.get(p);
            Station station = stations.get(mod.ist);
            for (int c = 0; c < mod.classes.size(); c++) {
                int r = mod.classes.get(c);
                double rate = mod.intensity(c, phases[p]);
                if (mod.arrival) {
                    // A silent phase (zero intensity) is an ON/OFF source: keep it as a
                    // rate rather than a disabled class, so that the class still exists
                    // in every stage and the rate-averaged limit averages a zero.
                    ((Source) station).setArrival(classes.get(r), new Exp(Math.max(rate, GlobalConstants.Zero)));
                } else {
                    if (rate <= GlobalConstants.Zero) {
                        throw new RuntimeException(String.format(
                                "Phase %d of the service process of class %d at station %d has zero completion "
                                        + "rate: the station never empties while the environment sits in that "
                                        + "stage, so the stage has no steady state and the random-environment "
                                        + "image is not defined. Model the stalled server as a breakdown stage "
                                        + "instead.", phases[p] + 1, r + 1, mod.ist + 1));
                    }
                    ((ServiceStation) station).setService(classes.get(r), new Exp(rate));
                }
            }
        }
        // The copy inherited the base model's cached NetworkStruct, so the edits
        // above are invisible until the struct is rebuilt.
        stageNet.refreshStruct(true);
        return stageNet;
    }
}
