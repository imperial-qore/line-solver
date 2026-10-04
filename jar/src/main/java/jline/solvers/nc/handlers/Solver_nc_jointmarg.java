/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.Pfqn_jointmarg;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.util.Utils;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Joint probability of the per-station TOTAL queue lengths.
 *
 * <p>This is NOT {@link Solver_nc_jointaggr}, which fixes the per-class
 * population of every station: each state here is the sum of jointaggr over the
 * whole fibre of per-class tables with these row sums, and that fibre grows
 * combinatorially. The permanent evaluates the sum in closed form
 * ({@link Pfqn_jointmarg}).</p>
 */
public final class Solver_nc_jointmarg {
    private Solver_nc_jointmarg() {}

    /** Result of a joint total-queue-length evaluation. */
    public static class Ret_jointmarg {
        /** Joint probability. */
        public final double Pr;
        /** Its logarithm. */
        public final double lPr;
        /** Log normalizing constant, whether supplied or computed here. */
        public final double lG;
        /** Wall-clock runtime in seconds. */
        public final double runtime;

        public Ret_jointmarg(double Pr, double lPr, double lG, double runtime) {
            this.Pr = Pr;
            this.lPr = lPr;
            this.lG = lG;
            this.runtime = runtime;
        }
    }

    /**
     * @param sn      network structure
     * @param options solver options
     * @param nvec    per-station total job counts
     * @param engine  permanent engine, "exact" by default
     * @param lG      log normalizing constant already known, or null
     * @return the probability, its logarithm, the constant and the runtime
     */
    public static Ret_jointmarg solver_nc_jointmarg(NetworkStruct sn, SolverOptions options,
                                                    Matrix nvec, String engine, Double lG) {
        long startTimeMillis = System.nanoTime();

        int M = sn.nstations;
        if (nvec.getNumRows() * nvec.getNumCols() != M) {
            line_error("Solver_nc_jointmarg", "The occupancy vector has "
                    + (nvec.getNumRows() * nvec.getNumCols()) + " entries but the model has " + M + " stations.");
        }

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain.copy();
        for (int i = 0; i < Lchain.getNumRows(); i++) {
            for (int j = 0; j < Lchain.getNumCols(); j++) {
                if (!Double.isFinite(Lchain.get(i, j))) {
                    Lchain.set(i, j, 0.0);
                }
            }
        }
        Matrix Nchain = dem.Nchain;

        checkSupported(sn, Nchain);

        List<Integer> inf = new ArrayList<Integer>();
        for (int ist = 0; ist < M; ist++) {
            if (Utils.isInf(sn.nservers.get(ist))) {
                inf.add(Integer.valueOf(ist));
            }
        }
        int[] infset = new int[inf.size()];
        for (int k = 0; k < inf.size(); k++) {
            infset[k] = inf.get(k).intValue();
        }

        Pfqn_jointmarg.Ret_jointmarg ret =
                Pfqn_jointmarg.pfqn_jointmarg(nvec, Lchain, Nchain, infset, lG, engine);

        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;
        return new Ret_jointmarg(ret.pjoint, ret.lpjoint, ret.lG, runtime);
    }

    /**
     * The permanent identity supplies one n_i! per queueing station and none per
     * infinite server. A multiserver or load-dependent station has neither, so
     * it is refused by name rather than approximated.
     */
    private static void checkSupported(NetworkStruct sn, Matrix Nchain) {
        for (int c = 0; c < Nchain.getNumRows() * Nchain.getNumCols(); c++) {
            if (!Double.isFinite(Nchain.get(c))) {
                line_error("Solver_nc_jointmarg", "getProbSysMarg requires a closed model: the joint law of the "
                        + "total queue lengths is not defined when a class has an infinite population.");
            }
        }
        for (int k = 0; k < sn.njobs.getNumRows() * sn.njobs.getNumCols(); k++) {
            if (Utils.isInf(sn.njobs.get(k))) {
                line_error("Solver_nc_jointmarg", "getProbSysMarg requires a closed model: the joint law of the "
                        + "total queue lengths is not defined when a class has an infinite population.");
            }
        }
        if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
            line_error("Solver_nc_jointmarg", "getProbSysMarg does not support load-dependent stations "
                    + "(sn.lldscaling is set): the permanent identity supplies exactly one n_i! per queueing station.");
        }
        if (sn.cdscaling != null && sn.cdscaling.size() > 0) {
            line_error("Solver_nc_jointmarg", "getProbSysMarg does not support class-dependent scaling "
                    + "(sn.cdscaling is set).");
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            double c = sn.nservers.get(ist);
            if (!Utils.isInf(c) && c > 1) {
                line_error("Solver_nc_jointmarg", "getProbSysMarg does not support the multiserver station "
                        + (ist + 1) + " (" + (int) c + " servers): the permanent identity supplies exactly one n_i! "
                        + "per queueing station.");
            }
        }
    }
}
