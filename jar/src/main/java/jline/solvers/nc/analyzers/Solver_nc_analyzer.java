/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.solvers.nc.SolverNC;
import jline.solvers.nc.handlers.Solver_nc;
import jline.util.matrix.Matrix;

public final class Solver_nc_analyzer {
    private Solver_nc_analyzer() {}

    public static NCResult solver_nc_analyzer(NetworkStruct sn, SolverOptions options) {
        long Tstart = System.nanoTime();
        Matrix nservers = sn.nservers;
        Matrix nserversFinite = nservers.copy();
        nserversFinite.removeInfinity();
        // The refusal is for an OPEN or MIXED multiserver model, as the message says:
        // a CLOSED one has an exact answer and MATLAB's solver_nc_analyzer.m:78 and
        // C++'s nc_dispatch.h:161 both gate on it. This copy omitted
        // `sn.njobs.hasInfinite()`, which its own ncld twin (Solver_ncld_analyzer:25)
        // carries, and so refused closed multiserver models too. Latent until
        // config.multiserver='seidmann' gave 'exact' a route to this analyzer with a
        // multiserver station still un-converted. See _kb/06-solver-catalog.md (NC)
        if (nserversFinite.elementMax() > 1 && sn.njobs.hasInfinite()
                && "exact".equals(options.method)) {
            throw new RuntimeException("NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the 'exact' option.");
        }
        NetworkStruct snfloor = sn.copy();
        NetworkStruct snceil = sn.copy();

        snfloor.njobs = sn.njobs.copy();
        snceil.njobs = sn.njobs.copy();

        Matrix eta = new Matrix(sn.njobs.getNumRows(), sn.njobs.getNumCols());
        boolean nonIntegerJob = false;
        for (int i = 0; i < eta.getNumRows(); i++) {
            for (int j = 0; j < eta.getNumCols(); j++) {
                snfloor.njobs.set(i, j, FastMath.floor(sn.njobs.get(i, j)));
                snceil.njobs.set(i, j, FastMath.ceil(sn.njobs.get(i, j)));
                eta.set(i, j, FastMath.abs(sn.njobs.get(i, j) - snfloor.njobs.get(i, j)));
                if (eta.get(i, j) > GlobalConstants.FineTol) {
                    nonIntegerJob = true;
                }
            }
        }

        NCResult res = new NCResult();
        if (nonIntegerJob) {
            SolverNC.SolverNCReturn retfloor = Solver_nc.solver_nc(snfloor, options);
            res.runtime = retfloor.runtime;
            SolverNC.SolverNCReturn retceil = Solver_nc.solver_nc(snceil, options);
            res.runtime += retceil.runtime;
            res.QN = retfloor.Q.add(1.0, eta.elementMult(retceil.Q.sub(retfloor.Q), null));
            res.UN = retfloor.U.add(1.0, eta.elementMult(retceil.U.sub(retfloor.U), null));
            res.RN = retfloor.R.add(1.0, eta.elementMult(retceil.R.sub(retfloor.R), null));
            res.TN = retfloor.T.add(1.0, eta.elementMult(retceil.T.sub(retfloor.T), null));
            res.XN = retfloor.X.add(1.0, eta.elementMult(retceil.X.sub(retfloor.X), null));

            res.CN = new Matrix(1, sn.nclasses);
            for (int k = 0; k < sn.nclasses; k++) {
                res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k));
            }

            res.lG = retfloor.lG + eta.elementSum() * (retfloor.lG - retceil.lG);
            res.it = retfloor.it + retceil.it;
            res.method = retceil.method;
        } else {
            SolverNC.SolverNCReturn ret = Solver_nc.solver_nc(sn, options);
            res.QN = ret.Q.copy();
            res.UN = ret.U.copy();
            res.RN = ret.R.copy();
            res.TN = ret.T.copy();
            res.XN = ret.X.copy();

            res.CN = new Matrix(1, sn.nclasses);
            for (int k = 0; k < sn.nclasses; k++) {
                res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k));
            }

            res.lG = ret.lG;
            res.it = ret.it;
            res.runtime = ret.runtime;
            res.method = ret.method;
        }
        res.runtime = (double) (System.nanoTime() - Tstart) / 1000000000.0;

        if (!Double.isNaN(res.lG) && !Double.isInfinite(res.lG)) {
            jline.io.LineConsole.step("normalizing constant obtained: log G = %.6g", res.lG);
        }
        return res;
    }
}
