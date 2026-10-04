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
import jline.solvers.nc.handlers.Solver_ncld;
import jline.util.matrix.Matrix;

public final class Solver_ncld_analyzer {
    private Solver_ncld_analyzer() {}

    public static NCResult solver_ncld_analyzer(NetworkStruct sn, SolverOptions options) {
        long Tstart = System.nanoTime();
        Matrix nservers = sn.nservers;
        Matrix nserversFinite = nservers.copy();
        nserversFinite.removeInfinity();
        if (nserversFinite.elementMax() > 1 && sn.njobs.hasInfinite() && "exact".equals(options.method)) {
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
            if ("exact".equals(options.method)) {
                throw new RuntimeException("NC load-dependent solver cannot provide exact solutions for fractional populations.");
            }
            SolverNC.SolverNCLDReturn retfloor = Solver_ncld.solver_ncld(snfloor, options);
            SolverNC.SolverNCLDReturn retceil = Solver_ncld.solver_ncld(snceil, options);
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
            SolverNC.SolverNCLDReturn ret = Solver_ncld.solver_ncld(sn, options);
            res.QN = ret.Q;
            res.UN = ret.U;
            res.RN = ret.R;
            res.TN = ret.T;
            res.XN = ret.X;

            res.CN = new Matrix(1, sn.nclasses);
            for (int k = 0; k < sn.nclasses; k++) {
                res.CN.set(0, k, sn.njobs.get(0, k) / res.XN.get(0, k));
            }

            res.lG = ret.lG;
            res.it = ret.it;
            res.method = ret.method;
        }
        res.runtime = (double) (System.nanoTime() - Tstart) / 1000000000.0;
        return res;
    }
}
