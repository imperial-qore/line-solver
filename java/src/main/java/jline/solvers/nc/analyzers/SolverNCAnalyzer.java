package jline.solvers.nc.analyzers;

import jline.lang.NetworkStruct;
import jline.lang.constant.GlobalConstants;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.solvers.nc.SolverNCResult;
import jline.util.Matrix;

public class SolverNCAnalyzer implements NCAnalyzer {
  @Override
  public void analyze(NetworkStruct sn, SolverOptions options, SolverNCResult res) {

    long Tstart = System.currentTimeMillis();
    Matrix nservers = sn.nservers;
    Matrix nserversFinite = nservers.clone();
    nserversFinite.removeInfinity();
    if (nserversFinite.elementMax() > 1 && options.method.equals("exact")) {
      throw new RuntimeException("NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the 'exact' option.");
    }
    NetworkStruct snfloor;
    NetworkStruct snceil;
    try {
       snfloor = (NetworkStruct) sn.clone();
       snceil = (NetworkStruct) sn.clone();
    } catch (CloneNotSupportedException e) {
      throw new RuntimeException("Cloning failed.");
    }

    snfloor.njobs = sn.njobs.clone();
    snceil.njobs = sn.njobs.clone();

    Matrix eta = new Matrix(sn.njobs.getNumRows(), sn.njobs.getNumCols());
    boolean nonIntegerJob = false;
    for (int i=0; i<eta.getNumRows(); i++) {
      for (int j=0; j<eta.getNumCols(); j++) {
        snfloor.njobs.set(i, j, Math.floor(sn.njobs.get(i, j)));
        snceil.njobs.set(i, j, Math.ceil(sn.njobs.get(i, j)));
        eta.set(i, j, Math.abs(sn.njobs.get(i, j) - snfloor.njobs.get(i, j)));
        if (eta.get(i, j) > GlobalConstants.FineTol) {
          nonIntegerJob = true;
        }
      }
    }
    if (nonIntegerJob) {
      SolverNC.SolverNCReturn retfloor = SolverNC.solver_nc(snfloor, options);
      res.runtime = retfloor.runtime;
      SolverNC.SolverNCReturn retceil = SolverNC.solver_nc(snceil, options);
      res.runtime += retceil.runtime;
      res.QN = retfloor.Q.add(1, eta.elementMult(retceil.Q.sub(1, retfloor.Q), null));
      res.UN = retfloor.U.add(1, eta.elementMult(retceil.U.sub(1, retfloor.U), null));
      res.RN = retfloor.R.add(1, eta.elementMult(retceil.R.sub(1, retfloor.R), null));
      res.TN = retfloor.T.add(1, eta.elementMult(retceil.T.sub(1, retfloor.T), null));
      res.CN = new Matrix(retfloor.C); // TODO: Decide how to handle C in this case.
      res.XN = retfloor.X.add(1, eta.elementMult(retceil.X.sub(1, retfloor.X), null));
      // TODO: ret.lG for more than one non-integer job quantity
      res.lG = retfloor.lG + eta.elementMult(new Matrix(retfloor.lG - retceil.lG), null).get(0, 0);
      res.it = retfloor.it + retceil.it;
    }
    else {
      SolverNC.SolverNCReturn ret = SolverNC.solver_nc(sn, options);
      res.QN = ret.Q;
      res.UN = ret.U;
      res.RN = ret.R;
      res.TN = ret.T;
      res.CN = new Matrix(ret.C);
      res.XN = ret.X;
      res.lG = ret.lG;
      res.it = ret.it;
      res.runtime = ret.runtime;
      res.method = ret.method;
    }
    res.runtime = System.currentTimeMillis() - Tstart;
  }
}
