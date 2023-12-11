package jline.solvers.fluid;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Erlang;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static java.lang.Double.NaN;
import static jline.solvers.fluid.FluidTestModels.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverFluidOtherTest {

  static double tol = 0.005; // Ideally should be 0.0001 but some results don't quite match MatLab

  @Test
  public void ex1ReturnsCorrectResultFromRunAnalyzer() {

    Network model = ex1();

    SolverOptions options = new SolverOptions(SolverType.Fluid);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.fluidResult;
    SolverResult result = solver.result;

    // method
    assertEquals("closing", result.method);

    // QN
    assertEquals(5, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(5, result.QN.getNumElements());
    assertEquals(0.4417, result.QN.get(0, 0), tol);
    assertEquals(0.15, result.QN.get(1, 0), tol);
    assertEquals(0.2, result.QN.get(2, 0), tol);
    assertEquals(0.075, result.QN.get(3, 0), tol);
    assertEquals(0.133, result.QN.get(4, 0), tol);

    // RN
    assertEquals(5, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(5, result.RN.getNumElements());
    assertEquals(0.2208, result.RN.get(0, 0), tol);
    assertEquals(0.25, result.RN.get(1, 0), tol);
    assertEquals(0.5, result.RN.get(2, 0), tol);
    assertEquals(0.125, result.RN.get(3, 0), tol);
    assertEquals(0.3333, result.RN.get(4, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(2, result.XN.get(0, 0), tol);

    // UN
    assertEquals(5, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(5, result.UN.getNumElements());
    assertEquals(0.4417, result.UN.get(0, 0), tol);
    assertEquals(0.15, result.UN.get(1, 0), tol);
    assertEquals(0.2, result.UN.get(2, 0), tol);
    assertEquals(0.075, result.UN.get(3, 0), tol);
    assertEquals(0.1333, result.UN.get(4, 0), tol);

    // TN
    assertEquals(5, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(5, result.TN.getNumElements());
    assertEquals(2, result.TN.get(0, 0), tol);
    assertEquals(0.6, result.TN.get(1, 0), tol);
    assertEquals(0.4, result.TN.get(2, 0), tol);
    assertEquals(0.5998, result.TN.get(3, 0), tol);
    assertEquals(0.4, result.TN.get(4, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(Double.POSITIVE_INFINITY, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(5, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(1, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    assertEquals(0, result.QNt[4][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.4417, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.15, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.2, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.075, result.QNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1333, result.QNt[4][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(5, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(1, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(0, result.UNt[4][0].get(0, 0), tol);
    assertEquals(0.4417, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.15, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.2, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.075, result.UNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1333, result.UNt[4][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(5, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(2, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(0, result.TNt[4][0].get(0, 0), tol);
    assertEquals(0.8834, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.6, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.4, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5998, result.TNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(0.4, result.TNt[4][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(1000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(0.4417, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.15, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0.1, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.1, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(0.075, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(0.0667, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(0.0667, fluidResult.odeStateVec.get(0, 6), tol);
  }

  @Test
  public void ex2ReturnsCorrectResultFromRunAnalyzer() {

    Network model = ex2();

    SolverOptions options = new SolverOptions(SolverType.Fluid);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.fluidResult;
    SolverResult result = solver.result;

    // method
    assertEquals("closing", result.method);

    // QN
    assertEquals(4, result.QN.getNumRows());
    assertEquals(3, result.QN.getNumCols());
    assertEquals(12, result.QN.getNumElements());
    assertEquals(0.6456, result.QN.get(0, 0), tol);
    assertEquals(0.6456, result.QN.get(0, 1), tol);
    assertEquals(0.8228, result.QN.get(0, 2), tol);
    assertEquals(0.0625, result.QN.get(1, 0), tol);
    assertEquals(0.0625, result.QN.get(1, 1), tol);
    assertEquals(0.0703, result.QN.get(1, 2), tol);
    assertEquals(0.0625, result.QN.get(2, 0), tol);
    assertEquals(0.0625, result.QN.get(2, 1), tol);
    assertEquals(0.1407, result.QN.get(2, 2), tol);
    assertEquals(0.1, result.QN.get(3, 0), tol);
    assertEquals(0.1, result.QN.get(3, 1), tol);
    assertEquals(0.2250, result.QN.get(3, 2), tol);

    // RN
    assertEquals(4, result.RN.getNumRows());
    assertEquals(3, result.RN.getNumCols());
    assertEquals(12, result.RN.getNumElements());
    assertEquals(0.3228, result.RN.get(0, 0), tol);
    assertEquals(0.3228, result.RN.get(0, 1), tol);
    assertEquals(0.8228, result.RN.get(0, 2), tol);
    assertEquals(0.0313, result.RN.get(1, 0), tol);
    assertEquals(0.0313, result.RN.get(1, 1), tol);
    assertEquals(0.0313, result.RN.get(1, 2), tol);
    assertEquals(0.0313, result.RN.get(2, 0), tol);
    assertEquals(0.0313, result.RN.get(2, 1), tol);
    assertEquals(0.0625, result.RN.get(2, 2), tol);
    assertEquals(0.0625, result.RN.get(3, 0), tol);
    assertEquals(0.0625, result.RN.get(3, 1), tol);
    assertEquals(0.125, result.RN.get(3, 2), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(3, result.XN.getNumCols());
    assertEquals(3, result.XN.getNumElements());
    assertEquals(2, result.XN.get(0, 0), tol);
    assertEquals(2, result.XN.get(0, 1), tol);
    assertEquals(1, result.XN.get(0, 2), tol);

    // UN
    assertEquals(4, result.UN.getNumRows());
    assertEquals(3, result.UN.getNumCols());
    assertEquals(12, result.UN.getNumElements());
    assertEquals(0.6456, result.UN.get(0, 0), tol);
    assertEquals(0.6456, result.UN.get(0, 1), tol);
    assertEquals(0.8228, result.UN.get(0, 2), tol);
    assertEquals(0.0625, result.UN.get(1, 0), tol);
    assertEquals(0.0625, result.UN.get(1, 1), tol);
    assertEquals(0.0703, result.UN.get(1, 2), tol);
    assertEquals(0.0625, result.UN.get(2, 0), tol);
    assertEquals(0.0625, result.UN.get(2, 1), tol);
    assertEquals(0.1407, result.UN.get(2, 2), tol);
    assertEquals(0.1, result.UN.get(3, 0), tol);
    assertEquals(0.1, result.UN.get(3, 1), tol);
    assertEquals(0.225, result.UN.get(3, 2), tol);

    // TN
    assertEquals(4, result.TN.getNumRows());
    assertEquals(3, result.TN.getNumCols());
    assertEquals(12, result.TN.getNumElements());
    assertEquals(2, result.TN.get(0, 0), tol);
    assertEquals(2, result.TN.get(0, 1), tol);
    assertEquals(1, result.TN.get(0, 2), tol);
    assertEquals(2, result.TN.get(1, 0), tol);
    assertEquals(2, result.TN.get(1, 1), tol);
    assertEquals(2.2487, result.TN.get(1, 2), tol);
    assertEquals(2, result.TN.get(2, 0), tol);
    assertEquals(2, result.TN.get(2, 1), tol);
    assertEquals(2.251, result.TN.get(2, 2), tol);
    assertEquals(1.6, result.TN.get(3, 0), tol);
    assertEquals(1.6, result.TN.get(3, 1), tol);
    assertEquals(1.7998, result.TN.get(3, 2), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(3, result.CN.getNumCols());
    assertEquals(3, result.CN.getNumElements());
    assertEquals(Double.POSITIVE_INFINITY, result.CN.get(0, 0), tol);
    assertEquals(Double.POSITIVE_INFINITY, result.CN.get(0, 1), tol);
    assertEquals(Double.POSITIVE_INFINITY, result.CN.get(0, 2), tol);

    // QNt
    assertEquals(4, result.QNt.length);
    assertEquals(3, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(1, result.QNt[0][0].get(0, 0), tol);
    assertEquals(1, result.QNt[0][1].get(0, 0), tol);
    assertEquals(1, result.QNt[0][2].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][2].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][1].get(0, 0), tol);
    assertEquals(0, result.QNt[2][2].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][1].get(0, 0), tol);
    assertEquals(0, result.QNt[3][2].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.6456, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.6456, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.8228, result.QNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.0703, result.QNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[2][1].get(Tmax - 1, 0), tol);
    assertEquals(0.1407, result.QNt[2][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1, result.QNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1, result.QNt[3][1].get(Tmax - 1, 0), tol);
    assertEquals(0.225, result.QNt[3][2].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(4, result.UNt.length);
    assertEquals(3, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(0.3333, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0.3333, result.UNt[0][1].get(0, 0), tol);
    assertEquals(0.3333, result.UNt[0][2].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][2].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][1].get(0, 0), tol);
    assertEquals(0, result.UNt[2][2].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][1].get(0, 0), tol);
    assertEquals(0, result.UNt[3][2].get(0, 0), tol);
    assertEquals(0.3054, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.3054, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.3892, result.UNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.0703, result.UNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[2][1].get(Tmax - 1, 0), tol);
    assertEquals(0.1407, result.UNt[2][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1, result.UNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1, result.UNt[3][1].get(Tmax - 1, 0), tol);
    assertEquals(0.225, result.UNt[3][2].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(4, result.TNt.length);
    assertEquals(3, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(0.6667, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0.6667, result.TNt[0][1].get(0, 0), tol);
    assertEquals(0.3333, result.TNt[0][2].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(0, result.TNt[1][2].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][1].get(0, 0), tol);
    assertEquals(0, result.TNt[2][2].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][1].get(0, 0), tol);
    assertEquals(0, result.TNt[3][2].get(0, 0), tol);
    assertEquals(0.6108, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.6108, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.3892, result.TNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(2, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.TNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(2.2487, result.TNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(2, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.TNt[2][1].get(Tmax - 1, 0), tol);
    assertEquals(2.251, result.TNt[2][2].get(Tmax - 1, 0), tol);
    assertEquals(1.6, result.TNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(1.6, result.TNt[3][1].get(Tmax - 1, 0), tol);
    assertEquals(1.7998, result.TNt[3][2].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(2000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(0.6456, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.6456, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0.8228, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(0.0703, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 6), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 7), tol);
    assertEquals(0.1407, fluidResult.odeStateVec.get(0, 8), tol);
    assertEquals(0.1, fluidResult.odeStateVec.get(0, 9), tol);
    assertEquals(0.1, fluidResult.odeStateVec.get(0, 10), tol);
    assertEquals(0.225, fluidResult.odeStateVec.get(0, 11), tol);
  }

  // Examples below this point were used as integration tests during the building of SolverFluid
  // They map to certain examples in "gettingstarted" or "examples" in LINE, some with tweaks i.e.
  // use of specific methods, stiff v. non-stiff, etc. Not used for performance evaluation

  @Test
  public void gettingStartedExample7ReturnsCorrectResultFromRunAnalyzer() {

    // Corresponds to "getting_started_ex7.m" in LINE
    String modelName = "getting_started_example_7";
    Network model = new Network(modelName);

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

    ClosedClass closedClass = new ClosedClass(model, "Class1", 5, delay, 0);
    delay.setService(closedClass, new Exp(1.0));
    queue.setService(closedClass, new Exp(0.5));

    model.link(model.serialRouting(delay, queue));

    SolverOptions options = new SolverOptions(SolverType.Fluid);
    options.iter_max = 200;
    SolverFluid solverFluid = new SolverFluid(model, options);

    solverFluid.options.stiff = false;
    solverFluid.runAnalyzer();
    SolverFluidResult fluidResult = solverFluid.fluidResult;
    SolverResult result = solverFluid.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(2, result.QN.getNumElements());
    assertEquals(0.5, result.QN.get(0, 0), tol);
    assertEquals(4.5, result.QN.get(1, 0), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(2, result.RN.getNumElements());
    assertEquals(1, result.RN.get(0, 0), tol);
    assertEquals(9, result.RN.get(1, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(0.5, result.XN.get(0, 0), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(2, result.UN.getNumElements());
    assertEquals(0.5, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);

    // TN
    assertEquals(2, result.TN.getNumRows(), 2);
    assertEquals(1, result.TN.getNumCols(), 1);
    assertEquals(2, result.TN.getNumElements(), 2);
    assertEquals(0.5, result.TN.get(0, 0), tol);
    assertEquals(0.5, result.TN.get(1, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(10, result.CN.get(0, 0), tol);

    // TODO: Do QNt et. al. need changing to add 't' matrix as a second column?
    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(5, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.5, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(4.5, result.QNt[1][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(5, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0.5, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(5, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0.5, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5, result.TNt[1][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(4000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(0.5, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(4.5, fluidResult.odeStateVec.get(0, 1), tol);
  }

  @Test
  public void exampleCdfRespT2ReturnsCorrectResultFromRunAnalyzerDefaultMethod() {

    // Corresponds to "example_cdfRespT_2.m" in LINE
    String modelName = "example_cdfRespT_2";
    Network model = new Network(modelName);

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);

    ClosedClass closedClass1 = new ClosedClass(model, "Class1", 1, delay, 0);
    ClosedClass closedClass2 = new ClosedClass(model, "Class2", 0, delay, 0);
    ClosedClass closedClass3 = new ClosedClass(model, "Class3", 0, delay, 0);

    delay.setService(closedClass1, new Exp(1.0));
    delay.setService(closedClass2, new Exp(1.0));
    delay.setService(closedClass3, new Exp(1.0));
    queue.setService(closedClass1, new Exp(1.0));
    queue.setService(closedClass2, new Erlang(0.5, 2));
    queue.setService(closedClass3, new Exp(1 / 0.01));

    RoutingMatrix routingMatrix =
            new RoutingMatrix(
                    model,
                    Arrays.asList(closedClass1, closedClass2, closedClass3),
                    Arrays.asList(delay, queue));
    routingMatrix.addConnection(closedClass1, closedClass1, delay, queue, 1.0);
    routingMatrix.addConnection(closedClass2, closedClass2, delay, queue, 1.0);
    routingMatrix.addConnection(closedClass1, closedClass2, queue, delay, 1.0);
    routingMatrix.addConnection(closedClass2, closedClass1, queue, delay, 1.0);
    routingMatrix.addConnection(delay, queue, closedClass3, 1.0);
    routingMatrix.addConnection(queue, delay, closedClass3, 1.0);
    model.link(routingMatrix);

    SolverOptions options = new SolverOptions(SolverType.Fluid);
    options.iter_max = 200;
    SolverFluid solverFluid = new SolverFluid(model, options);

    solverFluid.options.iter_max = 100;
    solverFluid.options.stiff = false;
    solverFluid.runAnalyzer();
    SolverFluidResult fluidResult = solverFluid.fluidResult;
    SolverResult result = solverFluid.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(3, result.QN.getNumCols());
    assertEquals(6, result.QN.getNumElements());
    assertEquals(0.1429, result.QN.get(0, 0), tol);
    assertEquals(0.1429, result.QN.get(0, 1), tol);
    assertEquals(0, result.QN.get(0, 2), tol);
    assertEquals(0.1429, result.QN.get(1, 0), tol);
    assertEquals(0.5714, result.QN.get(1, 1), tol);
    assertEquals(0, result.QN.get(1, 2), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(3, result.RN.getNumCols());
    assertEquals(6, result.RN.getNumElements());
    assertEquals(1, result.RN.get(0, 0), tol);
    assertEquals(1, result.RN.get(0, 1), tol);
    assertEquals(0, result.RN.get(0, 2), tol);
    assertEquals(1, result.RN.get(1, 0), tol);
    assertEquals(4, result.RN.get(1, 1), tol);
    assertEquals(0, result.RN.get(1, 2), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(3, result.XN.getNumCols());
    assertEquals(3, result.XN.getNumElements());
    assertEquals(0.1429, result.XN.get(0, 0), tol);
    assertEquals(0.1429, result.XN.get(0, 1), tol);
    assertEquals(0, result.XN.get(0, 2), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(3, result.UN.getNumCols());
    assertEquals(6, result.UN.getNumElements());
    assertEquals(0.1429, result.UN.get(0, 0), tol);
    assertEquals(0.1429, result.UN.get(0, 1), tol);
    assertEquals(0, result.UN.get(0, 2), tol);
    assertEquals(0.1429, result.UN.get(1, 0), tol);
    assertEquals(0.5714, result.UN.get(1, 1), tol);
    assertEquals(0, result.UN.get(1, 2), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(3, result.TN.getNumCols());
    assertEquals(6, result.TN.getNumElements());
    assertEquals(0.1429, result.TN.get(0, 0), tol);
    assertEquals(0.1429, result.TN.get(0, 1), tol);
    assertEquals(0, result.TN.get(0, 2), tol);
    assertEquals(0.1429, result.TN.get(1, 0), tol);
    assertEquals(0.1429, result.TN.get(1, 1), tol);
    assertEquals(0, result.TN.get(1, 2), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(3, result.CN.getNumCols());
    assertEquals(3, result.CN.getNumElements());
    assertEquals(7, result.CN.get(0, 0), tol);
    assertEquals(0, result.CN.get(0, 1), tol);
    assertEquals(NaN, result.CN.get(0, 2), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(3, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(1, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[0][1].get(0, 0), tol);
    assertEquals(0, result.QNt[0][2].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][2].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.1429, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.QNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5714, result.QNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.QNt[1][2].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(3, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(1, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[0][1].get(0, 0), tol);
    assertEquals(0, result.UNt[0][2].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][2].get(0, 0), tol);
    assertEquals(0.1429, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.UNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5714, result.UNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.UNt[1][2].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(3, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(1, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[0][1].get(0, 0), tol);
    assertEquals(0, result.TNt[0][2].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(0, result.TNt[1][2].get(0, 0), tol);
    assertEquals(0.1429, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.TNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.TNt[1][2].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(2000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(0.2857, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(0.2857, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(0, fluidResult.odeStateVec.get(0, 6), tol);
  }

  @Test
  public void exampleCdfRespT2ReturnsCorrectResultFromRunAnalyzerStatedepMethod() {

    // Corresponds to "example_cdfRespT_2.m" in LINE
    String modelName = "example_cdfRespT_2";
    Network model = new Network(modelName);

    Delay delay = new Delay(model, "Delay");
    Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);

    ClosedClass closedClass1 = new ClosedClass(model, "Class1", 1, delay, 0);
    ClosedClass closedClass2 = new ClosedClass(model, "Class2", 0, delay, 0);
    ClosedClass closedClass3 = new ClosedClass(model, "Class3", 0, delay, 0);

    delay.setService(closedClass1, new Exp(1.0));
    delay.setService(closedClass2, new Exp(1.0));
    delay.setService(closedClass3, new Exp(1.0));
    queue.setService(closedClass1, new Exp(1.0));
    queue.setService(closedClass2, new Erlang(0.5, 2));
    queue.setService(closedClass3, new Exp(1 / 0.01));

    RoutingMatrix routingMatrix =
            new RoutingMatrix(
                    model,
                    Arrays.asList(closedClass1, closedClass2, closedClass3),
                    Arrays.asList(delay, queue));
    routingMatrix.addConnection(closedClass1, closedClass1, delay, queue, 1.0);
    routingMatrix.addConnection(closedClass2, closedClass2, delay, queue, 1.0);
    routingMatrix.addConnection(closedClass1, closedClass2, queue, delay, 1.0);
    routingMatrix.addConnection(closedClass2, closedClass1, queue, delay, 1.0);
    routingMatrix.addConnection(delay, queue, closedClass3, 1.0);
    routingMatrix.addConnection(queue, delay, closedClass3, 1.0);
    model.link(routingMatrix);

    SolverOptions options = new SolverOptions(SolverType.Fluid);
    options.iter_max = 200;
    SolverFluid solverFluid = new SolverFluid(model, options);

    solverFluid.options.iter_max = 100;
    solverFluid.options.method = "statedep";
    solverFluid.options.stiff = false;
    solverFluid.runAnalyzer();
    SolverFluidResult fluidResult = solverFluid.fluidResult;
    SolverResult result = solverFluid.result;

    // method
    assertEquals("statedep", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(3, result.QN.getNumCols());
    assertEquals(6, result.QN.getNumElements());
    assertEquals(0.1429, result.QN.get(0, 0), tol);
    assertEquals(0.1429, result.QN.get(0, 1), tol);
    assertEquals(0, result.QN.get(0, 2), tol);
    assertEquals(0.1429, result.QN.get(1, 0), tol);
    assertEquals(0.5714, result.QN.get(1, 1), tol);
    assertEquals(0, result.QN.get(1, 2), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(3, result.RN.getNumCols());
    assertEquals(6, result.RN.getNumElements());
    assertEquals(1, result.RN.get(0, 0), tol);
    assertEquals(1, result.RN.get(0, 1), tol);
    assertEquals(0, result.RN.get(0, 2), tol);
    assertEquals(1, result.RN.get(1, 0), tol);
    assertEquals(4, result.RN.get(1, 1), tol);
    assertEquals(0, result.RN.get(1, 2), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(3, result.XN.getNumCols());
    assertEquals(3, result.XN.getNumElements());
    assertEquals(0.1429, result.XN.get(0, 0), tol);
    assertEquals(0.1429, result.XN.get(0, 1), tol);
    assertEquals(0, result.XN.get(0, 2), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(3, result.UN.getNumCols());
    assertEquals(6, result.UN.getNumElements());
    assertEquals(0.1429, result.UN.get(0, 0), tol);
    assertEquals(0.1429, result.UN.get(0, 1), tol);
    assertEquals(0, result.UN.get(0, 2), tol);
    assertEquals(0.1429, result.UN.get(1, 0), tol);
    assertEquals(0.5714, result.UN.get(1, 1), tol);
    assertEquals(0, result.UN.get(1, 2), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(3, result.TN.getNumCols());
    assertEquals(6, result.TN.getNumElements());
    assertEquals(0.1429, result.TN.get(0, 0), tol);
    assertEquals(0.1429, result.TN.get(0, 1), tol);
    assertEquals(0, result.TN.get(0, 2), tol);
    assertEquals(0.1429, result.TN.get(1, 0), tol);
    assertEquals(0.1429, result.TN.get(1, 1), tol);
    assertEquals(0, result.TN.get(1, 2), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(3, result.CN.getNumCols());
    assertEquals(3, result.CN.getNumElements());
    assertEquals(7, result.CN.get(0, 0), tol);
    assertEquals(0, result.CN.get(0, 1), tol);
    assertEquals(NaN, result.CN.get(0, 2), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(3, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(1, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[0][1].get(0, 0), tol);
    assertEquals(0, result.QNt[0][2].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][2].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.1429, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.QNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5714, result.QNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.QNt[1][2].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(3, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(1, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[0][1].get(0, 0), tol);
    assertEquals(0, result.UNt[0][2].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][2].get(0, 0), tol);
    assertEquals(0.1429, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.UNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5714, result.UNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.UNt[1][2].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(3, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(1, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[0][1].get(0, 0), tol);
    assertEquals(0, result.TNt[0][2].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(NaN, result.TNt[1][2].get(0, 0), tol);
    assertEquals(0.1429, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.TNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.1429, result.TNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0, result.TNt[1][2].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(2000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.1429, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(0.2857, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(0.2857, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(0, fluidResult.odeStateVec.get(0, 6), tol);
  }
}
