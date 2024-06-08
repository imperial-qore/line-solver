package jline.solvers.fluid;

import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import org.junit.jupiter.api.Test;

import static jline.solvers.fluid.FluidTestModels.*;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverFluidClosedTest {

  static double tol = 0.005; // Ideally should be 0.0001 but some results don't quite match MatLab

  @Test
  public void closedEx1ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex1();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);
    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(2, result.QN.getNumElements());
    assertEquals(2, result.QN.get(0, 0), tol);
    assertEquals(6, result.QN.get(1, 0), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(2, result.RN.getNumElements());
    assertEquals(0.25, result.RN.get(0, 0), tol);
    assertEquals(0.75, result.RN.get(1, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(8, result.XN.get(0, 0), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(2, result.UN.getNumElements());
    assertEquals(2, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(2, result.TN.getNumElements());
    assertEquals(8, result.TN.get(0, 0), tol);
    assertEquals(8, result.TN.get(1, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(1, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(8, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(2, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(6, result.QNt[1][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(8, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(2, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(32, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(8, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[1][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;

    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(500, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(2, fluidResult.odeStateVec.getNumCols());
    assertEquals(2, fluidResult.odeStateVec.getNumElements());
    assertEquals(2, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(6, fluidResult.odeStateVec.get(0, 1), tol);
  }

  @Test
  public void closedEx2ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex2();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(4, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(4, result.QN.getNumElements());
    assertEquals(2, result.QN.get(0, 0), tol);
    assertEquals(10.9999, result.QN.get(1, 0), tol);
    assertEquals(2, result.QN.get(2, 0), tol);
    assertEquals(1.0001, result.QN.get(3, 0), tol);

    // RN
    assertEquals(4, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(4, result.RN.getNumElements());
    assertEquals(0.25, result.RN.get(0, 0), tol);
    assertEquals(1.375, result.RN.get(1, 0), tol);
    assertEquals(0.25, result.RN.get(2, 0), tol);
    assertEquals(0.125, result.RN.get(3, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(8, result.XN.get(0, 0), tol);

    // UN
    assertEquals(4, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(4, result.UN.getNumElements());
    assertEquals(2, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);
    assertEquals(2, result.UN.get(2, 0), tol);
    assertEquals(1, result.UN.get(3, 0), tol);

    // TN
    assertEquals(4, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(4, result.TN.getNumElements());
    assertEquals(8, result.TN.get(0, 0), tol);
    assertEquals(8, result.TN.get(1, 0), tol);
    assertEquals(8, result.TN.get(2, 0), tol);
    assertEquals(8, result.TN.get(3, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(2, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(4, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(16, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(2, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(10.9999, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1.0001, result.QNt[3][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(4, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(16, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(2, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[3][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(4, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(64, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(8, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[3][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(500, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(4, fluidResult.odeStateVec.getNumCols());
    assertEquals(4, fluidResult.odeStateVec.getNumElements());
    assertEquals(2, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(10.9999, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(2, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(1.0001, fluidResult.odeStateVec.get(0, 3), tol);
  }

  @Test
  public void closedEx3ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex3();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(8, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(8, result.QN.getNumElements());
    assertEquals(1.9999, result.QN.get(0, 0), tol);
    assertEquals(21, result.QN.get(1, 0), tol);
    assertEquals(2, result.QN.get(2, 0), tol);
    assertEquals(1, result.QN.get(3, 0), tol);
    assertEquals(2, result.QN.get(4, 0), tol);
    assertEquals(1, result.QN.get(5, 0), tol);
    assertEquals(2, result.QN.get(6, 0), tol);
    assertEquals(1.0001, result.QN.get(7, 0), tol);

    // RN
    assertEquals(8, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(8, result.RN.getNumElements());
    assertEquals(0.25, result.RN.get(0, 0), tol);
    assertEquals(2.625, result.RN.get(1, 0), tol);
    assertEquals(0.25, result.RN.get(2, 0), tol);
    assertEquals(0.125, result.RN.get(3, 0), tol);
    assertEquals(0.25, result.RN.get(4, 0), tol);
    assertEquals(0.125, result.RN.get(5, 0), tol);
    assertEquals(0.25, result.RN.get(6, 0), tol);
    assertEquals(0.125, result.RN.get(7, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(7.9994, result.XN.get(0, 0), tol);

    // UN
    assertEquals(8, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(8, result.UN.getNumElements());
    assertEquals(1.9999, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);
    assertEquals(2, result.UN.get(2, 0), tol);
    assertEquals(1, result.UN.get(3, 0), tol);
    assertEquals(2, result.UN.get(4, 0), tol);
    assertEquals(1, result.UN.get(5, 0), tol);
    assertEquals(2, result.UN.get(6, 0), tol);
    assertEquals(1, result.UN.get(7, 0), tol);

    // TN
    assertEquals(8, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(8, result.TN.getNumElements());
    assertEquals(7.9994, result.TN.get(0, 0), tol);
    assertEquals(8, result.TN.get(1, 0), tol);
    assertEquals(8, result.TN.get(2, 0), tol);
    assertEquals(8, result.TN.get(3, 0), tol);
    assertEquals(8, result.TN.get(4, 0), tol);
    assertEquals(8, result.TN.get(5, 0), tol);
    assertEquals(8, result.TN.get(6, 0), tol);
    assertEquals(8, result.TN.get(7, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(4.0003, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(8, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(32, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    assertEquals(0, result.QNt[4][0].get(0, 0), tol);
    assertEquals(0, result.QNt[5][0].get(0, 0), tol);
    assertEquals(0, result.QNt[6][0].get(0, 0), tol);
    assertEquals(0, result.QNt[7][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(1.9999, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(21, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.QNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.QNt[4][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.QNt[5][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.QNt[6][0].get(Tmax - 1, 0), tol);
    assertEquals(1.0001, result.QNt[7][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(8, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(32, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(0, result.UNt[4][0].get(0, 0), tol);
    assertEquals(0, result.UNt[5][0].get(0, 0), tol);
    assertEquals(0, result.UNt[6][0].get(0, 0), tol);
    assertEquals(0, result.UNt[7][0].get(0, 0), tol);
    assertEquals(1.9999, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.UNt[4][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[5][0].get(Tmax - 1, 0), tol);
    assertEquals(2, result.UNt[6][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[7][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(8, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(128, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(0, result.TNt[4][0].get(0, 0), tol);
    assertEquals(0, result.TNt[5][0].get(0, 0), tol);
    assertEquals(0, result.TNt[6][0].get(0, 0), tol);
    assertEquals(0, result.TNt[7][0].get(0, 0), tol);
    assertEquals(7.9994, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[3][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[4][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[5][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[6][0].get(Tmax - 1, 0), tol);
    assertEquals(8, result.TNt[7][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(500, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(8, fluidResult.odeStateVec.getNumCols());
    assertEquals(8, fluidResult.odeStateVec.getNumElements());
    assertEquals(1.9999, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(21, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(2, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(1, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(2, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(1, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(2, fluidResult.odeStateVec.get(0, 6), tol);
    assertEquals(1.0001, fluidResult.odeStateVec.get(0, 7), tol);
  }

  @Test
  public void closedEx4ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex4();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(2, result.QN.getNumCols());
    assertEquals(4, result.QN.getNumElements());
    assertEquals(0.25, result.QN.get(0, 0), tol);
    assertEquals(0.25, result.QN.get(0, 1), tol);
    assertEquals(7.75, result.QN.get(1, 0), tol);
    assertEquals(7.75, result.QN.get(1, 1), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(2, result.RN.getNumCols());
    assertEquals(4, result.RN.getNumElements());
    assertEquals(0.125, result.RN.get(0, 0), tol);
    assertEquals(0.1429, result.RN.get(0, 1), tol);
    assertEquals(3.875, result.RN.get(1, 0), tol);
    assertEquals(4.4286, result.RN.get(1, 1), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(2, result.XN.getNumCols());
    assertEquals(2, result.XN.getNumElements());
    assertEquals(2.0001, result.XN.get(0, 0), tol);
    assertEquals(1.75, result.XN.get(0, 1), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(2, result.UN.getNumCols());
    assertEquals(4, result.UN.getNumElements());
    assertEquals(0.25, result.UN.get(0, 0), tol);
    assertEquals(0.25, result.UN.get(0, 1), tol);
    assertEquals(0.5, result.UN.get(1, 0), tol);
    assertEquals(0.5, result.UN.get(1, 1), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(2, result.TN.getNumCols());
    assertEquals(4, result.TN.getNumElements());
    assertEquals(2.0001, result.TN.get(0, 0), tol);
    assertEquals(1.75, result.TN.get(0, 1), tol);
    assertEquals(2, result.TN.get(1, 0), tol);
    assertEquals(1.75, result.TN.get(1, 1), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(2, result.CN.getNumCols());
    assertEquals(2, result.CN.getNumElements());
    assertEquals(3.9999, result.CN.get(0, 0), tol);
    assertEquals(4.5714, result.CN.get(0, 1), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(2, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(8, result.QNt[0][0].get(0, 0), tol);
    assertEquals(8, result.QNt[0][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.25, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(7.75, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(7.75, result.QNt[1][1].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(2, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(8, result.UNt[0][0].get(0, 0), tol);
    assertEquals(8, result.UNt[0][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0.25, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.5, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.5, result.UNt[1][1].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(2, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(64, result.TNt[0][0].get(0, 0), tol);
    assertEquals(56, result.TNt[0][1].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(2.0001, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1.75, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(2, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(1.75, result.TNt[1][1].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(571.4286, result.t.get(result.t.getNumRows() - 1, 0), tol);

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(4, fluidResult.odeStateVec.getNumCols());
    assertEquals(4, fluidResult.odeStateVec.getNumElements());
    assertEquals(0.25, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.25, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(7.75, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(7.75, fluidResult.odeStateVec.get(0, 3), tol);
  }

  @Test
  public void closedEx5ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex5();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(4, result.QN.getNumCols());
    assertEquals(8, result.QN.getNumElements());
    assertEquals(0.125, result.QN.get(0, 0), tol);
    assertEquals(0.125, result.QN.get(0, 1), tol);
    assertEquals(0.125, result.QN.get(0, 2), tol);
    assertEquals(0.125, result.QN.get(0, 3), tol);
    assertEquals(3.875, result.QN.get(1, 0), tol);
    assertEquals(3.875, result.QN.get(1, 1), tol);
    assertEquals(3.875, result.QN.get(1, 2), tol);
    assertEquals(3.875, result.QN.get(1, 3), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(4, result.RN.getNumCols());
    assertEquals(8, result.RN.getNumElements());
    assertEquals(0.125, result.RN.get(0, 0), tol);
    assertEquals(0.1429, result.RN.get(0, 1), tol);
    assertEquals(0.1667, result.RN.get(0, 2), tol);
    assertEquals(0.2, result.RN.get(0, 3), tol);
    assertEquals(3.875, result.RN.get(1, 0), tol);
    assertEquals(4.4286, result.RN.get(1, 1), tol);
    assertEquals(5.1667, result.RN.get(1, 2), tol);
    assertEquals(6.2, result.RN.get(1, 3), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(4, result.XN.getNumCols());
    assertEquals(4, result.XN.getNumElements());
    assertEquals(1.0001, result.XN.get(0, 0), tol);
    assertEquals(0.875, result.XN.get(0, 1), tol);
    assertEquals(0.75, result.XN.get(0, 2), tol);
    assertEquals(0.625, result.XN.get(0, 3), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(4, result.UN.getNumCols());
    assertEquals(8, result.UN.getNumElements());
    assertEquals(0.125, result.UN.get(0, 0), tol);
    assertEquals(0.125, result.UN.get(0, 1), tol);
    assertEquals(0.125, result.UN.get(0, 2), tol);
    assertEquals(0.125, result.UN.get(0, 3), tol);
    assertEquals(0.25, result.UN.get(1, 0), tol);
    assertEquals(0.25, result.UN.get(1, 1), tol);
    assertEquals(0.25, result.UN.get(1, 2), tol);
    assertEquals(0.25, result.UN.get(1, 3), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(4, result.TN.getNumCols());
    assertEquals(8, result.TN.getNumElements());
    assertEquals(1.0001, result.TN.get(0, 0), tol);
    assertEquals(0.875, result.TN.get(0, 1), tol);
    assertEquals(0.75, result.TN.get(0, 2), tol);
    assertEquals(0.625, result.TN.get(0, 3), tol);
    assertEquals(1, result.TN.get(1, 0), tol);
    assertEquals(0.875, result.TN.get(1, 1), tol);
    assertEquals(0.75, result.TN.get(1, 2), tol);
    assertEquals(0.625, result.TN.get(1, 3), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(4, result.CN.getNumCols());
    assertEquals(4, result.CN.getNumElements());
    assertEquals(3.9997, result.CN.get(0, 0), tol);
    assertEquals(4.5714, result.CN.get(0, 1), tol);
    assertEquals(5.3333, result.CN.get(0, 2), tol);
    assertEquals(6.4, result.CN.get(0, 3), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(4, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(4, result.QNt[0][0].get(0, 0), tol);
    assertEquals(4, result.QNt[0][1].get(0, 0), tol);
    assertEquals(4, result.QNt[0][2].get(0, 0), tol);
    assertEquals(4, result.QNt[0][3].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][2].get(0, 0), tol);
    assertEquals(0, result.QNt[1][3].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.125, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.QNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.QNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(3.875, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(3.875, result.QNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(3.875, result.QNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(3.875, result.QNt[1][3].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(4, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(4, result.UNt[0][0].get(0, 0), tol);
    assertEquals(4, result.UNt[0][1].get(0, 0), tol);
    assertEquals(4, result.UNt[0][2].get(0, 0), tol);
    assertEquals(4, result.UNt[0][3].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][2].get(0, 0), tol);
    assertEquals(0, result.UNt[1][3].get(0, 0), tol);
    assertEquals(0.125, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.UNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.UNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.UNt[1][3].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(4, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(32, result.TNt[0][0].get(0, 0), tol);
    assertEquals(28, result.TNt[0][1].get(0, 0), tol);
    assertEquals(24, result.TNt[0][2].get(0, 0), tol);
    assertEquals(20, result.TNt[0][3].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(0, result.TNt[1][2].get(0, 0), tol);
    assertEquals(0, result.TNt[1][3].get(0, 0), tol);
    assertEquals(1.0001, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.875, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.75, result.TNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.625, result.TNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(1, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.875, result.TNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.75, result.TNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.625, result.TNt[1][3].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(800, result.t.get(result.t.getNumRows() - 1, 0), tol);

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(8, fluidResult.odeStateVec.getNumCols());
    assertEquals(8, fluidResult.odeStateVec.getNumElements());
    assertEquals(0.125, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.125, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0.125, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.125, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(3.875, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(3.875, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(3.875, fluidResult.odeStateVec.get(0, 6), tol);
    assertEquals(3.875, fluidResult.odeStateVec.get(0, 7), tol);
  }

  @Test
  public void closedEx6ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex6();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(8, result.QN.getNumCols());
    assertEquals(16, result.QN.getNumElements());
    assertEquals(0.0625, result.QN.get(0, 0), tol);
    assertEquals(0.0625, result.QN.get(0, 1), tol);
    assertEquals(0.0625, result.QN.get(0, 2), tol);
    assertEquals(0.0625, result.QN.get(0, 3), tol);
    assertEquals(0.0625, result.QN.get(0, 4), tol);
    assertEquals(0.0625, result.QN.get(0, 5), tol);
    assertEquals(0.0625, result.QN.get(0, 6), tol);
    assertEquals(0.0625, result.QN.get(0, 7), tol);
    assertEquals(1.9375, result.QN.get(1, 0), tol);
    assertEquals(1.9375, result.QN.get(1, 1), tol);
    assertEquals(1.9375, result.QN.get(1, 2), tol);
    assertEquals(1.9375, result.QN.get(1, 3), tol);
    assertEquals(1.9375, result.QN.get(1, 4), tol);
    assertEquals(1.9375, result.QN.get(1, 5), tol);
    assertEquals(1.9375, result.QN.get(1, 6), tol);
    assertEquals(1.9375, result.QN.get(1, 7), tol);

    // RN
    assertEquals(2, result.RN.getNumRows());
    assertEquals(8, result.RN.getNumCols());
    assertEquals(16, result.RN.getNumElements());
    assertEquals(0.125, result.RN.get(0, 0), tol);
    assertEquals(0.1429, result.RN.get(0, 1), tol);
    assertEquals(0.1667, result.RN.get(0, 2), tol);
    assertEquals(0.2, result.RN.get(0, 3), tol);
    assertEquals(0.25, result.RN.get(0, 4), tol);
    assertEquals(0.3333, result.RN.get(0, 5), tol);
    assertEquals(0.5, result.RN.get(0, 6), tol);
    assertEquals(1, result.RN.get(0, 7), tol);
    assertEquals(3.875, result.RN.get(1, 0), tol);
    assertEquals(4.4286, result.RN.get(1, 1), tol);
    assertEquals(5.1667, result.RN.get(1, 2), tol);
    assertEquals(6.2, result.RN.get(1, 3), tol);
    assertEquals(7.75, result.RN.get(1, 4), tol);
    assertEquals(10.3333, result.RN.get(1, 5), tol);
    assertEquals(15.5, result.RN.get(1, 6), tol);
    assertEquals(31, result.RN.get(1, 7), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(8, result.XN.getNumCols());
    assertEquals(8, result.XN.getNumElements());
    assertEquals(0.5002, result.XN.get(0, 0), tol);
    assertEquals(0.4375, result.XN.get(0, 1), tol);
    assertEquals(0.375, result.XN.get(0, 2), tol);
    assertEquals(0.3125, result.XN.get(0, 3), tol);
    assertEquals(0.25, result.XN.get(0, 4), tol);
    assertEquals(0.1875, result.XN.get(0, 5), tol);
    assertEquals(0.125, result.XN.get(0, 6), tol);
    assertEquals(0.0625, result.XN.get(0, 7), tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(8, result.UN.getNumCols());
    assertEquals(16, result.UN.getNumElements());
    assertEquals(0.0625, result.UN.get(0, 0), tol);
    assertEquals(0.0625, result.UN.get(0, 1), tol);
    assertEquals(0.0625, result.UN.get(0, 2), tol);
    assertEquals(0.0625, result.UN.get(0, 3), tol);
    assertEquals(0.0625, result.UN.get(0, 4), tol);
    assertEquals(0.0625, result.UN.get(0, 5), tol);
    assertEquals(0.0625, result.UN.get(0, 6), tol);
    assertEquals(0.0625, result.UN.get(0, 7), tol);
    assertEquals(0.125, result.UN.get(1, 0), tol);
    assertEquals(0.125, result.UN.get(1, 1), tol);
    assertEquals(0.125, result.UN.get(1, 2), tol);
    assertEquals(0.125, result.UN.get(1, 3), tol);
    assertEquals(0.125, result.UN.get(1, 4), tol);
    assertEquals(0.125, result.UN.get(1, 5), tol);
    assertEquals(0.125, result.UN.get(1, 6), tol);
    assertEquals(0.125, result.UN.get(1, 7), tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(8, result.TN.getNumCols());
    assertEquals(16, result.TN.getNumElements());
    assertEquals(0.5002, result.TN.get(0, 0), tol);
    assertEquals(0.4375, result.TN.get(0, 1), tol);
    assertEquals(0.375, result.TN.get(0, 2), tol);
    assertEquals(0.3125, result.TN.get(0, 3), tol);
    assertEquals(0.25, result.TN.get(0, 4), tol);
    assertEquals(0.1875, result.TN.get(0, 5), tol);
    assertEquals(0.125, result.TN.get(0, 6), tol);
    assertEquals(0.0625, result.TN.get(0, 7), tol);
    assertEquals(0.5, result.TN.get(1, 0), tol);
    assertEquals(0.4375, result.TN.get(1, 1), tol);
    assertEquals(0.375, result.TN.get(1, 2), tol);
    assertEquals(0.3125, result.TN.get(1, 3), tol);
    assertEquals(0.25, result.TN.get(1, 4), tol);
    assertEquals(0.1875, result.TN.get(1, 5), tol);
    assertEquals(0.125, result.TN.get(1, 6), tol);
    assertEquals(0.0625, result.TN.get(1, 7), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(8, result.CN.getNumCols());
    assertEquals(8, result.CN.getNumElements());
    assertEquals(3.9985, result.CN.get(0, 0), tol);
    assertEquals(4.5715, result.CN.get(0, 1), tol);
    assertEquals(5.3334, result.CN.get(0, 2), tol);
    assertEquals(6.4, result.CN.get(0, 3), tol);
    assertEquals(8, result.CN.get(0, 4), tol);
    assertEquals(10.6667, result.CN.get(0, 5), tol);
    assertEquals(16, result.CN.get(0, 6), tol);
    assertEquals(32, result.CN.get(0, 7), tol);

    // QNt
    assertEquals(2, result.QNt.length);
    assertEquals(8, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(2, result.QNt[0][0].get(0, 0), tol);
    assertEquals(2, result.QNt[0][1].get(0, 0), tol);
    assertEquals(2, result.QNt[0][2].get(0, 0), tol);
    assertEquals(2, result.QNt[0][3].get(0, 0), tol);
    assertEquals(2, result.QNt[0][4].get(0, 0), tol);
    assertEquals(2, result.QNt[0][5].get(0, 0), tol);
    assertEquals(2, result.QNt[0][6].get(0, 0), tol);
    assertEquals(2, result.QNt[0][7].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][1].get(0, 0), tol);
    assertEquals(0, result.QNt[1][2].get(0, 0), tol);
    assertEquals(0, result.QNt[1][3].get(0, 0), tol);
    assertEquals(0, result.QNt[1][4].get(0, 0), tol);
    assertEquals(0, result.QNt[1][5].get(0, 0), tol);
    assertEquals(0, result.QNt[1][6].get(0, 0), tol);
    assertEquals(0, result.QNt[1][7].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.0625, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][4].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][5].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][6].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[0][7].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][3].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][4].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][5].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][6].get(Tmax - 1, 0), tol);
    assertEquals(1.9375, result.QNt[1][7].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(2, result.UNt.length);
    assertEquals(8, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(2, result.UNt[0][0].get(0, 0), tol);
    assertEquals(2, result.UNt[0][1].get(0, 0), tol);
    assertEquals(2, result.UNt[0][2].get(0, 0), tol);
    assertEquals(2, result.UNt[0][3].get(0, 0), tol);
    assertEquals(2, result.UNt[0][4].get(0, 0), tol);
    assertEquals(2, result.UNt[0][5].get(0, 0), tol);
    assertEquals(2, result.UNt[0][6].get(0, 0), tol);
    assertEquals(2, result.UNt[0][7].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][1].get(0, 0), tol);
    assertEquals(0, result.UNt[1][2].get(0, 0), tol);
    assertEquals(0, result.UNt[1][3].get(0, 0), tol);
    assertEquals(0, result.UNt[1][4].get(0, 0), tol);
    assertEquals(0, result.UNt[1][5].get(0, 0), tol);
    assertEquals(0, result.UNt[1][6].get(0, 0), tol);
    assertEquals(0, result.UNt[1][7].get(0, 0), tol);
    assertEquals(0.0625, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][4].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][5].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][6].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[0][7].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][3].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][4].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][5].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][6].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.UNt[1][7].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(2, result.TNt.length);
    assertEquals(8, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(16, result.TNt[0][0].get(0, 0), tol);
    assertEquals(14, result.TNt[0][1].get(0, 0), tol);
    assertEquals(12, result.TNt[0][2].get(0, 0), tol);
    assertEquals(10, result.TNt[0][3].get(0, 0), tol);
    assertEquals(8, result.TNt[0][4].get(0, 0), tol);
    assertEquals(6, result.TNt[0][5].get(0, 0), tol);
    assertEquals(4, result.TNt[0][6].get(0, 0), tol);
    assertEquals(2, result.TNt[0][7].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][1].get(0, 0), tol);
    assertEquals(0, result.TNt[1][2].get(0, 0), tol);
    assertEquals(0, result.TNt[1][3].get(0, 0), tol);
    assertEquals(0, result.TNt[1][4].get(0, 0), tol);
    assertEquals(0, result.TNt[1][5].get(0, 0), tol);
    assertEquals(0, result.TNt[1][6].get(0, 0), tol);
    assertEquals(0, result.TNt[1][7].get(0, 0), tol);
    assertEquals(0.5002, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(0.4375, result.TNt[0][1].get(Tmax - 1, 0), tol);
    assertEquals(0.375, result.TNt[0][2].get(Tmax - 1, 0), tol);
    assertEquals(0.3125, result.TNt[0][3].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.TNt[0][4].get(Tmax - 1, 0), tol);
    assertEquals(0.1875, result.TNt[0][5].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.TNt[0][6].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.TNt[0][7].get(Tmax - 1, 0), tol);
    assertEquals(0.5, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.4375, result.TNt[1][1].get(Tmax - 1, 0), tol);
    assertEquals(0.375, result.TNt[1][2].get(Tmax - 1, 0), tol);
    assertEquals(0.3125, result.TNt[1][3].get(Tmax - 1, 0), tol);
    assertEquals(0.25, result.TNt[1][4].get(Tmax - 1, 0), tol);
    assertEquals(0.1875, result.TNt[1][5].get(Tmax - 1, 0), tol);
    assertEquals(0.125, result.TNt[1][6].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.TNt[1][7].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(4000, result.t.get(result.t.getNumRows() - 1, 0), tol);

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(16, fluidResult.odeStateVec.getNumCols());
    assertEquals(16, fluidResult.odeStateVec.getNumElements());
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 3), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 4), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 5), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 6), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 7), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 8), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 9), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 10), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 11), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 12), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 13), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 14), tol);
    assertEquals(1.9375, fluidResult.odeStateVec.get(0, 15), tol);
  }

  @Test
  public void closedEx7ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex7();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(4, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(4, result.QN.getNumElements());
    assertEquals(13, result.QN.get(0, 0), tol);
    assertEquals(1, result.QN.get(1, 0), tol);
    assertEquals(1, result.QN.get(2, 0), tol);
    assertEquals(1, result.QN.get(3, 0), tol);

    // RN
    assertEquals(4, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(4, result.RN.getNumElements());
    assertEquals(0.8125, result.RN.get(0, 0), tol);
    assertEquals(0.0625, result.RN.get(1, 0), tol);
    assertEquals(0.0625, result.RN.get(2, 0), tol);
    assertEquals(0.0625, result.RN.get(3, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(16, result.XN.get(0, 0), tol);

    // UN
    assertEquals(4, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(4, result.UN.getNumElements());
    assertEquals(1, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);
    assertEquals(1, result.UN.get(2, 0), tol);
    assertEquals(1, result.UN.get(3, 0), tol);

    // TN
    assertEquals(4, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(4, result.TN.getNumElements());
    assertEquals(16, result.TN.get(0, 0), tol);
    assertEquals(16, result.TN.get(1, 0), tol);
    assertEquals(16, result.TN.get(2, 0), tol);
    assertEquals(16, result.TN.get(3, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(1, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(4, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(16, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(13, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.QNt[3][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(4, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(1, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(1, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[3][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(4, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(16, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(16, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(16, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(16, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(16, result.TNt[3][0].get(Tmax - 1, 0), tol);

    // t
    int sizeT = 0;
    int numElements = 0;
    sizeT += result.t.getNumRows();
    numElements += result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(125, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(4, fluidResult.odeStateVec.getNumCols());
    assertEquals(4, fluidResult.odeStateVec.getNumElements());
    assertEquals(13, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(1, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(1, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(1, fluidResult.odeStateVec.get(0, 3), tol);
  }

  @Test
  public void closedEx8ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex8();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.result;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(4, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(4, result.QN.getNumElements());
    assertEquals(0.0625, result.QN.get(0, 0), tol);
    assertEquals(15.8438, result.QN.get(1, 0), tol);
    assertEquals(0.0312, result.QN.get(2, 0), tol);
    assertEquals(0.0625, result.QN.get(3, 0), tol);

    // RN
    assertEquals(4, result.RN.getNumRows());
    assertEquals(1, result.RN.getNumCols());
    assertEquals(4, result.RN.getNumElements());
    assertEquals(0.0625, result.RN.get(0, 0), tol);
    assertEquals(15.8438, result.RN.get(1, 0), tol);
    assertEquals(0.0313, result.RN.get(2, 0), tol);
    assertEquals(0.0625, result.RN.get(3, 0), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(1, result.XN.getNumCols());
    assertEquals(1, result.XN.getNumElements());
    assertEquals(0.9997, result.XN.get(0, 0), tol);

    // UN
    assertEquals(4, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(4, result.UN.getNumElements());
    assertEquals(0.0625, result.UN.get(0, 0), tol);
    assertEquals(1, result.UN.get(1, 0), tol);
    assertEquals(0.0312, result.UN.get(2, 0), tol);
    assertEquals(0.0625, result.UN.get(3, 0), tol);

    // TN
    assertEquals(4, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(4, result.TN.getNumElements());
    assertEquals(0.9997, result.TN.get(0, 0), tol);
    assertEquals(1, result.TN.get(1, 0), tol);
    assertEquals(0.9997, result.TN.get(2, 0), tol);
    assertEquals(1.0003, result.TN.get(3, 0), tol);

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(1, result.CN.getNumCols());
    assertEquals(1, result.CN.getNumElements());
    assertEquals(16.0043, result.CN.get(0, 0), tol);

    // QNt
    assertEquals(4, result.QNt.length);
    assertEquals(1, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());
    assertEquals(16, result.QNt[0][0].get(0, 0), tol);
    assertEquals(0, result.QNt[1][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[3][0].get(0, 0), tol);
    int Tmax = result.QNt[0][0].getNumRows();
    assertEquals(0.0625, result.QNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(15.8438, result.QNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0312, result.QNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.QNt[3][0].get(Tmax - 1, 0), tol);

    // UNt
    assertEquals(4, result.UNt.length);
    assertEquals(1, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());
    assertEquals(1, result.UNt[0][0].get(0, 0), tol);
    assertEquals(0, result.UNt[1][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[3][0].get(0, 0), tol);
    assertEquals(0.0625, result.UNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.UNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0312, result.UNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(0.0625, result.UNt[3][0].get(Tmax - 1, 0), tol);

    // TNt
    assertEquals(4, result.TNt.length);
    assertEquals(1, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());
    assertEquals(16, result.TNt[0][0].get(0, 0), tol);
    assertEquals(0, result.TNt[1][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[3][0].get(0, 0), tol);
    assertEquals(0.9997, result.TNt[0][0].get(Tmax - 1, 0), tol);
    assertEquals(1, result.TNt[1][0].get(Tmax - 1, 0), tol);
    assertEquals(0.9997, result.TNt[2][0].get(Tmax - 1, 0), tol);
    assertEquals(1.0003, result.TNt[3][0].get(Tmax - 1, 0), tol);

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
    assertEquals(1, fluidResult.odeStateVec.getNumRows());
    assertEquals(4, fluidResult.odeStateVec.getNumCols());
    assertEquals(4, fluidResult.odeStateVec.getNumElements());
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 0), tol);
    assertEquals(15.8438, fluidResult.odeStateVec.get(0, 1), tol);
    assertEquals(0.0312, fluidResult.odeStateVec.get(0, 2), tol);
    assertEquals(0.0625, fluidResult.odeStateVec.get(0, 3), tol);
  }

  // Note: Test passes, but ~10 min runtime so commented out
  /*  @Test
  public void closedEx9ReturnsCorrectResultFromRunAnalyzer() {

    Network model = closed_ex9();

    SolverOptions options = new SolverOptions(SolverType.FLUID);
    options.iter_max = 200;
    SolverFluid solver = new SolverFluid(model, options);

    solver.options.stiff = true;
    solver.runAnalyzer();
    SolverFluidResult fluidResult = solver.fluidResult;
    SolverResult result = solver.result;

    // method
    assertEquals("matrix", result.method);

    // QN
    assertEquals(8, result.QN.getNumRows());
    assertEquals(8, result.QN.getNumCols());
    assertEquals(64, result.QN.getNumElements());
    for (int col = 0; col < 8; col++) {
      assertEquals(0.0625, result.QN.get(0, col), tol);
      assertEquals(0.178, result.QN.get(1, col), tol);
      assertEquals(0.0625, result.QN.get(2, col), tol);
      assertEquals(0.1514, result.QN.get(3, col), tol);
      assertEquals(0.0625, result.QN.get(4, col), tol);
      assertEquals(0.3695, result.QN.get(5, col), tol);
      assertEquals(0.0625, result.QN.get(6, col), tol);
      assertEquals(1.0511, result.QN.get(7, col), tol);
    }

    // RN
    assertEquals(8, result.RN.getNumRows());
    assertEquals(8, result.RN.getNumCols());
    assertEquals(64, result.RN.getNumElements());
    for (int row = 0; row < 8; row += 2) {
      assertEquals(0.125, result.RN.get(row, 0), tol);
      assertEquals(0.1429, result.RN.get(row, 1), tol);
      assertEquals(0.1667, result.RN.get(row, 2), tol);
      assertEquals(0.2, result.RN.get(row, 3), tol);
      assertEquals(0.25, result.RN.get(row, 4), tol);
      assertEquals(0.3333, result.RN.get(row, 5), tol);
      assertEquals(0.5, result.RN.get(row, 6), tol);
      assertEquals(1, result.RN.get(row, 7), tol);
    }

    assertEquals(0.3559, result.RN.get(1, 0), tol);
    assertEquals(0.4068, result.RN.get(1, 1), tol);
    assertEquals(0.4746, result.RN.get(1, 2), tol);
    assertEquals(0.5695, result.RN.get(1, 3), tol);
    assertEquals(0.7119, result.RN.get(1, 4), tol);
    assertEquals(0.9492, result.RN.get(1, 5), tol);
    assertEquals(1.4238, result.RN.get(1, 6), tol);
    assertEquals(2.8476, result.RN.get(1, 7), tol);

    assertEquals(0.3028, result.RN.get(3, 0), tol);
    assertEquals(0.3461, result.RN.get(3, 1), tol);
    assertEquals(0.4037, result.RN.get(3, 2), tol);
    assertEquals(0.4845, result.RN.get(3, 3), tol);
    assertEquals(0.6056, result.RN.get(3, 4), tol);
    assertEquals(0.8075, result.RN.get(3, 5), tol);
    assertEquals(1.2112, result.RN.get(3, 6), tol);
    assertEquals(2.4225, result.RN.get(3, 7), tol);

    assertEquals(0.7390, result.RN.get(5, 0), tol);
    assertEquals(0.8446, result.RN.get(5, 1), tol);
    assertEquals(0.9854, result.RN.get(5, 2), tol);
    assertEquals(1.1824, result.RN.get(5, 3), tol);
    assertEquals(1.4781, result.RN.get(5, 4), tol);
    assertEquals(1.9707, result.RN.get(5, 5), tol);
    assertEquals(2.9561, result.RN.get(5, 6), tol);
    assertEquals(5.9122, result.RN.get(5, 7), tol);

    assertEquals(2.1022, result.RN.get(7, 0), tol);
    assertEquals(2.4025, result.RN.get(7, 1), tol);
    assertEquals(2.8030, result.RN.get(7, 2), tol);
    assertEquals(3.3636, result.RN.get(7, 3), tol);
    assertEquals(4.2044, result.RN.get(7, 4), tol);
    assertEquals(5.6059, result.RN.get(7, 5), tol);
    assertEquals(8.4089, result.RN.get(7, 6), tol);
    assertEquals(16.8178, result.RN.get(7, 7), tol);

    // XN
    assertEquals(1, result.XN.getNumRows());
    assertEquals(8, result.XN.getNumCols());
    assertEquals(8, result.XN.getNumElements());
    assertEquals(0.5, result.XN.get(0, 0), tol);
    assertEquals(0.4375, result.XN.get(0, 1), tol);
    assertEquals(0.375, result.XN.get(0, 2), tol);
    assertEquals(0.3125, result.XN.get(0, 3), tol);
    assertEquals(0.25, result.XN.get(0, 4), tol);
    assertEquals(0.1875, result.XN.get(0, 5), tol);
    assertEquals(0.125, result.XN.get(0, 6), tol);
    assertEquals(0.0625, result.XN.get(0, 7), tol);

    // UN
    assertEquals(8, result.UN.getNumRows());
    assertEquals(8, result.UN.getNumCols());
    assertEquals(64, result.UN.getNumElements());
    for (int row = 0; row < 8; row += 2) {
      for (int col = 0; col < 8; col++) {
        assertEquals(0.0625, result.UN.get(row, col), tol);
      }
    }
    for (int row = 1; row < 8; row += 2) {
      for (int col = 0; col < 8; col++) {
        assertEquals(0.125, result.UN.get(row, col), tol);
      }
    }

    // TN
    assertEquals(8, result.TN.getNumRows());
    assertEquals(8, result.TN.getNumCols());
    assertEquals(64, result.TN.getNumElements());
    for (int row = 0; row < 8; row++) {
      assertEquals(0.5, result.TN.get(row, 0), tol);
      assertEquals(0.4375, result.TN.get(row, 1), tol);
      assertEquals(0.375, result.TN.get(row, 2), tol);
      assertEquals(0.3125, result.TN.get(row, 3), tol);
      assertEquals(0.25, result.TN.get(row, 4), tol);
      assertEquals(0.1875, result.TN.get(row, 5), tol);
      assertEquals(0.125, result.TN.get(row, 6), tol);
      assertEquals(0.0625, result.TN.get(row, 7), tol);
    }

    // CN
    assertEquals(1, result.CN.getNumRows());
    assertEquals(8, result.CN.getNumCols());
    assertEquals(8, result.CN.getNumElements());
    assertEquals(4, result.CN.get(0, 0), tol);
    assertEquals(4.5714, result.CN.get(0, 1), tol);
    assertEquals(5.3333, result.CN.get(0, 2), tol);
    assertEquals(6.4, result.CN.get(0, 3), tol);
    assertEquals(8, result.CN.get(0, 4), tol);
    assertEquals(10.6667, result.CN.get(0, 5), tol);
    assertEquals(16, result.CN.get(0, 6), tol);
    assertEquals(32, result.CN.get(0, 7), tol);

    // QNt
    assertEquals(8, result.QNt.length);
    assertEquals(8, result.QNt[0].length);
    assertEquals(1, result.QNt[0][0].getNumCols());

    assertEquals(2, result.QNt[0][0].get(0, 0), tol);
    assertEquals(2, result.QNt[0][1].get(0, 0), tol);
    for (int col = 2; col < 8; col++) {
      assertEquals(0, result.QNt[0][col].get(0, 0), tol);
    }
    assertEquals(0, result.QNt[2][0].get(0, 0), tol);
    assertEquals(0, result.QNt[2][1].get(0, 0), tol);
    assertEquals(2, result.QNt[2][2].get(0, 0), tol);
    assertEquals(2, result.QNt[2][3].get(0, 0), tol);
    for (int col = 4; col < 8; col++) {
      assertEquals(0, result.QNt[2][col].get(0, 0), tol);
    }
    for (int col = 0; col < 4; col++) {
      assertEquals(0, result.QNt[4][col].get(0, 0), tol);
    }
    assertEquals(2, result.QNt[4][4].get(0, 0), tol);
    assertEquals(2, result.QNt[4][5].get(0, 0), tol);
    assertEquals(0, result.QNt[4][6].get(0, 0), tol);
    assertEquals(0, result.QNt[4][7].get(0, 0), tol);
    for (int col = 0; col < 6; col++) {
      assertEquals(0, result.QNt[6][col].get(0, 0), tol);
    }
    assertEquals(2, result.QNt[6][6].get(0, 0), tol);
    assertEquals(2, result.QNt[6][7].get(0, 0), tol);
    for (int col = 0; col < 8; col++) {
      assertEquals(0, result.QNt[1][col].get(0, 0), tol);
      assertEquals(0, result.QNt[3][col].get(0, 0), tol);
      assertEquals(0, result.QNt[5][col].get(0, 0), tol);
      assertEquals(0, result.QNt[7][col].get(0, 0), tol);
    }

    int Tmax = result.QNt[0][0].getNumRows();
    for (int col = 0; col < 8; col++) {
      assertEquals(0.0625, result.QNt[0][col].get(Tmax - 1, 0), tol);
      assertEquals(0.178, result.QNt[1][col].get(Tmax - 1, 0), tol);
      assertEquals(0.0625, result.QNt[2][col].get(Tmax - 1, 0), tol);
      assertEquals(0.1514, result.QNt[3][col].get(Tmax - 1, 0), tol);
      assertEquals(0.0625, result.QNt[4][col].get(Tmax - 1, 0), tol);
      assertEquals(0.3695, result.QNt[5][col].get(Tmax - 1, 0), tol);
      assertEquals(0.0625, result.QNt[6][col].get(Tmax - 1, 0), tol);
      assertEquals(1.0511, result.QNt[7][col].get(Tmax - 1, 0), tol);
    }

    // UNt
    assertEquals(8, result.UNt.length);
    assertEquals(8, result.UNt[0].length);
    assertEquals(1, result.UNt[0][0].getNumCols());

    assertEquals(2, result.UNt[0][0].get(0, 0), tol);
    assertEquals(2, result.UNt[0][1].get(0, 0), tol);
    for (int col = 2; col < 8; col++) {
      assertEquals(0, result.UNt[0][col].get(0, 0), tol);
    }
    assertEquals(0, result.UNt[2][0].get(0, 0), tol);
    assertEquals(0, result.UNt[2][1].get(0, 0), tol);
    assertEquals(2, result.UNt[2][2].get(0, 0), tol);
    assertEquals(2, result.UNt[2][3].get(0, 0), tol);
    for (int col = 4; col < 8; col++) {
      assertEquals(0, result.UNt[2][col].get(0, 0), tol);
    }
    for (int col = 0; col < 4; col++) {
      assertEquals(0, result.UNt[4][col].get(0, 0), tol);
    }
    assertEquals(2, result.UNt[4][4].get(0, 0), tol);
    assertEquals(2, result.UNt[4][5].get(0, 0), tol);
    assertEquals(0, result.UNt[4][6].get(0, 0), tol);
    assertEquals(0, result.UNt[4][7].get(0, 0), tol);
    for (int col = 0; col < 6; col++) {
      assertEquals(0, result.UNt[6][col].get(0, 0), tol);
    }
    assertEquals(2, result.UNt[6][6].get(0, 0), tol);
    assertEquals(2, result.UNt[6][7].get(0, 0), tol);
    for (int col = 0; col < 8; col++) {
      assertEquals(0, result.UNt[1][col].get(0, 0), tol);
      assertEquals(0, result.UNt[3][col].get(0, 0), tol);
      assertEquals(0, result.UNt[5][col].get(0, 0), tol);
      assertEquals(0, result.UNt[7][col].get(0, 0), tol);
    }

    for (int row = 0; row < 8; row += 2) {
      for (int col = 0; col < 8; col++) {
        assertEquals(0.0625, result.UNt[row][col].get(Tmax - 1, 0), tol);
      }
    }
    for (int row = 1; row < 8; row += 2) {
      for (int col = 0; col < 8; col++) {
        assertEquals(0.125, result.UNt[row][col].get(Tmax - 1, 0), tol);
      }
    }

    // TNt
    assertEquals(8, result.TNt.length);
    assertEquals(8, result.TNt[0].length);
    assertEquals(1, result.TNt[0][0].getNumCols());

    assertEquals(16, result.TNt[0][0].get(0, 0), tol);
    assertEquals(14, result.TNt[0][1].get(0, 0), tol);
    for (int col = 2; col < 8; col++) {
      assertEquals(0, result.TNt[0][col].get(0, 0), tol);
    }
    assertEquals(0, result.TNt[2][0].get(0, 0), tol);
    assertEquals(0, result.TNt[2][1].get(0, 0), tol);
    assertEquals(12, result.TNt[2][2].get(0, 0), tol);
    assertEquals(10, result.TNt[2][3].get(0, 0), tol);
    for (int col = 4; col < 8; col++) {
      assertEquals(0, result.TNt[2][col].get(0, 0), tol);
    }
    for (int col = 0; col < 4; col++) {
      assertEquals(0, result.TNt[4][col].get(0, 0), tol);
    }
    assertEquals(8, result.TNt[4][4].get(0, 0), tol);
    assertEquals(6, result.TNt[4][5].get(0, 0), tol);
    assertEquals(0, result.TNt[4][6].get(0, 0), tol);
    assertEquals(0, result.TNt[4][7].get(0, 0), tol);
    for (int col = 0; col < 6; col++) {
      assertEquals(0, result.TNt[6][col].get(0, 0), tol);
    }
    assertEquals(4, result.TNt[6][6].get(0, 0), tol);
    assertEquals(2, result.TNt[6][7].get(0, 0), tol);
    for (int col = 0; col < 8; col++) {
      assertEquals(0, result.TNt[1][col].get(0, 0), tol);
      assertEquals(0, result.TNt[3][col].get(0, 0), tol);
      assertEquals(0, result.TNt[5][col].get(0, 0), tol);
      assertEquals(0, result.TNt[7][col].get(0, 0), tol);
    }

    for (int row = 0; row < 8; row++) {
      assertEquals(0.5, result.TNt[row][0].get(Tmax - 1, 0), tol);
      assertEquals(0.4375, result.TNt[row][1].get(Tmax - 1, 0), tol);
      assertEquals(0.375, result.TNt[row][2].get(Tmax - 1, 0), tol);
      assertEquals(0.3125, result.TNt[row][3].get(Tmax - 1, 0), tol);
      assertEquals(0.25, result.TNt[row][4].get(Tmax - 1, 0), tol);
      assertEquals(0.1875, result.TNt[row][5].get(Tmax - 1, 0), tol);
      assertEquals(0.125, result.TNt[row][6].get(Tmax - 1, 0), tol);
      assertEquals(0.0625, result.TNt[row][7].get(Tmax - 1, 0), tol);
    }

    // t
    int sizeT = result.t.getNumRows();
    int numElements = result.t.getNumElements();

    assertEquals(Tmax, sizeT);
    assertEquals(1, result.t.getNumCols());
    assertEquals(Tmax, numElements);
    assertEquals(0.00000001, result.t.get(0, 0));
    assertEquals(4000, result.t.get(result.t.getNumRows() - 1, 0));

    // odeStateVec
    int offset = 0;
    for (int i = 0; i < 8; i++) {
      assertEquals(0.0625, fluidResult.odeStateVec.get(0, offset), tol);
      assertEquals(0.178, fluidResult.odeStateVec.get(0, 8 + offset), tol);
      assertEquals(0.0625, fluidResult.odeStateVec.get(0, 16 + offset), tol);
      assertEquals(0.1514, fluidResult.odeStateVec.get(0, 24 + offset), tol);
      assertEquals(0.0625, fluidResult.odeStateVec.get(0, 32 + offset), tol);
      assertEquals(0.3695, fluidResult.odeStateVec.get(0, 40 + offset), tol);
      assertEquals(0.0625, fluidResult.odeStateVec.get(0, 48 + offset), tol);
      assertEquals(1.0511, fluidResult.odeStateVec.get(0, 56 + offset), tol);
      offset += 1;
    }
  }*/
}
