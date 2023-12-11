package jline.solvers.env;

import jline.examples.RandomEnvironment;
import jline.solvers.SolverResult;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverEnvExamplesTest {

  static double tol = 0.015;

  @Test
  public void exampleRandomEnvironment1ReturnsCorrectResult() {

    SolverEnv envSolver = jline.examples.RandomEnvironment.ex1();
    envSolver.getAvg();
    SolverResult result = envSolver.result;

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(2, result.QN.getNumElements());
    assertEquals(0.5590, result.QN.get(0, 0), 0.5590 * tol);
    assertEquals(0.4410, result.QN.get(1, 0), 0.4410 * tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(2, result.UN.getNumElements());
    assertEquals(0.5590, result.UN.get(0, 0), 0.5590 * tol);
    assertEquals(0.4410, result.UN.get(1, 0), 0.4410 * tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(2, result.TN.getNumElements());
    assertEquals(0.6946, result.TN.get(0, 0), 0.6946 * tol);
    assertEquals(0.6842, result.TN.get(1, 0), 0.6842 * tol);
  }

  @Test
  public void exampleRandomEnvironment2ReturnsCorrectResult() {

    SolverEnv envSolver = jline.examples.RandomEnvironment.ex2();
    envSolver.getAvg();
    SolverResult result = envSolver.result;

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(2, result.QN.getNumElements());
    assertEquals(0.4444, result.QN.get(0, 0), 0.4444 * tol);
    assertEquals(29.556, result.QN.get(1, 0), 29.556 * tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(2, result.UN.getNumElements());
    assertEquals(0.4444, result.UN.get(0, 0), 0.4444 * tol);
    assertEquals(1, result.UN.get(1, 0), 1 * tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(2, result.TN.getNumElements());
    assertEquals(0.9784, result.TN.get(0, 0), 0.9784 * tol);
    assertEquals(1, result.TN.get(1, 0), 1 * tol);
  }

  @Test
  public void exampleRandomEnvironment3ReturnsCorrectResult() {

    SolverEnv envSolver = jline.examples.RandomEnvironment.ex3();
    envSolver.getAvg();
    SolverResult result = envSolver.result;

    // QN
    assertEquals(2, result.QN.getNumRows());
    assertEquals(1, result.QN.getNumCols());
    assertEquals(2, result.QN.getNumElements());
    assertEquals(0.9446, result.QN.get(0, 0), 0.9446 * tol);
    assertEquals(1.0554, result.QN.get(1, 0), 1.0554 * tol);

    // UN
    assertEquals(2, result.UN.getNumRows());
    assertEquals(1, result.UN.getNumCols());
    assertEquals(2, result.UN.getNumElements());
    assertEquals(0.9446, result.UN.get(0, 0), 0.9446 * tol);
    assertEquals(0.8444, result.UN.get(1, 0), 0.8444 * tol);

    // TN
    assertEquals(2, result.TN.getNumRows());
    assertEquals(1, result.TN.getNumCols());
    assertEquals(2, result.TN.getNumElements());
    assertEquals(1.5542, result.TN.get(0, 0), 1.5542 * tol);
    assertEquals(1.5339, result.TN.get(1, 0), 1.5339 * tol);
  }
}
