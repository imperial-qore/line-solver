package jline.solvers.nc;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.SolverOptions;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static jline.solvers.nc.SolverNC.solver_nc_joint;
import static jline.solvers.nc.SolverNC.solver_nc_marg;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverNCTest {

  double tolerance = 1e-10;

  @Test
  public void SolverNCTestModel1() {
    Network model = new Network("example_closedModel");

    Delay node1 = new Delay(model, "Delay");
    Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

    ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

    node1.setService(jobclass1, Exp.fitMean(1.000000));
    node2.setService(jobclass1, Exp.fitMean(1.500000));

    RoutingMatrix routingMatrix = new RoutingMatrix(model,
            Arrays.asList(jobclass1),
            Arrays.asList(node1, node2));

    routingMatrix.addConnection(jobclass1, jobclass1, node1, node1, 0.700000); // (Delay,Class1) -> (Delay,Class1)
    routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
    routingMatrix.addConnection(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)

    model.link(routingMatrix);

    NetworkStruct sn = model.getStruct(true);

    SolverOptions options = new SolverOptions(SolverType.NC);


    SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, null);
    Matrix lPr1 = ret1.lPr;
    double G1 = ret1.G;
    double runtime1 = ret1.runtime;

    assertEquals(2, lPr1.numRows);
    assertEquals(1, lPr1.numCols);

    assertEquals(-9.341536168923934, lPr1.get(0), 9.341536168923934*tolerance);
    assertEquals(-9.341536168923934, lPr1.get(1), 9.341536168923934*tolerance);

    assertEquals(0.003142060751148, G1, tolerance);

    SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
    double Pr2 = ret2.Pr;
    double G2 = ret2.G;
    double lG2 = ret2.lG;
    double runtime2 = ret2.runtime;

    assertEquals(8.770460346419127e-05, Pr2, tolerance);
    assertEquals(0.003142060751148, G2, tolerance);
    assertEquals(-5.762876404151583, lG2, tolerance);

    System.out.println("solver_nc_marg runtime in ms: " + 1000*runtime1);
    System.out.println("solver_nc_joint runtime in ms: " + 1000*runtime2);
  }

  @Test
  public void SolverNCTestModel2() {
    Network model = new Network("example_loadDependentModel");

    Delay node1 = new Delay(model, "Delay");
    Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

    ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

    node1.setService(jobclass1, Exp.fitMean(1.000000));
    node2.setService(jobclass1, Exp.fitMean(1.000000));
    Matrix alpha = new Matrix(1, 10);
    alpha.set(0, 1);
    alpha.set(1, 2);
    alpha.set(2, 3);
    alpha.set(3, 4);
    alpha.set(4, 5);
    alpha.set(5, 5);
    alpha.set(6, 5);
    alpha.set(7, 5);
    alpha.set(8, 5);
    alpha.set(9, 5);
    node2.setLoadDependence(alpha);

    RoutingMatrix routingMatrix = new RoutingMatrix(model,
            Arrays.asList(jobclass1),
            Arrays.asList(node1, node2));

    routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node2, node1, 1.000000);

    model.link(routingMatrix);

    NetworkStruct sn = model.getStruct(true);

    SolverOptions options = new SolverOptions(SolverType.NC);

    SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, null);
    Matrix lPr1 = ret1.lPr;
    double G1 = ret1.G;
    double runtime1 = ret1.runtime;

    assertEquals(2, lPr1.numRows);
    assertEquals(1, lPr1.numCols);

    assertEquals(-16.104412563027751, lPr1.get(0), 16.104412563027751*tolerance);
    assertEquals(-16.104412563027751, lPr1.get(1), 16.104412563027751*tolerance);

    assertEquals(2.718281801146385, G1, tolerance);

    SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
    double Pr2 = ret2.Pr;
    double G2 = ret2.G;
    double lG2 = ret2.lG;
    double runtime2 = ret2.runtime;

    assertEquals(1.013777129816491e-07, Pr2, tolerance);
    assertEquals(2.718281801146385, G2, tolerance);
    assertEquals(0.999999989952234, lG2, tolerance);

    System.out.println("solver_nc_marg runtime in ms: " + 1000*runtime1);
    System.out.println("solver_nc_joint runtime in ms: " + 1000*runtime2);
  }

  @Test
  public void SolverNCTestModel3() {
    Network model = new Network("example_complexMultiClassModel");

    Queue node1 = new Queue(model, "Queue1", SchedStrategy.INF);
    Queue node2 = new Queue(model, "Queue2", SchedStrategy.INF);
    Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
    Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);

    ClosedClass jobclass1 = new ClosedClass(model, "Class1", 100, node1, 0);
    ClosedClass jobclass2 = new ClosedClass(model, "Class2", 33, node4, 1);

    node1.setService(jobclass1, new Exp(1.0));
    node1.setService(jobclass2, new Exp(Math.sqrt(1.0)));
    node2.setService(jobclass1, new Exp(2.0));
    node2.setService(jobclass2, new Exp(Math.sqrt(2.0)));
    node3.setService(jobclass1, new Exp(3.0));
    node3.setService(jobclass2, new Exp(Math.sqrt(3.0)));
    node4.setService(jobclass1, new Exp(4.0));
    node4.setService(jobclass2, new Exp(Math.sqrt(4.0)));

    RoutingMatrix routingMatrix = new RoutingMatrix(model,
            Arrays.asList(jobclass1, jobclass2),
            Arrays.asList(node1, node2, node3, node4));

    routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node2, node3, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node3, node4, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 1.000000);

    routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.000000);
    routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.000000);
    routingMatrix.addConnection(jobclass2, jobclass2, node4, node1, 1.000000);

    model.link(routingMatrix);

    NetworkStruct sn = model.getStruct(true);

    SolverOptions options = new SolverOptions(SolverType.NC);

    SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, null);
    Matrix lPr1 = ret1.lPr;
    double G1 = ret1.G;
    double runtime1 = ret1.runtime;

    assertEquals(4, lPr1.numRows);
    assertEquals(1, lPr1.numCols);

    assertEquals(-3.041632215808478e2, lPr1.get(0), 3.041632215808478*tolerance);
    assertEquals(-0.018344541300874e2, lPr1.get(1), tolerance);
    assertEquals(0.0, lPr1.get(2), tolerance);
    assertEquals(-0.471992884544684e2, lPr1.get(3), tolerance);

    assertEquals(6.406198113596468e-36, G1, 6.406198113596468e-36*tolerance);

    SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
    double Pr2 = ret2.Pr;
    double G2 = ret2.G;
    double lG2 = ret2.lG;
    double runtime2 = ret2.runtime;

    assertEquals(1.947180103098101e-133, Pr2, 1.947180103098101e-133*tolerance);
    assertEquals(6.406198113596468e-36, G2, 6.406198113596468e-36*tolerance);
    assertEquals(-81.035797370820802, lG2, 81.035797370820802*tolerance);

    System.out.println("solver_nc_marg runtime in ms: " + 1000*runtime1);
    System.out.println("solver_nc_joint runtime in ms: " + 1000*runtime2);
  }

  @Test
  public void SolverNCTestModel4() {
    Network model = new Network("example_complexMultiClassModel");

    Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
    Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
    Queue node3 = new Queue(model, "Queue3", SchedStrategy.INF);
    Queue node4 = new Queue(model, "Queue4", SchedStrategy.INF);

    ClosedClass jobclass1 = new ClosedClass(model, "Class1", 88, node1, 0);
    ClosedClass jobclass2 = new ClosedClass(model, "Class2", 33, node4, 1);

    node1.setService(jobclass1, new Exp(1.0));
    node1.setService(jobclass2, new Exp(Math.sqrt(1.0)));
    node2.setService(jobclass1, new Exp(2.0));
    node2.setService(jobclass2, new Exp(Math.sqrt(2.0)));
    node3.setService(jobclass1, new Exp(3.0));
    node3.setService(jobclass2, new Exp(Math.sqrt(3.0)));
    node4.setService(jobclass1, new Exp(4.0));
    node4.setService(jobclass2, new Exp(Math.sqrt(4.0)));

    RoutingMatrix routingMatrix = new RoutingMatrix(model,
            Arrays.asList(jobclass1, jobclass2),
            Arrays.asList(node1, node2, node3, node4));

    routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node2, node3, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node3, node4, 1.000000);
    routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 1.000000);

    routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.000000);
    routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.000000);
    routingMatrix.addConnection(jobclass2, jobclass2, node4, node1, 1.000000);

    model.link(routingMatrix);

    NetworkStruct sn = model.getStruct(true);

    SolverOptions options = new SolverOptions(SolverType.NC);

    SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, null);
    Matrix lPr1 = ret1.lPr;
    double G1 = ret1.G;
    double runtime1 = ret1.runtime;

    assertEquals(4, lPr1.numRows);
    assertEquals(1, lPr1.numCols);

    assertEquals(-0.804926658480257e2, lPr1.get(0), 0.804926658480257e2*tolerance);
    assertEquals(-0.008127259006488e2, lPr1.get(1), tolerance);
    assertEquals(0.0, lPr1.get(2), tolerance);
    assertEquals(-1.766646876121396e2, lPr1.get(3), 1.766646876121396e2*tolerance);

    assertEquals(1.984349911810360e+30, G1, 1.984349911810360e+30*tolerance);

    SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
    double Pr2 = ret2.Pr;
    double G2 = ret2.G;
    double lG2 = ret2.lG;
    double runtime2 = ret2.runtime;

    assertEquals(6.756257601443034e-78, Pr2, 6.756257601443034e-78*tolerance);
    assertEquals(1.984349911810360e+30, G2, 1.984349911810360e+30*tolerance);
    assertEquals(69.762844149973148, lG2, 69.762844149973148*tolerance);

    System.out.println("solver_nc_marg runtime in ms: " + 1000*runtime1);
    System.out.println("solver_nc_joint runtime in ms: " + 1000*runtime2);
  }
}
