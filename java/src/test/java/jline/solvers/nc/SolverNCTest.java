package jline.solvers.nc;

import jline.examples.GettingStarted;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.SolverOptions;
import jline.util.Matrix;
import jline.util.SerializableFunction;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static jline.solvers.nc.SolverNC.*;
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

    assertEquals(2, lPr1.getNumRows());
    assertEquals(1, lPr1.getNumCols());

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

    assertEquals(2, lPr1.getNumRows());
    assertEquals(1, lPr1.getNumCols());

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

    assertEquals(4, lPr1.getNumRows());
    assertEquals(1, lPr1.getNumCols());

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

    assertEquals(4, lPr1.getNumRows());
    assertEquals(1, lPr1.getNumCols());

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

  @Test
  public void SolverNCLDTestModel1() {
    int N = 16;
    int c = 2;
    Network model = new Network("example_closedModel");

    Delay node1 = new Delay(model, "Delay");
    Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

    ClosedClass jobclass1 = new ClosedClass(model, "Class1", 16, node1, 0);

    node1.setService(jobclass1, Exp.fitMean(1.000000));
    node2.setService(jobclass1, Exp.fitMean(1.500000));
    Matrix LD = new Matrix(1, N);
    for (int i=0; i<N; i++) {
      LD.set(i, Math.min(i + 1, c));
    }
    node2.setLoadDependence(LD);

    model.link(model.serialRouting(node1, node2));

    NetworkStruct sn = model.getStruct(true);

    SolverOptions options = new SolverOptions(SolverType.NC);
    SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
    options.method = "nr.logit";
    SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
    options.method = "nr.probit";
    SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
    options.method = "rd";
    SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);
    assertEquals(-2.576432645335951e+00, ret0.lG, tolerance);
    assertEquals(-2.350641375431591e+00, ret1.lG, 3.303193953352507e+00*tolerance*10000);
    assertEquals(-2.576432728076320e+00, ret2.lG, tolerance);
    assertEquals(-2.633413431541626e+00, ret3.lG, tolerance);

  }


  @Test
  public void SolverNCLDTestModel2() {
    int N = 4;
    int c = 2;
    Network ldmodel = new Network("ldmodel");
    Delay node1 = new Delay(ldmodel, "Delay");
    Queue node2 = new Queue(ldmodel, "Queue1", SchedStrategy.PS);
    ClosedClass jobclass1 = new ClosedClass(ldmodel, "Class1", N, node1, 0);
    ClosedClass jobclass2 = new ClosedClass(ldmodel, "Class2", N/2, node1, 0);
    node1.setService(jobclass1, Exp.fitMean(1.0));
    node1.setService(jobclass2, Exp.fitMean(2.0));
    node2.setService(jobclass1, Exp.fitMean(1.5));
    node2.setService(jobclass2, Exp.fitMean(2.5));
    Matrix LD = new Matrix(1, N + N/2);
    for (int i=0; i<N+N/2; i++) {
      LD.set(i, Math.min(i + 1, c));
    }
    node2.setLoadDependence(LD);
    RoutingMatrix P = ldmodel.initRoutingMatrix();
    P.set(jobclass1, jobclass1, ldmodel.serialRouting(node1, node2));
    P.set(jobclass2, jobclass2, ldmodel.serialRouting(node1, node2));
    ldmodel.link(P);
    NetworkStruct sn = ldmodel.getStruct(false);
    SolverOptions options = new SolverOptions(SolverType.NC);
    SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
    options.method = "nr.logit";
    SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
    options.method = "nr.probit";
    SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
    options.method = "rd";
    SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);

    assertEquals(4.115402209963468e+00, ret0.lG, tolerance);
    assertEquals(3.303193953352507e+00, ret1.lG, 3.303193953352507e+00*tolerance*10000);
    assertEquals(2.081621924892065e+00, ret2.lG, tolerance);
    assertEquals(3.980843543955730e+00, ret3.lG, tolerance);
  }

  @Test
  public void SolverNCLDTestModel3() {
    int N = 4;
    int c = 3;
    Network ldmodel = new Network("ldmodel");
    Delay node1 = new Delay(ldmodel, "Delay");
    Queue node2 = new Queue(ldmodel, "Queue1", SchedStrategy.PS);
    Queue node3 = new Queue(ldmodel, "Queue2", SchedStrategy.PS);
    ClosedClass jobclass1 = new ClosedClass(ldmodel, "Class1", N, node1, 0);
    ClosedClass jobclass2 = new ClosedClass(ldmodel, "Class2", N/2, node1, 0);
    node1.setService(jobclass1, Exp.fitMean(1.0));
    node1.setService(jobclass2, Exp.fitMean(2.0));

    node2.setService(jobclass1, Exp.fitMean(1.5));
    node2.setService(jobclass2, Exp.fitMean(2.5));
    Matrix LD = new Matrix(1, N+N/2);
    for (int i=0; i<N+N/2; i++) {
      LD.set(i, Math.min(i + 1, c));
    }
    node2.setLoadDependence(LD);

    node3.setService(jobclass1, Exp.fitMean(3.5));
    node3.setService(jobclass2, Exp.fitMean(4.5));
    node3.setLoadDependence(LD.clone());

    RoutingMatrix P = ldmodel.initRoutingMatrix();
    P.set(jobclass1, jobclass1, ldmodel.serialRouting(node1, node2, node3));
    P.set(jobclass2, jobclass2, ldmodel.serialRouting(node1, node2, node3));
    ldmodel.link(P);
    NetworkStruct sn = ldmodel.getStruct(true);
    SolverOptions options = new SolverOptions(SolverType.NC);
    SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
    options.method = "nr.logit";
    SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
    options.method = "nr.probit";
    SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
    options.method = "rd";
    SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);

    assertEquals(8.016437786118589e+00, ret0.lG, tolerance);
    assertEquals(7.305375561817382e+00, ret1.lG, tolerance*7.305375561817382e+00*10000);
    assertEquals(5.979290244081035e+00, ret2.lG, 5.979290244081035e+00*tolerance*100);
    assertEquals(7.755109452853959e+00, ret3.lG, tolerance);
  }

  @Test
  public void SolverNCLDTestModel4() {
    int N = 10;
    int c = 5;
    Network ldmodel = new Network("ldmodel");
    Delay node1 = new Delay(ldmodel, "Delay");
    Queue node2 = new Queue(ldmodel, "Queue1", SchedStrategy.PS);
    ClosedClass jobclass1 = new ClosedClass(ldmodel, "Class1", N, node1, 0);
    node1.setService(jobclass1, Exp.fitMean(1.0));
    node2.setService(jobclass1, Exp.fitMean(1.0));
    Matrix LD = new Matrix(1, N );
    for (int i=0; i<N; i++) {
      LD.set(i, Math.min(i + 1, c));
    }
    node2.setLoadDependence(LD);
    RoutingMatrix P = ldmodel.initRoutingMatrix();
    P.set(jobclass1, jobclass1, ldmodel.serialRouting(node1, node2));
    ldmodel.link(P);
    NetworkStruct sn = ldmodel.getStruct(false);
    SolverOptions options = new SolverOptions(SolverType.NC);
    SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
    options.method = "nr.logit";
    SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
    options.method = "nr.probit";
    SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
    options.method = "rd";
    SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);

    assertEquals(-7.957151694158279e+00, ret0.lG, tolerance);
    assertEquals(-7.731360424253920e+00, ret1.lG, tolerance);
    assertEquals(-7.957151776898646e+00, ret2.lG, tolerance);
    assertEquals(-8.635479492422949e+00, ret3.lG, tolerance);
  }

  @Test
  public void SolverNCLDTestModel5() {
    int N = 16;
    int c = 2;
    Network cdmodel = new Network("cdmodel");
    Delay node1 = new Delay(cdmodel, "Delay");
    Queue node2 = new Queue(cdmodel, "Queue1", SchedStrategy.PS);
    ClosedClass jobclass1 = new ClosedClass(cdmodel, "Class1", N, node1, 0);
    ClosedClass jobclass2 = new ClosedClass(cdmodel, "Class2", N/2, node1, 0);

    node1.setService(jobclass1, Exp.fitMean(1.0));
    node1.setService(jobclass2, Exp.fitMean(2.0));
    node2.setService(jobclass1, Exp.fitMean(1.5));
    node2.setService(jobclass2, Exp.fitMean(2.5));
    SerializableFunction<Matrix, Double> lcd = matrix -> Math.min(matrix.get(0), c);
    node2.setClassDependence(lcd);
    Matrix LD = new Matrix(1, N );
    for (int i=0; i<N; i++) {
      LD.set(i, Math.min(i + 1, c));
    }
    RoutingMatrix P = cdmodel.initRoutingMatrix();
    P.set(jobclass1, jobclass1, cdmodel.serialRouting(node1, node2));
    P.set(jobclass2, jobclass2, cdmodel.serialRouting(node1, node2));
    cdmodel.link(P);
    NetworkStruct sn = cdmodel.getStruct(false);
    SolverOptions options = new SolverOptions(SolverType.NC);
    options.method = "nr.logit";
    SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
    options.method = "nr.probit";
    SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
    options.method = "rd";
    SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);
    assertEquals(2.754023211796183e+01, ret1.lG, tolerance);
    assertEquals(2.554890090416877e+01, ret2.lG, tolerance);
    assertEquals(2.802932951958689e+01, ret3.lG, tolerance);
  }


}
