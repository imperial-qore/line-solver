// Copyright (c) 2012-2024, Imperial College London
// All rights reserved.

package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.lang.distributions.*;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.SolverEnv;
import jline.solvers.fluid.SolverFluid;
//import jline.solvers.fluid.smoothing.PStarSearcher;
import jline.util.Maths;
import jline.util.Matrix;

/**
 * Examples of models evolving in a random environment
 */
public class RandomEnvironment {

  // For System-Testing and Performance Evaluation
  // Corresponds to example_randomEnvironment_1.m in LINE
  public static SolverEnv example_randomEnvironment_1() {

    int N = 1;
    int M = 2;
    int E = 2;

    Env envModel = new Env("MyEnv", E);
    String[] envName = {"Stage1", "Stage2"};
    String[] envType = {"UP", "DOWN"};

    Matrix rate = new Matrix(M, E);
    rate.set(0, 0, 2);
    rate.set(0, 1, 1);
    rate.set(1, 0, 1);
    rate.set(1, 1, 2);

    Network[] envSubModel = new Network[E];
    envSubModel[0] = exGenModel(Matrix.extractColumn(rate, 0, null), N);
    envSubModel[1] = exGenModel(Matrix.extractColumn(rate, 1, null), N);

    for (int e = 0; e < E; e++) {
      envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
    }

    Matrix envRates = new Matrix(2, 2);
    envRates.set(0, 0, 0);
    envRates.set(0, 1, 1);
    envRates.set(1, 0, 0.5);
    envRates.set(1, 1, 0.5);

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        if (envRates.get(e, h) > 0) {
          envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
        }
      }
    }

    System.out.println(
        "The metasolver considers an environment with 2 stages and a queueing network with 2 stations.");
    System.out.println(
        "Every time the stage changes, the queueing network will modify the service rates of the stations.\n");

    // envModel.printStageTable();

    SolverOptions options = new SolverOptions(SolverType.ENV);
    options.iter_tol = 0.01;
    options.timespan[0] = 0;
    options.verbose = VerboseLevel.STD;

    SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
    fluidOptions.timespan[1] = 1000;
    fluidOptions.stiff = false;
    fluidOptions.setODEMaxStep(0.25);

    NetworkSolver[] solvers = new NetworkSolver[E];
    for (int e = 0; e < E; e++) {
      solvers[e] = new SolverFluid(envModel.getModel(e));
      solvers[e].options = fluidOptions;
    }

    return new SolverEnv(envModel, solvers, options);
  }

  // For System-Testing and Performance Evaluation
  // Corresponds to example_randomEnvironment_2.m in LINE
  public static SolverEnv example_randomEnvironment_2() {

    int N = 30;
    int M = 3;
    int E = 4;

    Env envModel = new Env("MyEnv", E);
    String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
    String[] envType = {"UP", "DOWN", "FAST", "SLOW"};

    Matrix rate = new Matrix(M, E);
    rate.set(0, 0, 4);
    rate.set(0, 1, 3);
    rate.set(0, 2, 2);
    rate.set(0, 3, 1);
    rate.set(1, 0, 1);
    rate.set(1, 1, 1);
    rate.set(1, 2, 1);
    rate.set(1, 3, 1);
    rate.set(2, 0, 1);
    rate.set(2, 1, 2);
    rate.set(2, 2, 3);
    rate.set(2, 3, 4);

    Network[] envSubModel = new Network[E];
    envSubModel[0] = exGenModel(Matrix.extractColumn(rate, 0, null), N);
    envSubModel[1] = exGenModel(Matrix.extractColumn(rate, 1, null), N);
    envSubModel[2] = exGenModel(Matrix.extractColumn(rate, 2, null), N);
    envSubModel[3] = exGenModel(Matrix.extractColumn(rate, 3, null), N);

    for (int e = 0; e < E; e++) {
      envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
    }

    Matrix envRates = new Matrix(4, 4);
    envRates.set(0, 0, 0);
    envRates.set(0, 1, 0.5);
    envRates.set(0, 2, 0);
    envRates.set(0, 3, 0);
    envRates.set(1, 0, 0);
    envRates.set(1, 1, 0);
    envRates.set(1, 2, 0.5);
    envRates.set(1, 3, 0.5);
    envRates.set(2, 0, 0.5);
    envRates.set(2, 1, 0);
    envRates.set(2, 2, 0);
    envRates.set(2, 3, 0.5);
    envRates.set(3, 0, 0.5);
    envRates.set(3, 1, 0.5);
    envRates.set(3, 2, 0);
    envRates.set(3, 3, 0);

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        if (envRates.get(e, h) > 0) {
          envModel.addTransition(e, h, Coxian.fitMeanAndSCV(1 / envRates.get(e, h), 0.5));
        }
      }
    }

    System.out.println(
        "The metasolver considers an environment with 4 stages and a queueing network with 3 stations.");
    System.out.println(
        "Every time the stage changes, the queueing network will modify the service rates of the stations.\n");

    // envModel.printStageTable();

    SolverOptions options = new SolverOptions(SolverType.ENV);
    options.iter_tol = 0.05;
    options.timespan[0] = 0;
    options.verbose = VerboseLevel.STD;

    SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
    fluidOptions.stiff = false;
    fluidOptions.setODEMaxStep(0.25);
    fluidOptions.verbose = VerboseLevel.SILENT;

    NetworkSolver[] solvers = new NetworkSolver[E];
    for (int e = 0; e < E; e++) {
      solvers[e] = new SolverFluid(envModel.getModel(e));
      solvers[e].options = fluidOptions;
    }

    return new SolverEnv(envModel, solvers, options);
  }

  // For System-Testing and Performance Evaluation
  // Corresponds to example_randomEnvironment_3.m in LINE
  public static SolverEnv example_randomEnvironment_3() {

    int N = 2;
    int M = 2;
    int E = 3;

    Env envModel = new Env("MyEnv", E);
    String[] envName = {"Stage1", "Stage2", "Stage3"};
    String[] envType = {"UP", "DOWN", "FAST"};

    Matrix rate = new Matrix(M, E);
    rate.set(0, 0, 3);
    rate.set(0, 1, 2);
    rate.set(0, 2, 1);
    rate.set(1, 0, 1);
    rate.set(1, 1, 2);
    rate.set(1, 2, 3);

    Network[] envSubModel = new Network[E];
    envSubModel[0] = exGenModel(Matrix.extractColumn(rate, 0, null), N);
    envSubModel[1] = exGenModel(Matrix.extractColumn(rate, 1, null), N);
    envSubModel[2] = exGenModel(Matrix.extractColumn(rate, 2, null), N);

    for (int e = 0; e < E; e++) {
      envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
    }

    Matrix envRates = Maths.circul(3);

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        if (envRates.get(e, h) > 0) {
          envModel.addTransition(e, h, Erlang.fitMeanAndOrder(1 / envRates.get(e, h), e + h + 2));
        }
      }
    }

    System.out.println(
        "The metasolver considers an environment with 3 stages and a queueing network with 2 stations.\n");

    // envModel.printStageTable();

    SolverOptions envOptions = new SolverOptions(SolverType.ENV);
    envOptions.iter_tol = 0.05;
    envOptions.timespan[0] = 0;

    SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
    fluidOptions.stiff = false;
    fluidOptions.setODEMaxStep(0.25);
    fluidOptions.verbose = VerboseLevel.SILENT;

    NetworkSolver[] solvers = new NetworkSolver[E];
    for (int e = 0; e < E; e++) {
      solvers[e] = new SolverFluid(envModel.getModel(e));
      solvers[e].options = fluidOptions;
    }

    return new SolverEnv(envModel, solvers, envOptions);
  }

  private static Network exGenModel(Matrix rate, int N) {

    Network model = new Network("qn1");

    Delay delay = new Delay(model, "Queue1");
    Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);

    ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
    delay.setService(cClass, new Exp(rate.get(0, 0)));
    queue.setService(cClass, new Exp(rate.get(1, 0)));

    RoutingMatrix routingMatrix = model.initRoutingMatrix();
    int numNodes = model.getNumberOfNodes();
    Matrix circulantMatrix = Maths.circul(numNodes);
    for (int row = 0; row < numNodes; row++) {
      for (int col = 0; col < numNodes; col++) {
        if (circulantMatrix.get(row, col) == 1) {
          routingMatrix.set(model.getNodes().get(row), model.getNodes().get(col));
        }
      }
    }
    model.link(routingMatrix);

    return model;
  }


  public static void main(String[] args) {

    //SolverEnv envSolver = EnvModel.ex4();
    //envSolver.printAvgTable();
    // envSolver.printEnsembleAvgTables();
  }
}
