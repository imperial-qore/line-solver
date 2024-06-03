// Copyright (c) 2012-2024, Imperial College London
// All rights reserved.

package jline.examples;

import jline.lang.*;
import jline.lang.constant.GlobalConstants;
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

import java.util.Arrays;
import java.util.Collections;

/**
 * Examples of models evolving in a random environment
 */
public class EnvironmentModel {

  // For System-Testing and Performance Evaluation
  // Corresponds to example_randomEnvironment_1.m in LINE
  public static SolverEnv ex1() {

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

    SolverOptions options = new SolverOptions(SolverType.Env);
    options.iter_tol = 0.01;
    options.timespan[0] = 0;
    options.verbose = VerboseLevel.STD;

    SolverOptions fluidOptions = new SolverOptions(SolverType.Fluid);
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
  public static SolverEnv ex2() {

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

    SolverOptions options = new SolverOptions(SolverType.Env);
    options.iter_tol = 0.05;
    options.timespan[0] = 0;
    options.verbose = VerboseLevel.STD;

    SolverOptions fluidOptions = new SolverOptions(SolverType.Fluid);
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
  public static SolverEnv ex3() {

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

    SolverOptions options = new SolverOptions(SolverType.Env);
    options.iter_tol = 0.05;
    options.timespan[0] = 0;

    SolverOptions fluidOptions = new SolverOptions(SolverType.Fluid);
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

  // For demonstration of State-Dependent Random Environments
  public static SolverEnv ex4() {

    int E = 2;
    Env envModel = new Env("MyEnv", E);
    String[] envName = {"Stage1", "Stage2"};
    String[] envType = {"FAST", "SLOW"};

    // Model in Stage 1
    Network modelStage1 = new Network("model");
    Queue queue1Stage1 = new Queue(modelStage1, "Queue1", SchedStrategy.PS);
    Queue queue2Stage1 = new Queue(modelStage1, "Queue2", SchedStrategy.PS);
    ClosedClass class1Stage1 = new ClosedClass(modelStage1, "Class1", 8, queue1Stage1, 0);
    queue1Stage1.setService(class1Stage1, new Exp(100));
    queue2Stage1.setService(class1Stage1, new Exp(10));
    modelStage1.link(modelStage1.serialRouting(queue1Stage1, queue2Stage1));
    envModel.addStage(0, envName[0], envType[0], modelStage1);

    // Model in Stage 2
    Network modelStage2 = new Network("model");
    Queue queue1Stage2 = new Queue(modelStage2, "Queue1", SchedStrategy.FCFS);
    Queue queue2Stage2 = new Queue(modelStage2, "Queue2", SchedStrategy.FCFS);
    ClosedClass class1Stage2 = new ClosedClass(modelStage2, "Class1", 8, queue1Stage2, 0);
    queue1Stage2.setService(class1Stage2, new Exp(1));
    queue2Stage2.setService(class1Stage2, new Exp(10));
    modelStage2.link(modelStage2.serialRouting(queue1Stage2, queue2Stage2));
    envModel.addStage(1, envName[1], envType[1], modelStage2);

    Matrix envRates = new Matrix(2, 2);
    envRates.set(0, 1, 100);
    envRates.set(1, 0, 0.01);

    Env.ResetEnvRatesFunction resetEnvRatesFunction =
        (originalDist, QExit, UExit, TExit) -> {
          double lambda = originalDist.getRate();
          lambda *= UExit.sumRows(0); // Time-averaged utilisation at Queue1
          return new Exp(Math.max(lambda, GlobalConstants.Zero));
        };

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        if (envRates.get(e, h) > 0) {
          envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
        }
      }
    }
    envModel.resetEnvRatesFun[0][1] = resetEnvRatesFunction;

    SolverOptions options = new SolverOptions(SolverType.Env);
    options.iter_tol = 0.01;
    options.timespan[0] = 0;
    options.method = "statedep";

    SolverOptions fluidOptions = new SolverOptions(SolverType.Fluid);
    fluidOptions.method = "matrix";
    fluidOptions.stiff = false;
    fluidOptions.setODEMaxStep(0.1);
    fluidOptions.verbose = VerboseLevel.SILENT;

    NetworkSolver[] solvers = new NetworkSolver[E];
    for (int e = 0; e < E; e++) {
      solvers[e] = new SolverFluid(envModel.getModel(e));
      solvers[e].options = fluidOptions;
    }

    return new SolverEnv(envModel, solvers, options);
  }

  // For demonstration of p-Norm Smoothing combined with SolverEnv
  public static SolverEnv ex5() {

    int E = 2;
    Env envModel = new Env("MyEnv", E);
    String[] envName = {"Stage1", "Stage2"};
    String[] envType = {"SLOW", "FAST"};

    // Model in Stage 1
    Network modelStage1 = new Network("model");
    Delay delayStage1 = new Delay(modelStage1, "Delay");
    Queue queueStage1 = new Queue(modelStage1, "Queue", SchedStrategy.PS);
    ClosedClass class1Stage1 = new ClosedClass(modelStage1, "Class1", 8, delayStage1, 0);
    delayStage1.setService(class1Stage1, new Exp(1));
    queueStage1.setService(class1Stage1, new Exp(8));
    modelStage1.link(modelStage1.serialRouting(delayStage1, queueStage1));
    envModel.addStage(0, envName[0], envType[0], modelStage1);

    // Model in Stage 2
    Network modelStage2 = new Network("model");
    Delay delayStage2 = new Delay(modelStage2, "Delay");
    Queue queueStage2 = new Queue(modelStage2, "Queue", SchedStrategy.PS);
    ClosedClass class1Stage2 = new ClosedClass(modelStage2, "Class1", 8, delayStage2, 0);
    delayStage2.setService(class1Stage2, new Exp(4));
    queueStage2.setService(class1Stage2, new Exp(8));
    modelStage2.link(modelStage2.serialRouting(delayStage2, queueStage2));
    envModel.addStage(1, envName[1], envType[1], modelStage2);

    Matrix envRates = new Matrix(2, 2);
    envRates.set(0, 1, 0.001); // Replace with 1000 to test Figure 6.11 results
    envRates.set(1, 0, 0.001); // Replace with 1000 to test Figure 6.11 results

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        if (envRates.get(e, h) > 0) {
          envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
        }
      }
    }

    SolverOptions options = new SolverOptions(SolverType.Env);
    options.iter_tol = 0.01;
    options.timespan[0] = 0;

    SolverOptions fluidOptions = new SolverOptions(SolverType.Fluid);
    fluidOptions.method = "matrix";
    fluidOptions.stiff = false;
    fluidOptions.setODEMaxStep(0.1);
    fluidOptions.verbose = VerboseLevel.SILENT;

    NetworkSolver[] solvers = new NetworkSolver[E];
    for (int e = 0; e < E; e++) {
      solvers[e] = new SolverFluid(envModel.getModel(e));
      solvers[e].options = fluidOptions;
      // Activating p-Norm Smoothing by providing each solver with an initial set of pStar values.
      // Going forwards, the pStar values are refreshed within SolverEnv's 'analyze' method
//      PStarSearcher searcher = new PStarSearcher();
//      Matrix targetQueueLengths = searcher.generateTargetQueueLengths(solvers[e].model);
//      PointValuePair pStarValues = searcher.findPStarValues(solvers[e].model, targetQueueLengths);
//      for (int i = 0; i < solvers[e].model.getNumberOfNodes(); i++) {
//        solvers[e].options.config.pstar.add(i, pStarValues.getPoint()[i]);
//      }
    }

    return new SolverEnv(envModel, solvers, options);
  }

  private static Network exGenModel(Matrix rate, int N) {

    Network model = new Network("qn1");

    Delay delay = new Delay(model, "Queue1");
    Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);

    ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
    delay.setService(cClass, new Exp(rate.get(0, 0)));
    queue.setService(cClass, new Exp(rate.get(1, 0)));

    RoutingMatrix routingMatrix =
        new RoutingMatrix(model, Collections.singletonList(cClass), Arrays.asList(delay, queue));
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

    SolverEnv envSolver = ex4();
    envSolver.printAvgTable();
    // envSolver.printEnsembleAvgTables();
  }
}
