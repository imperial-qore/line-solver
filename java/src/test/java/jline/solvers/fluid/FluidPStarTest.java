//// Copyright (c) 2012-2022, Imperial College London
//// All rights reserved.
//
//package jline.solvers.fluid;
//
//import jline.lang.*;
//import jline.lang.constant.SchedStrategy;
//import jline.lang.distributions.Coxian;
//import jline.lang.distributions.Exp;
//import jline.lang.nodes.Delay;
//import jline.lang.nodes.Queue;
//import jline.solvers.fluid.smoothing.PStarSearcher;
//import jline.util.Matrix;
//import org.apache.commons.math3.optim.*;
//
//import static jline.solvers.fluid.FluidTestModels.closed_ex1;
//
//public class FluidPStarTest {
//
//  // NOTE TO USER: the search-based p-norm smoothing implementation has been developed primarily for
//  // research and experimentation purposes. As such, currently, the CMA-ES parameters (such as
//  // the number of evaluations, population size, accuracy threshold etc.) are only editable via the
//  // PStarSearcher class. More information can be found in the PStarSearcher class.
//
//  public static void main(String[] args) {
//
//    // Step 1: Instantiate a model of choice
//    Network model = closed_ex1();
//
//    // Step 2: Determine target (accurate) queue lengths
//    PStarSearcher searcher = new PStarSearcher();
//    Matrix targetQueueLengths = searcher.generateTargetQueueLengths(model);
//
//    // Step 3: Search for sufficiently good set of pStar values
//    PointValuePair pStarValues = searcher.findPStarValues(model, targetQueueLengths);
//
//    // Step 4: Solve model using SolverFluid WITHOUT p-Norm Smoothing
//    SolverFluid solverFluidWithoutSmoothing = new SolverFluid(model);
//    solverFluidWithoutSmoothing.options.stiff = false;
//    solverFluidWithoutSmoothing.runAnalyzer();
//    Matrix QNFluidWithoutSmoothing = solverFluidWithoutSmoothing.result.QN;
//    double errorValue = 0;
//    for (int i = 0; i < model.getNumberOfNodes(); i++) {
//      for (int j = 0; j < model.getNumberOfClasses(); j++) {
//        errorValue +=
//            Math.pow(QNFluidWithoutSmoothing.get(i, j) - targetQueueLengths.get(i, j), 2);
//      }
//    }
//
//    // Step 5: Solve model using SolverFluid WITH p-Norm Smoothing
//    SolverFluid solverFluidWithSmoothing = new SolverFluid(model);
//    solverFluidWithSmoothing.options.stiff = false;
//    for (int i = 0; i < model.getNumberOfNodes(); i++) {
//      solverFluidWithSmoothing.options.config.pstar.add(i, pStarValues.getPoint()[i]);
//    }
//    solverFluidWithSmoothing.runAnalyzer();
//    Matrix QNFluidWithSmoothing = solverFluidWithSmoothing.result.QN;
//
//    // Print Outputs
//    System.out.println("\n%%%%%%%% Final Results %%%%%%%%\n");
//    System.out.println("Results: Target Queue Lengths");
//    for (int i = 0; i < QNFluidWithoutSmoothing.getNumRows(); i++) {
//      for (int j = 0; j < QNFluidWithoutSmoothing.getNumCols(); j++) {
//        System.out.format("Station %d   Class %d: %f\n", i, j, targetQueueLengths.get(i, 0));
//      }
//    }
//
//    System.out.println("\nResults: Q for SolverFluid WITHOUT p-Norm Smoothing");
//    for (int i = 0; i < QNFluidWithoutSmoothing.getNumRows(); i++) {
//      for (int j = 0; j < QNFluidWithoutSmoothing.getNumCols(); j++) {
//        System.out.format(
//            "Station %d   Class %d: %f\n", i, j, QNFluidWithoutSmoothing.get(i, 0));
//      }
//    }
//    System.out.format("Final Error Value: %f\n", errorValue);
//
//    System.out.println("\nResults: Q for SolverFluid WITH p-Norm Smoothing");
//    for (int i = 0; i < QNFluidWithSmoothing.getNumRows(); i++) {
//      for (int j = 0; j < QNFluidWithSmoothing.getNumCols(); j++) {
//        System.out.format("Station %d   Class %d: %f\n", i, j, QNFluidWithSmoothing.get(i, 0));
//      }
//    }
//
//    System.out.format("Final Error Value: %f\n", pStarValues.getValue());
//    for (int i = 0; i < QNFluidWithSmoothing.getNumRows(); i++) {
//      System.out.format("PStar for Station %d: %f\n", i, pStarValues.getPoint()[i]);
//    }
//    System.out.format("\nRuntime for CMA-ES: %f seconds", searcher.runTime / 1000000000.0);
//  }
//
//  public static Network ruuskanenExample2() {
//
//    Network model = new Network("RuuskanenExample2");
//
//    Delay delay = new Delay(model, "Delay");
//    Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
//    Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
//    queue1.setNumberOfServers(4);
//    queue2.setNumberOfServers(8);
//
//    ClosedClass class1 = new ClosedClass(model, "Job1", 50, delay, 0);
//    delay.setService(class1, new Exp(0.6));
//    queue1.setService(class1, Coxian.fitMeanAndSCV(0.5, 0.5));
//    queue2.setService(class1, Coxian.fitMeanAndSCV(1, 10));
//
//    model.link(model.serialRouting(delay, queue1, queue2));
//
//    return model;
//  }
//}
