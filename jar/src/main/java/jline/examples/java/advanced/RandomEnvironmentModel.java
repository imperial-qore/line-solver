/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Coxian;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.ENV;
import jline.solvers.fluid.FLD;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Examples of models evolving in a random environment
 */
public class RandomEnvironmentModel {

    /**
     * Helper method to generate queueing network models for environment stages.
     * <p>
     * Features:
     * - Creates a simple two-station closed network
     * - Delay station and PS queue in circular routing
     * - Service rates specified by the rate matrix parameter
     * - Used internally by all environment examples
     *
     * @param rate service rates for the stations
     * @param N    number of jobs in the closed class
     * @return configured queueing network model
     */
    private static Network exGenModel(Matrix rate, int N) {

        Network model = new Network("qn1");

        Delay delay = new Delay(model, "Queue1");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cclass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(cclass, new Exp(rate.value()));
        queue.setService(cclass, new Exp(rate.get(1, 0)));

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

    /**
     * Basic random environment model with 2 stages and 2 stations.
     * <p>
     * Features:
     * - Environment with 2 stages: Stage1 (UP), Stage2 (DOWN)
     * - Each stage has different service rates for the queueing network
     * - Exponential transitions between environment stages
     * - Single closed class with 1 job
     * - Fluid solver for each environment stage
     * - Corresponds to renv_twostages_repairmen.m in LINE
     *
     * @return configured environment model
     */
    public static Environment renv_twostages_repairmen() {

        int N = 1;
        int M = 2;
        int E = 2;

        Environment envModel = new Environment("MyEnv", E);
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

        //System.out.println(
        //        "The metasolver considers an environment with 2 stages and a queueing network with 2 stations.");
        //System.out.println(
        //        "Every time the stage changes, the queueing network will modify the service rates of the stations.\n");

        // envModel.printStageTable();

        return envModel;
    }

    /**
     * Complex random environment model with 4 stages and 3 stations.
     * <p>
     * Features:
     * - Environment with 4 stages: UP, DOWN, FAST, SLOW
     * - Varying service rates across stages and stations
     * - Coxian transitions between environment stages (SCV=0.5)
     * - Single closed class with 30 jobs
     * - Higher iteration tolerance and more complex dynamics
     * - Corresponds to renv_fourstages_repairmen.m in LINE
     *
     * @return configured complex environment model
     */
    public static Environment renv_fourstages_repairmen() {

        int N = 30;
        int M = 3;
        int E = 4;

        Environment envModel = new Environment("MyEnv", E);
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

        //System.out.println(
        //        "The metasolver considers an environment with 4 stages and a queueing network with 3 stations.");
        //System.out.println(
        //        "Every time the stage changes, the queueing network will modify the service rates of the stations.\n");

        // envModel.printStageTable();

        return envModel;
    }

    /**
     * Random environment model with circular transition structure.
     * <p>
     * Features:
     * - Environment with 3 stages: UP, DOWN, FAST
     * - Circular transition pattern between stages
     * - Erlang transitions with varying orders based on stage indices
     * - Single closed class with 2 jobs
     * - Demonstrates circular environment dynamics
     * - Corresponds to renv_threestages_repairmen.m in LINE
     *
     * @return configured environment model
     */
    public static Environment renv_threestages_repairmen() {

        int N = 2;
        int M = 2;
        int E = 3;

        Environment envModel = new Environment("MyEnv", E);
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

        //System.out.println(
        //        "The metasolver considers an environment with 3 stages and a queueing network with 2 stations.\n");

        // envModel.printStageTable();

        return envModel;
    }

    /**
     * Main method for testing and demonstrating random environment examples.
     *
     * <p>Currently contains commented code for running environment solvers
     * and printing average performance tables.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {

        //ENV envSolver = EnvModel.ex4();
        //envSolver.printAvgTable();
        // envSolver.printEnsembleAvgTables();
    }
}
