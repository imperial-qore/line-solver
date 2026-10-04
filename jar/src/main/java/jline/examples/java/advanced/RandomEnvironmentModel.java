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
import jline.lang.processes.MMPP2;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.ENV;
import jline.solvers.fluid.FLD;
import jline.solvers.mva.MVA;
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
     * `renv_genqn.m`: a Delay and a PS Queue in a cycle, one closed class of N
     * jobs. The reference exposes it as a FUNCTION rather than a script, and it
     * is the stage generator the repairmen environments above are built from,
     * so it is published here under its own name as well.
     *
     * @param rateDelay service rate of the Delay station
     * @param rateQueue service rate of the PS Queue
     * @param N         closed population
     * @return the configured two-station closed network
     */
    public static Network renv_genqn(double rateDelay, double rateQueue, int N) {
        Matrix rate = new Matrix(2, 1);
        rate.set(0, 0, rateDelay);
        rate.set(1, 0, rateQueue);
        return exGenModel(rate, N);
    }

    /**
     * `renv_map_fallback.m` / `example_mapqn2renv.m`: Think -> FCFS queue with
     * MMPP(2) service, one closed class of five.
     *
     * <p>The MMPP is the whole point of the model: MVA and NC have no native
     * non-renewal service, so the solver replaces the modulating chain by a
     * random environment and solves the stages instead.
     *
     * @param name the model name the reference gives it
     * @return the configured closed network with MMPP(2) service
     */
    public static Network mmppClosed(String name) {
        Network model = new Network(name);
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C1", 5, think, 0);
        think.setService(cclass, new Exp(1.0));
        // Slow phase 1, fast phase 10; phase switch rates 0.2 and 0.3.
        queue.setService(cclass, new MMPP2(1.0, 10.0, 0.2, 0.3));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /**
     * `renv_container_terminal.m`: one hour-stage network, N carriers cycling
     * yard-Delay -> quay-crane Queue.
     *
     * @param stageRate per-token yard completion rate of this hour
     * @param craneRate moves per hour served by one crane
     * @param nCranes   quay cranes (servers of the FCFS station)
     * @param N         straddle carriers circulating (the closed tokens)
     * @return the configured hour-stage network
     */
    public static Network terminalModel(double stageRate, double craneRate, int nCranes, int N) {
        Network qn = new Network("Terminal");
        Delay yard = new Delay(qn, "Yard");
        Queue cranes = new Queue(qn, "QuayCranes", SchedStrategy.FCFS);
        cranes.setNumberOfServers(nCranes);
        ClosedClass containers = new ClosedClass(qn, "Containers", N, yard, 0);
        yard.setService(containers, new Exp(stageRate));
        cranes.setService(containers, new Exp(craneRate));
        qn.link(Network.serialRouting(yard, cranes));
        return qn;
    }

    /**
     * `renv_container_terminal.m`: the same hour-stage with the internal
     * terminal handling collapsed into one closed load-dependent FES, serving
     * at mu(n) when n containers are inside it.
     *
     * @param stageRate per-token yard completion rate of this hour
     * @param fesRate   the Norton rate curve mu(1..N); its length is the population
     * @return the configured hour-stage network with a load-dependent FES
     */
    public static Network terminalFESModel(double stageRate, Matrix fesRate) {
        int N = (int) fesRate.length();
        Network qn = new Network("TerminalFES");
        Delay yard = new Delay(qn, "Yard");
        Queue fesq = new Queue(qn, "TerminalFES", SchedStrategy.PS);
        ClosedClass containers = new ClosedClass(qn, "Containers", N, yard, 0);
        yard.setService(containers, new Exp(stageRate));
        fesq.setService(containers, new Exp(fesRate.get(0)));   // the base rate mu(1)
        Matrix alpha = new Matrix(1, N);
        for (int i = 0; i < N; i++) {
            alpha.set(0, i, fesRate.get(i) / fesRate.get(0));   // alpha(n) = mu(n)/mu(1)
        }
        fesq.setLoadDependence(alpha);
        qn.link(Network.serialRouting(yard, fesq));
        return qn;
    }

    /**
     * `renv_container_terminal.m`: Norton flow-equivalent rates, the throughput
     * of the isolated internal subnetwork (quay cranes -> stacking cranes, both
     * processor sharing) at n = 1..N containers, with the rest of the terminal
     * short-circuited by a near-instantaneous delay.
     *
     * @param muQuay  quay-crane service rate
     * @param muStack stacking-crane service rate
     * @param N       largest population the curve is needed at
     * @return the row vector mu(1..N)
     */
    public static Matrix fesRateCurve(double muQuay, double muStack, int N) {
        Matrix mu = new Matrix(1, N);
        for (int n = 1; n <= N; n++) {
            Network sub = new Network("TerminalInternals");
            Delay ref = new Delay(sub, "ShortCircuit");
            Queue quay = new Queue(sub, "QuayCranes", SchedStrategy.PS);
            Queue stack = new Queue(sub, "StackingCranes", SchedStrategy.PS);
            ClosedClass cls = new ClosedClass(sub, "Containers", n, ref, 0);
            ref.setService(cls, new Exp(1e6));               // ~instantaneous short-circuit
            quay.setService(cls, new Exp(muQuay));
            stack.setService(cls, new Exp(muStack));
            sub.link(Network.serialRouting(ref, quay, stack));
            SolverOptions mvaOpt = new SolverOptions(SolverType.MVA);
            mvaOpt.method = "exact";
            mvaOpt.verbose = VerboseLevel.SILENT;
            // The short-circuit's throughput IS the subnetwork completion rate
            // at this occupancy: every token passes it once per cycle.
            mu.set(0, n - 1, new MVA(sub, mvaOpt).getAvgTput().get(0, 0));
        }
        return mu;
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
