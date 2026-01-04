/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.solvers.env;

import jline.lang.ClosedClass;
import jline.lang.Env;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Coxian;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.RoutingMatrix;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.SolverFluid;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;

/**
 * Examples of models evolving in a random environment
 */
public class SolverEnvTestFixtures {

    // Corresponds to renv_twostages_repairmen.m in LINE
    public static SolverENV renv_twostages_repairmen() {

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

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";
        fluidOptions.setODEMaxStep(0.25);

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }


    // Corresponds to renv_fourstages_repairmen.m in LINE
    public static SolverENV renv_fourstages_repairmen() {

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
        // envModel.printStageTable();

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 0.05;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    // Corresponds to renv_threestages_repairmen.m in LINE
    public static SolverENV renv_threestages_repairmen() {

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

        SolverOptions envOptions = new SolverOptions(SolverType.ENV);
        envOptions.iter_tol = 0.05;
        envOptions.timespan[0] = 0;
        envOptions.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.stiff = false;
        fluidOptions.setODEMaxStep(0.25);
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, envOptions);
    }

    private static Network exGenModel(Matrix rate, int N) {

        Network model = new Network("qn1");

        Delay delay = new Delay(model, "Queue1");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.FCFS);

        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(cClass, new Exp(rate.value()));
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

    public static SolverENV example_randomEnvironment_4() {
        Matrix sq1SR = new Matrix(new double[]{1, 2});
        Matrix sq2SR = new Matrix(new double[]{4, 5});
        Matrix[] sqSRArray = {sq1SR, sq2SR};
        MatrixCell sqSR = new MatrixCell(sqSRArray) ;

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1.0, -1.0 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});

        SolverENV envSolver = delayAndSwitchingQueues(100, 7, 0.2, sqSR, transition);
        envSolver.setStateDepMethod("stateindep");
        return envSolver;
    }

    public static SolverENV example_randomEnvironment_5() {
        // two classes, each following a different MMPP
        int N = 50;
        int M = 2;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};
        Matrix rate = new Matrix(3, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel_RE5(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel_RE5(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModel_RE5(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModel_RE5(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(2, 2);
        Matrix Q2 = new Matrix(2, 2);
        Q1.set(0, 0, -0.5);
        Q1.set(0, 1, 0.5);
        Q1.set(1, 0, 2);
        Q1.set(1, 1, -2);
        Q2.set(0, 0, -0.2);
        Q2.set(0, 1, 0.2);
        Q2.set(1, 0, 1);
        Q2.set(1, 1, -1);
        Matrix envRates = Q1.krons(Q2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModel_RE5(Matrix rate, int N) {
        Network model = new Network("qn12");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N- N/3, delay, 0);
        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.value()));
        queue.setService(cClass1, new Exp(rate.get(1, 0)));
        queue.setService(cClass2, new Exp(rate.get(2, 0)));

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

    public static SolverENV example_randomEnvironment_6() {
        Matrix sq1SR = new Matrix(new double[]{1, 2});
        Matrix sq2SR = new Matrix(new double[]{4, 5});
        Matrix[] sqSRArray = {sq1SR, sq2SR};
        MatrixCell sqSR = new MatrixCell(sqSRArray) ;

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1.0, -1.0 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});

        SolverENV envSolver = SolverEnvTestFixtures.delayAndSwitchingQueues(100, 4, 0.2, sqSR, transition);
        envSolver.setStateDepMethod("statedep");
        return envSolver;
    }

    public static SolverENV example_randomEnvironment_7() {
        Matrix sq1SR = new Matrix(new double[]{1, 2});
        Matrix sq2SR = new Matrix(new double[]{4, 5});
        Matrix sq3SR = new Matrix(new double[]{3, 6});
        Matrix[] sqSRArray = {sq1SR, sq2SR, sq3SR};
        MatrixCell sqSR = new MatrixCell(sqSRArray) ;

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1.0, -1.0 }});
        Matrix Q3 = new Matrix(new double[][] {
                { -0.25,  0.25 },
                {  0.5, -0.5 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2, Q3});

        SolverENV envSolver = SolverEnvTestFixtures.delayAndSwitchingQueues(100, 4, 0.2, sqSR, transition);
        envSolver.setStateDepMethod("statedep");
        return envSolver;
    }

    public static SolverENV example_randomEnvironment_8() {
        Matrix sq1SR = new Matrix(new double[]{1, 2});
        Matrix sq2SR = new Matrix(new double[]{4, 5});
        Matrix[] sqSRArray = {sq1SR, sq2SR};
        MatrixCell sqSR = new MatrixCell(sqSRArray) ;

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1.0, -1.0 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});

        Matrix routing = new Matrix(new double[][] {
                { 0, 0.5, 0.5 },
                { 1, 0, 0 },
                { 1, 0, 0 }});

        SolverENV envSolver = SolverEnvTestFixtures.delayAndSwitchingQueues(100, 4, 0.2, sqSR, transition, routing);
        envSolver.setStateDepMethod("statedep");
        return envSolver;
    }

    public static SolverENV example_randomEnvironment_9() {

        int N = 100;
        int M = 3;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};

        Matrix rate = new Matrix(M, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel1d2sq(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel1d2sq(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModel1d2sq(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModel1d2sq(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1, -1 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});
        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        Matrix shape1 = new Matrix(new double[][] {
                { -4,  4 },
                {  1, -1 }});

        Matrix shape2 = new Matrix(new double[][] {
                { -25,  25 },
                {  4, -4 }});
        Matrix Shape = shape1.krons(shape2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0 && Shape.get(e, h) > 0) {
                    int k = (int) Shape.get(e, h);
                    envModel.addTransition(e, h, new Erlang(k * envRates.get(e, h), k));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-8;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "matrix";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    public static SolverENV renv_twostages_repairmen0() {
        Matrix sq1SR = new Matrix(new double[]{1, 2});
        Matrix sq2SR = new Matrix(new double[]{4, 5});
        Matrix[] sqSRArray = {sq1SR, sq2SR};
        MatrixCell sqSR = new MatrixCell(sqSRArray) ;

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1.0, -1.0 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});

        SolverENV envSolver = delayAndSwitchingQueues(100, 7, 0.2, sqSR, transition);
        envSolver.setStateDepMethod("stateindep");
        return envSolver;
    }

    /**
     * Builds a SolverENV for a closed network with one delay station and multiple switching queues.
     *
     *
     * @param N           the total population of customers circulating in the closed network
     * @param s           the number of servers allocated to each switching queue (finite-server queues)
     * @param delaySR     the service rate of the infinite-server delay station
     * @param sqSR        a MatrixCell where each entry  is a column vector of service rates for the
     *                    i-th switching queue across its phases
     * @param transition  a MatrixCell of infinitesimal-generator matrices representing phase-transition
     *                    rates for each switching queue; the overall environment transition generator
     *                    is obtained by kronecker-summing these matrices
     * @return            a SolverENV instance configured with the multi-stage environment and per-stage
     *                    SolverFluid solvers
     */
    public static SolverENV delayAndSwitchingQueues(int N, int s, double delaySR, MatrixCell sqSR, MatrixCell transition) {
        int numQueues = sqSR.size();
        int M = numQueues + 1;

        int[] stageCounts = new int[numQueues];
        for (int i = 0; i < numQueues; i++) {
            stageCounts[i] = sqSR.get(i).getNumElements();
        }

        int E = 1;
        for (int c : stageCounts) E *= c;

        Matrix index = generateCombinations(stageCounts, E);

        Env envModel = new Env("MyEnv", E);

        String[] envName = new String[E];
        String[] envType = new String[E];
        for (int i = 0; i < E; i++) {
            envName[i] = "Stage" + (i+1);
            envType[i] = "Type" + (i+1);
        }

        Matrix rate = new Matrix(M, E);
        for (int e = 0; e < E; e++) {
            rate.set(0, e, delaySR);
            for (int i = 0; i < numQueues; i++) {
                rate.set(i + 1, e, sqSR.get(i).get((int) index.get(e, i) - 1,0));
            }
        }

        Network[] envSubModel = new Network[E];
        for (int e = 0; e < E; e++) {
            envSubModel[e] = envGenD1SQN(Matrix.extractColumn(rate, e, null), N, s);
        }

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-6;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = true;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method ="closing";


        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Matrix generateCombinations(int[] stageCounts, int E) {
        int sqNum = stageCounts.length;
        Matrix combinations = new Matrix(sqNum, E);
        List<int[]> combList = new ArrayList<>();

        int[] currentCombination = new int[stageCounts.length];

        generateRecursive(0, stageCounts, currentCombination, combList);

        Matrix result = new Matrix(combList.size(), stageCounts.length);
        for (int i = 0; i < combList.size(); i++) {
            for (int j = 0; j < stageCounts.length; j++) {
                result.set(i, j, combList.get(i)[j] + 1);
            }
        }
        return result;
    }

    private static void generateRecursive(int index, int[] stageCounts, int[] current, List<int[]> output) {
        if (index == stageCounts.length) {
            output.add(current.clone());
        } else {
            for (int i = 0; i < stageCounts[index]; i++) {
                current[index] = i;
                generateRecursive(index + 1, stageCounts, current, output);
            }
        }
    }



    private static Network envGenD1SQN(Matrix rate, int N, int s) {
        Network model = new Network("qn2");
        Delay delay = new Delay(model, "Delay");
        List<Queue> queues = new ArrayList<>();
        for (int i = 1; i < rate.getNumRows(); i++) {
            Queue queue = new Queue(model, "Queue" + i, SchedStrategy.FCFS);
            queue.setNumberOfServers(s);
            queues.add(queue);
        }
        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);


        delay.setService(cClass, new Exp(rate.get(0, 0)));
        for (int i = 0; i < queues.size(); i++) {
            queues.get(i).setService(cClass, new Exp(rate.get(i + 1, 0)));
        }

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

    public static SolverENV delayAndSwitchingQueues(int N, int s, double delaySR, MatrixCell sqSR, MatrixCell transition, Matrix circulantMatrix) {
        int numQueues = sqSR.size();
        int M = numQueues + 1;

        int[] stageCounts = new int[numQueues];
        for (int i = 0; i < numQueues; i++) {
            stageCounts[i] = sqSR.get(i).getNumElements();
        }

        int E = 1;
        for (int c : stageCounts) E *= c;

        Matrix index = generateCombinations(stageCounts, E);

        Env envModel = new Env("MyEnv", E);

        String[] envName = new String[E];
        String[] envType = new String[E];
        for (int i = 0; i < E; i++) {
            envName[i] = "Stage" + (i+1);
            envType[i] = "Type" + (i+1);
        }

        Matrix rate = new Matrix(M, E);
        for (int e = 0; e < E; e++) {
            rate.set(0, e, delaySR);
            for (int i = 0; i < numQueues; i++) {
                rate.set(i + 1, e, sqSR.get(i).get((int) index.get(e, i) - 1,0));
            }
        }

        if (circulantMatrix.getNumRows() != M || circulantMatrix.getNumCols() != M) {
            throw new IllegalArgumentException("Circulant matrix must be square with size equal to stations.");
        }
        Matrix circulantSum = circulantMatrix.sumRows();
        for (int i = 0; i < circulantSum.getNumRows(); i++) {
            if (Math.abs(circulantSum.get(i, 0) - 1) > 1E-6) {
                throw new IllegalArgumentException("Circulant matrix must sum to 1 for each row.");
            }
        }

        Network[] envSubModel = new Network[E];
        for (int e = 0; e < E; e++) {
            envSubModel[e] = envGenD1SQN(Matrix.extractColumn(rate, e, null), N, s, circulantMatrix);
        }

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-6;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network envGenD1SQN(Matrix rate, int N, int s, Matrix circulantMatrix) {
        Network model = new Network("qn2");
        Delay delay = new Delay(model, "Delay");
        List<Queue> queues = new ArrayList<>();
        for (int i = 1; i < rate.getNumRows(); i++) {
            Queue queue = new Queue(model, "Queue" + i, SchedStrategy.FCFS);
            queue.setNumberOfServers(s);
            queues.add(queue);
        }
        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);


        delay.setService(cClass, new Exp(rate.get(0, 0)));
        for (int i = 0; i < queues.size(); i++) {
            queues.get(i).setService(cClass, new Exp(rate.get(i + 1, 0)));
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        int numNodes = model.getNumberOfNodes();
        for (int row = 0; row < numNodes; row++) {
            for (int col = 0; col < numNodes; col++) {
                if (circulantMatrix.get(row, col) > 0 ) {
                    routingMatrix.set(model.getNodes().get(row), model.getNodes().get(col), circulantMatrix.get(row, col));
                }
            }
        }
        model.link(routingMatrix);
        return model;
    }

    private static Network exGenModel1d2sq(Matrix rate, int N) {

        Network model = new Network("qn1");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(cClass, new Exp(rate.value()));
        queue.setService(cClass, new Exp(rate.get(1, 0)));
        queue2.setService(cClass, new Exp(rate.get(2, 0)));

        queue.setNumberOfServers(7);
        queue2.setNumberOfServers(7);

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


    public static SolverENV renv_twostages_repairmen1() {
        // MMPP is placed at the delay station, initial test for CTMC multiclass solver
        int N = 10;
        int M = 1;
        int E = 2;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2"};
        String[] envType = {"UP", "DOWN"};

        Matrix rate = new Matrix(2, E);
        rate.set(0, 0, 2);
        rate.set(0, 1, 0.5);
        rate.set(1, 0, 1);
        rate.set(1, 1, 4);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel_RE11(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel_RE11(Matrix.extractColumn(rate, 1, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix envRates = new Matrix(2, 2);
        envRates.set(0, 0, 0);
        envRates.set(0, 1, 1);
        envRates.set(1, 0, 0.5);
        envRates.set(1, 1, 0);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "statedep";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModel_RE11(Matrix rate, int N) {
        Network model = new Network("qn11");

        Delay delay = new Delay(model, "Queue1");

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N- N/3, delay, 0);
        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.get(1, 0)));

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

    public static SolverENV renv_twostages_repairmen2() {
        // two classes, one follows a MMPP, the other follows an exponential distribution
        int N = 50;
        int M = 2;
        int E = 2;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2"};
        String[] envType = {"UP", "DOWN"};
        Matrix rate = new Matrix(3, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(1, 0, 2);
        rate.set(1, 1, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel_RE12(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel_RE12(Matrix.extractColumn(rate, 1, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix envRates = new Matrix(2, 2);
        envRates.set(0, 0, 0);
        envRates.set(0, 1, 2);
        envRates.set(1, 0, 0.5);
        envRates.set(1, 1, 0);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "statedep";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModel_RE12(Matrix rate, int N) {
        Network model = new Network("qn12");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N- N/3, delay, 0);
        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.value()));
        queue.setService(cClass1, new Exp(rate.get(1, 0)));
        queue.setService(cClass2, new Exp(rate.get(2, 0)));


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

    public static SolverENV renv_twostages_repairmen3() {
        // two classes, each following a different MMPP at two queues
        int N = 12;
        int M = 3;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};
        Matrix rate = new Matrix(3, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel_RE13(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel_RE13(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModel_RE13(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModel_RE13(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(2, 2);
        Matrix Q2 = new Matrix(2, 2);
        Q1.set(0, 0, -0.5);
        Q1.set(0, 1, 0.5);
        Q1.set(1, 0, 2);
        Q1.set(1, 1, -2);
        Q2.set(0, 0, -0.2);
        Q2.set(0, 1, 0.2);
        Q2.set(1, 0, 1);
        Q2.set(1, 1, -1);
        Matrix envRates = Q1.krons(Q2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModel_RE13(Matrix rate, int N) {
        Network model = new Network("qn13");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N- N/3, delay, 0);
        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.value()));
        queue.setService(cClass1, new Exp(rate.get(1, 0)));
        queue.setService(cClass2, new Exp(rate.get(1, 0)));
        queue2.setService(cClass1, new Exp(rate.get(2, 0)));
        queue2.setService(cClass2, new Exp(rate.get(2, 0)));



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

    public static SolverENV renv_twostages_repairmen4() {
        int N = 50;
        int M = 3;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};

        Matrix rate = new Matrix(M, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModelRE14(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModelRE14(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModelRE14(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModelRE14(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1, -1 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});
        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        Matrix shape1 = new Matrix(new double[][] {
                { -4,  4 },
                {  1, -1 }});

        Matrix shape2 = new Matrix(new double[][] {
                { -25,  25 },
                {  4, -4 }});
        Matrix Shape = shape1.krons(shape2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0 && Shape.get(e, h) > 0) {
                    int k = (int) Shape.get(e, h);
                    envModel.addTransition(e, h, new Erlang(k * envRates.get(e, h), k));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-6;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "matrix";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModelRE14(Matrix rate, int N) {

        Network model = new Network("qn14");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(cClass, new Exp(rate.value()));
        queue.setService(cClass, new Exp(rate.get(1, 0)));
        queue2.setService(cClass, new Exp(rate.get(2, 0)));

        queue.setNumberOfServers(4);
        queue2.setNumberOfServers(4);

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


    public static SolverENV renv_twostages_repairmen5() {
        int N = 70;
        int M = 2;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};

        Matrix rate = new Matrix(3, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModelRE15(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModelRE15(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModelRE15(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModelRE15(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1, -1 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});
        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        Matrix shape1 = new Matrix(new double[][] {
                { -4,  4 },
                {  1, -1 }});

        Matrix shape2 = new Matrix(new double[][] {
                { -25,  25 },
                {  4, -4 }});
        Matrix Shape = shape1.krons(shape2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0 && Shape.get(e, h) > 0) {
                    int k = (int) Shape.get(e, h);
                    envModel.addTransition(e, h, new Erlang(k * envRates.get(e, h), k));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-6;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "matrix";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModelRE15(Matrix rate, int N) {
        Network model = new Network("qn15");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N - N/3, delay, 0);

        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.value()));
        queue.setService(cClass1, new Exp(rate.get(1, 0)));
        queue.setService(cClass2, new Exp(rate.get(2, 0)));

        queue.setNumberOfServers(3);

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


    public static SolverENV renv_twostages_repairmen6() {
        int N = 50;
        int M = 3;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};

        Matrix rate = new Matrix(M, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 1);
        rate.set(1, 1, 1);
        rate.set(1, 2, 2);
        rate.set(1, 3, 2);
        rate.set(2, 0, 4);
        rate.set(2, 1, 5);
        rate.set(2, 2, 4);
        rate.set(2, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModelRE16(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModelRE16(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModelRE16(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModelRE16(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(new double[][] {
                { -0.5,  0.5 },
                {  2.0, -2.0 }});
        Matrix Q2 = new Matrix(new double[][] {
                { -0.2,  0.2 },
                {  1, -1 }});
        MatrixCell transition = new MatrixCell(new Matrix[]{Q1, Q2});
        Matrix envRates = transition.get(0);
        for (int i = 1; i < transition.size(); i++) {
            envRates = envRates.krons(transition.get(i));
        }

        Matrix shape1 = new Matrix(new double[][] {
                { -4,  4 },
                {  1, -1 }});

        Matrix shape2 = new Matrix(new double[][] {
                { -25,  25 },
                {  4, -4 }});
        Matrix Shape = shape1.krons(shape2);

        for(int e = 0; e < E; e++){
            for(int h = 0; h < E; h++){
                double lambda = envRates.get(e,h);
                double k  = Shape.get(e,h);
                if( lambda > 0 && k > 0){

                    double mean = k / lambda;
                    double scv  = 1.0/k;
                    Coxian cox = Coxian.fitMeanAndSCV(mean, mean * mean * scv, 2.0);
                    envModel.addTransition(e, h, cox);
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-6;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;
        options.iter_max = 1000;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "matrix";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModelRE16(Matrix rate, int N) {

        Network model = new Network("qn16");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cClass = new ClosedClass(model, "Class1", N, delay, 0);
        delay.setService(cClass, new Exp(rate.value()));
        queue.setService(cClass, new Exp(rate.get(1, 0)));
        queue2.setService(cClass, new Exp(rate.get(2, 0)));

        queue.setNumberOfServers(4);
        queue2.setNumberOfServers(4);

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


    public static SolverENV renv_twostages_repairmen7() {
        // two classes, each following a different MMPP at two queues
        int N = 12;
        int M = 3;
        int E = 4;

        Env envModel = new Env("MyEnv", E);
        String[] envName = {"Stage1", "Stage2", "Stage3", "Stage4"};
        String[] envType = {"UP", "DOWN", "FAST", "SLOW"};
        Matrix rate = new Matrix(7, E);
        rate.set(0, 0, 0.2);
        rate.set(0, 1, 0.2);
        rate.set(0, 2, 0.2);
        rate.set(0, 3, 0.2);
        rate.set(1, 0, 0.5);
        rate.set(1, 1, 0.5);
        rate.set(1, 2, 0.5);
        rate.set(1, 3, 0.5);
        rate.set(2, 0, 1);
        rate.set(2, 1, 1);
        rate.set(2, 2, 1);
        rate.set(2, 3, 1);
        rate.set(3, 0, 1);
        rate.set(3, 1, 1);
        rate.set(3, 2, 2);
        rate.set(3, 3, 2);
        rate.set(4, 0, 3);
        rate.set(4, 1, 3);
        rate.set(4, 2, 3);
        rate.set(4, 3, 3);
        rate.set(5, 0, 2);
        rate.set(5, 1, 2);
        rate.set(5, 2, 2);
        rate.set(5, 3, 2);
        rate.set(6, 0, 4);
        rate.set(6, 1, 5);
        rate.set(6, 2, 4);
        rate.set(6, 3, 5);

        Network[] envSubModel = new Network[E];
        envSubModel[0] = exGenModel_RE17(Matrix.extractColumn(rate, 0, null), N);
        envSubModel[1] = exGenModel_RE17(Matrix.extractColumn(rate, 1, null), N);
        envSubModel[2] = exGenModel_RE17(Matrix.extractColumn(rate, 2, null), N);
        envSubModel[3] = exGenModel_RE17(Matrix.extractColumn(rate, 3, null), N);

        for (int e = 0; e < E; e++) {
            envModel.addStage(e, envName[e], envType[e], envSubModel[e]);
        }

        Matrix Q1 = new Matrix(2, 2);
        Matrix Q2 = new Matrix(2, 2);
        Q1.set(0, 0, -0.5);
        Q1.set(0, 1, 0.5);
        Q1.set(1, 0, 2);
        Q1.set(1, 1, -2);
        Q2.set(0, 0, -0.2);
        Q2.set(0, 1, 0.2);
        Q2.set(1, 0, 1);
        Q2.set(1, 1, -1);
        Matrix envRates = Q1.krons(Q2);

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (envRates.get(e, h) > 0) {
                    envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
                }
            }
        }

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 1E-5;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.SILENT;

        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.verbose = VerboseLevel.SILENT;
        fluidOptions.method = "closing";

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverFluid(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }

        return new SolverENV(envModel, solvers, options);
    }

    private static Network exGenModel_RE17(Matrix rate, int N) {
        Network model = new Network("qn13");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass cClass1 = new ClosedClass(model, "Class1", N/3, delay, 0);
        ClosedClass cClass2 = new ClosedClass(model, "Class2", N/2, delay, 0);
        ClosedClass cClass3 = new ClosedClass(model, "Class3", N - N/3 - N/2, delay, 0);

        // one class: Exp+Exp, second class: MMPP+Exp, third class: Exp+MMPP
        delay.setService(cClass1, new Exp(rate.value()));
        delay.setService(cClass2, new Exp(rate.value()));
        delay.setService(cClass3, new Exp(rate.value()));
        queue.setService(cClass1, new Exp(rate.get(1, 0)));
        queue.setService(cClass2, new Exp(rate.get(2, 0)));
        queue.setService(cClass3, new Exp(rate.get(3, 0)));
        queue2.setService(cClass1, new Exp(rate.get(4, 0)));
        queue2.setService(cClass2, new Exp(rate.get(5, 0)));
        queue2.setService(cClass3, new Exp(rate.get(6, 0)));

//        queue.setNumberOfServers(2);
//        queue2.setNumberOfServers(2);

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
}




