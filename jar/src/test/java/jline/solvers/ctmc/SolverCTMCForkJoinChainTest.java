/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ctmc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * getGenerator and getStateSpace on a CLOSED fork-join model.
 *
 * <p>Both used to enumerate on the un-augmented struct. A fork firing does not
 * conserve the per-chain population, so the population lattice returns an empty
 * local space for the Join: getGenerator handed back a 0x0 matrix WITHOUT
 * raising, and getStateSpace then ran off the end of an empty space. A caller
 * that trusted the generator received a vacuous chain rather than an error --
 * every loop over its states simply did nothing. See BUGS.md BUG-88.
 *
 * <p>The reference sizes are the ones MATLAB and native Python produce for the
 * same model: 16 states after the vanishing fork-occupied and join-firable
 * states are removed by the stochastic complement.
 */
public class SolverCTMCForkJoinChainTest {

    private static VerboseLevel originalVerboseLevel;

    @BeforeAll
    public static void setUpClass() {
        originalVerboseLevel = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    @AfterAll
    public static void tearDownClass() {
        GlobalConstants.setVerbose(originalVerboseLevel);
    }

    private static Network closedForkJoin(int population) {
        Network model = new Network("fjclosed");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        ClosedClass jobclass = new ClosedClass(model, "class1", population, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue1.setService(jobclass, new Exp(1.0));
        queue2.setService(jobclass, new Exp(1.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(jobclass, jobclass, delay, fork, 1.0);
        routing.set(jobclass, jobclass, fork, queue1, 1.0);
        routing.set(jobclass, jobclass, fork, queue2, 1.0);
        routing.set(jobclass, jobclass, queue1, join, 1.0);
        routing.set(jobclass, jobclass, queue2, join, 1.0);
        routing.set(jobclass, jobclass, join, delay, 1.0);
        model.link(routing);
        return model;
    }

    @Test
    public void closedForkJoinGeneratorIsNotEmpty() {
        Matrix infGen = new SolverCTMC(closedForkJoin(2)).getGenerator().infGen;
        assertEquals(16, infGen.getNumRows(), "closed fork-join generator size");
        assertEquals(16, infGen.getNumCols(), "closed fork-join generator size");
        // A generator's rows sum to zero; an all-zero matrix would pass a size
        // check alone, so the conservation is asserted too.
        for (int i = 0; i < infGen.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < infGen.getNumCols(); j++) {
                rowSum += infGen.get(i, j);
            }
            assertEquals(0.0, rowSum, 1e-6, "generator row " + i + " does not sum to zero");
        }
    }

    @Test
    public void closedForkJoinStateSpaceMatchesTheGenerator() {
        SolverCTMC solver = new SolverCTMC(closedForkJoin(2));
        Matrix space = solver.getStateSpace().stateSpace;
        // The state space labels the generator, so the two must have the same
        // number of rows: enumerating them separately gave 35 against 16.
        assertEquals(16, space.getNumRows(), "closed fork-join state space size");
        assertTrue(space.getNumCols() > 0, "closed fork-join state space has no columns");
    }

    @Test
    public void nonForkJoinModelIsUnaffected() {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "c1");
        source.setArrival(jobclass, new Exp(0.8));
        queue.setService(jobclass, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));
        SolverCTMC solver = new SolverCTMC(model);
        solver.options.cutoff = Matrix.singleton(3);
        assertEquals(4, solver.getGenerator().infGen.getNumRows(), "M/M/1 generator size");
        assertEquals(4, solver.getStateSpace().stateSpace.getNumRows(), "M/M/1 state space size");
    }
}
