package jline.lang.state;

import jline.lang.*;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.util.matrix.Matrix;

/**
 * Test model fixtures for State event testing, replicating spn_basic_closed.m
 * and test_SPN_events.m from the MATLAB test suite.
 */
public class PetriStateTestFixtures {

    /**
     * Creates a test SPN model with 3 modes at T1 (all 1 server):
     * - Mode 1: Exponential distribution (mean=1, rate=1)
     * - Mode 2: Erlang distribution (mean=1, order=2, gives 2 phases with rate=2)
     * - Mode 3: HyperExp distribution (mean=1, SCV=4)
     *
     * Used for afterEvent tests (checks 1-5 in test_SPN_events.m).
     * Matches spn_basic_closed.m model structure.
     */
    public static Network createTestModelWithThreeModes() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 1);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1);

        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1.0, 2));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1);

        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 1);
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobClass.getPopulation()));

        return model;
    }

    /**
     * Creates a test SPN model with mode1 having 2 servers.
     * Matches MATLAB: model.reset(); T{1}.setNumberOfServers(mode1,2);
     * - Mode 1: Exp(1), 2 servers
     * - Mode 2: Erlang(1,2), 1 server
     * - Mode 3: HyperExp(1,4), 1 server
     *
     * Used for afterGlobalEvent tests (checks 6-12 in test_SPN_events.m).
     */
    public static Network createTestModelMode1TwoServers() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 2);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1);

        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1.0, 2));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1);

        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 1);
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobClass.getPopulation()));

        return model;
    }

    /**
     * Creates a test SPN model with mode1 and mode3 having 2 servers.
     * Matches MATLAB: model.reset(); T{1}.setNumberOfServers(mode3,2);
     * (after mode1 was already set to 2 servers)
     * - Mode 1: Exp(1), 2 servers
     * - Mode 2: Erlang(1,2), 1 server
     * - Mode 3: HyperExp(1,4), 2 servers
     *
     * Used for afterGlobalEvent check 13 in test_SPN_events.m.
     */
    public static Network createTestModelMode1Mode3TwoServers() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 2);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1);

        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1.0, 2));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1);

        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 2);
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobClass.getPopulation()));

        return model;
    }
}
