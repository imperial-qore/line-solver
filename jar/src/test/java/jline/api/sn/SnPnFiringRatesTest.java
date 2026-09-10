package jline.api.sn;

import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.RoutingMatrix;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

/**
 * Covers the decline paths of the Place firing-rate recovery.
 *
 * <p>The departure equation may only reference TIMED modes, because an
 * immediate firing takes zero time and is never counted in TN. A Place drained
 * only by immediate modes therefore contributes no departure row. When NO Place
 * contributes one the system is homogeneous, the pseudoinverse returns the zero
 * vector, and reporting it would mark every Place idle; the recovery must
 * decline instead so the caller keeps the throughputs it already had.
 */
public class SnPnFiringRatesTest {

    /** Two Places in a cycle, both drained by an IMMEDIATE mode. */
    private static Network allImmediate() {
        Network model = new Network("spn_all_immediate");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        ClosedClass jobclass = new ClosedClass(model, "Class1", 1, P1, 0);

        Mode m1 = T1.addMode("M1");
        T1.setDistribution(m1, new Immediate());
        T1.setTimingStrategy(m1, TimingStrategy.IMMEDIATE);
        T1.setFiringPriorities(m1, 1);
        T1.setFiringWeights(m1, 1.0);
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);

        Mode m2 = T2.addMode("M2");
        T2.setDistribution(m2, new Immediate());
        T2.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
        T2.setFiringPriorities(m2, 1);
        T2.setFiringWeights(m2, 1.0);
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, P1, 1);

        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(jobclass, jobclass, P1, T1, 1.0);
        rm.set(jobclass, jobclass, P2, T2, 1.0);
        rm.set(jobclass, jobclass, T1, P2, 1.0);
        rm.set(jobclass, jobclass, T2, P1, 1.0);
        model.link(rm);

        P1.setState(Matrix.singleton(1));
        P2.setState(Matrix.singleton(0));
        return model;
    }

    /**
     * The same cycle with T1 made TIMED. One departure row now survives, so the
     * recovery proceeds; this is the control that shows the test above fails for
     * the stated reason and not because the fixture is rejected earlier.
     */
    private static Network oneTimed() {
        Network model = new Network("spn_one_timed");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        ClosedClass jobclass = new ClosedClass(model, "Class1", 1, P1, 0);

        Mode m1 = T1.addMode("M1");
        T1.setDistribution(m1, new Exp(3));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);

        Mode m2 = T2.addMode("M2");
        T2.setDistribution(m2, new Immediate());
        T2.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
        T2.setFiringPriorities(m2, 1);
        T2.setFiringWeights(m2, 1.0);
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, P1, 1);

        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(jobclass, jobclass, P1, T1, 1.0);
        rm.set(jobclass, jobclass, P2, T2, 1.0);
        rm.set(jobclass, jobclass, T1, P2, 1.0);
        rm.set(jobclass, jobclass, T2, P1, 1.0);
        model.link(rm);

        P1.setState(Matrix.singleton(1));
        P2.setState(Matrix.singleton(0));
        return model;
    }

    @Test
    public void declinesWhenEveryPlaceIsDrainedOnlyByImmediateModes() {
        NetworkStruct sn = allImmediate().getStruct();
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        // Whatever the analyzer measured is uninformative here: an immediate
        // firing is not a timed event, so these entries carry no firing rate.
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                TN.set(i, k, 0.0);
            }
        }
        SnPnFiringRates.Ret ret = SnPnFiringRates.snPnFiringRates(sn, TN, false);
        assertNull(ret.rates,
                "recovery must decline when no timed mode drains any Place, "
                        + "otherwise the homogeneous system returns all-zero rates");
    }

    @Test
    public void recoversWhenAtLeastOneTimedModeDrainsAPlace() {
        NetworkStruct sn = oneTimed().getStruct();
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                TN.set(i, k, 1.5);
            }
        }
        SnPnFiringRates.Ret ret = SnPnFiringRates.snPnFiringRates(sn, TN, false);
        assertNotNull(ret.rates, "one timed departure row is enough to determine the rates");
        assertEquals(2, ret.rates.getNumRows(), "one rate per (transition, mode) pair");
        // T1 is the only measured mode and drains P1 at 1.5; the balance
        // equations carry that through to the immediate T2.
        assertEquals(1.5, ret.rates.get(0, 0), 1e-9, "timed rate must equal the measured throughput");
        assertEquals(1.5, ret.rates.get(1, 0), 1e-9, "immediate rate is pinned by token balance");
    }
}
