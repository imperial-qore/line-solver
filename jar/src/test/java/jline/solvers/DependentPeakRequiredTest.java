package jline.solvers;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import jline.util.SerializableFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The peak rate that normalizes utilization is a MODEL INPUT, and a model that
 * omits it is refused rather than guessed at.
 *
 * <p>Utilization at a rate-dependent station is reported as U = T*E[S]/peak,
 * the same fraction-of-capacity the T*S/c of an ordinary multiserver station
 * is. For load dependence the peak is max(c, max alpha) and the user has
 * already supplied every alpha, so it is always known. For CLASS and JOINT
 * dependence the scaling is a HANDLE: recovering max_n beta(n) would mean
 * sweeping the population lattice, which needs a bound the handle does not
 * carry -- and an open class has no bound at all. So the peak must be
 * declared, and MATLAB's setLimitedClassDependence and
 * getLimitedClassDependencePeak.m both error when it is not.
 *
 * <p>The JAR used to be the outlier: the one-argument setters accepted the
 * declaration silently, and getLimitedClassDependencePeak then SWEPT a peak out
 * of the handle over getNumberOfJobs() -- which is 0 for every open class, so
 * the "peak" of an unbounded beta was beta(0). That is a number that reads as a
 * utilization and is not one, which is worse than an error.
 */
public class DependentPeakRequiredTest {

    private static final SerializableFunction<Matrix, Matrix> BETA =
            (Matrix n) -> {
                Matrix out = new Matrix(1, 1);
                out.set(0, 0, Math.min(n.get(0), 2.0));
                return out;
            };

    /** Closed Delay -> Queue with one class of N jobs; nothing dependent yet. */
    private static Network closedModel() {
        Network model = new Network("cdpeak");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "Class1", 3, delay, 0);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(1.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Matrix peak(double v) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, v);
        return m;
    }

    @Test
    public void testQueueSetClassDependenceWithoutPeakIsRefused() {
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        RuntimeException e =
                assertThrows(RuntimeException.class, () -> queue.setClassDependence(BETA));
        assertTrue(e.getMessage().contains("peakRatePerClass"),
                "the refusal must name the remedy, got: " + e.getMessage());
    }

    @Test
    public void testStationSetLimitedClassDependenceWithoutPeakIsRefused() {
        // the "don't expose" setter beneath setClassDependence, which the JSON
        // reader and SolverLN both call, and which used to accept silently
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        assertThrows(RuntimeException.class, () -> queue.setLimitedClassDependence(BETA));
    }

    @Test
    public void testStationSetLimitedJointDependenceWithoutPeakIsRefused() {
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        assertThrows(RuntimeException.class, () -> queue.setLimitedJointDependence(BETA));
    }

    @Test
    public void testNonPositivePeakIsRefused() {
        // a peak of zero would divide, and a negative one would flip the sign;
        // MATLAB rejects both in the same setter
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        assertThrows(RuntimeException.class, () -> queue.setClassDependence(BETA, peak(0.0)));
        assertThrows(RuntimeException.class, () -> queue.setClassDependence(BETA, peak(-2.0)));
        assertThrows(RuntimeException.class, () -> queue.setJointDependence(BETA, peak(0.0)));
    }

    @Test
    public void testDeclaredPeakIsAccepted() {
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        assertDoesNotThrow(() -> queue.setClassDependence(BETA, peak(2.0)));
        assertEquals(2.0, model.getLimitedClassDependencePeak().get(queue).get(0), 1e-12,
                "the declared peak must reach the struct unchanged");
    }

    @Test
    public void testPeakIsNotSweptOutOfTheHandle() {
        // beta maxes at 2 over the lattice but is DECLARED as 5, and the struct
        // must carry 5: the declaration is the normalizer, not the sweep. A
        // solver that swept instead would silently report 2.5x this station's
        // utilization.
        Network model = closedModel();
        Queue queue = (Queue) model.getStations().get(1);
        queue.setClassDependence(BETA, peak(5.0));
        assertEquals(5.0, model.getLimitedClassDependencePeak().get(queue).get(0), 1e-12);
    }

    @Test
    public void testOpenClassHandleHasNoLatticeToSweep() {
        // THE CASE THAT MOTIVATES THE ERROR. beta(n) = 1 + n is unbounded, and
        // getNumberOfJobs() is infinite for an open class, so there is no
        // lattice: the old sweep evaluated beta at n = 0 and called 1 the peak.
        Network model = new Network("openCd");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.5));
        queue.setService(oclass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));

        SerializableFunction<Matrix, Matrix> growing = (Matrix n) -> {
            Matrix out = new Matrix(1, 1);
            out.set(0, 0, 1.0 + n.get(0));
            return out;
        };
        assertThrows(RuntimeException.class, () -> queue.setClassDependence(growing));
        queue.setClassDependence(growing, peak(11.0));
        assertEquals(11.0, model.getLimitedClassDependencePeak().get(queue).get(0), 1e-12);
    }
}
