package jline.solvers.ctmc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * The CTMC transient path taken by fast adaptive uniformization,
 * {@code options.config.transient_method = "fau"}.
 *
 * <p>The oracle is the ODE trajectory the same analyzer produces on the same
 * generator, sampled on the same grid: the two integrate the same forward equation
 * and must agree to the integrator's own tolerance. Tightening {@code fau_epsilon}
 * by six orders must then move the FAU answer by far less than that gap, which is
 * what shows the residual is the integrator's and not the uniformization's.
 *
 * <p>An M/M/1 truncated at six jobs is used rather than a closed model because the
 * truncation makes the state space small while leaving a genuine transient.
 */
public class SolverCTMCFauTransientTest {

    private static final double TSPAN = 10.0;
    private static final double STEP = 0.25;

    private static Network buildModel() {
        Network model = new Network("mm1_transient");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.5));
        queue.setService(oclass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static Matrix queueTrajectory(String transientMethod, double epsilon) {
        Network model = buildModel();
        SolverOptions o = new SolverCTMC(model).defaultOptions();
        o.cutoff = new Matrix(1, 1);
        o.cutoff.set(0, 0, 6);
        o.timespan = new double[]{0.0, TSPAN};
        o.timestep = STEP;
        if (transientMethod != null) {
            o.config.transient_method = transientMethod;
        }
        if (epsilon > 0) {
            o.config.fau_epsilon = epsilon;
        }
        SolverCTMC s = new SolverCTMC(model, o);
        s.getTranAvg();
        return s.getResults().QNt[1][0];
    }

    private static double maxDeviation(Matrix a, Matrix b) {
        assertEquals(a.getNumRows(), b.getNumRows(), "the two trajectories use different grids");
        double dev = 0.0;
        for (int i = 0; i < a.getNumRows(); i++) {
            dev = Math.max(dev, Math.abs(a.get(i, 0) - b.get(i, 0)));
        }
        return dev;
    }

    @Test
    public void fauAnswersWhatTheIntegratorAnswers() {
        Matrix ode = queueTrajectory(null, -1);
        Matrix fau = queueTrajectory("fau", -1);
        double dev = maxDeviation(ode, fau);
        assertTrue(dev < 1e-3, "fau deviates from the integrator by " + dev);
        // Guard against the silent fallback: a flat curve would satisfy the
        // comparison above vacuously if neither path ran a transient.
        assertTrue(fau.elementMax() - fau.elementMin() > 1e-3, "the fau curve is flat");
    }

    @Test
    public void tighteningEpsilonMovesFauLessThanTheGapToTheIntegrator() {
        Matrix ode = queueTrajectory(null, -1);
        Matrix loose = queueTrajectory("fau", -1);
        Matrix tight = queueTrajectory("fau", 1e-12);
        double residual = maxDeviation(loose, tight);
        double gap = maxDeviation(loose, ode);
        assertTrue(residual < 1e-5, "the fau answer moved by " + residual + " under a tighter tolerance");
        assertTrue(residual < gap, "the deviation from the integrator is not the integrator's");
    }

    @Test
    public void anUnknownTransientMethodIsRefused() {
        try {
            queueTrajectory("uniformization", -1);
            throw new AssertionError("an unknown transient_method was accepted");
        } catch (RuntimeException e) {
            assertTrue(String.valueOf(e.getMessage()).contains("transient_method"),
                    "unexpected message: " + e.getMessage());
        }
    }
}
