package jline.lang.processes;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.util.SerializableFunction;

/**
 * The Laplace-Stieltjes transform off the real axis.
 *
 * A transform is evaluated at a COMPLEX argument by everything that inverts it
 * or locates its roots, so sn.lst carries the complex overload. The oracles here
 * are the closed forms themselves, continued to the complex plane: mu/(mu+s) for
 * the exponential, (r/(r+s))^k for the Erlang, exp(-sd) for the deterministic
 * law, and (e^-sa - e^-sb)/(s(b-a)) for the uniform. The generic CDF-sum
 * fallback on Distribution is checked against the one law it must reproduce
 * without a closed form of its own.
 */
public class ComplexLstTest {

    private static void assertComplexEquals(Complex expected, Complex actual, double tol) {
        assertEquals(expected.getReal(), actual.getReal(), tol);
        assertEquals(expected.getImaginary(), actual.getImaginary(), tol);
    }

    @Test
    public void exponentialMatchesItsClosedFormOffTheRealAxis() {
        final double mu = 1.4;
        final Exp d = new Exp(mu);
        final Complex s = new Complex(0.5, 1.0);
        // mu / (mu + s)
        assertComplexEquals(new Complex(mu, 0.0).divide(s.add(mu)), d.evalLST(s), 1e-12);
        // the real overload stays consistent with the complex one
        assertEquals(d.evalLST(0.5), d.evalLST(new Complex(0.5, 0.0)).getReal(), 1e-12);
    }

    @Test
    public void erlangMatchesItsClosedFormOffTheRealAxis() {
        final Erlang d = Erlang.fitMeanAndSCV(1.0, 0.5);   // two phases of rate 2
        final Complex s = new Complex(0.3, -0.7);
        final Complex one = new Complex(2.0, 0.0).divide(s.add(2.0));
        assertComplexEquals(one.multiply(one), d.evalLST(s), 1e-12);
    }

    @Test
    public void deterministicIsExactWhereTheGenericSumCannotSeeIt() {
        final double v = 0.75;
        final Det d = new Det(v);
        final Complex s = new Complex(0.4, 2.0);
        assertComplexEquals(s.multiply(-v).exp(), d.evalLST(s), 1e-12);
    }

    @Test
    public void uniformMatchesItsClosedForm() {
        final double a = 0.1, b = 0.9;
        final Uniform d = new Uniform(a, b);
        final Complex s = new Complex(0.6, 0.8);
        final Complex ref = s.multiply(-a).exp().subtract(s.multiply(-b).exp())
                .divide(s.multiply(b - a));
        assertComplexEquals(ref, d.evalLST(s), 1e-12);
        assertEquals(1.0, d.evalLST(new Complex(0.0, 0.0)).getReal(), 1e-12);
    }

    @Test
    public void genericCdfSumReproducesAKnownTransform() {
        // Gamma has no complex override, so it exercises the Distribution
        // fallback; its transform is (1 + s/beta)^-alpha with alpha = 1/scv.
        final Gamma d = (Gamma) Gamma.fitMeanAndSCV(1.0, 0.5);
        final Complex s = new Complex(0.5, 0.0);
        assertEquals(d.evalLST(0.5), d.evalLST(s).getReal(), 1e-3);
    }

    @Test
    public void networkStructCarriesTheComplexTransform() {
        Network model = new Network("lst");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "K");
        OpenClass cls = new OpenClass(model, "C");
        src.setArrival(cls, new Exp(0.5));
        q.setService(cls, Erlang.fitMeanAndSCV(1.0, 0.5));
        model.link(model.serialRouting(src, q, snk));
        jline.lang.NetworkStruct sn = model.getStruct(true);
        assertNotNull(sn.lst);
        SerializableFunction<Complex, Complex> f = sn.lst.get(sn.stations.get(1)).get(cls);
        assertNotNull(f);
        final Complex s = new Complex(0.3, -0.7);
        final Complex one = new Complex(2.0, 0.0).divide(s.add(2.0));
        assertComplexEquals(one.multiply(one), f.apply(s), 1e-12);
    }
}
