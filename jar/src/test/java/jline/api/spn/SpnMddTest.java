package jline.api.spn;

import jline.api.mdd.MDD;
import jline.api.mdd.MddDescriptor;
import jline.api.mdd.MddMcdResult;
import jline.api.mdd.Mdd_descriptor;
import jline.api.mdd.Mdd_mcd;
import jline.api.mdd.Mdd_reachset;

import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for the stochastic-Petri-net front end of the decision-diagram
 * aggregation, Spn_mdd, and for the exactness certificate Mdd_mcd returns.
 *
 * A cyclic net of single-server timed transitions IS the closed cyclic queueing
 * network of the same rates, so the SPN and queueing descriptors can be checked
 * against each other without any external reference: they must agree exactly.
 * The MATLAB twin is spn_mdd.m and the python twin is api/mdd/spn.py; all three
 * are compared numerically outside the suite.
 */
public class SpnMddTest {

    private static Network cyclicSpn(int nplaces, int ntokens, double[] rates,
                                     Distribution[] procs) {
        Network model = new Network("spn");
        Place[] pl = new Place[nplaces];
        Transition[] tr = new Transition[nplaces];
        for (int i = 0; i < nplaces; i++) {
            pl[i] = new Place(model, "P" + i);
            tr[i] = new Transition(model, "T" + i);
        }
        ClosedClass jc = new ClosedClass(model, "Class1", ntokens, pl[0]);
        for (int i = 0; i < nplaces; i++) {
            Mode m = tr[i].addMode("M" + i);
            tr[i].setDistribution(m, procs != null && procs[i] != null
                    ? procs[i] : new Exp(rates[i]));
            tr[i].setNumberOfServers(m, 1);
            tr[i].setEnablingConditions(m, jc, pl[i], 1);
            tr[i].setFiringOutcome(m, jc, pl[(i + 1) % nplaces], 1);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int i = 0; i < nplaces; i++) {
            P.set(jc, jc, pl[i], tr[i], 1.0);
            P.set(jc, jc, tr[i], pl[(i + 1) % nplaces], 1.0);
        }
        model.link(P);
        for (int i = 0; i < nplaces; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? ntokens : 0));
        }
        return model;
    }

    private static double[] placeTokens(Spn_mdd.SpnResult r, MddMcdResult out) {
        int np = r.info.nplacelevels;
        double[] q = new double[np];
        for (int l = 0; l < np; l++) {
            q[l] = out.QLen[l];
        }
        return q;
    }

    @Test
    public void cyclicNetReproducesTheClosedQueueingNetwork() {
        // a place with a single-server timed transition IS a station, so the two
        // descriptors must agree to numerical noise on the same rates
        double[] rates = {1.0, 1.5, 2.0};
        int N = 4;
        Spn_mdd.SpnResult r = Spn_mdd.spn_mdd(cyclicSpn(3, N, rates, null));
        MddMcdResult spnOut = Mdd_mcd.mdd_mcd(r.mdds, r.desc);

        double[] servers = {1.0, 1.0, 1.0};
        double[][] P = new double[3][3];
        for (int i = 0; i < 3; i++) {
            P[i][(i + 1) % 3] = 1.0;
        }
        MddDescriptor qn = Mdd_descriptor.mdd_descriptor(rates, P, servers, N, null, null);
        MDD qmdd = Mdd_reachset.mdd_reachset(qn.domain, qn.init, qn.nextfun);
        MddMcdResult qnOut = Mdd_mcd.mdd_mcd(qmdd.toStruct(), qn);

        assertEquals(qmdd.cardinality(), r.info.mdd.cardinality(),
                "the SPN and the queueing network must have the same reachable set");
        double[] spnQ = placeTokens(r, spnOut);
        for (int i = 0; i < 3; i++) {
            assertEquals(qnOut.QLen[i], spnQ[i], 1e-9,
                    "place " + i + " must hold the queue length of station " + i);
        }
    }

    @Test
    public void theMarkingInvariantIsCarriedAndChecked() {
        double[] rates = {1.0, 1.5, 2.0, 2.5};
        Spn_mdd.SpnResult r = Spn_mdd.spn_mdd(cyclicSpn(4, 4, rates, null));
        MddMcdResult out = Mdd_mcd.mdd_mcd(r.mdds, r.desc);
        double total = 0;
        for (int l = 0; l < r.info.nplacelevels; l++) {
            total += out.QLen[l];
        }
        assertEquals(4.0, total, 1e-6, "the token count must be conserved");
        assertEquals(4.0, r.desc.invariantValue, 1e-12,
                "the place invariant value must be the initial token count");
    }

    @Test
    public void preemptiveRepeatIsRefusedUnlessResumeIsAsked() {
        // a mode is disabled whenever its input place is empty, so a phase-type
        // firing time makes the repeat-vs-resume difference observable and the
        // descriptor must not silently pick resume
        double[] rates = {1.0, 1.5, 2.0};
        Distribution[] erl = new Distribution[3];
        for (int i = 0; i < 3; i++) {
            erl[i] = Erlang.fitMeanAndOrder(1.0 / rates[i], 2);
        }
        assertThrows(RuntimeException.class,
                () -> Spn_mdd.spn_mdd(cyclicSpn(3, 3, rates, erl)),
                "preemptive repeat has no Kronecker form and must be refused");

        Spn_mdd.SpnOptions resume = new Spn_mdd.SpnOptions();
        resume.phmemory = "resume";
        Spn_mdd.SpnResult r = Spn_mdd.spn_mdd(cyclicSpn(3, 3, rates, erl), resume);
        MddMcdResult out = Mdd_mcd.mdd_mcd(r.mdds, r.desc);
        double total = 0;
        for (int l = 0; l < r.info.nplacelevels; l++) {
            total += out.QLen[l];
        }
        assertEquals(3.0, total, 1e-6,
                "the resume semantics must still conserve the token count");
        assertTrue(r.info.everDisabled[0], "the guard must have observed a disabled mode");
    }

    @Test
    public void phaseLevelsAreAddedOnlyForMultiPhaseModes() {
        double[] rates = {1.0, 1.5, 2.0};
        Spn_mdd.SpnResult plain = Spn_mdd.spn_mdd(cyclicSpn(3, 3, rates, null));
        assertEquals(3, plain.desc.K, "an all-exponential net needs only place levels");
        for (int l = 0; l < plain.desc.K; l++) {
            assertEquals(1, plain.info.levelkind[l], "every level must be a place level");
        }

        Distribution[] erl = new Distribution[3];
        for (int i = 0; i < 3; i++) {
            erl[i] = Erlang.fitMeanAndOrder(1.0 / rates[i], 2);
        }
        Spn_mdd.SpnOptions resume = new Spn_mdd.SpnOptions();
        resume.phmemory = "resume";
        Spn_mdd.SpnResult ph = Spn_mdd.spn_mdd(cyclicSpn(3, 3, rates, erl), resume);
        assertEquals(6, ph.desc.K, "each 2-phase mode adds one phase level");
        assertEquals(3, ph.info.nplacelevels);
    }

    @Test
    public void exactnessCertificateFiresExactlyWhenNothingIsShared() {
        // K=2 exponential: each level-2 node is reached by one path, so
        // conditioning on the node equals conditioning on the path
        double[][] P2 = {{0.0, 1.0}, {1.0, 0.0}};
        MddDescriptor d2 = Mdd_descriptor.mdd_descriptor(new double[]{1.0, 1.5}, P2,
                new double[]{1.0, 1.0}, 4, null, null);
        MDD m2 = Mdd_reachset.mdd_reachset(d2.domain, d2.init, d2.nextfun);
        MddMcdResult o2 = Mdd_mcd.mdd_mcd(m2.toStruct(), d2);
        assertTrue(o2.noAggregation, "a tree-shaped diagram must certify as exact");
        for (int k = 0; k < o2.pathsPerLevel.length; k++) {
            assertEquals(1.0, o2.pathsPerLevel[k], 1e-12,
                    "no node may be shared when the certificate fires");
        }

        // K=3 shares nodes, so the certificate must NOT fire even though the
        // model is product form and the answer happens to be exact anyway
        double[][] P3 = new double[3][3];
        for (int i = 0; i < 3; i++) {
            P3[i][(i + 1) % 3] = 1.0;
        }
        MddDescriptor d3 = Mdd_descriptor.mdd_descriptor(new double[]{1.0, 1.5, 2.0}, P3,
                new double[]{1.0, 1.0, 1.0}, 4, null, null);
        MDD m3 = Mdd_reachset.mdd_reachset(d3.domain, d3.init, d3.nextfun);
        MddMcdResult o3 = Mdd_mcd.mdd_mcd(m3.toStruct(), d3);
        assertFalse(o3.noAggregation,
                "a shared diagram must not be certified; false means 'not certified', "
                        + "never 'approximate'");
        double mx = 0;
        for (int k = 0; k < o3.pathsPerLevel.length; k++) {
            mx = Math.max(mx, o3.pathsPerLevel[k]);
        }
        assertTrue(mx > 1.0, "sharing must show up as more than one path per node");
    }
}
