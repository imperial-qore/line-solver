package jline.api.spn;

import jline.api.mdd.MDD;
import jline.api.mdd.MddDescriptor;
import jline.api.mdd.MddMcdResult;
import jline.api.mdd.MddStruct;
import jline.api.mdd.Mdd_descriptor;
import jline.api.mdd.Mdd_mcd;
import jline.api.mdd.Mdd_reachset;
import jline.api.mdd.Mdd_rec;
import jline.api.spn.Spn_metrics.SpnMetricsResult;
import jline.api.spn.Spn_rec_enabled.SpnEnabling;
import jline.api.spn.Spn_sinvariants.SpnInvariants;

import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for MDD-rec and the SPN measures built on it: Mdd_rec,
 * Mdd_rec.mdd_rec_marginal, Spn_rec_enabled, Spn_metrics, Spn_sinvariants,
 * Spn_conv.
 *
 * <p>S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant
 * for product-form models of distributed systems with synchronisation", FGCS 111
 * (2020) 475-490.</p>
 *
 * <p>THREE INDEPENDENT ORACLES, because a normalising constant is a single
 * number that a wrong recursion can still produce plausibly:</p>
 *
 * <ol>
 * <li>EXPLICIT SUM. G is also the sum over the enumerated reachable set of the
 * product of the g_l. The diagram walk and the enumeration share no code path,
 * so agreement pins the recursion itself.</li>
 * <li>THE CONVOLUTION. Spn_conv decomposes {m : S m = V} instead of walking the
 * diagram; on this net the two must agree exactly (FGCS Sec. 5.2).</li>
 * <li>THE OTHER CODEBASES. The marginals must reproduce MATLAB's, the C++
 * port's and native python's queue lengths on the same net, and the mode
 * throughputs must match the level aggregation's X, which comes from a
 * different algorithm on a different descriptor.</li>
 * </ol>
 *
 * <p>THE MODEL is the 3-place cyclic net at N = 4 with firing rates
 * {1, 1.5, 2}: a Gordon-Newell network, so its product form is known in closed
 * form, g_l(n) = (1/mu_l)^n with unit visit ratios, and no product-form TEST is
 * needed to obtain the g_l -- which is the part the paper itself declares out of
 * scope.</p>
 */
public class SpnRecTest {

    private static final double[] RATES = {1.0, 1.5, 2.0};
    private static final int NJOBS = 4;
    /** MATLAB, the C++ port and native python on this net, at %.12f. */
    private static final double[] QLEN_REF = {2.249874392899, 1.069837548149, 0.680288058952};

    private static Network cyclicSpn(int ntokens) {
        Network model = new Network("spn");
        Place[] pl = new Place[3];
        Transition[] tr = new Transition[3];
        for (int i = 0; i < 3; i++) {
            pl[i] = new Place(model, "P" + i);
            tr[i] = new Transition(model, "T" + i);
        }
        ClosedClass jc = new ClosedClass(model, "Class1", ntokens, pl[0]);
        for (int i = 0; i < 3; i++) {
            Mode m = tr[i].addMode("fire");
            tr[i].setDistribution(m, new Exp(RATES[i]));
            tr[i].setNumberOfServers(m, 1);
            tr[i].setEnablingConditions(m, jc, pl[i], 1);
            tr[i].setFiringOutcome(m, jc, pl[(i + 1) % 3], 1);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int i = 0; i < 3; i++) {
            P.set(jc, jc, pl[i], tr[i], 1.0);
            P.set(jc, jc, tr[i], pl[(i + 1) % 3], 1.0);
        }
        model.link(P);
        for (int i = 0; i < 3; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? ntokens : 0));
        }
        return model;
    }

    /** Gordon-Newell factors: g_l(n) = (v_l/mu_l)^n, unit visit ratios on a cycle. */
    private static double[][] gordonNewellG(MddStruct mdds) {
        double[][] g = new double[mdds.K][];
        for (int l = 0; l < mdds.K; l++) {
            g[l] = new double[mdds.domain[l]];
            double x = 1.0;
            for (int n = 0; n < mdds.domain[l]; n++) {
                g[l][n] = x;
                x /= RATES[l];
            }
        }
        return g;
    }

    private static Spn_mdd.SpnResult solved() {
        return Spn_mdd.spn_mdd(cyclicSpn(NJOBS));
    }

    @Test
    public void mddRecMatchesTheExplicitSum() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        int[][] states = r.info.mdd.enumerate();
        assertEquals(15, states.length);              // C(4+2,2), the closed lattice
        double expected = 0.0;
        for (int s = 0; s < states.length; s++) {
            double p = 1.0;
            for (int l = 0; l < r.mdds.K; l++) {
                p *= g[l][states[s][l]];
            }
            expected += p;
        }
        assertEquals(expected, Mdd_rec.mdd_rec(r.mdds, g), 1e-12 * Math.max(1.0, expected));
    }

    @Test
    public void anEmptyMaskIsTheUnmaskedRecursion() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        boolean[][] mask = new boolean[r.mdds.K][];
        for (int l = 0; l < r.mdds.K; l++) {
            mask[l] = new boolean[r.mdds.domain[l]];
            for (int v = 0; v < r.mdds.domain[l]; v++) {
                mask[l][v] = true;
            }
        }
        assertEquals(Mdd_rec.mdd_rec(r.mdds, g), Mdd_rec.mdd_rec_masked(r.mdds, g, mask), 1e-12);
    }

    @Test
    public void theMarginalsPartitionTheNormalisingConstant() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        double G = Mdd_rec.mdd_rec(r.mdds, g);
        for (int l = 0; l < r.info.nplacelevels; l++) {
            double[] mass = Mdd_rec.mdd_rec_marginal(r.mdds, g, l);
            double sum = 0.0;
            for (int k = 0; k < mass.length; k++) {
                sum += mass[k];
            }
            assertEquals(G, sum, 1e-12 * G);
        }
    }

    @Test
    public void sinvariantsFindTheTokenConservationLaw() {
        SpnInvariants inv = Spn_sinvariants.spn_sinvariants(cyclicSpn(NJOBS).getStruct());
        // the only minimal-support invariant of a cycle is "the tokens are conserved"
        assertEquals(1, inv.S.length);
        for (int p = 0; p < 3; p++) {
            assertEquals(1L, inv.S[0][p]);
        }
        assertEquals(NJOBS, inv.V[0]);
        assertEquals(NJOBS, inv.m0[0]);
        assertEquals(0L, inv.m0[1]);
        assertEquals(0L, inv.m0[2]);
    }

    @Test
    public void theConvolutionAgreesWithTheDiagramWalk() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        SpnInvariants inv = Spn_sinvariants.spn_sinvariants(cyclicSpn(NJOBS).getStruct());
        double G = Mdd_rec.mdd_rec(r.mdds, g);
        assertEquals(G, Spn_conv.spn_conv(inv, g), 1e-12 * G);
    }

    @Test
    public void metricsReproduceTheReferenceQueueLengths() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        SpnMetricsResult met = Spn_metrics.spn_metrics(r.mdds, g, r.info);
        double total = 0.0;
        for (int l = 0; l < r.info.nplacelevels; l++) {
            assertEquals(QLEN_REF[l], met.tokens[l], 1e-9);
            assertEquals(1.0 - met.marginal[l][0], met.placeUtil[l], 1e-12);
            total += met.tokens[l];
        }
        assertEquals(NJOBS, total, 1e-12);
    }

    @Test
    public void modeThroughputMatchesTheLevelAggregation() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        SpnMetricsResult met = Spn_metrics.spn_metrics(r.mdds, g, r.info);
        double[][] P = new double[3][3];
        double[] mu = new double[3];
        double[] servers = new double[3];
        for (int i = 0; i < 3; i++) {
            P[i][(i + 1) % 3] = 1.0;
            mu[i] = RATES[i];
            servers[i] = 1.0;
        }
        MddDescriptor qdesc = Mdd_descriptor.mdd_descriptor(mu, P, servers, NJOBS, null, null);
        MDD qdiag = Mdd_reachset.mdd_reachset(qdesc.domain, qdesc.init, qdesc.nextfun);
        MddMcdResult qout = Mdd_mcd.mdd_mcd(qdiag.toStruct(), qdesc);
        // a cycle carries one flow, so every transition sees the same throughput
        for (int e = 0; e < met.modeTput.length; e++) {
            assertEquals(qout.X[0], met.modeTput[e], 1e-6 * qout.X[0]);
        }
        for (int l = 0; l < r.info.nplacelevels; l++) {
            assertEquals(qout.X[0], met.placeTput[l], 1e-6 * qout.X[0]);
        }
    }

    @Test
    public void enablingDegreesAreADistribution() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        double G = Mdd_rec.mdd_rec(r.mdds, g);
        SpnMetricsResult met = Spn_metrics.spn_metrics(r.mdds, g, r.info);
        for (int e = 0; e < r.info.modes.size(); e++) {
            SpnEnabling en = Spn_rec_enabled.spn_rec_enabled(r.mdds, g, r.info.modes.get(e),
                    r.info.nplacelevels);
            assertEquals(NJOBS, en.maxDegree);        // one token per firing set
            assertEquals(G, en.ge[0], 1e-12 * G);
            double sum = 0.0;
            for (int k = 0; k < en.eq.length; k++) {
                sum += en.eq[k];
            }
            assertEquals(G, sum, 1e-12 * G);
            assertEquals(met.modeUtil[e], en.ge[1] / G, 1e-12);
            // P(e >= k) is non-increasing by construction; a mask that grew the
            // set would be a silently wrong measure rather than an error
            for (int k = 0; k + 1 < en.ge.length; k++) {
                assertTrue(en.ge[k + 1] <= en.ge[k] + 1e-15);
            }
        }
    }

    @Test
    public void aModeWithNoInputPlaceIsRefused() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        Spn_mdd.SpnMode mde = r.info.modes.get(0);
        double[] saved = mde.enab.clone();
        for (int l = 0; l < mde.enab.length; l++) {
            mde.enab[l] = 0.0;
        }
        assertThrows(RuntimeException.class,
                () -> Spn_rec_enabled.spn_rec_enabled(r.mdds, g, mde, r.info.nplacelevels));
        System.arraycopy(saved, 0, mde.enab, 0, saved.length);
    }

    @Test
    public void aShapeMismatchInGIsRefused() {
        Spn_mdd.SpnResult r = solved();
        double[][] g = gordonNewellG(r.mdds);
        double[][] shorter = new double[g.length - 1][];
        System.arraycopy(g, 0, shorter, 0, g.length - 1);
        assertThrows(RuntimeException.class, () -> Mdd_rec.mdd_rec(r.mdds, shorter));
    }
}
