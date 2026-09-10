package jline.api.spn;

import jline.api.spn.Spn_metrics.SpnMetricsResult;
import jline.api.spn.Spn_pf.SpnPfResult;
import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for Spn_pf: derive the product form of a stochastic Petri net.
 *
 * <p>The MDD-rec paper (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490) takes the
 * g_l as GIVEN -- deriving them is declared out of scope in its Sec. 3.2 -- so
 * everything the api/spn functions do downstream was, until Spn_pf, unreachable
 * from a solver. These tests pin the derivation itself.</p>
 *
 * <p>FOUR INDEPENDENT ORACLES, because a product form that is merely plausible
 * is worse than none:</p>
 *
 * <ol>
 * <li>GLOBAL BALANCE. pi built from the derived g_l must satisfy pi Q = 0 on the
 * generator assembled independently from the net's own rate law. This is the
 * definition of stationarity and shares no code path with the derivation.</li>
 * <li>THE CLOSED FORM. On the cyclic net the Gordon-Newell factors are known,
 * so y must come out proportional to 1/mu.</li>
 * <li>THE CERTIFICATE. Deficiency, linkage classes, stoichiometric rank and weak
 * reversibility are structural facts that can be read off by hand here.</li>
 * <li>THE OTHER CODEBASES. y, G and the token counts are pinned at 12
 * decimals.</li>
 * </ol>
 *
 * <p>THE MODELS are the 3-place cyclic net at N = 4 with rates {1, 1.5, 2},
 * whose token counts are the value the four codebases already agree on for
 * Spn_rec; the same net under INFINITE-SERVER firing, which is mass action
 * rather than a constant rate and so selects the other psi; and a FORK-JOIN net,
 * whose marking is not a conserved job population at all -- a mode consumes one
 * token and produces two -- which is the case the MDD-rec paper exists to serve
 * and which SolverCTMC cannot solve.</p>
 */
public class SpnPfTest {

    private static final double[] RATES = {1.0, 1.5, 2.0};
    private static final int NJOBS = 4;

    /** MATLAB, native python and the C++ port on these nets, at %.12f. */
    private static final double[] CYCLIC_Y =
            {1.442249570307, 0.961499713538, 0.721124785154};
    private static final double CYCLIC_G = 19.934426352559;
    private static final double[] CYCLIC_TOKENS =
            {2.249874392899, 1.069837548149, 0.680288058952};
    private static final double[] CYCLIC_MODEX =
            {0.869201138838, 0.869201138838, 0.869201138838};
    private static final double[] FJ_Y =
            {1.028383350947, 1.381974961646, 1.381974961646, 0.703630713806};
    private static final double FJ_G = 20.320471787308;
    private static final double[] FJ_TOKENS =
            {0.710262429604, 1.864198441415, 1.864198441415, 0.425539128981};
    private static final double[] FJ_MODEX =
            {0.607360882508, 0.607360882508, 0.607360882508};

    /**
     * P0 -&gt; T0 -&gt; P1 -&gt; T1 -&gt; P2 -&gt; T2 -&gt; P0, one token class.
     *
     * <p>servers = 1 makes every mode fire at its rate constant, which is
     * psi = 1; a server count AT OR ABOVE the token population makes it fire at
     * rate lambda*min(m_p, servers) = lambda*m_p, which is mass action. The
     * MATLAB, python and C++ twins pass an infinite count for the same effect;
     * setNumberOfServers takes an Integer here, and on this net the two are the
     * same function because m_p never exceeds the population.</p>
     */
    private static Network cyclicSpn(int ntokens, int servers) {
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
            tr[i].setNumberOfServers(m, Integer.valueOf(servers));
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

    /**
     * P0 -(Tf)-&gt; P1 + P2 -(Tj)-&gt; P3 -(Tb)-&gt; P0.
     *
     * <p>Tf consumes ONE token and produces TWO, Tj the reverse, so the marking
     * is not a conserved population and the net has no queueing-network
     * counterpart. Its place invariant is (2, 1, 1, 2).</p>
     */
    private static Network forkJoinSpn(int ntokens) {
        Network model = new Network("fj");
        Place[] pl = new Place[4];
        for (int i = 0; i < 4; i++) {
            pl[i] = new Place(model, "P" + i);
        }
        Transition tf = new Transition(model, "Tf");
        Transition tj = new Transition(model, "Tj");
        Transition tb = new Transition(model, "Tb");
        ClosedClass jc = new ClosedClass(model, "C", ntokens, pl[0]);

        Mode mf = tf.addMode("f");
        tf.setDistribution(mf, new Exp(1.3));
        tf.setNumberOfServers(mf, 1);
        tf.setEnablingConditions(mf, jc, pl[0], 1);
        tf.setFiringOutcome(mf, jc, pl[1], 1);
        tf.setFiringOutcome(mf, jc, pl[2], 1);

        Mode mj = tj.addMode("j");
        tj.setDistribution(mj, new Exp(0.7));
        tj.setNumberOfServers(mj, 1);
        tj.setEnablingConditions(mj, jc, pl[1], 1);
        tj.setEnablingConditions(mj, jc, pl[2], 1);
        tj.setFiringOutcome(mj, jc, pl[3], 1);

        Mode mb = tb.addMode("b");
        tb.setDistribution(mb, new Exp(1.9));
        tb.setNumberOfServers(mb, 1);
        tb.setEnablingConditions(mb, jc, pl[3], 1);
        tb.setFiringOutcome(mb, jc, pl[0], 1);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, pl[0], tf, 1.0);
        P.set(jc, jc, tf, pl[1], 1.0);
        P.set(jc, jc, tf, pl[2], 1.0);
        P.set(jc, jc, pl[1], tj, 1.0);
        P.set(jc, jc, pl[2], tj, 1.0);
        P.set(jc, jc, tj, pl[3], 1.0);
        P.set(jc, jc, pl[3], tb, 1.0);
        P.set(jc, jc, tb, pl[0], 1.0);
        model.link(P);
        for (int i = 0; i < 4; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? ntokens : 0));
        }
        return model;
    }

    /**
     * max |pi Q| over the enumerated reachable set.
     *
     * <p>The generator is assembled here from the net's own rate law,
     * independently of the derivation, so agreement is evidence and not a
     * tautology.</p>
     */
    private static double balanceResidual(SpnPfResult pf) {
        Spn_mdd.SpnInfo info = pf.spn.info;
        int L = info.nplacelevels;
        int[][] st = info.mdd.enumerate();
        int n = st.length;
        Map<String, Integer> idx = new HashMap<String, Integer>();
        for (int i = 0; i < n; i++) {
            idx.put(key(st[i], L), Integer.valueOf(i));
        }
        double[][] Q = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int e = 0; e < info.modes.size(); e++) {
                Spn_mdd.SpnMode mde = info.modes.get(e);
                boolean enabled = true;
                for (int l = 0; l < L; l++) {
                    if (st[i][l] < mde.enab[l]) {
                        enabled = false;
                    }
                }
                if (!enabled) {
                    continue;
                }
                double rate = mde.D1[0][0];
                if ("massaction".equals(pf.kind)) {
                    for (int l = 0; l < L; l++) {
                        for (int j = 0; j < (int) mde.enab[l]; j++) {
                            rate *= (st[i][l] - j);
                        }
                    }
                } else {
                    double deg = Double.POSITIVE_INFINITY;
                    for (int l = 0; l < L; l++) {
                        if (mde.enab[l] > 0) {
                            deg = Math.min(deg, Math.floor(st[i][l] / mde.enab[l]));
                        }
                    }
                    if (Double.isInfinite(deg)) {
                        deg = 1;
                    }
                    rate *= Math.min(deg, mde.srv);
                }
                int[] t = new int[L];
                for (int l = 0; l < L; l++) {
                    t[l] = (int) (st[i][l] - mde.enab[l] + mde.fire[l]);
                }
                Integer at = idx.get(key(t, L));
                assertTrue(at != null, "a firing leaves the reachable set");
                Q[i][at.intValue()] += rate;
            }
        }
        double[] pi = new double[n];
        double tot = 0;
        for (int i = 0; i < n; i++) {
            double p = 1;
            for (int l = 0; l < L; l++) {
                p *= pf.g[l][st[i][l]];
            }
            pi[i] = p;
            tot += p;
        }
        for (int i = 0; i < n; i++) {
            pi[i] /= tot;
            Q[i][i] -= sum(Q[i]);
        }
        double res = 0;
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int i = 0; i < n; i++) {
                s += pi[i] * Q[i][j];
            }
            res = Math.max(res, Math.abs(s));
        }
        return res;
    }

    private static double sum(double[] v) {
        double s = 0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
        }
        return s;
    }

    private static String key(int[] v, int L) {
        StringBuilder sb = new StringBuilder();
        for (int l = 0; l < L; l++) {
            sb.append(v[l]).append(',');
        }
        return sb.toString();
    }

    @Test
    public void cyclicProductFormIsTheGordonNewellOne() {
        SpnPfResult pf = Spn_pf.spn_pf(cyclicSpn(NJOBS, 1));
        assertEquals("geometric", pf.kind);
        // y is proportional to 1/mu: the closed-form Gordon-Newell factors, up
        // to the gauge the minimum-norm solution fixes.
        for (int l = 0; l < 3; l++) {
            assertEquals(RATES[0] / RATES[l], pf.y[l] / pf.y[0], 1e-11);
            assertEquals(CYCLIC_Y[l], pf.y[l], 1e-11);
        }
    }

    @Test
    public void cyclicCertificateIsTheStructuralOne() {
        SpnPfResult pf = Spn_pf.spn_pf(cyclicSpn(NJOBS, 1));
        // three complexes {e0, e1, e2}, one linkage class (the cycle is
        // strongly connected), stoichiometric rank 2, deficiency 3 - 1 - 2 = 0
        assertEquals(3, pf.complexes.length);
        assertEquals(1, pf.linkage);
        assertEquals(2, pf.srank);
        assertEquals(0, pf.deficiency);
        assertTrue(pf.weaklyReversible);
        assertTrue(pf.residual < 1e-12);
    }

    @Test
    public void theDerivedLawSatisfiesGlobalBalance() {
        assertTrue(balanceResidual(Spn_pf.spn_pf(cyclicSpn(NJOBS, 1))) < 1e-12);
        assertTrue(balanceResidual(Spn_pf.spn_pf(forkJoinSpn(3))) < 1e-12);
        assertTrue(balanceResidual(
                Spn_pf.spn_pf(cyclicSpn(NJOBS, NJOBS))) < 1e-12);
    }

    @Test
    public void infiniteServerFiringSelectsMassAction() {
        SpnPfResult pf = Spn_pf.spn_pf(cyclicSpn(NJOBS, NJOBS));
        assertEquals("massaction", pf.kind);
        // the factors carry the 1/k! that psi = prod 1/m! puts there
        for (int l = 0; l < 3; l++) {
            double fact = 1;
            for (int k = 0; k < pf.g[l].length; k++) {
                if (k > 0) {
                    fact *= k;
                }
                assertEquals(Math.pow(pf.y[l], k) / fact, pf.g[l][k], 1e-11);
            }
        }
    }

    @Test
    public void cyclicValuesArePinnedAcrossTheCodebases() {
        SpnPfResult pf = Spn_pf.spn_pf(cyclicSpn(NJOBS, 1));
        SpnMetricsResult met = Spn_metrics.spn_metrics(pf.spn.mdds, pf.g, pf.spn.info);
        assertEquals(CYCLIC_G, met.G, 1e-9);
        for (int l = 0; l < 3; l++) {
            assertEquals(CYCLIC_TOKENS[l], met.tokens[l], 1e-9);
            assertEquals(CYCLIC_MODEX[l], met.modeTput[l], 1e-9);
        }
    }

    @Test
    public void forkJoinValuesArePinnedAcrossTheCodebases() {
        SpnPfResult pf = Spn_pf.spn_pf(forkJoinSpn(3));
        SpnMetricsResult met = Spn_metrics.spn_metrics(pf.spn.mdds, pf.g, pf.spn.info);
        assertEquals(0, pf.deficiency);
        assertTrue(pf.weaklyReversible);
        assertEquals(FJ_G, met.G, 1e-9);
        for (int l = 0; l < 4; l++) {
            assertEquals(FJ_Y[l], pf.y[l], 1e-11);
            assertEquals(FJ_TOKENS[l], met.tokens[l], 1e-9);
        }
        // every token that forks must later join and return, so on this cycle
        // the three modes share one throughput -- a conservation law the
        // derivation was never told about
        for (int e = 0; e < 3; e++) {
            assertEquals(FJ_MODEX[e], met.modeTput[e], 1e-9);
            assertEquals(met.modeTput[0], met.modeTput[e], 1e-11);
        }
    }

    @Test
    public void anInhibitorArcIsRefusedByName() {
        Network model = cyclicSpn(NJOBS, 1);
        Transition t0 = null;
        Place p2 = null;
        for (int i = 0; i < model.getNodes().size(); i++) {
            jline.lang.nodes.Node nd = model.getNodes().get(i);
            if (t0 == null && nd instanceof Transition) {
                t0 = (Transition) nd;
            }
            if (nd instanceof Place && "P2".equals(nd.getName())) {
                p2 = (Place) nd;
            }
        }
        t0.setInhibitingConditions(t0.getModes().get(0), model.getClasses().get(0), p2, 2);
        RuntimeException ex = assertThrows(RuntimeException.class,
                () -> Spn_pf.spn_pf(model));
        assertTrue(ex.getMessage().contains("inhibitor"), ex.getMessage());
    }
}
