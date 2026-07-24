/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fes;

import jline.VerboseLevel;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static jline.TestTools.LOOSE_FINE_TOL;

/**
 * Exactness contract for FES (Flow-Equivalent Server) aggregation on
 * PRODUCT-FORM MULTICLASS closed networks.
 *
 * Theory (Kritzinger, van Wyk and Krzesinski 1982, "A generalisation of
 * Norton's theorem for multiclass queueing networks"): replacing a subset of
 * centres in a network satisfying local balance by a single composite centre
 * with the state-dependent per-class throughput tau_r(n1,...,nR) leaves the
 * marginal distribution of the remaining (complement) centres UNCHANGED.
 * Hence, for a product-form multiclass model, every per-class metric at the
 * complement stations must match the exact full-model solution.
 *
 * These models are product-form (PS queues + Delay, BCMP), so SolverNC 'exact'
 * is the ground truth on both the full and the reduced network.
 *
 * Findings from the MATLAB reference validation (matlab/tests/test_fes_exactness.m):
 *   - SINGLE-EXIT subsets (tandem): FES is exact to machine precision (~1e-16).
 *   - MULTI-EXIT subsets: currently NOT exact (~33% error) because
 *     ModelAdapter.aggregateFES weights the FES exit routing UNIFORMLY
 *     (visitRatios = ones(1,nSub)/nSub) instead of using the subnetwork's true
 *     relative arrival rates xi_cr (Kritzinger Eq. 4.6). See testMultiExit_*.
 *
* Exact for closed product-form multiclass networks after FESAggregator was fixed
 * (stale-population struct refresh, original-routing complement paths, cdscaling stores
 * the composite rate). Exit routing was fixed in FESAggregator to weight the FES -> complement split
 * by the per-class subnetwork visit ratios (was uniform); the MATLAB pipeline is
 * now exact for all cases below (matlab/tests/test_fes_exactness.m passes to
 * ~1e-16). This test remains disabled because the jar solver still does not
 * consume the class-dependence handle correctly (see FESAggregatorTest): even the single-exit
 * tandem cases come out 5-96% off. Enable once the jar SolverNC/AMVA-LD path
 * applies the per-class cdscaling rates. The exit-routing fix is a prerequisite, done.
 */
public class FESExactnessTest {

    private static SolverOptions exactNc() {
        SolverOptions o = Solver.defaultOptions();
        o.verbose = VerboseLevel.SILENT;
        o.method = "exact";
        return o;
    }

    private static NetworkAvgTable solveExact(Network m) {
        return new SolverNC(m, exactNc()).getAvgTable();
    }

    /**
     * Assert that every per-class metric (QLen, Util, Tput, RespT) at the named
     * complement stations matches between the full and reduced models.
     */
    private static void assertComplementMatches(NetworkAvgTable full, NetworkAvgTable reduced,
                                                List<String> complementStations, double tol) {
        List<String> fSt = full.getStationNames(),  fCl = full.getClassNames();
        List<String> rSt = reduced.getStationNames(), rCl = reduced.getClassNames();
        List<List<Double>> fM = Arrays.asList(full.getQLen(), full.getUtil(), full.getTput(), full.getRespT());
        List<List<Double>> rM = Arrays.asList(reduced.getQLen(), reduced.getUtil(), reduced.getTput(), reduced.getRespT());
        String[] metricName = {"QLen", "Util", "Tput", "RespT"};

        for (int i = 0; i < fSt.size(); i++) {
            if (!complementStations.contains(fSt.get(i))) continue;
            // locate matching row in reduced table by (station, class)
            int j = -1;
            for (int k = 0; k < rSt.size(); k++) {
                if (rSt.get(k).equals(fSt.get(i)) && rCl.get(k).equals(fCl.get(i))) { j = k; break; }
            }
            assertTrue(j >= 0, "Complement row missing in reduced model: " + fSt.get(i) + "/" + fCl.get(i));
            for (int mIdx = 0; mIdx < metricName.length; mIdx++) {
                double a = fM.get(mIdx).get(i), b = rM.get(mIdx).get(j);
                double relErr = Math.abs(a - b) / Math.max(Math.abs(a), 1e-9);
                assertTrue(relErr < tol, String.format(
                    "%s/%s %s: full=%.10g reduced=%.10g relErr=%.3e (tol=%.1e)",
                    fSt.get(i), fCl.get(i), metricName[mIdx], a, b, relErr, tol));
            }
        }
    }

    /** Closed PF multiclass tandem: Delay -> Q1 -> Q2 -> Q3 -> Delay (PS queues). */
    private static Network buildTandem() {
        Network model = new Network("PF_tandem");
        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 3, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, delay, 0);
        delay.setService(c1, Exp.fitMean(1.0)); delay.setService(c2, Exp.fitMean(1.5));
        q1.setService(c1, Exp.fitMean(0.5));    q1.setService(c2, Exp.fitMean(0.8));
        q2.setService(c1, Exp.fitMean(0.3));    q2.setService(c2, Exp.fitMean(0.6));
        q3.setService(c1, Exp.fitMean(0.4));    q3.setService(c2, Exp.fitMean(0.7));
        RoutingMatrix P = model.initRoutingMatrix();
        for (ClosedClass c : new ClosedClass[]{c1, c2}) {
            P.set(c, c, delay, q1, 1.0);
            P.set(c, c, q1, q2, 1.0);
            P.set(c, c, q2, q3, 1.0);
            P.set(c, c, q3, delay, 1.0);
        }
        model.link(P);
        return model;
    }

    private void runTandemCase(List<String> subsetNames, List<String> complementNames) {
        Network full = buildTandem();
        NetworkAvgTable fullT = solveExact(full);

        // resolve subset Station objects by name
        List<Station> subset = new java.util.ArrayList<>();
        for (Node n : full.getNodes()) {
            if (n instanceof Station && subsetNames.contains(n.getName())) {
                subset.add((Station) n);
            }
        }
        FESResult fes = ModelAdapter.aggregateFES(full, subset);
        NetworkAvgTable redT = solveExact(fes.getFesModel());

        assertComplementMatches(fullT, redT, complementNames, LOOSE_FINE_TOL);
    }

    @Test
    public void testExactTandem_aggQ1Q2() {
        runTandemCase(Arrays.asList("Q1", "Q2"), Arrays.asList("Delay", "Q3"));
    }

    @Test
    public void testExactTandem_aggQ2Q3() {
        runTandemCase(Arrays.asList("Q2", "Q3"), Arrays.asList("Delay", "Q1"));
    }

    @Test
    public void testExactTandem_aggQ1Q2Q3() {
        runTandemCase(Arrays.asList("Q1", "Q2", "Q3"), Arrays.asList("Delay"));
    }

    /**
     * Multi-exit subset: Delay -> Q1 -> {Q2 (0.6), Q3 (0.4)}; Q2,Q3 -> Delay.
     * Aggregating {Q1, Q2} yields a subset with two exits (Q2->Delay and the
     * direct Q1->Q3 branch). By Norton's theorem this is exact once the FES exit
     * routing is weighted by the subnetwork visit ratios and the composite rate
     * is the subnetwork departure rate.
     */
    @Test
    public void testMultiExit_shouldBeExact() {
        Network model = new Network("PF_prob");
        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, delay, 0);
        delay.setService(c1, Exp.fitMean(1.0)); delay.setService(c2, Exp.fitMean(1.2));
        q1.setService(c1, Exp.fitMean(0.5));    q1.setService(c2, Exp.fitMean(0.7));
        q2.setService(c1, Exp.fitMean(0.4));    q2.setService(c2, Exp.fitMean(0.5));
        q3.setService(c1, Exp.fitMean(0.6));    q3.setService(c2, Exp.fitMean(0.3));
        RoutingMatrix P = model.initRoutingMatrix();
        for (ClosedClass c : new ClosedClass[]{c1, c2}) {
            P.set(c, c, delay, q1, 1.0);
            P.set(c, c, q1, q2, 0.6);
            P.set(c, c, q1, q3, 0.4);
            P.set(c, c, q2, delay, 1.0);
            P.set(c, c, q3, delay, 1.0);
        }
        model.link(P);

        NetworkAvgTable fullT = solveExact(model);
        List<Station> subset = Arrays.asList(q1, q2);
        FESResult fes = ModelAdapter.aggregateFES(model, subset);
        NetworkAvgTable redT = solveExact(fes.getFesModel());

        assertComplementMatches(fullT, redT, Arrays.asList("Delay", "Q3"), LOOSE_FINE_TOL);
    }
}
