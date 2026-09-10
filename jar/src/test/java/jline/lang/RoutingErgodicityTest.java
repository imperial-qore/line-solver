/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lang;

import java.util.List;

import org.junit.jupiter.api.Test;

import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Routing reducibility and its repair, against MATLAB
 * {@code @MNetwork/isRoutingErgodic}, {@code getReducibilityInfo},
 * {@code getAbsorbingStations} and {@code makeErgodic}.
 *
 * <p>The reducible fixture is D -&gt; Q1 -&gt; Q2 -&gt; Q2: Q2 absorbs, D and Q1
 * are transient, and there are three strongly connected components.
 */
public class RoutingErgodicityTest {

    private static Network reducible() {
        Network model = new Network("red");
        Delay d = new Delay(model, "D");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 2, d, 0);
        d.setService(c, new Exp(1.0));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        Matrix M = new Matrix(3, 3);
        M.set(0, 1, 1.0);
        M.set(1, 2, 1.0);
        M.set(2, 2, 1.0);
        P.set(c, c, M);
        model.link(P);
        return model;
    }

    private static Network ergodic() {
        Network model = new Network("erg");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 2, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(d, q));
        return model;
    }

    @Test
    public void reducibleStructureMatchesMatlab() {
        Network model = reducible();
        Network.RoutingErgodicityResult erg = model.isRoutingErgodic();
        assertFalse(erg.isErgodic);
        assertTrue(erg.isReducible);
        assertEquals(3, erg.numSCCs);
        assertEquals(1, erg.absorbingStations.size());
        assertEquals("Q2", erg.absorbingStations.get(0));
        assertTrue(erg.transientStations.contains("D"));
        assertTrue(erg.transientStations.contains("Q1"));
    }

    @Test
    public void reducibilityInfoSuggestsOneFixPerAbsorbingStation() {
        RoutingErgodicity.ReducibilityInfo info = reducible().getReducibilityInfo();
        assertFalse(info.isRoutingErgodic);
        assertEquals(1, info.suggestedFixes.size());
        assertEquals("Route jobs from Q2 back to D (e.g., P{class}(Q2, D) = 1.0)",
                info.suggestedFixes.get(0));
    }

    @Test
    public void absorbingStationsAreReturnedAsObjects() {
        List<Station> abs = reducible().getAbsorbingStations();
        assertEquals(1, abs.size());
        assertEquals("Q2", abs.get(0).getName());
    }

    @Test
    public void makeErgodicRedirectsTheAbsorbingRow() {
        Network model = reducible();
        RoutingMatrix P = model.makeErgodic();
        JobClass c = model.getClasses().get(0);
        Matrix M = P.get(c, c);
        // MATLAB returns [0 1 0; 0 0 1; 1 0 0]
        double[][] want = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(want[i][j], M.get(i, j), 0.0, "P(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void anErgodicRoutingIsLeftAlone() {
        Network model = ergodic();
        assertTrue(model.isRoutingErgodic().isErgodic);
        RoutingErgodicity.ReducibilityInfo info = model.getReducibilityInfo();
        assertTrue(info.isRoutingErgodic);
        assertTrue(info.suggestedFixes.isEmpty());
        assertTrue(model.getAbsorbingStations().isEmpty());
    }
}
