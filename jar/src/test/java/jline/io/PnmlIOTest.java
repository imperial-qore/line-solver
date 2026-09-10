package jline.io;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Path;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * PNML place/transition import and export.
 *
 * The oracle is not a recorded file. A round trip through {@link PnmlIO} must leave every
 * CTMC metric of the model unchanged, which a consistent error on both sides of the
 * interchange would not survive; and the untimed net read from another tool is checked
 * against a marginal computed by hand from the chain it defines.
 *
 * Mirrors matlab/scratch-level checks of pnml_save.m / pnml_load.m and the native-Python
 * python/tests/test_pnml_io.py.
 */
public class PnmlIOTest {

    /** Two places, one timed mode each, arc weight 2 on one side. */
    private static Network buildTwoModes() {
        Network model = new Network("twomodes");
        Place p1 = new Place(model, "P1");
        Place p2 = new Place(model, "P2");
        Transition t1 = new Transition(model, "T1");
        Transition t2 = new Transition(model, "T2");
        ClosedClass jc = new ClosedClass(model, "Class1", 4, p1, 0);

        Mode m1 = t1.addMode("Mode1");
        t1.setDistribution(m1, new Exp(2));
        t1.setEnablingConditions(m1, jc, p1, 2);
        t1.setFiringOutcome(m1, jc, p2, 2);

        Mode m2 = t2.addMode("Mode2");
        t2.setDistribution(m2, new Erlang(1.5, 2));
        t2.setEnablingConditions(m2, jc, p2, 1);
        t2.setFiringOutcome(m2, jc, p1, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, p1, t1, 1.0);
        R.set(jc, jc, p2, t2, 1.0);
        R.set(jc, jc, t1, p2, 1.0);
        R.set(jc, jc, t2, p1, 1.0);
        model.link(R);
        p1.setMarking(4);
        p2.setMarking(0);
        return model;
    }

    /** One transition carrying TWO modes, plus an immediate mode with an inhibitor arc. */
    private static Network buildMultiMode() {
        Network model = new Network("multimode");
        Place q1 = new Place(model, "Q1");
        Place q2 = new Place(model, "Q2");
        Transition u1 = new Transition(model, "U1");
        Transition u2 = new Transition(model, "U2");
        ClosedClass jc = new ClosedClass(model, "Class1", 3, q1, 0);

        Mode fast = u1.addMode("Fast");
        u1.setDistribution(fast, new Exp(3));
        u1.setEnablingConditions(fast, jc, q1, 1);
        u1.setFiringOutcome(fast, jc, q2, 1);
        Mode slow = u1.addMode("Slow");
        u1.setDistribution(slow, new Exp(1));
        u1.setEnablingConditions(slow, jc, q1, 2);
        u1.setFiringOutcome(slow, jc, q2, 2);

        Mode back = u2.addMode("Back");
        u2.setTimingStrategy(back, TimingStrategy.IMMEDIATE);
        u2.setFiringWeights(back, 2.5);
        u2.setFiringPriorities(back, 1);
        u2.setEnablingConditions(back, jc, q2, 1);
        u2.setInhibitingConditions(back, jc, q1, 3);
        u2.setFiringOutcome(back, jc, q1, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, q1, u1, 1.0);
        R.set(jc, jc, u1, q2, 1.0);
        R.set(jc, jc, q2, u2, 1.0);
        R.set(jc, jc, q1, u2, 1.0);
        R.set(jc, jc, u2, q1, 1.0);
        model.link(R);
        q1.setMarking(3);
        q2.setMarking(0);
        return model;
    }

    private static void assertSameMetrics(Network a, Network b, int cutoff) {
        NetworkAvgTable ta = new SolverCTMC(a, "cutoff", cutoff, "verbose", VerboseLevel.SILENT).getAvgTable();
        NetworkAvgTable tb = new SolverCTMC(b, "cutoff", cutoff, "verbose", VerboseLevel.SILENT).getAvgTable();
        List<String> na = ta.getStationNames();
        List<String> nb = tb.getStationNames();
        assertEquals(na.size(), nb.size(), "station count");
        for (int i = 0; i < na.size(); i++) {
            int j = nb.indexOf(na.get(i));
            assertTrue(j >= 0, "station " + na.get(i) + " missing after the round trip");
            assertEquals(ta.getQLen().get(i), tb.getQLen().get(j), 1e-12, "QLen at " + na.get(i));
            assertEquals(ta.getUtil().get(i), tb.getUtil().get(j), 1e-12, "Util at " + na.get(i));
            assertEquals(ta.getRespT().get(i), tb.getRespT().get(j), 1e-12, "RespT at " + na.get(i));
            assertEquals(ta.getTput().get(i), tb.getTput().get(j), 1e-12, "Tput at " + na.get(i));
        }
    }

    @Test
    public void testRoundTripKeepsEveryMetric(@TempDir Path dir) throws IOException {
        Network model = buildTwoModes();
        String f = dir.resolve("twomodes.pnml").toString();
        PnmlIO.save(model, f);
        Network back = PnmlIO.load(f);
        assertSameMetrics(model, back, 4);
    }

    @Test
    public void testRoundTripRegroupsModesAndKeepsInhibitors(@TempDir Path dir) throws IOException {
        Network model = buildMultiMode();
        String f = dir.resolve("multimode.pnml").toString();
        PnmlIO.save(model, f);

        // Two modes of one transition are written as two PNML transitions and must be
        // regrouped into ONE LINE transition by the reader, not left as two.
        String text = new String(Files.readAllBytes(new File(f).toPath()), Charset.forName("UTF-8"));
        assertTrue(text.contains("id=\"U1.Fast\""), "the first mode of U1 should be its own PNML transition");
        assertTrue(text.contains("id=\"U1.Slow\""), "the second mode of U1 should be its own PNML transition");
        assertTrue(text.contains("<type value=\"inhibitor\"/>"), "the inhibitor arc should be written");

        Network back = PnmlIO.load(f);
        int transitions = 0;
        Transition u1 = null;
        Transition u2 = null;
        Place q1 = null;
        for (int i = 0; i < back.getNodes().size(); i++) {
            jline.lang.nodes.Node nd = back.getNodes().get(i);
            if (nd instanceof Transition) {
                transitions++;
                Transition tr = (Transition) nd;
                if ("U1".equals(tr.getName())) {
                    u1 = tr;
                }
                if ("U2".equals(tr.getName())) {
                    u2 = tr;
                }
            } else if (nd instanceof Place && "Q1".equals(nd.getName())) {
                q1 = (Place) nd;
            }
        }
        assertEquals(2, transitions, "the two modes of U1 must regroup into one transition");
        assertEquals(2, u1.getNumberOfModes(), "U1 keeps both of its modes");
        Mode back0 = u2.getModes().get(0);
        assertEquals(TimingStrategy.IMMEDIATE, u2.timingStrategies.get(back0), "U2 stays immediate");
        assertEquals(2.5, u2.firingWeights.get(0, 0), 1e-12, "the firing weight survives the file");
        assertEquals(3.0, u2.inhibitingConditions.get(back0).get(back.getNodeIndex(q1), 0), 1e-12,
                "the inhibitor threshold survives the file");

        assertSameMetrics(model, back, 3);
    }

    /**
     * A field lost in the round trip changes the second file. Comparing metrics can only
     * see what the solver reads; comparing the FILE sees every attribute the writer emits,
     * so this is the stronger oracle of the two and needs no solver at all.
     */
    private static void assertWriteReadWriteIsIdentical(Network model, Path dir, String stem) throws IOException {
        String first = dir.resolve(stem + "-1.pnml").toString();
        String second = dir.resolve(stem + "-2.pnml").toString();
        PnmlIO.save(model, first);
        PnmlIO.save(PnmlIO.load(first), second);
        String a = new String(Files.readAllBytes(new File(first).toPath()), Charset.forName("UTF-8"));
        String b = new String(Files.readAllBytes(new File(second).toPath()), Charset.forName("UTF-8"));
        assertEquals(a, b, "the second write of " + stem + " differs, so the read lost something");
    }

    @Test
    public void testWriteReadWriteIsByteIdentical(@TempDir Path dir) throws IOException {
        assertWriteReadWriteIsIdentical(buildTwoModes(), dir, "twomodes");
        assertWriteReadWriteIsIdentical(buildMultiMode(), dir, "multimode");
    }

    @Test
    public void testUntimedNetFromAnotherToolReadsAndSolves(@TempDir Path dir) throws IOException {
        // No toolspecific block and no timing: every transition is read as Exp(1) with one
        // server. Two indistinguishable tokens then cycle between two places, so the
        // number in `busy` is a 3-state birth-death chain with equal rates and the
        // stationary distribution is uniform: mean queue length 1 at each place.
        String xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                + "<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">\n"
                + "  <net id=\"foreign\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">\n"
                + "    <page id=\"p0\">\n"
                + "      <place id=\"ready\"><initialMarking><text>2</text></initialMarking></place>\n"
                + "      <place id=\"busy\"><initialMarking><text>0</text></initialMarking></place>\n"
                + "      <transition id=\"start\"/>\n"
                + "      <transition id=\"finish\"/>\n"
                + "      <arc id=\"e1\" source=\"ready\" target=\"start\"/>\n"
                + "      <arc id=\"e2\" source=\"start\" target=\"busy\"/>\n"
                + "      <arc id=\"e3\" source=\"busy\" target=\"finish\"/>\n"
                + "      <arc id=\"e4\" source=\"finish\" target=\"ready\"/>\n"
                + "    </page>\n"
                + "  </net>\n"
                + "</pnml>\n";
        File f = dir.resolve("foreign.pnml").toFile();
        Files.write(f.toPath(), xml.getBytes(Charset.forName("UTF-8")));

        Network model = PnmlIO.load(f.toString());
        NetworkAvgTable t = new SolverCTMC(model, "cutoff", 2, "verbose", VerboseLevel.SILENT).getAvgTable();
        List<String> names = t.getStationNames();
        assertEquals(1.0, t.getQLen().get(names.indexOf("ready")), 1e-9, "mean tokens in ready");
        assertEquals(1.0, t.getQLen().get(names.indexOf("busy")), 1e-9, "mean tokens in busy");
        assertEquals(2.0 / 3.0, t.getTput().get(names.indexOf("busy")), 1e-9, "firing rate of the cycle");
    }

    @Test
    public void testOpenNetIsRefused(@TempDir Path dir) {
        Network model = new Network("open");
        Source source = new Source(model, "Source");
        Place p = new Place(model, "P1");
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "Class1", 0);
        String f = dir.resolve("bad.pnml").toString();
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() throws Throwable {
                PnmlIO.save(model, f);
            }
        });
        assertTrue(e.getMessage().contains("unbounded token source"), "the refusal should name the reason: " + e.getMessage());
    }
}
