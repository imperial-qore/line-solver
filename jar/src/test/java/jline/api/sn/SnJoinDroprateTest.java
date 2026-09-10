package jline.api.sn;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jline.examples.java.models.Gallery;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Node;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

/**
 * A Join is the one station where LossRate = ArvR - Tput does not hold: ArvR counts the SIBLINGS
 * offered and Tput the PARENT jobs released. Reading the identity there charges (N-1)/N of the
 * offered traffic as lost at every join, standard joins included; these tests pin the rule that
 * replaces it.
 */
public class SnJoinDroprateTest {

    private static final double TOL = 1e-9;

    private static Node joinNode(NetworkStruct sn) {
        for (int i = 0; i < sn.nodes.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Join) {
                return sn.nodes.get(i);
            }
        }
        throw new IllegalStateException("no Join node");
    }

    private static int joinStation(NetworkStruct sn) {
        for (int i = 0; i < sn.nodes.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Join) {
                return (int) sn.nodeToStation.get(i);  // 1 x nnodes row vector
            }
        }
        throw new IllegalStateException("no Join node");
    }

    @Test
    public void siblingsAreCountedAtTheFork() {
        NetworkStruct snq = Gallery.gallery_fj_quorum().getStruct();
        assertEquals(3, SnJoinSiblings.snJoinSiblings(snq, joinNode(snq)));
        NetworkStruct sns = Gallery.gallery_fj_closed().getStruct();
        assertEquals(2, SnJoinSiblings.snJoinSiblings(sns, joinNode(sns)));
    }

    @Test
    public void standardJoinLosesNothing() {
        NetworkStruct sn = Gallery.gallery_fj_closed().getStruct();
        int ist = joinStation(sn);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        TN.fill(0.0);
        AN.fill(0.0);
        TN.set(ist, 0, 0.75);
        AN.set(ist, 0, 2 * 0.75);
        Matrix d = SnJoinDroprate.snJoinDroprate(sn, TN, AN);
        assertEquals(0.0, d.get(ist, 0), TOL);
    }

    @Test
    public void quorumDiscardsTheStragglers() {
        NetworkStruct sn = Gallery.gallery_fj_quorum().getStruct();
        int ist = joinStation(sn);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        TN.fill(0.0);
        AN.fill(0.0);
        double x = 1.6;
        TN.set(ist, 0, x);
        AN.set(ist, 0, 3 * x);
        Matrix d = SnJoinDroprate.snJoinDroprate(sn, TN, AN);
        assertEquals((3 - 2) * x, d.get(ist, 0), TOL);
    }

    @Test
    public void droprateIsZeroAwayFromTheJoin() {
        NetworkStruct sn = Gallery.gallery_fj_quorum().getStruct();
        int ist = joinStation(sn);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        TN.fill(1.6);
        AN.fill(4.8);
        Matrix d = SnJoinDroprate.snJoinDroprate(sn, TN, AN);
        for (int i = 0; i < sn.nstations; i++) {
            if (i != ist) {
                assertEquals(0.0, d.get(i, 0), TOL);
            }
        }
    }

    @Test
    public void droprateNeverNegative() {
        NetworkStruct sn = Gallery.gallery_fj_quorum().getStruct();
        int ist = joinStation(sn);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        TN.fill(0.0);
        AN.fill(0.0);
        TN.set(ist, 0, 5.0);  // inconsistent with AN on purpose
        AN.set(ist, 0, 1.0);
        Matrix d = SnJoinDroprate.snJoinDroprate(sn, TN, AN);
        assertTrue(d.get(ist, 0) >= 0.0);
    }

    @Test
    public void quorumChainCapacityIsUnbounded() {
        // A quorum join leaves stragglers in flight past their parent's next fork, so
        // a branch station is NOT bounded by the class population.
        // The JAR spells "unbounded" as Integer.MAX_VALUE, which is the sentinel
        // refreshCapacity itself tests against (chainCap >= Integer.MAX_VALUE).
        NetworkStruct snq = Gallery.gallery_fj_quorum().getStruct();
        for (int i = 0; i < snq.classcap.getNumRows(); i++) {
            double c = snq.classcap.get(i, 0);
            assertTrue(Double.isInfinite(c) || c >= Integer.MAX_VALUE,
                    "classcap row " + i + " should be unbounded under a quorum join, got " + c);
        }
        NetworkStruct sns = Gallery.gallery_fj_closed().getStruct();
        for (int i = 0; i < sns.classcap.getNumRows(); i++) {
            double c = sns.classcap.get(i, 0);
            assertTrue(Double.isFinite(c) && c < Integer.MAX_VALUE,
                    "classcap row " + i + " should be the class population under a standard join,"
                            + " got " + c);
        }
    }

    @Test
    public void lossTableJoinRowIsLosslessUnderAStandardJoin() {
        // The old ArvR - Tput rule reported LossRatio = 1/2 here. Nothing is lost.
        Network model = Gallery.gallery_fj_closed();
        NetworkStruct sn = model.getStruct();
        int ist = joinStation(sn);
        SolverMVA solver = new SolverMVA(model);
        solver.getAvg();
        Matrix d = solver.result.DropRateJoin;
        assertTrue(d != null, "DropRateJoin should be populated for a fork-join model");
        assertEquals(0.0, d.get(ist, 0), 1e-6);
    }

    @Test
    public void lossTableJoinRowUnderAQuorum() {
        // One of the three siblings is discarded, so the drop rate is exactly X.
        Network model = Gallery.gallery_fj_quorum();
        NetworkStruct sn = model.getStruct();
        int ist = joinStation(sn);
        SolverNC solver = new SolverNC(model);
        solver.getAvg();
        Matrix d = solver.result.DropRateJoin;
        assertTrue(d != null, "DropRateJoin should be populated for a fork-join model");
        assertEquals(solver.result.TN.get(ist, 0), d.get(ist, 0), 1e-6);
    }
}
