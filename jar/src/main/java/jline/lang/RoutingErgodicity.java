/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lang;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;

/**
 * The repair side of routing reducibility: what to do once
 * {@link Network#isRoutingErgodic()} has said the routing is reducible.
 *
 * <p>Port of the MATLAB {@code @MNetwork} methods {@code getReducibilityInfo},
 * {@code getAbsorbingStations} and {@code makeErgodic}. The DETECTION half,
 * {@code isRoutingErgodic}, was already in {@link Network} and is reused as is
 * -- there is one adjacency build and one SCC decomposition, not two.
 *
 * <p>The whole family inspects the ROUTING STRUCTURE ONLY, not the state space.
 * Ergodic routing is NECESSARY but not sufficient for the chain to be ergodic:
 * state-dependent routing, finite buffers and blocking can still make the CTMC
 * reducible. The converse is the useful direction -- reducible routing is a
 * defect the modeller can see and fix before any solver runs, and it is what
 * otherwise surfaces as "the generator has no recurrent state" inside
 * ctmc_solve.
 */
public final class RoutingErgodicity {

    private RoutingErgodicity() {
    }

    /**
     * {@link Network.RoutingErgodicityResult} with the repairs attached, the
     * MATLAB {@code getReducibilityInfo} struct.
     */
    public static final class ReducibilityInfo {
        /** True when the routing is ergodic (irreducible). */
        public boolean isRoutingErgodic = true;
        /** True when the routing creates a reducible structure. */
        public boolean isReducible = false;
        /** Names of the absorbing stations, the Sink excluded. */
        public List<String> absorbingStations = new ArrayList<String>();
        /** Names of the stations lying in a transient component. */
        public List<String> transientStations = new ArrayList<String>();
        /** Number of strongly connected components of the routing graph. */
        public int numSCCs = 1;
        /** One suggested repair per absorbing station; empty when ergodic. */
        public List<String> suggestedFixes = new ArrayList<String>();
    }

    /** The reducibility structure plus a suggested repair per absorbing station. */
    public static ReducibilityInfo getReducibilityInfo(Network model) {
        Network.RoutingErgodicityResult erg = model.isRoutingErgodic();
        ReducibilityInfo info = new ReducibilityInfo();
        info.isRoutingErgodic = erg.isErgodic;
        info.isReducible = erg.isReducible;
        info.absorbingStations = new ArrayList<String>(erg.absorbingStations);
        info.transientStations = new ArrayList<String>(erg.transientStations);
        info.numSCCs = erg.numSCCs;
        if (erg.isErgodic) {
            return info;
        }
        String returnTarget = defaultTargetName(model, info.absorbingStations);
        for (int i = 0; i < info.absorbingStations.size(); i++) {
            String abs = info.absorbingStations.get(i);
            if (returnTarget != null && !returnTarget.equals(abs)) {
                info.suggestedFixes.add("Route jobs from " + abs + " back to " + returnTarget
                        + " (e.g., P{class}(" + abs + ", " + returnTarget + ") = 1.0)");
            }
        }
        return info;
    }

    /** The absorbing stations as model objects, the Sink excluded. */
    public static List<Station> getAbsorbingStations(Network model) {
        Network.RoutingErgodicityResult erg = model.isRoutingErgodic();
        List<Station> out = new ArrayList<Station>();
        for (int i = 0; i < erg.absorbingStations.size(); i++) {
            Station st = model.getStationByName(erg.absorbingStations.get(i));
            if (st != null) {
                out.add(st);
            }
        }
        return out;
    }

    /**
     * A routing matrix that makes the network ergodic, by redirecting every
     * absorbing station to TARGETNAME.
     *
     * <p>Does NOT relink: the caller applies it with {@code model.link(P)}. That
     * is deliberate -- relinking rebuilds the struct, and the modeller should
     * see the repair before it is applied.
     *
     * @param targetName node to route absorbing stations to, or null to take the
     *                   first Delay, else the first non-absorbing station
     */
    public static RoutingMatrix makeErgodic(Network model, String targetName) {
        Network.RoutingErgodicityResult erg = model.isRoutingErgodic();
        RoutingMatrix P = seedFromCurrentRouting(model);
        if (erg.isErgodic) {
            line_warning("makeErgodic", "Routing is already ergodic. No changes needed.");
            return P;
        }
        if (erg.absorbingStations.isEmpty()) {
            line_warning("makeErgodic", "No absorbing stations found, but the network is not "
                    + "ergodic. Manual intervention is required.");
            return P;
        }

        String target = (targetName == null || targetName.isEmpty())
                ? defaultTargetName(model, erg.absorbingStations) : targetName;
        if (target == null) {
            line_error("makeErgodic", "Cannot find a suitable target node for routing.");
        }
        Node targetNode = model.getNodeByName(target);
        if (targetNode == null) {
            line_error("makeErgodic", "Target node \"" + target + "\" not found in the network.");
        }
        int targetIdx = model.getNodeIndex(targetNode);

        List<JobClass> classes = model.getClasses();
        for (int i = 0; i < erg.absorbingStations.size(); i++) {
            String absName = erg.absorbingStations.get(i);
            Node absNode = model.getNodeByName(absName);
            if (absNode == null) {
                continue;
            }
            int absIdx = model.getNodeIndex(absNode);
            if (absIdx == targetIdx) {
                line_warning("makeErgodic", "Cannot route %s to itself. Skipping.", absName);
                continue;
            }
            for (int r = 0; r < classes.size(); r++) {
                for (int s = 0; s < classes.size(); s++) {
                    Matrix blk = P.get(classes.get(r), classes.get(s));
                    if (blk == null || absIdx >= blk.getNumRows()) {
                        continue;
                    }
                    for (int j = 0; j < blk.getNumCols(); j++) {
                        blk.set(absIdx, j, 0.0);
                    }
                }
                Matrix diag = P.get(classes.get(r), classes.get(r));
                if (diag != null && absIdx < diag.getNumRows() && targetIdx < diag.getNumCols()) {
                    diag.set(absIdx, targetIdx, 1.0);
                }
            }
        }
        return P;
    }

    /** A RoutingMatrix carrying the model's current linked routing. */
    private static RoutingMatrix seedFromCurrentRouting(Network model) {
        RoutingMatrix P = model.initRoutingMatrix();
        Map<JobClass, Map<JobClass, Matrix>> cur = model.getLinkedRoutingMatrix();
        if (cur == null) {
            return P;
        }
        List<JobClass> classes = model.getClasses();
        for (int r = 0; r < classes.size(); r++) {
            Map<JobClass, Matrix> row = cur.get(classes.get(r));
            if (row == null) {
                continue;
            }
            for (int s = 0; s < classes.size(); s++) {
                Matrix blk = row.get(classes.get(s));
                if (blk != null) {
                    P.set(classes.get(r), classes.get(s), blk.copy());
                }
            }
        }
        return P;
    }

    /** First Delay, else the first station that is not itself absorbing. */
    private static String defaultTargetName(Network model, List<String> absorbing) {
        List<Node> nodes = model.getNodes();
        for (int i = 0; i < nodes.size(); i++) {
            if (nodes.get(i) instanceof Delay) {
                return nodes.get(i).getName();
            }
        }
        String first = null;
        for (int i = 0; i < nodes.size(); i++) {
            if (!(nodes.get(i) instanceof Station)) {
                continue;
            }
            String nm = nodes.get(i).getName();
            if (first == null) {
                first = nm;
            }
            if (!absorbing.contains(nm)) {
                return nm;
            }
        }
        return first;
    }
}
