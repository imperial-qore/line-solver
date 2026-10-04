/**
 * Fork-Join Visit Ratio Calculator via Auxiliary SPN Models
 *
 * Computes fork-join node visit ratios by building, for each class that passes
 * through a fork-join pair, an auxiliary closed Stochastic Petri Net capturing
 * the fork/join synchronization semantics.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.io.InputOutput;
import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.RoutingMatrix;
import jline.lang.constant.NodeType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

public final class SnFjVisitsSpn {
    private SnFjVisitsSpn() {}

    /**
     * Compute fork-join node visit ratios via auxiliary SPN models.
     *
     * The net is solved through SolverCTMC, as in MATLAB. Its throughput ratios
     * are NOT uniform across the station Places: the pre-fork transition
     * consumes all B tokens at once and the Join returns them, so a station
     * inside a fork-join region fires once per B firings of the cycle. After the
     * normalisation on the reference station, a station outside the region
     * carries 1, a station inside it carries 1/B, and a Fork or a Join, which
     * holds no Place, carries 0. The Python and C++ ports evaluate that rule
     * analytically instead of enumerating the state space.
     */
    public static List<Matrix> snFjVisitsSpn(NetworkStruct sn) {
        int I = sn.nnodes;
        int K = sn.nclasses;
        int nchains = sn.nchains;
        Map<Integer, Matrix> inchain = sn.inchain;
        Matrix refstat = sn.refstat;

        List<Matrix> nodevisits = new ArrayList<Matrix>();
        for (int c = 0; c < nchains; c++) {
            Matrix m = new Matrix(I, K);
            m.zero();
            nodevisits.add(m);
        }

        // Quick check: if no fork-join pairs, return zeros
        boolean hasFJ = false;
        for (int i = 0; i < sn.fj.getNumRows() && !hasFJ; i++) {
            for (int j = 0; j < sn.fj.getNumCols(); j++) {
                if (sn.fj.get(i, j) > 0) {
                    hasFJ = true;
                    break;
                }
            }
        }
        if (!hasFJ) {
            return nodevisits;
        }

        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            if (inchain_c == null) continue;
            for (int kidx = 0; kidx < inchain_c.getNumCols(); kidx++) {
                int r = (int) inchain_c.get(0, kidx);

                double[][] P_r = new double[I][I];
                for (int i = 0; i < I; i++) {
                    for (int j = 0; j < I; j++) {
                        P_r[i][j] = sn.rtnodes.get(i * K + r, j * K + r);
                    }
                }

                int refnode = (int) sn.stationToNode.get((int) refstat.get(r));
                boolean[] visited = new boolean[I];
                visited[refnode] = true;
                boolean changed = true;
                while (changed) {
                    changed = false;
                    for (int i = 0; i < I; i++) {
                        if (visited[i]) {
                            for (int j = 0; j < I; j++) {
                                if (P_r[i][j] > 0 && !visited[j]) {
                                    visited[j] = true;
                                    changed = true;
                                }
                            }
                        }
                    }
                }

                boolean hasFork = false;
                for (int i = 0; i < I; i++) {
                    if (visited[i] && sn.nodetype.get(i) == NodeType.Fork) {
                        hasFork = true;
                        break;
                    }
                }
                if (!hasFork) {
                    for (int i = 0; i < I; i++) {
                        if (visited[i]) {
                            nodevisits.get(c).set(i, r, 1.0);
                        }
                    }
                    continue;
                }

                double[] visits_r = buildAndSolveSpn(sn, P_r, visited, r, refnode);
                for (int i = 0; i < I; i++) {
                    nodevisits.get(c).set(i, r, visits_r[i]);
                }
            }

            int refclass0 = (int) inchain_c.get(0, 0);
            int refnode_c = (int) sn.stationToNode.get((int) refstat.get(refclass0));
            for (int kidx = 0; kidx < inchain_c.getNumCols(); kidx++) {
                int r = (int) inchain_c.get(0, kidx);
                double normVal = nodevisits.get(c).get(refnode_c, r);
                if (normVal > GlobalConstants.FineTol) {
                    for (int i = 0; i < I; i++) {
                        nodevisits.get(c).set(i, r, nodevisits.get(c).get(i, r) / normVal);
                    }
                }
            }
        }

        return nodevisits;
    }

    private static List<Integer> resolveForkDests(NetworkStruct sn, double[][] P_r, boolean[] visited, int forkNd) {
        return resolveForkDests(sn, P_r, visited, forkNd, 0, 1.0, null);
    }

    /**
     * Recursively resolve a Fork's destinations down to station nodes.
     *
     * <p>When {@code weights} is non-null it collects the EXPECTED number of
     * tasks each leaf receives per firing of the outermost fork: P(branch fires)
     * times E[tasks on that link], multiplied down through any nesting. A plain
     * fork gives every link exactly 1, so the weighted leaf count collapses to
     * the leaf count it always was.</p>
     */
    private static List<Integer> resolveForkDests(NetworkStruct sn, double[][] P_r, boolean[] visited,
                                                  int forkNd, int r, double w, List<Double> weights) {
        List<Integer> stDests = new ArrayList<Integer>();
        ForkNodeParam fp = null;
        Object po = (sn.nodeparam != null) ? sn.nodeparam.get(sn.nodes.get(forkNd)) : null;
        if (po instanceof ForkNodeParam) fp = (ForkNodeParam) po;
        boolean hasFan = fp != null && fp.fanOutLink != null;
        for (int bd = 0; bd < sn.nnodes; bd++) {
            if (P_r[forkNd][bd] > 0 && visited[bd]) {
                double wbd = w;
                if (hasFan) wbd = w * fp.fanOutProb.get(bd, r) * fp.fanOutLink.get(bd, r);
                if (sn.nodetype.get(bd) == NodeType.Fork) {
                    stDests.addAll(resolveForkDests(sn, P_r, visited, bd, r, wbd, weights));
                } else if (sn.isstation.get(bd, 0) > 0) {
                    stDests.add(bd);
                    if (weights != null) weights.add(Double.valueOf(wbd));
                }
            }
        }
        return stDests;
    }

    private static double[] buildAndSolveSpn(NetworkStruct sn, double[][] P_r, boolean[] visited, int r, int refnode) {
        int I = sn.nnodes;
        double[] visits_r = new double[I];

        List<Integer> visitedNodes = new ArrayList<Integer>();
        for (int i = 0; i < I; i++) {
            if (visited[i]) visitedNodes.add(i);
        }
        if (visitedNodes.isEmpty()) return visits_r;

        List<Integer> stationNodes = new ArrayList<Integer>();
        List<Integer> forkNodes = new ArrayList<Integer>();
        List<Integer> joinNodes = new ArrayList<Integer>();
        for (Integer nd : visitedNodes) {
            NodeType nt = sn.nodetype.get(nd);
            if (nt == NodeType.Fork) {
                forkNodes.add(nd);
            } else if (nt == NodeType.Join) {
                joinNodes.add(nd);
            } else if (sn.isstation.get(nd, 0) > 0 && nt != NodeType.Source && nt != NodeType.Sink) {
                stationNodes.add(nd);
            }
        }

        Map<Integer, Integer> joinLeaves = new HashMap<Integer, Integer>();
        for (int pass = 0; pass < joinNodes.size(); pass++) {
            for (Integer jnd : joinNodes) {
                int lc = 0;
                for (int srcNd = 0; srcNd < I; srcNd++) {
                    if (P_r[srcNd][jnd] > 0 && visited[srcNd]) {
                        if (sn.nodetype.get(srcNd) == NodeType.Join && joinLeaves.containsKey(srcNd)) {
                            lc += joinLeaves.get(srcNd);
                        } else {
                            lc += 1;
                        }
                    }
                }
                joinLeaves.put(jnd, lc);
            }
        }

        double B = 0.0;
        for (Integer fnd : forkNodes) {
            boolean isOutermost = true;
            for (int srcNd = 0; srcNd < I; srcNd++) {
                if (P_r[srcNd][fnd] > 0 && visited[srcNd]) {
                    if (sn.nodetype.get(srcNd) == NodeType.Join) {
                        isOutermost = false;
                        break;
                    }
                }
            }
            if (isOutermost) {
                List<Double> leafw = new ArrayList<Double>();
                List<Integer> leafs = resolveForkDests(sn, P_r, visited, fnd, r, 1.0, leafw);
                if (leafs.isEmpty()) {
                    line_error(mfilename(new Object() {}),
                            "A Fork reaches no station on any branch, so the auxiliary net has "
                            + "nothing to synchronize.");
                }
                // The EXPECTED sibling count, which is the leaf count exactly
                // when every link is certain and carries one task.
                double expected = 0.0;
                for (Double d : leafw) expected += d.doubleValue();
                if (expected <= 0.0) {
                    line_error(mfilename(new Object() {}),
                            "A Fork emits no task in expectation, so its Join can never fire; at "
                            + "least one branch must be certain to emit at least one task.");
                }
                if (expected > B) B = expected;
            }
        }
        if (B == 0.0) B = 1.0;

        // THE NET BELOW IS THE PLAIN FORK'S NET. It puts B tokens in a closed
        // SPN whose pre-fork transition emits ONE token per branch, so it is
        // balanced only when every branch carries exactly one task: give it a
        // fork whose links carry three each, and it deadlocks and reports zero
        // throughput everywhere. Under a variable forking level B is the
        // expected sibling count, which is generally not an integer either --
        // and a fractional token population is not a net anyone can enumerate.
        // The rule that solve produces is available in closed form and is what
        // answers whenever the fork varies. The plain fork keeps taking the
        // solve, which is what holds the two in step.
        boolean isVariableFork = false;
        for (Integer fnd : forkNodes) {
            Object po = (sn.nodeparam != null) ? sn.nodeparam.get(sn.nodes.get(fnd)) : null;
            if (!(po instanceof ForkNodeParam)) continue;
            ForkNodeParam fp = (ForkNodeParam) po;
            if (fp.fanOutLink == null) continue;
            for (int k = 0; k < fp.fanOutLink.getNumRows() && !isVariableFork; k++)
                for (int rr = 0; rr < fp.fanOutLink.getNumCols() && !isVariableFork; rr++) {
                    if (fp.fanOutProb.get(k, rr) == 0.0) continue;  // link not taken
                    if (fp.fanOutProb.get(k, rr) != 1.0) isVariableFork = true;
                    if (fp.fanOutLink.get(k, rr) != fp.fanOut) isVariableFork = true;
                    if (fp.fanOutDist != null && fp.fanOutDist[k][rr] != null) isVariableFork = true;
                }
        }
        if (isVariableFork || Math.abs(B - Math.round(B)) > GlobalConstants.FineTol) {
            return closedFormVisits(sn, P_r, visited, forkNodes, stationNodes, B);
        }
        int Bint = (int) Math.round(B);

        Network model = new Network("fj_spn_aux");

        Place[] places = new Place[I];
        for (Integer nd : stationNodes) {
            places[nd] = new Place(model, "P_" + sn.nodenames.get(nd));
        }

        Place[] preJoin = new Place[I];
        for (Integer nd : stationNodes) {
            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[nd][dstNd] > 0 && visited[dstNd] && sn.nodetype.get(dstNd) == NodeType.Join) {
                    preJoin[nd] = new Place(model, "P_" + sn.nodenames.get(nd) + "_done");
                    break;
                }
            }
        }

        for (Integer jnd : joinNodes) {
            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[jnd][dstNd] > 0 && visited[dstNd] && sn.nodetype.get(dstNd) == NodeType.Join) {
                    preJoin[jnd] = new Place(model, "P_" + sn.nodenames.get(jnd) + "_done");
                    break;
                }
            }
        }

        Place[] interJF = new Place[I];
        for (Integer jnd : joinNodes) {
            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[jnd][dstNd] > 0 && visited[dstNd] && sn.nodetype.get(dstNd) == NodeType.Fork) {
                    interJF[jnd] = new Place(model, "P_" + sn.nodenames.get(jnd) + "_to_" + sn.nodenames.get(dstNd));
                    break;
                }
            }
        }

        ClosedClass jobclass = new ClosedClass(model, "Token", Bint, places[refnode]);

        List<Transition> transitions = new ArrayList<Transition>();
        List<List<Place>> transInfoInPlaces = new ArrayList<List<Place>>();
        List<List<Place>> transInfoOutPlaces = new ArrayList<List<Place>>();

        for (Integer nd : stationNodes) {
            int ist = (int) sn.nodeToStation.get(nd);

            boolean isTimed = false;
            double rate = 0.0;
            if (ist >= 0 && !sn.rates.isEmpty() && ist < sn.rates.getNumRows() && r < sn.rates.getNumCols()) {
                double rateVal = sn.rates.get(ist, r);
                if (rateVal > 0 && !Double.isInfinite(rateVal)) {
                    rate = rateVal;
                    isTimed = true;
                }
            }

            List<Place> outPlaces = new ArrayList<Place>();
            int enableCount = 1;
            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[nd][dstNd] <= 0 || !visited[dstNd]) continue;

                if (sn.nodetype.get(dstNd) == NodeType.Fork) {
                    List<Integer> forkDests = resolveForkDests(sn, P_r, visited, dstNd);
                    enableCount = Bint;
                    for (Integer fd : forkDests) {
                        if (places[fd] != null) outPlaces.add(places[fd]);
                    }
                } else if (sn.nodetype.get(dstNd) == NodeType.Join) {
                    if (preJoin[nd] != null) outPlaces.add(preJoin[nd]);
                } else if (sn.isstation.get(dstNd, 0) > 0) {
                    if (places[dstNd] != null) outPlaces.add(places[dstNd]);
                }
            }

            if (outPlaces.isEmpty()) continue;

            String tName = "T_svc_" + sn.nodenames.get(nd);
            Transition T = new Transition(model, tName);
            Mode mode = T.addMode("serve");

            if (isTimed) {
                T.setDistribution(mode, new Exp(rate));
            } else {
                T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
                T.setFiringWeights(mode, 1.0);
            }

            T.setEnablingConditions(mode, jobclass, places[nd], enableCount);
            T.setFiringOutcome(mode, jobclass, places[nd], -enableCount);
            for (Place op : outPlaces) {
                T.setFiringOutcome(mode, jobclass, op, 1);
            }

            transitions.add(T);
            List<Place> inList = new ArrayList<Place>();
            inList.add(places[nd]);
            transInfoInPlaces.add(inList);
            transInfoOutPlaces.add(outPlaces);
        }

        for (Integer jnd : joinNodes) {
            List<Place> inPlaces = new ArrayList<Place>();
            List<Integer> inSrcNodes = new ArrayList<Integer>();
            for (int srcNd = 0; srcNd < I; srcNd++) {
                if (P_r[srcNd][jnd] > 0 && visited[srcNd]) {
                    if (preJoin[srcNd] != null) {
                        inPlaces.add(preJoin[srcNd]);
                        inSrcNodes.add(srcNd);
                    }
                }
            }
            if (inPlaces.isEmpty()) continue;

            int[] enCounts = new int[inPlaces.size()];
            for (int ii = 0; ii < enCounts.length; ii++) enCounts[ii] = 1;
            for (int ii = 0; ii < inSrcNodes.size(); ii++) {
                int srcNd = inSrcNodes.get(ii);
                if (sn.nodetype.get(srcNd) == NodeType.Join && joinLeaves.containsKey(srcNd)) {
                    enCounts[ii] = joinLeaves.get(srcNd);
                }
            }

            int produceCount;
            if (joinLeaves.containsKey(jnd)) {
                produceCount = joinLeaves.get(jnd);
            } else {
                int sum = 0;
                for (int ec : enCounts) sum += ec;
                produceCount = sum;
            }

            List<Place> outPlaces = new ArrayList<Place>();
            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[jnd][dstNd] <= 0 || !visited[dstNd]) continue;

                if (sn.nodetype.get(dstNd) == NodeType.Fork) {
                    if (interJF[jnd] != null) outPlaces.add(interJF[jnd]);
                } else if (sn.nodetype.get(dstNd) == NodeType.Join) {
                    if (preJoin[jnd] != null) outPlaces.add(preJoin[jnd]);
                } else if (sn.isstation.get(dstNd, 0) > 0) {
                    if (places[dstNd] != null) outPlaces.add(places[dstNd]);
                }
            }
            if (outPlaces.isEmpty()) continue;

            String tName = "T_join_" + sn.nodenames.get(jnd);
            Transition T = new Transition(model, tName);
            Mode mode = T.addMode("sync");
            T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
            T.setFiringWeights(mode, 1.0);

            for (int ii = 0; ii < inPlaces.size(); ii++) {
                T.setEnablingConditions(mode, jobclass, inPlaces.get(ii), enCounts[ii]);
                T.setFiringOutcome(mode, jobclass, inPlaces.get(ii), -enCounts[ii]);
            }
            for (Place op : outPlaces) {
                T.setFiringOutcome(mode, jobclass, op, produceCount);
            }

            transitions.add(T);
            transInfoInPlaces.add(inPlaces);
            transInfoOutPlaces.add(outPlaces);
        }

        for (Integer jnd : joinNodes) {
            if (interJF[jnd] == null) continue;

            for (int dstNd = 0; dstNd < I; dstNd++) {
                if (P_r[jnd][dstNd] <= 0 || !visited[dstNd]) continue;
                if (sn.nodetype.get(dstNd) != NodeType.Fork) continue;

                List<Integer> forkDests = resolveForkDests(sn, P_r, visited, dstNd);
                if (forkDests.isEmpty()) continue;

                String tName = "T_fork_" + sn.nodenames.get(jnd) + "_" + sn.nodenames.get(dstNd);
                Transition T = new Transition(model, tName);
                Mode mode = T.addMode("fork");
                T.setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
                T.setFiringWeights(mode, 1.0);

                T.setEnablingConditions(mode, jobclass, interJF[jnd], Bint);
                T.setFiringOutcome(mode, jobclass, interJF[jnd], -Bint);
                List<Place> forkOut = new ArrayList<Place>();
                for (Integer fd : forkDests) {
                    if (places[fd] != null) {
                        T.setFiringOutcome(mode, jobclass, places[fd], 1);
                        forkOut.add(places[fd]);
                    }
                }

                transitions.add(T);
                List<Place> inList = new ArrayList<Place>();
                inList.add(interJF[jnd]);
                transInfoInPlaces.add(inList);
                transInfoOutPlaces.add(forkOut);
            }
        }

        RoutingMatrix R = model.initRoutingMatrix();
        for (int tidx = 0; tidx < transitions.size(); tidx++) {
            Transition T = transitions.get(tidx);
            List<Place> inPl = transInfoInPlaces.get(tidx);
            List<Place> outPl = transInfoOutPlaces.get(tidx);
            for (Place ip : inPl) {
                R.set(jobclass, jobclass, ip, T, 1.0);
            }
            for (Place op : outPl) {
                R.set(jobclass, jobclass, T, op, 1.0);
            }
        }
        model.link(R);

        places[refnode].setState(Bint);
        for (int i = 0; i < I; i++) {
            if (places[i] != null && i != refnode) {
                places[i].setState(0);
            }
            if (preJoin[i] != null) preJoin[i].setState(0);
            if (interJF[i] != null) interJF[i].setState(0);
        }

        try {
            SolverCTMC solver = new SolverCTMC(model);
            NetworkAvgTable avg = (NetworkAvgTable) solver.getAvgTable();

            List<String> stationNames = avg.getStationNames();
            List<Double> tputVals = avg.getTput();

            for (Integer nd : visitedNodes) {
                if (places[nd] != null) {
                    String placeName = places[nd].getName();
                    for (int row = 0; row < stationNames.size(); row++) {
                        if (stationNames.get(row).equals(placeName)) {
                            visits_r[nd] = tputVals.get(row);
                            break;
                        }
                    }
                }
            }
        } catch (Exception e) {
            InputOutput.line_warning("snFjVisitsSpn", "SPN CTMC solve failed for class " + r + ": " + e.getMessage());
        }

        return visits_r;
    }

    /**
     * The visit rule the SPN solve produces, evaluated directly.
     *
     * <p>Used when the circulating population B is FRACTIONAL, which a variable
     * forking level makes it: there is no net with a fractional token count to
     * build, so there is nothing to solve. A station outside the fork-join
     * region carries 1, one inside it carries 1/B, and a Fork or a Join, which
     * holds no Place, carries 0.</p>
     */
    private static double[] closedFormVisits(NetworkStruct sn, double[][] P_r, boolean[] visited,
                                             List<Integer> forkNodes, List<Integer> stationNodes,
                                             double B) {
        int I = sn.nnodes;
        double[] visits_r = new double[I];
        boolean[] inRegion = new boolean[I];
        for (Integer fnd : forkNodes) {
            boolean outermost = true;
            for (int srcNd = 0; srcNd < I && outermost; srcNd++)
                if (P_r[srcNd][fnd] > 0 && visited[srcNd]
                        && sn.nodetype.get(srcNd) == NodeType.Join) outermost = false;
            if (!outermost) continue;  // a serial stage, inside a region already fixed
            List<Integer> frontier = new ArrayList<Integer>();
            frontier.add(fnd);
            for (int h = 0; h < frontier.size(); h++) {
                int nd = frontier.get(h).intValue();
                for (int j = 0; j < I; j++) {
                    if (!(P_r[nd][j] > 0) || !visited[j]) continue;
                    if (sn.nodetype.get(j) == NodeType.Join || inRegion[j]) continue;
                    inRegion[j] = true;
                    frontier.add(Integer.valueOf(j));
                }
            }
        }
        for (Integer nd : stationNodes) {
            visits_r[nd.intValue()] = inRegion[nd.intValue()] ? (1.0 / B) : 1.0;
        }
        return visits_r;
    }
}
