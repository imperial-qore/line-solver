/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.GlobalConstants;
import jline.io.CodeGenSupport.Lang;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.Replayer;
import jline.lang.processes.Trace;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;

import static jline.io.CodeGenSupport.fmtD;
import static jline.io.CodeGenSupport.fmtF;
import static jline.io.CodeGenSupport.fmtInt;
import static jline.io.CodeGenSupport.mstr;
import static jline.io.InputOutput.line_error;

/**
 * Generates a MATLAB script that rebuilds a {@link Network} with the LINE MATLAB API.
 *
 * <p>Port of {@code matlab/src/io/QN2MATLAB.m}: the script defines {@code model}, the cell arrays
 * {@code node} and {@code jobclass}, and a routing cell array {@code P} read from {@code sn.rtnodes},
 * in the same three blocks and with the same comments as the MATLAB generator. Processes are written
 * from their first two moments by the same rules as {@link QN2JAVA}.</p>
 *
 * <p>Differences from the MATLAB text, each made so the script rebuilds the same model: a number that
 * MATLAB's {@code %f} or {@code %d} would round is printed in full; a scheduling weight other than 1
 * (DPS, GPS and their priority variants) is passed to {@code setService}, as {@code QN2JAVA.m} already
 * does; a Fork's tasks-per-link other than 1 is restored through {@code setTasksPerLink}. Node types
 * MATLAB would silently omit (Cache, Logger, Place, Transition) are refused by name instead.</p>
 */
public final class QN2MATLAB {

    private QN2MATLAB() {
    }

    /** Returns the generated script for a model named after {@code model.getName()}. */
    public static String generate(Network model) {
        return generate(model, model.getName());
    }

    /** Returns the generated script. */
    public static String generate(final Network model, final String modelName) {
        return CodeGenSupport.capture(new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                QN2MATLAB.emit(model, modelName, out);
            }
        });
    }

    /** Writes the generated script to a stream (e.g. {@code System.out}), which is flushed but not closed. */
    public static void write(Network model, String modelName, PrintStream out) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw);
        pw.flush();
    }

    /** Writes the generated script to a writer, which is flushed but not closed. */
    public static void write(Network model, String modelName, Writer out) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw);
        pw.flush();
    }

    /** Writes the generated script to {@code filename}, overwriting it. */
    public static void write(final Network model, final String modelName, String filename) {
        CodeGenSupport.toFile(filename, new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                QN2MATLAB.emit(model, modelName, out);
            }
        });
    }

    static void emit(Network model, String modelName, PrintWriter out) {
        final Lang M = Lang.MATLAB;
        if (modelName == null) {
            modelName = "myModel";
        }
        NetworkStruct sn = model.getStruct();
        int K = sn.nclasses;

        out.printf("model = Network('%s');\n", mstr(modelName));
        out.print("\n%% Block 1: nodes");
        out.print("\n");
        for (int i = 0; i < sn.nnodes; i++) {
            int n = i + 1;
            String name = mstr(sn.nodenames.get(i));
            NodeType nt = sn.nodetype.get(i);
            switch (nt) {
                case Source:
                    out.printf("node{%d} = Source(model, '%s');\n", n, name);
                    break;
                case Delay:
                    out.printf("node{%d} = DelayStation(model, '%s');\n", n, name);
                    break;
                case Queue: {
                    int ist = (int) sn.nodeToStation.get(i);
                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    out.printf("node{%d} = Queue(model, '%s', SchedStrategy.%s);\n", n, name, sched.name());
                    double c = sn.nservers.get(ist);
                    if (c > 1) {
                        out.printf("node{%d}.setNumServers(%s);\n", n, fmtInt(c, M));
                    }
                    break;
                }
                case Router:
                    out.printf("node{%d} = Router(model, '%s');\n", n, name);
                    break;
                case Fork: {
                    out.printf("node{%d} = Fork(model, '%s');\n", n, name);
                    double fanOut = QN2JAVA.forkFanOut(sn, i);
                    if (fanOut != 1.0) {
                        out.printf("node{%d}.setTasksPerLink(%s);\n", n, fmtInt(fanOut, M));
                    }
                    break;
                }
                case Join:
                    out.printf("node{%d} = Join(model, '%s', node{%d});\n", n, name, QN2JAVA.forkOf(sn, i) + 1);
                    break;
                case Sink:
                    out.printf("node{%d} = Sink(model, '%s');\n", n, name);
                    break;
                case ClassSwitch:
                    out.printf("node{%d} = Router(model, '%s'); %% Class switching is embedded in the routing matrix \n", n, name);
                    break;
                default:
                    line_error("QN2MATLAB", "QN2MATLAB cannot generate code for node " + sn.nodenames.get(i) + " of type " + nt + ".");
            }
        }

        out.print("\n%% Block 2: classes\n");
        for (int k = 0; k < K; k++) {
            String cname = mstr(sn.classnames.get(k));
            double njobs = sn.njobs.get(k);
            int prio = (int) sn.classprio.get(k);
            if (Double.isInfinite(njobs)) {
                out.printf("jobclass{%d} = OpenClass(model, '%s', %d);\n", k + 1, cname, prio);
            } else {
                out.printf("jobclass{%d} = ClosedClass(model, '%s', %s, node{%d}, %d);\n",
                        k + 1, cname, fmtD(njobs, M), QN2JAVA.referenceNode(sn, k) + 1, prio);
            }
        }
        out.print("\n");

        for (int ist = 0; ist < sn.nstations; ist++) {
            Station st = sn.stations.get(ist);
            int inode = (int) sn.stationToNode.get(ist);
            if (sn.nodetype.get(inode) == NodeType.Join) {
                continue;
            }
            boolean ext = sn.sched.get(st) == SchedStrategy.EXT;
            String setter = ext ? "setArrival" : "setService";
            for (int k = 0; k < K; k++) {
                JobClass jc = sn.jobclasses.get(k);
                String tag = " % (" + sn.nodenames.get(inode) + "," + sn.classnames.get(k) + ")\n";
                String head = "node{" + (inode + 1) + "}." + setter + "(jobclass{" + (k + 1) + "}, ";
                double w = ext ? 1.0 : sn.schedparam.get(ist, k);
                String wtxt = (w != 1.0) ? ", " + fmtF(w, M) : "";
                Distribution d = QN2JAVA.userDistribution(st, jc);
                if (d instanceof Replayer) {
                    String kind = (d instanceof Trace) ? "Trace" : "Replayer";
                    out.print(head + kind + "('" + mstr(((Replayer) d).getFileName()) + "')" + wtxt + ");" + tag);
                    continue;
                }
                QN2JAVA.Moments mo = QN2JAVA.moments(sn, st, jc);
                if (mo.disabled) {
                    out.print(head + "Disabled.getInstance());" + tag);
                } else if (mo.scv >= 0.5) {
                    if (Math.abs(mo.scv - 1.0) < 1e-12) { // JAR map_scv can return 1+eps where MATLAB's returns exactly 1
                        if (mo.mean < GlobalConstants.CoarseTol) {
                            out.print(head + "Immediate());" + tag);
                        } else {
                            out.print(head + "Exp.fitMean(" + fmtF(mo.mean, M) + ")" + wtxt + ");" + tag);
                        }
                    } else {
                        out.print(head + "APH.fitMeanAndSCV(" + fmtF(mo.mean, M) + "," + fmtF(mo.scv, M) + ")" + wtxt + ");" + tag);
                    }
                } else {
                    int nPhases = QN2JAVA.erlangPhases(mo.scv, sn.nodenames.get(inode), sn.classnames.get(k));
                    out.print(head + "Erlang(" + fmtF(nPhases / mo.mean, M) + "," + fmtF(nPhases, M) + ")" + wtxt + ");" + tag);
                }
            }
        }

        out.print("\n%% Block 3: topology");
        out.print("\n");
        out.print("P = model.initRoutingMatrix(); % initialize routing matrix \n");
        for (int k = 0; k < K; k++) {
            for (int c = 0; c < K; c++) {
                for (int i = 0; i < sn.nnodes; i++) {
                    NodeType nti = sn.nodetype.get(i);
                    // the JAR rtnodes keeps a Source row for closed classes, which MATLAB and Python do not; a Source emits no closed jobs
                    if (nti == NodeType.Sink || (nti == NodeType.Source && !Double.isInfinite(sn.njobs.get(k)))) {
                        continue;
                    }
                    for (int m = 0; m < sn.nnodes; m++) {
                        double p = sn.rtnodes.get(i * K + k, m * K + c);
                        if (p > 0) {
                            String v = (nti == NodeType.Fork) ? "1.0" : fmtD(p, M);
                            out.printf("P{%d,%d}(%d,%d) = %s; %% (%s,%s) -> (%s,%s)\n", k + 1, c + 1, i + 1, m + 1, v,
                                    sn.nodenames.get(i), sn.classnames.get(k), sn.nodenames.get(m), sn.classnames.get(c));
                        }
                    }
                }
            }
        }
        out.print("model.link(P);\n");
    }
}
