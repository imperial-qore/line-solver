/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.GlobalConstants;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_scv;
import jline.io.CodeGenSupport.Lang;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.NodeParam;
import jline.lang.nodes.Node;
import jline.lang.nodes.ServiceStation;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.Replayer;
import jline.lang.processes.Trace;
import jline.util.matrix.MatrixCell;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;

import static jline.io.CodeGenSupport.fmtF;
import static jline.io.CodeGenSupport.fmtInt;
import static jline.io.CodeGenSupport.jstr;
import static jline.io.InputOutput.line_error;

/**
 * Generates the Java (JLINE) source of a method that rebuilds a {@link Network}.
 *
 * <p>Port of {@code matlab/src/io/QN2JAVA.m}. The output is the body of a static method
 * {@code public static Network ex()} (or only its statements when {@code headers} is false) laid out in
 * the same three blocks as the MATLAB generator: nodes, classes with their arrival and service
 * processes, and the routing matrix read from {@code sn.rtnodes}. It compiles against
 * {@code jline.lang.*}, {@code jline.lang.nodes.*}, {@code jline.lang.processes.*} and
 * {@code jline.lang.constant.*}.</p>
 *
 * <p>Processes are written as MATLAB writes them, from their first two moments: an SCV of 1 becomes
 * {@code Exp.fitMean}, an SCV of at least 0.5 becomes {@code APH.fitMeanAndSCV}, a lower SCV becomes an
 * Erlang with {@code round(1/SCV)} phases, and a mean below {@link GlobalConstants#CoarseTol} becomes
 * {@code Immediate}. A trace-driven process keeps its file. So mean-value measures of the rebuilt model
 * equal the original's, while a non-exponential process is represented by its two-moment fit. Class
 * switching nodes become Routers whose switching is carried by the routing matrix, exactly as in MATLAB.</p>
 *
 * <p>Differences from the MATLAB text, each made so the output compiles and rebuilds the same model:
 * distributions are instantiated with {@code new}; a number MATLAB's {@code %f} would round is printed
 * in full; an Erlang phase count is an integer literal; a Fork row is written with probability 1 and a
 * tasks-per-link other than 1 is restored through {@code setTasksPerLink}. Node types MATLAB would
 * silently omit (Cache, Logger, Place, Transition) are refused by name instead.</p>
 */
public final class QN2JAVA {

    private QN2JAVA() {
    }

    /** Returns the generated method, with headers, for a model named after {@code model.getName()}. */
    public static String generate(Network model) {
        return generate(model, model.getName(), true);
    }

    /** Returns the generated method, with headers. */
    public static String generate(Network model, String modelName) {
        return generate(model, modelName, true);
    }

    /** Returns the generated code; {@code headers} adds the {@code ex()} method signature and return. */
    public static String generate(final Network model, final String modelName, final boolean headers) {
        return CodeGenSupport.capture(new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                QN2JAVA.emit(model, modelName, out, headers);
            }
        });
    }

    /** Writes the generated method, with headers, to a stream (e.g. {@code System.out}). */
    public static void write(Network model, String modelName, PrintStream out) {
        write(model, modelName, out, true);
    }

    /** Writes the generated code to a stream, which is flushed but not closed. */
    public static void write(Network model, String modelName, PrintStream out, boolean headers) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw, headers);
        pw.flush();
    }

    /** Writes the generated method, with headers, to a writer. */
    public static void write(Network model, String modelName, Writer out) {
        write(model, modelName, out, true);
    }

    /** Writes the generated code to a writer, which is flushed but not closed. */
    public static void write(Network model, String modelName, Writer out, boolean headers) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw, headers);
        pw.flush();
    }

    /** Writes the generated method, with headers, to {@code filename}. */
    public static void write(Network model, String modelName, String filename) {
        write(model, modelName, filename, true);
    }

    /** Writes the generated code to {@code filename}, overwriting it. */
    public static void write(final Network model, final String modelName, String filename, final boolean headers) {
        CodeGenSupport.toFile(filename, new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                QN2JAVA.emit(model, modelName, out, headers);
            }
        });
    }

    static void emit(Network model, String modelName, PrintWriter out, boolean headers) {
        final Lang J = Lang.JAVA;
        if (modelName == null) {
            modelName = "myModel";
        }
        NetworkStruct sn = model.getStruct();
        int K = sn.nclasses;

        if (headers) {
            out.print("\tpublic static Network ex() {\n");
        }
        out.printf("\t\tNetwork model = new Network(\"%s\");\n", jstr(modelName));
        out.print("\n\t\t// Block 1: nodes");
        out.print("\t\t\t\n");

        // Block 1: nodes
        for (int i = 0; i < sn.nnodes; i++) {
            int n = i + 1;
            String name = jstr(sn.nodenames.get(i));
            NodeType nt = sn.nodetype.get(i);
            switch (nt) {
                case Source:
                    out.printf("\t\tSource node%d = new Source(model, \"%s\");\n", n, name);
                    break;
                case Delay:
                    out.printf("\t\tDelay node%d = new Delay(model, \"%s\");\n", n, name);
                    break;
                case Queue: {
                    int ist = (int) sn.nodeToStation.get(i);
                    SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                    out.printf("\t\tQueue node%d = new Queue(model, \"%s\", SchedStrategy.%s);\n", n, name, sched.name());
                    double c = sn.nservers.get(ist);
                    if (c > 1) {
                        out.printf("\t\tnode%d.setNumberOfServers(%s);\n", n, fmtInt(c, J));
                    }
                    break;
                }
                case Router:
                    out.printf("\t\tRouter node%d = new Router(model, \"%s\");\n", n, name);
                    break;
                case Fork: {
                    out.printf("\t\tFork node%d = new Fork(model, \"%s\");\n", n, name);
                    double fanOut = forkFanOut(sn, i);
                    if (fanOut != 1.0) {
                        out.printf("\t\tnode%d.setTasksPerLink(%s);\n", n, fmtInt(fanOut, J));
                    }
                    break;
                }
                case Join:
                    out.printf("\t\tJoin node%d = new Join(model, \"%s\", node%d);\n", n, name, forkOf(sn, i) + 1);
                    break;
                case Sink:
                    out.printf("\t\tSink node%d = new Sink(model, \"%s\");\n", n, name);
                    break;
                case ClassSwitch:
                    out.printf("\t\tRouter node%d = new Router(model, \"%s\"); // Dummy node, class switching is embedded in the routing matrix P \n", n, name);
                    break;
                default:
                    line_error("QN2JAVA", "QN2JAVA cannot generate code for node " + sn.nodenames.get(i) + " of type " + nt + ".");
            }
        }

        // Block 2: classes
        out.print("\n\t\t// Block 2: classes\n");
        for (int k = 0; k < K; k++) {
            String cname = jstr(sn.classnames.get(k));
            double njobs = sn.njobs.get(k);
            int prio = (int) sn.classprio.get(k);
            if (Double.isInfinite(njobs)) {
                out.printf("\t\tOpenClass jobclass%d = new OpenClass(model, \"%s\", %d);\n", k + 1, cname, prio);
            } else {
                int refNode = referenceNode(sn, k);
                out.printf("\t\tClosedClass jobclass%d = new ClosedClass(model, \"%s\", %d, node%d, %d);\n",
                        k + 1, cname, (long) njobs, refNode + 1, prio);
            }
        }
        out.print("\t\t\n");

        // arrival and service processes
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
                String tag = " // (" + sn.nodenames.get(inode) + "," + sn.classnames.get(k) + ")\n";
                String head = "\t\tnode" + (inode + 1) + "." + setter + "(jobclass" + (k + 1) + ", ";
                double w = ext ? 1.0 : sn.schedparam.get(ist, k);
                String wtxt = (w != 1.0) ? ", " + fmtF(w, J) : "";
                Distribution d = userDistribution(st, jc);
                if (d instanceof Replayer) {
                    String kind = (d instanceof Trace) ? "Trace" : "Replayer";
                    out.print(head + "new " + kind + "(\"" + jstr(((Replayer) d).getFileName()) + "\")" + wtxt + ");" + tag);
                    continue;
                }
                Moments mo = moments(sn, st, jc);
                if (mo.disabled) {
                    out.print(head + "Disabled.getInstance());" + tag);
                } else if (mo.scv >= 0.5) {
                    if (Math.abs(mo.scv - 1.0) < 1e-12) { // JAR map_scv can return 1+eps where MATLAB's returns exactly 1
                        if (mo.mean < GlobalConstants.CoarseTol) {
                            out.print(head + "Immediate.getInstance());" + tag);
                        } else {
                            out.print(head + "Exp.fitMean(" + fmtF(mo.mean, J) + ")" + wtxt + ");" + tag);
                        }
                    } else {
                        out.print(head + "APH.fitMeanAndSCV(" + fmtF(mo.mean, J) + "," + fmtF(mo.scv, J) + ")" + wtxt + ");" + tag);
                    }
                } else {
                    int nPhases = erlangPhases(mo.scv, sn.nodenames.get(inode), sn.classnames.get(k));
                    out.print(head + "new Erlang(" + fmtF(nPhases / mo.mean, J) + "," + nPhases + ")" + wtxt + ");" + tag);
                }
            }
        }

        // Block 3: topology
        out.print("\n\t\t// Block 3: topology");
        out.print("\t\n");
        out.print("\t\tRoutingMatrix routingMatrix = model.initRoutingMatrix(); \n");
        out.print("\t\n");
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
                            double v = (nti == NodeType.Fork) ? 1.0 : p;
                            out.printf("\t\troutingMatrix.set(jobclass%d, jobclass%d, node%d, node%d, %s); // (%s,%s) -> (%s,%s)\n",
                                    k + 1, c + 1, i + 1, m + 1, fmtF(v, J), sn.nodenames.get(i), sn.classnames.get(k),
                                    sn.nodenames.get(m), sn.classnames.get(c));
                        }
                    }
                }
            }
        }
        out.print("\n\t\tmodel.link(routingMatrix);\n\n");
        if (headers) {
            out.print("\t\treturn model;\n");
            out.print("\t}\n");
        }
    }

    // ---- helpers shared with QN2MATLAB ------------------------------------------------------------

    /** First two moments of a (station, class) process read from {@code sn.proc}, as MATLAB reads them. */
    static final class Moments {
        final double mean;
        final double scv;
        final boolean disabled;

        Moments(double mean, double scv, boolean disabled) {
            this.mean = mean;
            this.scv = scv;
            this.disabled = disabled;
        }
    }

    static Moments moments(NetworkStruct sn, Station st, JobClass jc) {
        ProcessType pt = (sn.procid != null && sn.procid.get(st) != null) ? sn.procid.get(st).get(jc) : null;
        MatrixCell ph = (sn.proc != null && sn.proc.get(st) != null) ? sn.proc.get(st).get(jc) : null;
        if (pt == ProcessType.DISABLED || ph == null || ph.size() < 2 || Double.isNaN(ph.get(0).get(0, 0))) {
            return new Moments(Double.NaN, Double.NaN, true);
        }
        return new Moments(Map_mean.map_mean(ph), Map_scv.map_scv(ph), false);
    }

    /** MATLAB's {@code max(1, round(1/SCV))}, refusing the SCV-0 case that MATLAB would print as Inf. */
    static int erlangPhases(double scv, String node, String cls) {
        double r = Math.rint(1.0 / scv);
        if (Double.isInfinite(r) || Double.isNaN(r) || r > Integer.MAX_VALUE) {
            line_error("QN2JAVA", "The process of class " + cls + " at " + node + " has SCV " + scv
                    + ", which has no Erlang representation.");
        }
        // Math.round rounds halves up as MATLAB's round does for positives; rint would round them to even
        return (int) Math.max(1, Math.round(1.0 / scv));
    }

    /** The user-level distribution of a station, used only to recognise trace-driven processes. */
    static Distribution userDistribution(Station st, JobClass jc) {
        try {
            if (st instanceof Source) {
                return ((Source) st).getArrivalProcess(jc);
            }
            if (st instanceof ServiceStation) {
                return ((ServiceStation) st).getServiceProcess(jc);
            }
        } catch (RuntimeException ignored) {
            // a class the station does not serve has no process; it is written from sn.proc as Disabled
        }
        return null;
    }

    /** Node index of the Fork a Join closes, read from {@code sn.fj}; a Join that closes none is refused. */
    static int forkOf(NetworkStruct sn, int joinIdx) {
        if (sn.fj != null) {
            for (int f = 0; f < sn.fj.getNumRows(); f++) {
                if (sn.fj.get(f, joinIdx) != 0) {
                    return f;
                }
            }
        }
        line_error("QN2JAVA", "Join '" + sn.nodenames.get(joinIdx) + "' closes no Fork: the model cannot be written as source.");
        return -1;
    }

    /** Tasks per link of a Fork node, 1 when not configured. */
    static double forkFanOut(NetworkStruct sn, int nodeIdx) {
        if (sn.nodeparam == null) {
            return 1.0;
        }
        Node node = sn.nodes.get(nodeIdx);
        NodeParam np = sn.nodeparam.get(node);
        if (np instanceof ForkNodeParam) {
            double f = ((ForkNodeParam) np).fanOut;
            if (!Double.isNaN(f) && f > 0) {
                return f;
            }
        }
        return 1.0;
    }

    /**
     * Reference node of a closed class, as MATLAB's zeroPopRefNode picks it. A class with jobs uses its reference
     * station. A class with no jobs takes the chain's choice (the reference station of a populated closed class of its
     * chain, else its own sn.refstat), since refreshRoutingMatrix rejects a chain whose classes name different
     * reference stations; only when that is not a valid station does it fall back to the first station serving it.
     */
    static int referenceNode(NetworkStruct sn, int k) {
        if (sn.njobs.get(k) > 0) {
            return (int) sn.stationToNode.get((int) sn.refstat.get(k));
        }
        double cand = Double.NaN;
        int chain = -1;
        for (int c = 0; c < sn.nchains && chain < 0; c++) {
            if (sn.chains.get(c, k) > 0) {
                chain = c;
            }
        }
        if (chain >= 0) {
            for (int r = 0; r < sn.nclasses; r++) {
                double nj = sn.njobs.get(r);
                if (sn.chains.get(chain, r) > 0 && nj > 0 && !Double.isInfinite(nj)) {
                    cand = sn.refstat.get(r);
                    break;
                }
            }
            if (Double.isNaN(cand)) {
                cand = sn.refstat.get(k);
            }
        }
        int ref = -1;
        if (cand >= 0 && cand < sn.nstations && cand == Math.rint(cand)) {
            ref = (int) cand;
        } else {
            JobClass jc = sn.jobclasses.get(k);
            for (int ist = 0; ist < sn.nstations; ist++) {
                if (!moments(sn, sn.stations.get(ist), jc).disabled) {
                    ref = ist;
                    break;
                }
            }
        }
        if (ref < 0) {
            line_error("QN2JAVA", "Class '" + sn.classnames.get(k) + "' has no reference station and no station serves it.");
        }
        return (int) sn.stationToNode.get(ref);
    }
}
