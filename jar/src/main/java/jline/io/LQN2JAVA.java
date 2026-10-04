/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.io.CodeGenSupport.Lang;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.CallType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.layered.Task;
import jline.lang.processes.Distribution;
import jline.lang.processes.Replayer;
import jline.util.matrix.Matrix;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import static jline.io.CodeGenSupport.fmtG;
import static jline.io.CodeGenSupport.jarr;
import static jline.io.CodeGenSupport.jmat;
import static jline.io.CodeGenSupport.jnum;
import static jline.io.CodeGenSupport.fmtInt;
import static jline.io.CodeGenSupport.jstr;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;

/**
 * Generates a Java (JLINE) program that rebuilds a {@link LayeredNetwork} and solves it with SolverLN.
 *
 * <p>Port of {@code matlab/src/io/LQN2JAVA.m}. The output is a class {@code TestSolver<NAME>} in package
 * {@code jline.examples} whose {@code main} declares the processors, tasks, entries and activities
 * (with bindings, synchronous and asynchronous calls and replies), then the activity precedences grouped
 * as MATLAB groups them (serial, loop, or-fork, and-fork, or-join, and-join), and finally runs
 * {@code SolverLN.getEnsembleAvg()}.</p>
 *
 * <p>Think times and host demands are written as MATLAB's {@code javaDist} writes them: each law by its own
 * constructor and parameters (Weibull as (shape, scale), Erlang as (rate, phases), HyperExp in its scalar or
 * vector form, PH/APH/MAP/ME/RAP as matrices), numbers through {@code jnum} ({@code 4.0}, shortest round-trip
 * spelling otherwise), and only an unrecognised class as {@code APH.fitMeanAndSCV} with a warning.</p>
 *
 * <p>Differences from the MATLAB text, each made so the output compiles and rebuilds the same model:
 * the import block lists the same packages in another order; a call mean, loop count or or-fork
 * probability that MATLAB's {@code %g} would round is printed in full; the class name drops every
 * character a Java identifier cannot hold, not only spaces.
 * Precedences are read from each task's declared {@link ActivityPrecedence} list rather than
 * reconstructed from {@code sn.graph}: the reconstruction loses a loop whose count is at most 1 and,
 * after the first and-fork, carries the previous fork's branches into the next one. Output order within
 * each group is the MATLAB order (by activity index; joins by descending index of the joined activity).</p>
 */
public final class LQN2JAVA {

    private LQN2JAVA() {
    }

    /** Returns the generated program for a model named after {@code model.getName()}. */
    public static String generate(LayeredNetwork model) {
        return generate(model, model.getName());
    }

    /** Returns the generated program. */
    public static String generate(final LayeredNetwork model, final String modelName) {
        return CodeGenSupport.capture(new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                LQN2JAVA.emit(model, modelName, out);
            }
        });
    }

    /** Writes the generated program to a stream (e.g. {@code System.out}), which is flushed but not closed. */
    public static void write(LayeredNetwork model, String modelName, PrintStream out) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw);
        pw.flush();
    }

    /** Writes the generated program to a writer, which is flushed but not closed. */
    public static void write(LayeredNetwork model, String modelName, Writer out) {
        PrintWriter pw = CodeGenSupport.wrap(out);
        emit(model, modelName, pw);
        pw.flush();
    }

    /** Writes the generated program to {@code filename}, overwriting it. */
    public static void write(final LayeredNetwork model, final String modelName, String filename) {
        CodeGenSupport.toFile(filename, new CodeGenSupport.Emitter() {
            public void emit(PrintWriter out) {
                LQN2JAVA.emit(model, modelName, out);
            }
        });
    }

    /** The class name MATLAB derives from the model name, reduced to a valid Java identifier. */
    public static String className(String modelName) {
        StringBuilder b = new StringBuilder("TestSolver");
        String up = modelName.replace(" ", "").toUpperCase(java.util.Locale.ROOT);
        for (int i = 0; i < up.length(); i++) {
            char ch = up.charAt(i);
            if (Character.isJavaIdentifierPart(ch)) {
                b.append(ch);
            }
        }
        return b.toString();
    }

    static void emit(LayeredNetwork model, String modelName, PrintWriter out) {
        final Lang J = Lang.JAVA;
        if (modelName == null) {
            modelName = "myLayeredModel";
        }
        LayeredNetworkStruct sn = model.getStruct();
        LqnChecks.refuseUnsupported(model, sn, "LQN2JAVA");

        out.print("package jline.examples;\n\n");
        out.print("import java.util.ArrayList;\n");
        out.print("import jline.lang.*;\n");
        out.print("import jline.lang.layered.*;\n");
        out.print("import jline.lang.constant.*;\n");
        out.print("import jline.lang.processes.*;\n");
        out.print("import jline.util.matrix.Matrix;\n");
        out.print("import jline.solvers.ln.SolverLN;\n\n");
        out.printf("public class %s {\n\n", className(modelName));
        out.print("\tpublic static void main(String[] args) throws Exception{\n\n");
        out.printf("\tLayeredNetwork model = new LayeredNetwork(\"%s\");\n", jstr(modelName));
        out.print("\n");

        // host processors
        for (int h = 0; h < sn.nhosts; h++) {
            out.printf("\tProcessor P%d = new Processor(model, \"%s\", %s, SchedStrategy.%s);\n", h + 1,
                    jstr(sn.names.get(h)), fmtInt(sn.mult.get(h), J), schedFeature(sn.sched.get(h)));
            if (sn.repl.get(h) != 1) {
                out.printf("P%d.setReplication(%d);\n", h + 1, (long) sn.repl.get(h));
            }
        }
        out.print("\n");

        // tasks
        for (int t = 0; t < sn.ntasks; t++) {
            int tidx = sn.tshift + t;
            out.printf("\tTask T%d = new Task(model, \"%s\", %s, SchedStrategy.%s).on(P%d);\n", t + 1,
                    jstr(sn.names.get(tidx)), fmtInt(sn.mult.get(tidx), J), schedFeature(sn.sched.get(tidx)),
                    (int) sn.parent.get(tidx) + 1);
            if (sn.repl.get(tidx) != 1) {
                out.printf("\tT%d.setReplication(%d);\n", t + 1, (long) sn.repl.get(tidx));
            }
            Distribution th = sn.think.get(tidx);
            if (th != null && sn.think_type.get(tidx) != ProcessType.DISABLED) {
                out.printf("\tT%d.setThinkTime(%s);\n", t + 1, javaDist(th, "LQN2JAVA"));
            }
        }
        out.print("\n");

        // entries
        for (int e = 0; e < sn.nentries; e++) {
            int eidx = sn.eshift + e;
            out.printf("\tEntry E%d = new Entry(model, \"%s\").on(T%d);\n", e + 1, jstr(sn.names.get(eidx)),
                    (int) sn.parent.get(eidx) - sn.tshift + 1);
        }
        out.print("\n");

        // activities
        for (int a = 0; a < sn.nacts; a++) {
            int aidx = sn.ashift + a;
            int tidx = (int) sn.parent.get(aidx);
            StringBuilder boundTo = new StringBuilder();
            for (int e = 0; e < sn.nentries; e++) {
                if (sn.graph.get(sn.eshift + e, aidx) != 0) {
                    boundTo.append(".boundTo(E").append(e + 1).append(')');
                }
            }
            StringBuilder repliesTo = new StringBuilder();
            if (sn.sched.get(tidx) != SchedStrategy.REF) {
                for (int e = 0; e < sn.nentries; e++) {
                    if (sn.replygraph.get(a, e) != 0 && sn.sched.get((int) sn.parent.get(sn.eshift + e)) != SchedStrategy.REF) {
                        repliesTo.append(".repliesTo(E").append(e + 1).append(')');
                    }
                }
            }
            StringBuilder calls = new StringBuilder();
            for (int c = 0; c < sn.ncalls; c++) {
                if ((int) sn.callpair.get(c, 0) != aidx) {
                    continue;
                }
                int target = (int) sn.callpair.get(c, 1) - sn.eshift + 1;
                CallType ct = sn.calltype.get(c);
                if (ct == CallType.SYNC) {
                    calls.append(".synchCall(E").append(target).append(',').append(fmtG(sn.callproc_mean.get(c), J)).append(')');
                } else if (ct == CallType.ASYNC) {
                    calls.append(".asynchCall(E").append(target).append(',').append(fmtG(sn.callproc_mean.get(c), J)).append(')');
                }
            }
            out.printf("\tActivity A%d = new Activity(model, \"%s\", %s).on(T%d);", a + 1, jstr(sn.names.get(aidx)),
                    javaDist(sn.hostdem.get(aidx), "LQN2JAVA"),
                    tidx - sn.tshift + 1);
            if (boundTo.length() > 0) {
                out.printf(" A%d%s;", a + 1, boundTo);
            }
            if (calls.length() > 0) {
                out.printf(" A%d%s;", a + 1, calls);
            }
            if (repliesTo.length() > 0) {
                out.printf(" A%d%s;", a + 1, repliesTo);
            }
            out.print("\n");
            // activity think time; the Activity default is Immediate, so only a non-default one is written
            Distribution ath = sn.actthink == null ? null : sn.actthink.get(aidx);
            ProcessType at = sn.actthink_type == null ? null : sn.actthink_type.get(aidx);
            if (ath != null && at != ProcessType.DISABLED && at != ProcessType.IMMEDIATE) {
                out.printf("\tA%d.setThinkTime(%s);\n", a + 1, javaDist(ath, "LQN2JAVA"));
            }
        }
        out.print("\n");

        emitPrecedences(model, sn, out);

        out.print("\n\t// Model solution \n");
        out.print("\tSolverLN solver = new SolverLN(model);\n");
        out.print("\tsolver.getEnsembleAvg();\n");
        out.print("\t}\n}\n");
    }

    /** The constant MATLAB's {@code SchedStrategy.toFeature} names: FCFSPRIO is spelled HOL, every other policy by its name. */
    static String schedFeature(SchedStrategy s) {
        return s == SchedStrategy.FCFSPRIO ? "HOL" : s.name();
    }

    /**
     * Java constructor call rebuilding {@code d} with the same parameters, as MATLAB's {@code javaDist}: each law is
     * written with its own parameters in constructor order, and only an unrecognised class falls back to an APH
     * fitted to its mean and SCV, with a warning. No distribution at all is the Immediate default.
     */
    static String javaDist(Distribution d, String caller) {
        if (d == null) {
            return "new Immediate()";
        }
        String cls = d.getClass().getSimpleName();
        switch (cls) {
            case "Immediate":
            case "Disabled":
                return "new " + cls + "()";
            case "Exp":
            case "Det":
            case "Geometric":
            case "Poisson":
            case "Bernoulli":
                return "new " + cls + "(" + jnum(scalar(d, 1)) + ")";
            case "Erlang":
                return "new Erlang(" + jnum(scalar(d, 1)) + ", " + Math.round(scalar(d, 2)) + ")";
            case "HyperExp":
                if (d.getParam(1).getValue() instanceof Number) {
                    return "new HyperExp(" + jnum(scalar(d, 1)) + ", " + jnum(scalar(d, 2)) + ", " + jnum(scalar(d, 3)) + ")";
                }
                return "new HyperExp(" + jarr(matrix(d, 1).toArray1D()) + ", " + jarr(matrix(d, 2).toArray1D()) + ")";
            case "Gamma":
            case "Lognormal":
            case "Uniform":
            case "Pareto":
            case "Normal":
            case "DiscreteUniform":
                return "new " + cls + "(" + jnum(scalar(d, 1)) + ", " + jnum(scalar(d, 2)) + ")";
            case "Binomial":
                return "new Binomial(" + Math.round(scalar(d, 1)) + ", " + jnum(scalar(d, 2)) + ")";
            case "Weibull":
                // params are (scale alpha, shape r), the constructor takes (shape, scale)
                return "new Weibull(" + jnum(scalar(d, 2)) + ", " + jnum(scalar(d, 1)) + ")";
            case "Coxian":
                return "new Coxian(" + jmat(matrix(d, 1)) + ", " + jmat(matrix(d, 2)) + ")";
            case "PH":
            case "APH":
                // the JAR stores (n, alpha, T), the constructor takes (alpha, T)
                return "new " + cls + "(" + jmat(matrix(d, 2)) + ", " + jmat(matrix(d, 3)) + ")";
            case "MAP":
            case "ME":
            case "RAP":
                return "new " + cls + "(" + jmat(matrix(d, 1)) + ", " + jmat(matrix(d, 2)) + ")";
            case "MMPP2":
                return "new MMPP2(" + jnum(scalar(d, 1)) + ", " + jnum(scalar(d, 2)) + ", " + jnum(scalar(d, 3)) + ", "
                        + jnum(scalar(d, 4)) + ")";
            case "Zipf":
                return "new Zipf(" + jnum(scalar(d, 3)) + ", " + Math.round(scalar(d, 4)) + ")";
            case "DiscreteSampler":
                return "new DiscreteSampler(" + jmat(matrix(d, 1)) + ", " + jmat(matrix(d, 2)) + ")";
            case "Replayer":
            case "Trace": {
                String f = ((Replayer) d).getFileName();
                if (f == null) {
                    line_error(caller, caller + " cannot write a " + cls + " built from in-memory samples: it names no trace file.");
                }
                return "new " + cls + "(\"" + f.replace("\\", "\\\\").replace("\"", "\\\"") + "\")";
            }
            default:
                line_warning(caller, "%s writes the %s distribution as an APH fitted to its mean and SCV.", caller, cls);
                return "APH.fitMeanAndSCV(" + jnum(d.getMean()) + ", " + jnum(d.getSCV()) + ")";
        }
    }

    /** Parameter {@code k} of {@code d} as a double. */
    private static double scalar(Distribution d, int k) {
        Object v = d.getParam(k).getValue();
        if (v instanceof Number) {
            return ((Number) v).doubleValue();
        }
        return matrix(d, k).get(0, 0);
    }

    /** Parameter {@code k} of {@code d} as a Matrix; a list or array is a row. */
    private static Matrix matrix(Distribution d, int k) {
        Object v = d.getParam(k).getValue();
        if (v instanceof Matrix) {
            return (Matrix) v;
        }
        double[] row;
        if (v instanceof double[]) {
            row = (double[]) v;
        } else if (v instanceof List) {
            List<?> l = (List<?>) v;
            row = new double[l.size()];
            for (int i = 0; i < row.length; i++) {
                row[i] = ((Number) l.get(i)).doubleValue();
            }
        } else if (v instanceof Number) {
            row = new double[]{((Number) v).doubleValue()};
        } else {
            line_error("LQN2JAVA", "LQN2JAVA cannot read parameter " + k + " of a " + d.getClass().getSimpleName() + " holding a "
                    + (v == null ? "null" : v.getClass().getSimpleName()) + ".");
            return null;
        }
        Matrix m = new Matrix(1, row.length);
        for (int i = 0; i < row.length; i++) {
            m.set(0, i, row[i]);
        }
        return m;
    }

    // ---- precedences ---------------------------------------------------------------------------

    /** A declared precedence with its owning task, keyed for MATLAB's output order. */
    private static final class Prec {
        final int task;          // 1-based local task index
        final ActivityPrecedence p;
        final int key;           // activity index that orders it within its group

        Prec(int task, ActivityPrecedence p, int key) {
            this.task = task;
            this.p = p;
            this.key = key;
        }
    }

    private static int actIndex(LayeredNetworkStruct sn, String name) {
        for (int a = 0; a < sn.nacts; a++) {
            if (name.equals(sn.names.get(sn.ashift + a))) {
                return a;
            }
        }
        line_error("LQN2JAVA", "Precedence names an unknown activity '" + name + "'.");
        return -1;
    }

    /** Names sorted by activity index, as MATLAB's find() returns them; {@code perm} receives the permutation. */
    static List<String> sortedByIndex(final LayeredNetworkStruct sn, List<String> names, List<Integer> perm) {
        final List<Integer> order = new ArrayList<Integer>();
        for (int i = 0; i < names.size(); i++) {
            order.add(i);
        }
        final List<Integer> idx = new ArrayList<Integer>();
        for (String n : names) {
            idx.add(actIndex(sn, n));
        }
        Collections.sort(order, new Comparator<Integer>() {
            public int compare(Integer x, Integer y) {
                return Integer.compare(idx.get(x), idx.get(y));
            }
        });
        List<String> out = new ArrayList<String>();
        for (Integer o : order) {
            out.add(names.get(o));
            if (perm != null) {
                perm.add(o);
            }
        }
        return out;
    }

    /** Groups every declared precedence by the MATLAB section it belongs to. */
    static List<List<Prec>> groupPrecedences(LayeredNetwork model, LayeredNetworkStruct sn) {
        List<List<Prec>> g = new ArrayList<List<Prec>>();
        for (int i = 0; i < 7; i++) {
            g.add(new ArrayList<Prec>());
        }
        for (int t = 0; t < sn.ntasks; t++) {
            Task task = model.getTasks().get(t);
            for (ActivityPrecedence p : task.getPrecedences()) {
                String pre = p.getPreType();
                String post = p.getPostType();
                int firstPre = actIndex(sn, p.getPreActs().get(0));
                if (ActivityPrecedenceType.PRE_SEQ.equals(pre) && ActivityPrecedenceType.POST_SEQ.equals(post)) {
                    g.get(0).add(new Prec(t + 1, p, firstPre));
                } else if (ActivityPrecedenceType.POST_LOOP.equals(post)) {
                    g.get(1).add(new Prec(t + 1, p, firstPre));
                } else if (ActivityPrecedenceType.PRE_SEQ.equals(pre) && ActivityPrecedenceType.POST_OR.equals(post)) {
                    g.get(2).add(new Prec(t + 1, p, firstPre));
                } else if (ActivityPrecedenceType.PRE_SEQ.equals(pre) && ActivityPrecedenceType.POST_AND.equals(post)) {
                    g.get(3).add(new Prec(t + 1, p, firstPre));
                } else if (ActivityPrecedenceType.PRE_OR.equals(pre) && ActivityPrecedenceType.POST_SEQ.equals(post)) {
                    g.get(4).add(new Prec(t + 1, p, -actIndex(sn, p.getPostActs().get(0))));
                } else if (ActivityPrecedenceType.PRE_AND.equals(pre) && ActivityPrecedenceType.POST_SEQ.equals(post)) {
                    g.get(5).add(new Prec(t + 1, p, -actIndex(sn, p.getPostActs().get(0))));
                } else {
                    g.get(6).add(new Prec(t + 1, p, firstPre));
                }
            }
        }
        for (List<Prec> l : g) {
            Collections.sort(l, new Comparator<Prec>() {
                public int compare(Prec x, Prec y) {
                    return Integer.compare(x.key, y.key);
                }
            });
        }
        return g;
    }

    private static void emitPrecedences(LayeredNetwork model, LayeredNetworkStruct sn, PrintWriter out) {
        final Lang J = Lang.JAVA;
        List<List<Prec>> g = groupPrecedences(model, sn);
        boolean hasPreActs = false;
        boolean hasPostActs = false;
        boolean hasProbs = false;

        for (Prec q : g.get(0)) {
            for (String a : q.p.getPreActs()) {
                for (String b : q.p.getPostActs()) {
                    out.printf("\tT%d.addPrecedence(ActivityPrecedence.Serial(\"%s\", \"%s\"));\n", q.task, jstr(a), jstr(b));
                }
            }
        }
        for (Prec q : g.get(1)) {
            out.print("\n\t// Loop Activity Precedence \n");
            hasPreActs = declareList(out, "precActs", hasPreActs, false);
            for (String b : q.p.getPostActs()) {
                out.printf("\tprecActs.add(\"%s\");\n", jstr(b));
            }
            out.printf("\tT%d.addPrecedence(ActivityPrecedence.Loop(\"%s\", precActs, Matrix.singleton(%s)));\n", q.task,
                    jstr(q.p.getPreActs().get(0)), fmtG(q.p.getPostParams().value(), J));
        }
        for (Prec q : g.get(2)) {
            List<Integer> perm = new ArrayList<Integer>();
            List<String> posts = sortedByIndex(sn, q.p.getPostActs(), perm);
            out.print("\n\t// OrFork Activity Precedence \n");
            hasPreActs = declareList(out, "precActs", hasPreActs, false);
            if (!hasProbs) {
                out.printf("\tMatrix probs = new Matrix(1,%d);\n", posts.size());
                hasProbs = true;
            } else {
                out.printf("\tprobs = new Matrix(1,%d);\n", posts.size());
            }
            for (String b : posts) {
                out.printf("\tprecActs.add(\"%s\");\n", jstr(b));
            }
            Matrix pr = q.p.getPostParams();
            for (int j = 0; j < posts.size(); j++) {
                out.printf("\tprobs.set(0,%d,%s);\n", j, fmtG(pr.get(perm.get(j)), J));
            }
            out.printf("\tT%d.addPrecedence(ActivityPrecedence.OrFork(\"%s\", precActs, probs));\n", q.task,
                    jstr(q.p.getPreActs().get(0)));
        }
        for (Prec q : g.get(3)) {
            out.print("\n\t// AndFork Activity Precedence \n");
            hasPostActs = declareList(out, "postActs", hasPostActs, true);
            for (String b : sortedByIndex(sn, q.p.getPostActs(), null)) {
                out.printf("\tpostActs.add(\"%s\");\n", jstr(b));
            }
            out.printf("\tT%d.addPrecedence(ActivityPrecedence.AndFork(\"%s\", postActs));\n", q.task,
                    jstr(q.p.getPreActs().get(0)));
        }
        for (Prec q : g.get(4)) {
            out.print("\n\t// OrJoin Activity Precedence \n");
            hasPreActs = declareList(out, "precActs", hasPreActs, false);
            for (String a : sortedByIndex(sn, q.p.getPreActs(), null)) {
                out.printf("\tprecActs.add(\"%s\");\n", jstr(a));
            }
            out.printf("\tT%d.addPrecedence(ActivityPrecedence.OrJoin(precActs, \"%s\"));\n", q.task,
                    jstr(q.p.getPostActs().get(0)));
        }
        for (Prec q : g.get(5)) {
            out.print("\n\t// AndJoin Activity Precedence \n");
            hasPreActs = declareList(out, "precActs", hasPreActs, false);
            for (String a : sortedByIndex(sn, q.p.getPreActs(), null)) {
                out.printf("\tprecActs.add(\"%s\");\n", jstr(a));
            }
            Matrix quorum = q.p.getPreParams();
            if (quorum == null || quorum.isEmpty()) {
                out.printf("\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, \"%s\"));\n", q.task,
                        jstr(q.p.getPostActs().get(0)));
            } else {
                out.printf("\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, \"%s\", Matrix.singleton(%s)));\n",
                        q.task, jstr(q.p.getPostActs().get(0)), fmtG(quorum.value(), J));
            }
        }
        for (Prec q : g.get(6)) {
            // a precedence that is a join and a fork at once has no factory method; use the general constructor
            out.print("\n\t// Activity Precedence \n");
            out.printf("\tT%d.addPrecedence(new ActivityPrecedence(%s, %s, \"%s\", \"%s\", %s, %s));\n", q.task,
                    javaList(q.p.getPreActs()), javaList(q.p.getPostActs()), jstr(q.p.getPreType()),
                    jstr(q.p.getPostType()), javaRow(q.p.getPreParams()), javaRow(q.p.getPostParams()));
        }
    }

    private static boolean declareList(PrintWriter out, String var, boolean declared, boolean matlabSpace) {
        if (!declared) {
            out.printf("\tArrayList<String> %s = new ArrayList<String>();\n", var);
        } else {
            // MATLAB's LQN2JAVA writes the and-fork reassignment with a stray space; kept for identical text
            out.printf(matlabSpace ? "\t %s = new ArrayList<String>();\n" : "\t%s = new ArrayList<String>();\n", var);
        }
        return true;
    }

    private static String javaList(List<String> names) {
        StringBuilder b = new StringBuilder("java.util.Arrays.asList(");
        for (int i = 0; i < names.size(); i++) {
            if (i > 0) {
                b.append(", ");
            }
            b.append('"').append(jstr(names.get(i))).append('"');
        }
        return b.append(')').toString();
    }

    private static String javaRow(Matrix m) {
        if (m == null) {
            return "(Matrix) null";
        }
        StringBuilder b = new StringBuilder("new Matrix(new double[][]{");
        for (int i = 0; i < m.getNumRows(); i++) {
            b.append(i > 0 ? ", {" : "{");
            for (int j = 0; j < m.getNumCols(); j++) {
                if (j > 0) {
                    b.append(", ");
                }
                b.append(CodeGenSupport.exact(m.get(i, j), Lang.JAVA));
            }
            b.append('}');
        }
        return b.append("})").toString();
    }
}
