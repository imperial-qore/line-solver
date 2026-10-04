/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.lqn;

import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.LayeredNetworkStruct;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Conservation laws of a layered queueing network, enumerated from its structure.
 *
 * <p>A layered model is not free to report any tuple of throughputs, think times
 * and utilizations: five families of relations tie them together, and every one of
 * them is fixed by the STRUCTURE of the model alone. This class walks a
 * {@link LayeredNetworkStruct} and emits them, one record per relation, with the
 * index sets to aggregate over, the constant coefficients and a printable form.
 * Nothing is solved here.</p>
 *
 * <p>The families, with {@code kind} as emitted:</p>
 *
 * <p><b>little</b> -- Little's law on a task's THREAD POOL. The threads of task t
 * form a closed cycle of one delay stage (the surrogate think time SolverLN imputes
 * to the task, plus the declared think time of a reference task) and one service
 * stage (holding a request from above). With B(t,k) the mean number of threads of t
 * busy serving caller class k -- the per-class utilization in JOB units --</p>
 *
 * <pre>    X(t)*(Z(t) + z(t)) + sum_k B(t,k) = N(t)</pre>
 *
 * <p>which is the update {@code SolverLN.updateThinkTimes} iterates on. In the
 * [0,1]-normalized utilization LINE reports for a queueing station,
 * B(t,k) = N(t)*U(t,k), giving X(t)*(Z(t)+z(t)) = N(t)*(1 - sum_k U(t,k)); at an
 * infinite server the utilization is already a job count, so B = U. The caller
 * classes k are the CALLS targeting an entry of t -- the in-edges of t in the call
 * graph -- plus each entry of t carrying an OPEN ARRIVAL, a stream that holds a
 * thread exactly as a call does and that a task can have alongside its callers.
 * Both are structural neighbours of t, which makes the relation node-local;
 * {@code termisentry} says which of the two a term is.</p>
 *
 * <p><b>callflow</b> -- throughput conservation across one call,
 * X(c) = X(src(c))*y(c), with y(c) the mean number of calls and src(c) the
 * dispatching activity (the dispatching ENTRY for a forwarding call).</p>
 *
 * <p><b>entryflow</b> -- the requests an entry serves are the calls reaching it
 * plus its open-arrival stream, X(e) = sum_c X(c) + lambda(e).</p>
 *
 * <p><b>actflow</b> -- an activity executes v(a) times per invocation of its entry,
 * X(a) = X(e)*v(a), with v from the activity precedence graph. An AND-JOIN is the
 * one place where flow does not add up -- its target executes once per fork, not
 * once per branch -- so the arcs into a join are scaled by 1/(number of joined
 * branches).</p>
 *
 * <p><b>hostutil</b> -- the utilization law at a processor,
 * sum_a X(a)*D(a) = m(h)*U(h); m(h)*U(h) is a job count and the factor m(h) drops
 * at an infinite server.</p>
 *
 * <p>Together these close the system: {@code little} alone is one equation per task
 * and admits the all-zero solution, so a physics-informed loss built on it should
 * carry the flow and utilization families as well.</p>
 *
 * <p><b>TWO INDEX SPACES, and they differ by one.</b> The JAR
 * {@link LayeredNetworkStruct} is 0-BASED: an element is {@code 0..nidx-1}, a call is
 * {@code 0..ncalls-1}, {@code callpair} is {@code (ncalls,2)} with columns 0/1, and a
 * not-found index is -1 because 0 is the first host. Everything this class EMITS is
 * 1-BASED instead -- {@code Relation.target}, {@code Relation.terms},
 * {@code Result.visits} and the rows and columns of the incidence matrices are the
 * struct index plus one -- which is also how a {@link Solution} is indexed, so that a
 * padded slot 0 can carry "unset" and the records read like the MATLAB and C++ twins.
 * Internally the class works in the STRUCT's space throughout: the conversion happens
 * exactly twice, at {@code readSolution} on the way in and at the {@code r.target} /
 * {@code r.terms} assignments on the way out. Reading a struct field with a 1-based
 * index is the standing trap here -- most of them are shorter than {@code nidx}
 * ({@code mult} and {@code sched} stop after the tasks), so the slip surfaces as a
 * silent NaN rather than as an error; see {@code _kb/11-conventions-and-gotchas.md}.</p>
 *
 * <p>Twin of the MATLAB {@code lqn_balance_equations.m}, the Python
 * {@code line_solver.api.lqn.balance_equations} and the C++
 * {@code line/api/lqn/lqn_balance_equations.h}. MATLAB and C++ are 1-based in their
 * struct as well, so they need no conversion; Python is 0-based on both sides and
 * emits 0-based indices.</p>
 *
 * @since LINE 3.0
 */
public final class LqnBalanceEquations {

    private LqnBalanceEquations() {
    }

    /** ActivityPrecedenceType.PRE_AND, carried on the JOINED predecessors of a join. */
    private static final int ID_PRE_AND = ActivityPrecedenceType.ID_PRE_AND;

    /**
     * Column of {@code lqn.callpair} holding the DISPATCHING element of a call (the
     * issuing activity, or the issuing entry of a forwarding call). The JAR matrix is
     * {@code (ncalls, 2)} with 0-based columns; MATLAB and C++ spell the same column
     * {@code callpair(cidx,1)}.
     */
    private static final int CALL_SRC = 0;
    /** Column of {@code lqn.callpair} holding the CALLED entry; MATLAB {@code (cidx,2)}. */
    private static final int CALL_DST = 1;

    /**
     * One conservation law, as an aggregation over a node's structural neighbours.
     */
    public static final class Relation {
        /** little | callflow | entryflow | actflow | hostutil. */
        public String kind = "";
        /**
         * Which case produced a {@code little} record (ref, inf, queueing, fwd,
         * arrival), the call type of a {@code callflow}, or the server kind of a
         * {@code hostutil}.
         */
        public String branch = "";
        /**
         * Absolute index the relation is anchored on, 1-BASED: the struct index plus
         * one, and {@code lqn.cshift + cidx + 1} for a call. See the class javadoc.
         */
        public int target;
        public String targetname = "";
        /**
         * Indices to aggregate over, 1-BASED. A {@code little} record mixes the two
         * families {@code termisentry} distinguishes: a call term is a 1-based CALL
         * index, an entry term a 1-based ELEMENT index.
         */
        public int[] terms = new int[0];
        /**
         * Per term, true when it is an ENTRY class (an open-arrival stream, or the
         * self-driven cycle of a reference task) rather than a CALL class; the two
         * read their rate from different vectors.
         */
        public boolean[] termisentry = new boolean[0];
        /** Constant coefficient of each term. */
        public double[] coeff = new double[0];
        /** Right-hand-side constant. */
        public double rhsconst;
        public double mult = Double.NaN;
        public double maxmult = Double.NaN;
        public double repl = 1.0;
        /**
         * The per-class utilization must be multiplied by {@code mult} to reach job
         * units (a queueing server); false at an infinite server.
         */
        public boolean scaled;
        public boolean phase2;
        public boolean setup;
        /**
         * The relation cannot serve as a residual: infinite multiplicity, empty term
         * set, or a zero call count.
         */
        public boolean degenerate;
        /**
         * Set once instantiated: the equality is unattainable because the task is
         * saturated and its think time is pinned at zero. Use a one-sided (hinge)
         * form there.
         */
        public boolean clamped;
        public String text = "";
        public double lhs = Double.NaN;
        public double rhs = Double.NaN;
        public double residual = Double.NaN;
        public double relresidual = Double.NaN;
        public double[] perclassutil = new double[0];

        @Override
        public String toString() {
            return kind + " " + targetname + ": " + text;
        }
    }

    /**
     * The iterates of a solved layered model.
     *
     * <p>The five vectors map onto the public fields of the same name on
     * {@code jline.solvers.ln.SolverLN}; this class does not reference the solver so
     * that the API layer keeps no dependency on it. All are indexed 1-BASED by absolute
     * element index, matching {@link Relation#target} -- which is NOT how SolverLN
     * stores them, nor how the 0-based {@link LayeredNetworkStruct} does, so use
     * {@link #fromLayeredIterates} rather than passing them straight through.</p>
     *
     * <p>{@code un} is the REPORTED utilization, which is not the same quantity as
     * the {@code util} iterate: the iterate holds a task's utilization as a SERVER in
     * its own task layer -- the U that closes the thread-pool cycle -- and it is left
     * at zero on a host. The {@code hostutil} family needs the reported one, so it is
     * supplied separately rather than read off the solver: {@code getEnsembleAvg}
     * re-enters {@code iterate()} in every codebase, and a diagnostic must not re-run
     * a fixed point as a side effect of being asked a question. A {@code hostutil}
     * record whose host has no {@code un} entry is emitted symbolically with a NaN
     * residual.</p>
     */
    public static final class Solution {
        public Matrix tput;
        public Matrix util;
        public Matrix thinkt;
        public Matrix servt;
        public Matrix residt;
        public Matrix un;

        public Solution() {
        }

        public Solution(Matrix tput, Matrix util, Matrix thinkt, Matrix servt,
                        Matrix residt, Matrix un) {
            this.tput = tput;
            this.util = util;
            this.thinkt = thinkt;
            this.servt = servt;
            this.residt = residt;
            this.un = un;
        }

        /**
         * A Solution built from the iterate vectors of {@code jline.solvers.ln.SolverLN}.
         *
         * <p>THE SHIFT IS THE POINT. SolverLN allocates its iterates as
         * {@code new Matrix(1, lqn.nidx)}, so element {@code idx} of the struct lives at
         * COLUMN {@code idx} and there is no padded slot. A {@link Solution} is 1-based,
         * matching {@link Relation#target} and the MATLAB and C++ twins, so a vector
         * coming off the solver has to be shifted UP by one on the way in; {@code compute}
         * shifts it back DOWN before indexing the struct with it. Reading it unshifted
         * silently attributes every element's iterate to its neighbour and leaves one end
         * unset, which no assertion on a single metric would catch.</p>
         *
         * <p>{@code un} is not among them: the reported utilization comes from the
         * ensemble average table as a {@code List<Double>}, which the caller indexes
         * itself, so assign that field directly.</p>
         */
        public static Solution fromLayeredIterates(Matrix tput, Matrix util, Matrix thinkt,
                                                   Matrix servt, Matrix residt) {
            return new Solution(shiftUp(tput), shiftUp(util), shiftUp(thinkt),
                    shiftUp(servt), shiftUp(residt), null);
        }

        /** A 0-based-over-elements row vector as a 1-based one, with slot 0 unused. */
        private static Matrix shiftUp(Matrix v) {
            if (v == null) {
                return null;
            }
            int n = (v.getNumRows() == 1) ? v.getNumCols() : v.getNumRows();
            Matrix out = new Matrix(1, n + 1);
            out.set(0, 0, Double.NaN);
            for (int i = 0; i < n; i++) {
                out.set(0, i + 1, (v.getNumRows() == 1) ? v.get(0, i) : v.get(i, 0));
            }
            return out;
        }
    }

    /** The relation set of one layered model. */
    public static final class Result {
        public List<Relation> eqs = new ArrayList<Relation>();
        /**
         * (nidx+1) expected executions of each activity per invocation of its entry,
         * indexed 1-BASED by absolute element index; slot 0 is the unused pad.
         */
        public double[] visits = new double[0];
        public List<String> convention = new ArrayList<String>();
        /**
         * (ntasks+1, ncalls+1) 1 where call c is a caller class of task t, both axes
         * 1-BASED over the local task and call index, with row and column 0 unused.
         */
        public Matrix A_little;
        /** (ncalls+1, nidx+1) call-to-source incidence weighted by y(c), 1-based. */
        public Matrix A_flow;
        /** (nhosts+1, nidx+1) host-to-activity incidence weighted by D(a), 1-based. */
        public Matrix A_host;
        /** Largest absolute residual over the instantiated relations, NaN without a solution. */
        public double maxresidual = Double.NaN;
        public List<String> text = new ArrayList<String>();

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < text.size(); i++) {
                if (i > 0) {
                    sb.append('\n');
                }
                sb.append(text.get(i));
            }
            return sb.toString();
        }

        public void print() {
            System.out.println(toString());
        }
    }

    /** Enumerate the relations of LQN symbolically. */
    public static Result compute(LayeredNetworkStruct lqn) {
        return compute(lqn, null);
    }

    /**
     * Enumerate the relations of LQN and, when SOL is given, instantiate each one and
     * report its residual.
     *
     * <p>Conventions: rates and populations in a {@code little} record are PER
     * REPLICA, matching {@code updateThinkTimes} -- X is tput/repl and N is the
     * multiplicity of one copy. Elsewhere throughputs and utilizations are as the
     * solver reports them, totalled over replicas, which is why the server count in
     * {@code hostutil} is {@code mult} and not {@code mult*repl}. N(t) is
     * {@code lqn.mult}; SolverLN iterates on {@code njobs}, which carries the
     * interlocking corrections and may be {@code maxmult} under replication, so both
     * are returned per record. S(k) is the entry SERVICE time (phase 1 plus phase 2),
     * the time a thread is held, not the residence time the caller waits for; the
     * difference is the phase-2 tail, flagged by {@code phase2}.</p>
     */
    public static Result compute(LayeredNetworkStruct lqn, Solution sol) {
        int nidx = lqn.nidx;
        double[] visits = actVisits(lqn);
        Map<Integer, List<Integer>> callsInto = callsInto(lqn);
        Sol s = (sol == null) ? null : readSolution(lqn, sol);

        Result out = new Result();
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            Pool pool = threadPoolBranch(lqn, tidx, callsInto.get(tidx));
            if (pool == null) {
                continue;
            }
            out.eqs.add(littleRecord(lqn, tidx, pool, s));
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            out.eqs.add(callflowRecord(lqn, cidx, s));
        }
        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            int[] inc = incomingCalls(lqn, eidx);
            double lam = arrivalRate(lqn, eidx);
            if (inc.length == 0 && lam == 0.0) {
                continue;   // a reference entry is driven by its own cycle, not by flow
            }
            out.eqs.add(entryflowRecord(lqn, eidx, inc, lam, s));
        }
        for (int a = 0; a < lqn.nacts; a++) {
            int aidx = lqn.ashift + a;
            int eidx = entryOfActivity(lqn, aidx);
            if (eidx < 0) {
                continue;
            }
            out.eqs.add(actflowRecord(lqn, aidx, eidx, visits[aidx + 1], s));
        }
        for (int h = 0; h < lqn.nhosts; h++) {
            out.eqs.add(hostutilRecord(lqn, lqn.hshift + h, s));
        }

        out.visits = visits;
        out.convention = conventionText();
        incidences(lqn, out, nidx);
        if (s != null) {
            double worst = Double.NaN;
            for (Relation r : out.eqs) {
                if (r.degenerate || Double.isNaN(r.residual)) {
                    continue;
                }
                double v = Math.abs(r.residual);
                if (Double.isNaN(worst) || v > worst) {
                    worst = v;
                }
            }
            out.maxresidual = worst;
        }
        out.text = report(lqn, out, s != null);
        return out;
    }

    // ====================================================================
    // record construction
    // ====================================================================

    /** A thread pool's case label and its caller classes. */
    private static final class Pool {
        String branch;
        int[] terms;
        boolean[] isentry;
    }

    /**
     * Which thread-pool case task TIDX falls in, and the caller classes of its pool.
     *
     * <p>The term set is uniform across the cases: every request stream that can hold
     * a thread of the task contributes one class. That is each CALL targeting one of
     * its entries, plus each entry carrying an OPEN ARRIVAL, plus, for a reference
     * task, its own entries, since the cycle of a reference task closes on itself and
     * no layer above drives it. The branch label follows the case analysis of
     * {@code updateThinkTimes}, which is what decides whether the utilization is a job
     * count (infinite server) or is normalized to [0,1] (every other discipline).</p>
     */
    private static Pool threadPoolBranch(LayeredNetworkStruct lqn, int tidx, List<Integer> cin) {
        List<Integer> ents = entriesOf(lqn, tidx);
        List<Integer> arv = new ArrayList<Integer>();
        for (int eidx : ents) {
            if (arrivalRate(lqn, eidx) > 0) {
                arv.add(eidx);
            }
        }
        List<Integer> calls = (cin == null) ? new ArrayList<Integer>() : cin;
        List<Integer> terms = new ArrayList<Integer>();
        List<Boolean> isentry = new ArrayList<Boolean>();
        for (int c : calls) {
            terms.add(c);
            isentry.add(Boolean.FALSE);
        }
        for (int e : arv) {
            terms.add(e);
            isentry.add(Boolean.TRUE);
        }
        Pool p = new Pool();
        if (at(lqn.isref, tidx) != 0) {
            p.branch = "ref";
            for (int e : ents) {
                if (!arv.contains(e)) {
                    terms.add(e);
                    isentry.add(Boolean.TRUE);
                }
            }
            p.terms = toIntArray(terms);
            p.isentry = toBoolArray(isentry);
            return p;
        }
        boolean blocking = false;
        for (int c : calls) {
            if (lqn.calltype.get(c) != CallType.FWD) {
                blocking = true;
                break;
            }
        }
        if (blocking) {
            p.branch = (lqn.sched.get(tidx) == SchedStrategy.INF) ? "inf" : "queueing";
        } else if (!calls.isEmpty()) {
            p.branch = "fwd";       // a forwarded request holds a thread too
        } else if (!arv.isEmpty()) {
            p.branch = "arrival";
        } else {
            return null;            // no caller, no arrival: no cycle to close
        }
        p.terms = toIntArray(terms);
        p.isentry = toBoolArray(isentry);
        return p;
    }

    private static Relation littleRecord(LayeredNetworkStruct lqn, int tidx, Pool pool, Sol s) {
        Relation r = new Relation();
        r.kind = "little";
        r.branch = pool.branch;
        r.target = tidx + 1;
        r.targetname = elemName(lqn, tidx);
        r.terms = plusOne(pool.terms);
        r.termisentry = pool.isentry;
        r.coeff = ones(pool.terms.length);
        r.mult = at(lqn.mult, tidx);
        r.maxmult = at(lqn.maxmult, tidx);
        r.repl = Math.max(1.0, nanTo(at(lqn.repl, tidx), 1.0));
        r.rhsconst = r.mult;
        r.scaled = !"inf".equals(pool.branch);
        r.phase2 = hasPhase2(lqn, tidx);
        r.setup = at(lqn.hassetup, tidx) != 0;
        r.degenerate = !isFinite(r.rhsconst) || pool.terms.length == 0;

        double z = refThinkTime(lqn, tidx);
        String[] names = new String[pool.terms.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = termName(lqn, pool.isentry[i], pool.terms[i]);
        }
        r.text = littleText(r, z, names);

        if (s != null) {
            double X = nanTo(at(s.tput, tidx), 0.0) / r.repl;
            double[] B = new double[pool.terms.length];
            double sumB = 0;
            for (int i = 0; i < B.length; i++) {
                if (pool.isentry[i]) {
                    B[i] = nanTo(at(s.tput, pool.terms[i]), 0.0)
                            * nanTo(at(s.servt, pool.terms[i]), 0.0);
                } else {
                    int dst = (int) lqn.callpair.get(pool.terms[i], CALL_DST);
                    B[i] = nanTo(s.calltput[pool.terms[i]], 0.0) * nanTo(at(s.servt, dst), 0.0);
                }
                sumB += B[i];
            }
            double zt = nanTo(at(s.thinkt, tidx), 0.0);
            r.lhs = X * (zt + z) + sumB;
            r.rhs = r.rhsconst;
            r.residual = r.lhs - r.rhs;
            r.relresidual = r.residual / Math.max(Math.abs(r.rhs), 1e-12);
            r.clamped = isFinite(r.rhsconst) && (r.rhsconst - sumB - X * z) < 0;
            r.perclassutil = new double[B.length];
            boolean jobUnits = lqn.sched.get(tidx) == SchedStrategy.INF
                    || !isFinite(r.mult) || r.mult <= 0;
            for (int i = 0; i < B.length; i++) {
                r.perclassutil[i] = jobUnits ? B[i] : B[i] / r.mult;
            }
        }
        return r;
    }

    private static String littleText(Relation r, double z, String[] names) {
        String X = "X(" + r.targetname + ")";
        String Z = "Z(" + r.targetname + ")";
        StringBuilder lhs = new StringBuilder(X + "*(" + Z + (z > 0 ? " + " + num(z) : "") + ")");
        for (String nm : names) {
            lhs.append(" + B(").append(nm).append(")");
        }
        StringBuilder txt = new StringBuilder(lhs + " = " + num(r.rhsconst));
        if (r.scaled && names.length > 0) {
            StringBuilder us = new StringBuilder();
            for (String nm : names) {
                us.append(" - U(").append(r.targetname).append(",").append(nm).append(")");
            }
            txt.append("\n           equivalently  ").append(X).append("*").append(Z)
                    .append(" = ").append(num(r.rhsconst)).append("*(1").append(us).append(")");
        }
        return txt.toString();
    }

    private static Relation callflowRecord(LayeredNetworkStruct lqn, int cidx, Sol s) {
        Relation r = new Relation();
        r.kind = "callflow";
        r.branch = callTypeName(lqn.calltype.get(cidx));
        r.target = lqn.cshift + cidx + 1;
        r.targetname = callName(lqn, cidx);
        int src = (int) lqn.callpair.get(cidx, CALL_SRC);
        double y = mapAt(lqn.callproc_mean, cidx, 0.0);
        r.terms = new int[]{src + 1};
        r.coeff = new double[]{y};
        r.degenerate = (y == 0.0);
        r.text = "X(" + r.targetname + ") = X(" + elemName(lqn, src) + ") * " + num(y);
        if (s != null) {
            r.lhs = nanTo(s.calltput[cidx], 0.0);
            r.rhs = nanTo(at(s.tput, src), 0.0) * y;
            r.residual = r.lhs - r.rhs;
            r.relresidual = r.residual / Math.max(Math.abs(r.rhs), 1e-12);
        }
        return r;
    }

    private static Relation entryflowRecord(LayeredNetworkStruct lqn, int eidx, int[] inc,
                                            double lam, Sol s) {
        Relation r = new Relation();
        r.kind = "entryflow";
        r.target = eidx + 1;
        r.targetname = elemName(lqn, eidx);
        r.terms = plusOne(inc);
        r.coeff = ones(inc.length);
        r.rhsconst = lam;
        StringBuilder txt = new StringBuilder("X(" + r.targetname + ") =");
        for (int i = 0; i < inc.length; i++) {
            txt.append(i == 0 ? " X(" : " + X(").append(callName(lqn, inc[i])).append(")");
        }
        if (lam > 0) {
            txt.append(inc.length == 0 ? " " : " + ").append(num(lam)).append("   (open arrival)");
        }
        r.text = txt.toString();
        if (s != null) {
            r.lhs = nanTo(at(s.tput, eidx), 0.0);
            double rhs = lam;
            for (int c : inc) {
                rhs += nanTo(s.calltput[c], 0.0);
            }
            r.rhs = rhs;
            r.residual = r.lhs - r.rhs;
            r.relresidual = r.residual / Math.max(Math.abs(r.rhs), 1e-12);
        }
        return r;
    }

    private static Relation actflowRecord(LayeredNetworkStruct lqn, int aidx, int eidx,
                                          double v, Sol s) {
        Relation r = new Relation();
        r.kind = "actflow";
        r.target = aidx + 1;
        r.targetname = elemName(lqn, aidx);
        r.terms = new int[]{eidx + 1};
        r.coeff = new double[]{v};
        r.degenerate = !isFinite(v);
        r.text = "X(" + r.targetname + ") = X(" + elemName(lqn, eidx) + ") * " + num(v);
        if (s != null) {
            r.lhs = nanTo(at(s.tput, aidx), 0.0);
            r.rhs = nanTo(at(s.tput, eidx), 0.0) * v;
            r.residual = r.lhs - r.rhs;
            r.relresidual = r.residual / Math.max(Math.abs(r.rhs), 1e-12);
        }
        return r;
    }

    private static Relation hostutilRecord(LayeredNetworkStruct lqn, int hidx, Sol s) {
        Relation r = new Relation();
        r.kind = "hostutil";
        r.target = hidx + 1;
        r.targetname = elemName(lqn, hidx);
        r.mult = at(lqn.mult, hidx);
        r.repl = Math.max(1.0, nanTo(at(lqn.repl, hidx), 1.0));
        r.scaled = lqn.sched.get(hidx) != SchedStrategy.INF;
        List<Integer> acts = new ArrayList<Integer>();
        List<Double> dem = new ArrayList<Double>();
        List<Integer> tasks = lqn.tasksof.get(hidx);
        if (tasks != null) {
            for (int tidx : tasks) {
                List<Integer> as = lqn.actsof.get(tidx);
                if (as == null) {
                    continue;
                }
                for (int aidx : as) {
                    double d = mapAt(lqn.hostdem_mean, aidx, 0.0);
                    if (d == 0.0 || Double.isNaN(d)) {
                        continue;
                    }
                    acts.add(aidx);
                    dem.add(d);
                }
            }
        }
        int[] acts0 = toIntArray(acts);      // struct space, for the reads below
        r.terms = plusOne(acts0);
        r.coeff = toDoubleArray(dem);
        // The server count is the declared multiplicity ALONE, not mult*repl: a
        // replicated host reports its throughputs and its utilization as TOTALS over
        // the copies, so the extra factor would double-count the replication.
        double m = r.mult;
        if (!r.scaled || !isFinite(m)) {
            m = 1.0;
        }
        r.rhsconst = m;
        r.branch = r.scaled ? "queueing" : "inf";
        r.degenerate = acts.isEmpty();
        StringBuilder lhs = new StringBuilder();
        for (int i = 0; i < r.terms.length; i++) {
            if (i > 0) {
                lhs.append(" + ");
            }
            lhs.append("X(").append(elemName(lqn, acts0[i])).append(")*").append(num(r.coeff[i]));
        }
        r.text = (lhs.length() == 0 ? "0" : lhs.toString())
                + " = " + num(m) + "*U(" + r.targetname + ")";
        if (s != null && s.un != null && !Double.isNaN(at(s.un, hidx))) {
            double v = 0;
            for (int i = 0; i < r.terms.length; i++) {
                v += nanTo(at(s.tput, acts0[i]), 0.0) * r.coeff[i];
            }
            r.lhs = v;
            r.rhs = m * at(s.un, hidx);
            r.residual = r.lhs - r.rhs;
            r.relresidual = r.residual / Math.max(Math.abs(r.rhs), 1e-12);
        }
        return r;
    }

    // ====================================================================
    // structural helpers
    // ====================================================================

    /**
     * Expected executions of every activity per invocation of its entry.
     *
     * <p>The activity precedence arcs of a task are a transient Markov chain whose
     * absorbing state is the reply, so the visit counts of the block solve
     * v = e0*(I-P)^-1: a loop back-edge of weight 1-1/count returns count, and an
     * AND-fork row summing above one returns the branching expectation. Call arcs
     * leave the block and drop out.</p>
     */
    private static double[] actVisits(LayeredNetworkStruct lqn) {
        double[] v = new double[lqn.nidx + 1];
        Matrix G = joinScaledGraph(lqn);
        int n = G.getNumRows();
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            List<Integer> as = lqn.actsof.get(tidx);
            if (as == null || as.isEmpty()) {
                continue;
            }
            int[] A = uniqueSorted(as);
            Map<Integer, Integer> pos = new HashMap<Integer, Integer>();
            for (int i = 0; i < A.length; i++) {
                pos.put(A[i], i);
            }
            // (I-P)' as the coefficient matrix, so that solving gives v' from e0'
            Matrix M = new Matrix(A.length, A.length);
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A.length; j++) {
                    double p = (A[i] < n && A[j] < n) ? G.get(A[i], A[j]) : 0.0;
                    M.set(j, i, (i == j ? 1.0 : 0.0) - p);
                }
            }
            for (int eidx : entriesOf(lqn, tidx)) {
                Matrix e0 = new Matrix(A.length, 1);
                boolean any = false;
                for (int i = 0; i < A.length; i++) {
                    if (eidx < n && A[i] < n && G.get(eidx, A[i]) != 0) {
                        e0.set(i, 0, 1.0);
                        any = true;
                    }
                }
                if (!any) {
                    continue;
                }
                Matrix x = new Matrix(A.length, 1);
                if (!Matrix.solve(M, e0, x)) {
                    continue;
                }
                for (int i = 0; i < A.length; i++) {
                    v[A[i] + 1] += x.get(i, 0);
                }
            }
        }
        return v;
    }

    /**
     * The precedence graph with the arcs into every AND-join target divided by the
     * number of branches the join waits for.
     *
     * <p>An AND-JOIN is the one place where flow does not add up: its target executes
     * ONCE per fork, not once per branch, so summing the inbound arcs would count it
     * as many times as there are branches. Dividing recovers the rate of one branch
     * exactly when the branches carry equal rate, the case for a well-formed
     * fork/join block. {@code actpretype} marks the joined PREDECESSORS, so a join
     * target is any successor of one.</p>
     */
    private static Matrix joinScaledGraph(LayeredNetworkStruct lqn) {
        Matrix G = lqn.graph.copy();
        if (lqn.actpretype == null) {
            return G;
        }
        int n = G.getNumRows();
        List<Integer> andpre = new ArrayList<Integer>();
        for (int i = 0; i < lqn.nidx && i < n; i++) {
            if ((int) at(lqn.actpretype, i) == ID_PRE_AND) {
                andpre.add(i);
            }
        }
        if (andpre.isEmpty()) {
            return G;
        }
        for (int j = 0; j < G.getNumCols(); j++) {
            List<Integer> joined = new ArrayList<Integer>();
            for (int i : andpre) {
                if (G.get(i, j) != 0) {
                    joined.add(i);
                }
            }
            if (joined.size() > 1) {
                for (int i : joined) {
                    G.set(i, j, G.get(i, j) / joined.size());
                }
            }
        }
        return G;
    }

    /** Calls targeting each task, i.e. the caller classes of its thread pool. */
    private static Map<Integer, List<Integer>> callsInto(LayeredNetworkStruct lqn) {
        Map<Integer, List<Integer>> out = new HashMap<Integer, List<Integer>>();
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            int dst = (int) lqn.callpair.get(cidx, CALL_DST);
            if (dst < 0 || dst >= lqn.nidx) {
                continue;
            }
            int tidx = (int) at(lqn.parent, dst);
            List<Integer> l = out.get(tidx);
            if (l == null) {
                l = new ArrayList<Integer>();
                out.put(tidx, l);
            }
            l.add(cidx);
        }
        return out;
    }

    private static int[] incomingCalls(LayeredNetworkStruct lqn, int eidx) {
        List<Integer> inc = new ArrayList<Integer>();
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if ((int) lqn.callpair.get(cidx, CALL_DST) == eidx) {
                inc.add(cidx);
            }
        }
        return toIntArray(inc);
    }

    private static double arrivalRate(LayeredNetworkStruct lqn, int eidx) {
        double m = mapAt(lqn.arrival_mean, eidx, Double.NaN);
        if (!Double.isNaN(m) && isFinite(m) && m > 1e-12) {
            return 1.0 / m;
        }
        return 0.0;
    }

    /**
     * Entry an activity belongs to: the one whose activity block contains it. An
     * activity shared by two entries is attributed to the first, matching how
     * actVisits accumulates its visit counts.
     */
    private static int entryOfActivity(LayeredNetworkStruct lqn, int aidx) {
        int tidx = (int) at(lqn.parent, aidx);
        for (int eidx : entriesOf(lqn, tidx)) {
            List<Integer> as = lqn.actsof.get(eidx);
            if (as != null && as.contains(aidx)) {
                return eidx;
            }
        }
        return -1;
    }

    private static boolean hasPhase2(LayeredNetworkStruct lqn, int tidx) {
        if (lqn.actphase == null) {
            return false;
        }
        for (int eidx : entriesOf(lqn, tidx)) {
            List<Integer> as = lqn.actsof.get(eidx);
            if (as == null) {
                continue;
            }
            for (int aidx : as) {
                int a = aidx - lqn.ashift;
                if (a >= 0 && a < lqn.nacts && lqn.actphase.get(0, a) > 1) {
                    return true;
                }
            }
        }
        return false;
    }

    private static void incidences(LayeredNetworkStruct lqn, Result out, int nidx) {
        out.A_little = new Matrix(lqn.ntasks + 1, lqn.ncalls + 1);
        out.A_flow = new Matrix(lqn.ncalls + 1, nidx + 1);
        out.A_host = new Matrix(lqn.nhosts + 1, nidx + 1);
        for (Relation r : out.eqs) {
            if ("little".equals(r.kind)) {
                // the call classes only; an entry class (open arrival, or the
                // self-driven cycle of a reference task) is not a call
                for (int i = 0; i < r.terms.length; i++) {
                    if (!r.termisentry[i]) {
                        out.A_little.set(r.target - lqn.tshift, r.terms[i], 1.0);
                    }
                }
            } else if ("callflow".equals(r.kind)) {
                out.A_flow.set(r.target - lqn.cshift, r.terms[0], r.coeff[0]);
            } else if ("hostutil".equals(r.kind)) {
                for (int i = 0; i < r.terms.length; i++) {
                    out.A_host.set(r.target - lqn.hshift, r.terms[i], r.coeff[i]);
                }
            }
        }
    }

    // ====================================================================
    // solution adapter
    // ====================================================================

    /** The solution vectors plus the per-call rates derived from them. */
    private static final class Sol {
        Matrix tput;
        Matrix util;
        Matrix thinkt;
        Matrix servt;
        Matrix residt;
        Matrix un;
        double[] calltput;
    }

    /**
     * The caller's 1-based {@link Solution} read into the struct's own 0-based element
     * space, so that every internal lookup indexes {@code lqn} and the iterates with the
     * same number. This is the ONLY place the two spaces meet on the way in.
     */
    private static Sol readSolution(LayeredNetworkStruct lqn, Solution sol) {
        Sol s = new Sol();
        s.tput = shiftDown(sol.tput);
        s.util = shiftDown(sol.util);
        s.thinkt = shiftDown(sol.thinkt);
        s.servt = shiftDown(sol.servt);
        s.residt = shiftDown(sol.residt);
        s.un = shiftDown(sol.un);
        // A call throughput is not an iterate of SolverLN: a call inherits the rate of
        // its dispatching element scaled by the mean call count.
        s.calltput = new double[lqn.ncalls];
        Arrays.fill(s.calltput, Double.NaN);
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            int src = (int) lqn.callpair.get(cidx, CALL_SRC);
            double y = mapAt(lqn.callproc_mean, cidx, 0.0);
            if (src >= 0 && src < lqn.nidx) {
                s.calltput[cidx] = at(s.tput, src) * y;
            }
        }
        return s;
    }

    // ====================================================================
    // reporting
    // ====================================================================

    private static List<String> conventionText() {
        List<String> t = new ArrayList<String>();
        t.add("Conventions:");
        t.add("  little   : rates and populations PER REPLICA (X = tput/repl, N = mult of one copy).");
        t.add("             B(t,k) is the per-class utilization in JOB units; for a queueing task");
        t.add("             B = mult*U with U in [0,1], at an infinite server B = U directly.");
        t.add("             S(k) is the entry SERVICE time (phase 1 + phase 2), the thread hold time.");
        t.add("  other    : throughputs and utilizations as the solver reports them, totalled over replicas.");
        t.add("  N(t)     : lqn.mult. SolverLN iterates on njobs (interlocking corrections, maxmult under");
        t.add("             replication), reported per record as mult/maxmult.");
        return t;
    }

    private static final String[][] TITLES = {
            {"little", "thread-pool Little's law"},
            {"callflow", "call-flow balance"},
            {"entryflow", "entry-flow balance"},
            {"actflow", "activity-flow balance"},
            {"hostutil", "host utilization law"},
    };

    private static List<String> report(LayeredNetworkStruct lqn, Result out, boolean hasSol) {
        List<String> txt = new ArrayList<String>();
        txt.add(String.format("LQN balance equations: %d hosts, %d tasks, %d entries, %d activities, %d calls",
                lqn.nhosts, lqn.ntasks, lqn.nentries, lqn.nacts, lqn.ncalls));
        txt.addAll(out.convention);
        for (String[] kt : TITLES) {
            List<Integer> sel = new ArrayList<Integer>();
            for (int i = 0; i < out.eqs.size(); i++) {
                if (kt[0].equals(out.eqs.get(i).kind)) {
                    sel.add(i);
                }
            }
            if (sel.isEmpty()) {
                continue;
            }
            txt.add("");
            txt.add(String.format("--- %s  (kind='%s', %d relations) ---", kt[1], kt[0], sel.size()));
            for (int i : sel) {
                Relation r = out.eqs.get(i);
                StringBuilder head = new StringBuilder(String.format("[%3d] %-14s", i + 1, r.targetname));
                if ("little".equals(r.kind) || "hostutil".equals(r.kind)) {
                    head.append(String.format(" %-9s mult=%s repl=%s", r.branch, num(r.mult), num(r.repl)));
                } else if (!r.branch.isEmpty()) {
                    head.append(String.format(" %-9s", r.branch));
                }
                txt.add(head.toString());
                for (String line : r.text.split("\n")) {
                    txt.add("       " + line);
                }
                List<String> notes = new ArrayList<String>();
                if (r.phase2) {
                    notes.add("phase-2 tail on Z");
                }
                if (r.setup) {
                    notes.add("setup charge on Z");
                }
                if (r.degenerate) {
                    notes.add("DEGENERATE (not usable as a residual)");
                }
                if (r.clamped) {
                    notes.add("SATURATED (equality unattainable, use a hinge)");
                }
                if (hasSol && !r.degenerate && Double.isNaN(r.residual)) {
                    notes.add("not instantiated (pass un for the host utilization law)");
                }
                if (!notes.isEmpty()) {
                    StringBuilder nb = new StringBuilder("       note: ");
                    for (int k = 0; k < notes.size(); k++) {
                        if (k > 0) {
                            nb.append(", ");
                        }
                        nb.append(notes.get(k));
                    }
                    txt.add(nb.toString());
                }
                if (hasSol && !r.degenerate && !Double.isNaN(r.residual)) {
                    txt.add(String.format("       lhs=%.6g  rhs=%.6g  residual=%.3e  rel=%.3e",
                            r.lhs, r.rhs, r.residual, r.relresidual));
                }
            }
        }
        if (hasSol) {
            int nd = 0;
            for (Relation r : out.eqs) {
                if (!r.degenerate) {
                    nd++;
                }
            }
            txt.add("");
            txt.add(String.format("max |residual| over %d non-degenerate relations: %.3e",
                    nd, out.maxresidual));
        }
        return txt;
    }

    // ====================================================================
    // small utilities
    // ====================================================================

    /**
     * Element of a struct vector by absolute index, whatever orientation it is stored
     * in: the LayeredNetworkStruct keeps some of them as rows and some as columns.
     */
    private static double at(Matrix m, int i) {
        if (m == null || i < 0) {
            return Double.NaN;
        }
        if (m.getNumRows() == 1) {
            return i < m.getNumCols() ? m.get(0, i) : Double.NaN;
        }
        if (m.getNumCols() == 1) {
            return i < m.getNumRows() ? m.get(i, 0) : Double.NaN;
        }
        return Double.NaN;
    }

    private static double mapAt(Map<Integer, Double> m, int i, double dflt) {
        if (m == null) {
            return dflt;
        }
        Double v = m.get(i);
        return (v == null || Double.isNaN(v)) ? dflt : v;
    }

    private static List<Integer> entriesOf(LayeredNetworkStruct lqn, int tidx) {
        List<Integer> e = lqn.entriesof.get(tidx);
        return (e == null) ? new ArrayList<Integer>() : e;
    }

    private static double refThinkTime(LayeredNetworkStruct lqn, int tidx) {
        if (at(lqn.isref, tidx) == 0) {
            return 0.0;
        }
        double z = mapAt(lqn.think_mean, tidx, 0.0);
        return (isFinite(z) && z >= 0) ? z : 0.0;
    }

    private static String elemName(LayeredNetworkStruct lqn, int idx) {
        String nm = (lqn.hashnames == null) ? null : lqn.hashnames.get(idx);
        return (nm == null || nm.isEmpty()) ? ("#" + idx) : nm;
    }

    private static String callName(LayeredNetworkStruct lqn, int cidx) {
        String nm = (lqn.callhashnames == null) ? null : lqn.callhashnames.get(cidx);
        return (nm == null || nm.isEmpty()) ? ("C" + cidx) : nm;
    }

    private static String termName(LayeredNetworkStruct lqn, boolean isEntry, int idx) {
        return isEntry ? elemName(lqn, idx) : callName(lqn, idx);
    }

    private static String callTypeName(CallType ct) {
        if (ct == CallType.SYNC) {
            return "sync";
        }
        if (ct == CallType.ASYNC) {
            return "async";
        }
        if (ct == CallType.FWD) {
            return "fwd";
        }
        return "";
    }

    private static String num(double x) {
        if (Double.isInfinite(x)) {
            return x > 0 ? "Inf" : "-Inf";
        }
        if (Double.isNaN(x)) {
            return "NaN";
        }
        if (x == Math.rint(x) && Math.abs(x) < 1e15) {
            return Long.toString((long) x);
        }
        return trimZeros(String.format("%.6g", x));
    }

    private static String trimZeros(String s) {
        if (s.indexOf('.') < 0 || s.indexOf('e') >= 0 || s.indexOf('E') >= 0) {
            return s;
        }
        int end = s.length();
        while (end > 0 && s.charAt(end - 1) == '0') {
            end--;
        }
        if (end > 0 && s.charAt(end - 1) == '.') {
            end--;
        }
        return s.substring(0, end);
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static double nanTo(double x, double dflt) {
        return Double.isNaN(x) ? dflt : x;
    }

    /**
     * A 1-based-over-elements row vector read back in the struct's 0-based space: slot 0
     * of the input is the unused pad, so element {@code i} of the result is slot
     * {@code i+1} of the input.
     */
    private static Matrix shiftDown(Matrix v) {
        if (v == null) {
            return null;
        }
        boolean row = v.getNumRows() == 1;
        int n = row ? v.getNumCols() : v.getNumRows();
        if (n <= 1) {
            return null;   // nothing but the pad: the quantity was not supplied
        }
        Matrix out = new Matrix(1, n - 1);
        for (int i = 0; i + 1 < n; i++) {
            out.set(0, i, row ? v.get(0, i + 1) : v.get(i + 1, 0));
        }
        return out;
    }

    /** Struct indices as the 1-based indices a {@link Relation} exposes. */
    private static int[] plusOne(int[] a) {
        int[] b = new int[a.length];
        for (int i = 0; i < a.length; i++) {
            b[i] = a[i] + 1;
        }
        return b;
    }

    private static double[] ones(int n) {
        double[] a = new double[n];
        Arrays.fill(a, 1.0);
        return a;
    }

    private static int[] toIntArray(List<Integer> l) {
        int[] a = new int[l.size()];
        for (int i = 0; i < a.length; i++) {
            a[i] = l.get(i);
        }
        return a;
    }

    private static double[] toDoubleArray(List<Double> l) {
        double[] a = new double[l.size()];
        for (int i = 0; i < a.length; i++) {
            a[i] = l.get(i);
        }
        return a;
    }

    private static boolean[] toBoolArray(List<Boolean> l) {
        boolean[] a = new boolean[l.size()];
        for (int i = 0; i < a.length; i++) {
            a[i] = l.get(i);
        }
        return a;
    }

    private static int[] uniqueSorted(List<Integer> l) {
        int[] a = toIntArray(l);
        Arrays.sort(a);
        int n = 0;
        for (int i = 0; i < a.length; i++) {
            if (i == 0 || a[i] != a[i - 1]) {
                a[n++] = a[i];
            }
        }
        return Arrays.copyOf(a, n);
    }
}
