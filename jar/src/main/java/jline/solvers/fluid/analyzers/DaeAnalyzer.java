/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import jline.GlobalConstants;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.moments.FluidLyapunov;
import jline.solvers.fluid.moments.FluidMomentTerms;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.ode.Rodas;
import org.apache.commons.math3.util.FastMath;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * Differential-algebraic formulation of the min-normal closure, backing
 * {@code options.method='dae'}. Java twin of MATLAB's {@code solver_fluid_dae}.
 *
 * <p>{@link MinNormalAnalyzer} already solves a differential system (the mean)
 * coupled to an algebraic one (the covariance). It solves them by SUCCESSIVE
 * SUBSTITUTION: integrate the mean to its fixed point at a held variance, solve
 * the Lyapunov equation there, extract sigma2, repeat, up to 20 times and only
 * to CoarseTol. This analyzer states the same closure as one system and solves
 * it as one system. NOTHING ABOUT THE CLOSURE CHANGES -- the drift, the rate
 * factors and the Lyapunov equation are taken unmodified from
 * {@link FluidMomentTerms} and {@link FluidLyapunov} -- only the way the coupled
 * equations are discharged.</p>
 *
 * <p>Two modes, chosen by the horizon:</p>
 *
 * <p>STEADY STATE (unbounded horizon, the usual case) solves</p>
 * <pre>
 *     0 = D r(x, sigma2)               drift residual, one row per state
 *     0 = C x - Nchain                 population conservation, per closed chain
 *     0 = sigma2 - sigmaOf(x, sigma2)  closure consistency, per closable station
 * </pre>
 * <p>simultaneously by a damped projected Newton. THE FIXED POINT IS NOT FOUND
 * BY INTEGRATING TO IT: integrating a stable ODE until it stops moving is a poor
 * way to solve f(x)=0, because the cost is set by the slowest mode of the model
 * rather than by the accuracy wanted. One seed trajectory is still integrated,
 * cheaply and at the first-order closure, because Newton needs a point inside
 * the basin; everything after that is algebraic.</p>
 *
 * <p>TRANSIENT (finite horizon) integrates the index-1 DAE with a SINGULAR mass
 * matrix, the covariance advancing alongside the mean. This and {@code kp} are
 * the only fluid methods producing a time-varying second moment; {@code
 * minnormal} evaluates its whole transient at the single stationary variance.
 * The integrator is {@link Rodas}, because LSODA solves y' = f and cannot carry
 * a mass matrix at all.</p>
 *
 * <p>WHAT THE ALGEBRAIC CONSTRAINT BUYS. Population conservation otherwise holds
 * only to integrator tolerance: it is a CONSEQUENCE of the drift (the rows of D
 * sum to zero on a closed chain), never an equation. Writing it as a constraint
 * enforces it to solver tolerance, and it is also what makes the Newton system
 * solvable at all -- the drift Jacobian is singular along exactly the conserved
 * directions, the same singularity {@link FluidLyapunov} works around by
 * projecting onto range(D), so the constraint rows supply the missing rank
 * instead of a pseudo-inverse hiding it.</p>
 *
 * <p>WHY THE COVARIANCE IS NOT A NEWTON UNKNOWN. Sigma is nstate^2 entries, so a
 * Jacobian over it is quartic work -- strictly worse than the cubic Lyapunov
 * solves it would replace. Sigma is LINEAR in itself for a held x, so it is
 * eliminated by one Lyapunov solve per residual evaluation and only sigma2, M
 * numbers, joins x in the unknown vector.</p>
 *
 * <p>Java 8: no {@code var}, no {@code List.of}, no switch expressions.</p>
 *
 * @see FluidMomentTerms
 * @see Rodas
 */
public class DaeAnalyzer extends ClosingAndStateDepMethodsAnalyzer {

    /** State-level stationary covariance at the converged fixed point. */
    public Matrix sigmaMatrix;
    /** Per-(station,class) queue-length variance. */
    public Matrix qVar;
    /** State coordinates of each (station,class), flattened as {@code i*K+k}. */
    public int[][] classBlock;
    /** State coordinates of each station. */
    public int[][] stationBlock;
    /** Per-station population variance. */
    public double[] sigma2;
    /** The same variance as it enters the DRIFT: zero where nothing is closed. */
    public double[] sigma2Drift;
    /** Newton steps taken, summed over the active-set passes. */
    public int outerIters;
    /** Infinity norm of the residual the solve stopped at. */
    public double residual = Double.POSITIVE_INFINITY;
    /** Whether that residual reached {@code options.tol}. */
    public boolean converged;
    /** Largest violation of the conservation rows at the reported point. */
    public double conservation;

    /** Transient state covariance, one matrix per reported time. */
    public Matrix[] sigmat;
    /** Transient per-(station,class) variance, flattened as {@code i*K+k}. */
    public Matrix[] qVart;
    /** The time grid the transient covariance is reported on. */
    public Matrix tvar;

    /** Human-readable origin of each capacity constraint. */
    public String[] capacityLabel;
    /** Right-hand side of each capacity constraint. */
    public Matrix capacityB;
    /** Value each constrained quantity actually took. */
    public Matrix capacityValue;
    /** Which constraints bound at the reported point. */
    public int[] capacityActive;
    /** Mass held OUTSIDE a capped region, per staging coordinate. */
    public Matrix staging;
    /** Region each staging coordinate belongs to. */
    public int[] stagingRegion;
    /** Class each staging coordinate carries. */
    public int[] stagingClass;
    /** Total blocked mass; station queues sum to N minus this. */
    public double blocked;
    /** Rate each waiting queue drains at, one per binding cap. */
    public Matrix drain;

    /** True where a cap stages the blocked job rather than holding or losing it. */
    public boolean[] capacityStaged;

    /** Which region each cap came from, -1 for a station buffer. */
    public int[] capacityRegion;

    /** Which station each cap came from, -1 for a region cap. */
    public int[] capacityStation;

    /**
     * Every time the trajectory made a cap start or stop binding, as
     * {t, row, kind} with kind 0 = release, 1 = activate, 2 = a crossing the cap
     * could not hold. Empty for a steady-state solve.
     */
    public List<double[]> capacitySwitches = new ArrayList<double[]>();

    private static final double ZERO = 1e-14;
    private static final double FINE = 1e-8;

    // -----------------------------------------------------------------------
    // finite capacity regions, as linear admission constraints
    // -----------------------------------------------------------------------

    /**
     * {@code A x <= b}, one row per cap, with the region each row came from.
     *
     * <p>EVERY FORM THE REGION CAN TAKE IS ONE FAMILY. LINE stores four
     * different limits -- a region-global job cap, a per-class job cap, a memory
     * budget with per-class sizes, and an arbitrary linear pair (A,b) -- and they
     * are all the same object once written against the state, with {@code Arow}
     * the per-class weight on the region's coordinates and zero outside it.</p>
     */
    static final class CapacityConstraints {
        Matrix A;
        Matrix As;            // the same rows over the staging coordinates
        double[] b;
        int[] region;         // -1 for a station cap
        int[] station;        // -1 for a region cap
        int[] klassRow;       // -1 when the row limits several classes
        boolean[] staged;     // true where the blocked job waits in a room
        String[] label;
        boolean[][] member;   // (nregions x nstate)
        int[] klass;          // coordinate -> class
        int nregions;

        boolean isEmpty() {
            return b == null || b.length == 0;
        }
    }

    /** The waiting room outside each capped region, as fluid coordinates. */
    static final class CapacityStaging {
        int n;
        int[] region;
        int[] klass;
        boolean[] adm;
        int[] admRegion;
        int[] admStage;
        boolean[][] gatedBy;  // (ncon x n) which rows gate which room
        Matrix Dn;
        Matrix Dp;
    }

    /**
     * Which events each cap throttles, and what happens to the mass it stops.
     *
     * <p>THE THREE ANSWERS ARE THE MODEL, AND THEY ARE NOT INTERCHANGEABLE. LINE's
     * own semantics decides which one a cap gets, and the choice is visible in the
     * answer -- a job held upstream is still counted at that station, a lost job is
     * counted nowhere and breaks flow balance across the cap on purpose, a staged
     * job is counted in neither and is reported as blocked mass:</p>
     *
     * <ul>
     *   <li>STAGED: a finite capacity region under a waiting queue. The job
     *   COMPLETES upstream service and waits outside the region, which is what JMT
     *   and LDES simulate; the upstream station empties exactly as it would with no
     *   region.</li>
     *   <li>HELD: a station buffer reached by a CLOSED class. State.arrivalIsLost
     *   refuses to lose a closed job -- population conservation is a defining
     *   invariant -- and returns an empty successor instead, which DISABLES the
     *   upstream departure until room frees. The job is therefore still at the
     *   upstream station, so the fluid analogue scales the WHOLE event.</li>
     *   <li>LOSS: a station buffer reached by an OPEN class. The same predicate
     *   loses it: the arrival event still fires and only the carried flow is
     *   admitted, so the fluid analogue scales the ARRIVAL leg alone.</li>
     * </ul>
     */
    static final class CapacityGates {
        boolean[][] gate;     // (ncon x nevents)
        boolean[][] held;
        boolean[][] loss;
        Matrix Dn;            // removal at real stations
        Matrix DnExt;         // removal from the EXT source pool
        Matrix Dp;            // arrival
    }

    /** The two legs of every event under the active caps, and the rooms. */
    static final class Legs {
        double[] rup;         // the rate each event FIRES at
        double[] rin;         // the rate mass LANDS at
        double[] ds;          // the derivative of each room
        double[] drain;       // each room's total outflow
    }

    /**
     * Every capacity limit in the model as a linear constraint on the fluid state.
     *
     * <p>TWO DECLARATIONS, ONE FAMILY OF ROWS. A finite capacity region caps a SET
     * of stations jointly (addRegion); a station cap limits one station's own
     * buffer (setCapacity/setClassCapacity). LINE stores the region form four ways
     * -- a region-global job cap, a per-class job cap, a memory budget with
     * per-class sizes, and an arbitrary linear pair (A,b) -- and the station form
     * two, a total and a per-class buffer. All six are the same object once written
     * against the state, {@code Arow x <= b} with Arow the per-class weight on the
     * covered coordinates and zero elsewhere, so the solver carries one mechanism
     * rather than six.</p>
     *
     * <p>WHAT DIFFERS IS NOT THE ROW BUT WHERE THE BLOCKED JOB GOES, which is
     * decided per admission event in {@link #capacityGates}, not here.</p>
     *
     * <p>WHAT IS REFUSED, AND WHY IT IS NOT A CONSTRAINT. Only a waiting queue or a
     * drop is a constraint on this drift. BAS/BBS/RSRD give the upstream station a
     * blocked-server state and the retrial rules add an orbit; each changes the
     * event set itself, so each needs a different drift rather than a constraint on
     * this one, and is refused by name.</p>
     */
    static CapacityConstraints capacityConstraints(NetworkStruct sn, FluidMomentTerms terms) {
        int M = terms.M;
        int K = terms.K;
        int nstate = terms.nstate;
        CapacityConstraints con = new CapacityConstraints();
        con.A = new Matrix(0, nstate);
        con.As = new Matrix(0, 0);
        con.b = new double[0];
        con.region = new int[0];
        con.station = new int[0];
        con.klassRow = new int[0];
        con.staged = new boolean[0];
        con.label = new String[0];
        con.klass = new int[nstate];
        con.nregions = FastMath.max(sn.nregions, 0);
        con.member = new boolean[FastMath.max(con.nregions, 1)][nstate];

        int[] coordStation = new int[nstate];
        int[] coordClass = new int[nstate];
        Arrays.fill(coordStation, -1);
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                int kic = (int) terms.Kic.get(i, c);
                if (kic > 0) {
                    int lo = (int) terms.qIndices.get(i, c);
                    for (int s = lo; s < lo + kic; s++) {
                        coordStation[s] = i;
                        coordClass[s] = c;
                    }
                }
            }
        }
        con.klass = coordClass;

        List<double[]> rows = new ArrayList<double[]>();
        List<Double> bs = new ArrayList<Double>();
        List<Integer> rgs = new ArrayList<Integer>();
        List<Integer> sts = new ArrayList<Integer>();
        List<Integer> kls = new ArrayList<Integer>();
        List<Boolean> sdg = new ArrayList<Boolean>();
        List<String> labels = new ArrayList<String>();

        for (int f = 0; f < con.nregions; f++) {
            boolean[] members = new boolean[M];
            boolean any = false;
            if (sn.regionmembers != null && sn.regionmembers.size() > f
                    && sn.regionmembers.get(f) != null) {
                Matrix mm = sn.regionmembers.get(f);
                for (int i = 0; i < M && i < mm.length(); i++) {
                    members[i] = mm.get(i) != 0;
                    any |= members[i];
                }
            }
            if (!any) {
                continue;
            }
            boolean[] inRegion = new boolean[nstate];
            for (int s = 0; s < nstate; s++) {
                inRegion[s] = coordStation[s] >= 0 && members[coordStation[s]];
            }
            con.member[f] = inRegion;

            // -- refusals, named before anything is built --------------------
            if (sn.regionrule != null && sn.regionrule.getNumRows() > f) {
                for (int r = 0; r < K && r < sn.regionrule.getNumCols(); r++) {
                    if (sn.regionrule.get(f, r) != DropStrategy.WaitingQueue.getID()) {
                        line_error(mfilename(new Object() {
                        }), String.format("Region %d applies a drop rule other than a waiting queue to "
                                + "class %d. Only a waiting queue is a constraint on the fluid drift: it "
                                + "conserves the population and throttles the admission flow. The other "
                                + "rules change the event set instead, so they need a different drift "
                                + "rather than an algebraic equation on this one.", f + 1, r + 1));
                    }
                }
            }
            if (sn.regionweight != null && sn.regionweight.getNumRows() > f) {
                for (int r = 0; r < K && r < sn.regionweight.getNumCols(); r++) {
                    if (FastMath.abs(sn.regionweight.get(f, r) - 1.0) > ZERO) {
                        line_error(mfilename(new Object() {
                        }), String.format("Region %d sets per-class admission weights, which decide WHICH "
                                + "blocked class enters when capacity frees up. The constraint form "
                                + "throttles the admission flow in proportion to its own rate and carries "
                                + "no such priority, so the weights would be ignored silently. Use "
                                + "SolverCTMC, SolverJMT or SolverSSA.", f + 1));
                    }
                }
            }

            Matrix rmat = (sn.region != null && sn.region.size() > f) ? sn.region.get(f) : null;
            int memberRow = 0;
            for (int i = 0; i < M; i++) {
                if (members[i]) {
                    memberRow = i;
                    break;
                }
            }

            // -- 1. region-global job cap ------------------------------------
            if (rmat != null && rmat.getNumCols() >= K + 1) {
                double gcap = rmat.get(memberRow, K);
                if (isFiniteCap(gcap)) {
                    double[] row = new double[nstate];
                    for (int s = 0; s < nstate; s++) {
                        if (inRegion[s]) {
                            row[s] = 1.0;
                        }
                    }
                    addRow(rows, bs, rgs, sts, kls, sdg, labels, row, gcap, f, -1, -1, true,
                            "region " + (f + 1) + " global job cap");
                }
            }
            // -- 2. per-class job caps ---------------------------------------
            if (rmat != null) {
                for (int r = 0; r < K && r < rmat.getNumCols(); r++) {
                    double ccap = rmat.get(memberRow, r);
                    if (!isFiniteCap(ccap)) {
                        continue;
                    }
                    double[] row = new double[nstate];
                    boolean anyc = false;
                    for (int s = 0; s < nstate; s++) {
                        if (inRegion[s] && coordClass[s] == r) {
                            row[s] = 1.0;
                            anyc = true;
                        }
                    }
                    if (anyc) {
                        addRow(rows, bs, rgs, sts, kls, sdg, labels, row, ccap, f, -1, r, true,
                                "region " + (f + 1) + " class " + (r + 1) + " job cap");
                    }
                }
            }
            // -- 3. region-global memory budget ------------------------------
            // A per-class memory limit has already been folded into the per-class
            // job cap by the struct refresh, so only the region-global budget is
            // left: the job-count row with each class weighted by its size.
            if (sn.regionmaxmem != null && sn.regionmaxmem.size() > f
                    && sn.regionmaxmem.get(f) != null) {
                Matrix mem = sn.regionmaxmem.get(f);
                if (memberRow < mem.length()) {
                    double gmem = mem.get(memberRow);
                    if (isFiniteCap(gmem)) {
                        double[] row = new double[nstate];
                        boolean anym = false;
                        for (int s = 0; s < nstate; s++) {
                            if (!inRegion[s]) {
                                continue;
                            }
                            double sz = 1.0;
                            if (sn.regionsz != null && sn.regionsz.getNumRows() > f
                                    && coordClass[s] < sn.regionsz.getNumCols()) {
                                sz = sn.regionsz.get(f, coordClass[s]);
                            }
                            row[s] = sz;
                            anym |= sz != 0;
                        }
                        if (anym) {
                            addRow(rows, bs, rgs, sts, kls, sdg, labels, row, gmem, f, -1, -1, true,
                                    "region " + (f + 1) + " memory budget");
                        }
                    }
                }
            }
            // -- 4. explicit linear constraints ------------------------------
            if (sn.regionlincon != null && sn.regionlincon.containsKey(f)
                    && sn.regionlincon.get(f) != null) {
                MatrixCell pair = sn.regionlincon.get(f);
                if (pair.size() >= 2 && pair.get(0) != null && pair.get(1) != null) {
                    Matrix alin = pair.get(0);
                    Matrix blin = pair.get(1);
                    for (int c = 0; c < alin.getNumRows(); c++) {
                        double[] row = new double[nstate];
                        boolean anyl = false;
                        for (int r = 0; r < K && r < alin.getNumCols(); r++) {
                            double a = alin.get(c, r);
                            if (a == 0) {
                                continue;
                            }
                            for (int s = 0; s < nstate; s++) {
                                if (inRegion[s] && coordClass[s] == r) {
                                    row[s] = a;
                                    anyl = true;
                                }
                            }
                        }
                        if (anyl) {
                            addRow(rows, bs, rgs, sts, kls, sdg, labels, row, blin.get(c, 0), f, -1,
                                    -1, true, "region " + (f + 1) + " linear constraint " + (c + 1));
                        }
                    }
                }
            }
        }

        // -- 5. per-station buffers ------------------------------------------
        // A STATION CAP IS THE ONE-STATION CASE OF THE SAME ROW, and every fluid
        // method other than this one ignores it outright: nothing in the fluid tree
        // reads sn.cap or sn.classcap, so a capped station was integrated as an
        // unbounded one and the table reported more jobs in the buffer than the
        // buffer holds.
        for (int i = 0; i < M; i++) {
            if (terms.isExt[i]) {
                continue;   // a source holds no jobs, so its cap caps nothing
            }
            boolean anyAt = false;
            for (int s = 0; s < nstate; s++) {
                anyAt |= coordStation[s] == i;
            }
            if (!anyAt) {
                continue;
            }
            if (sn.cap != null && i < sn.cap.length()) {
                double gcap = sn.cap.get(i);
                // Integer.MAX_VALUE is how the JAR spells "no buffer at all": a row
                // built from it can never bind, and it would still cost an
                // active-set test on every pass and appear in the report.
                if (Double.isFinite(gcap) && gcap >= 0 && gcap < Integer.MAX_VALUE) {
                    double[] row = new double[nstate];
                    for (int s = 0; s < nstate; s++) {
                        if (coordStation[s] == i) {
                            row[s] = 1.0;
                        }
                    }
                    addRow(rows, bs, rgs, sts, kls, sdg, labels, row, gcap, -1, i, -1, false,
                            "station " + (i + 1) + " buffer");
                }
            }
            if (sn.classcap != null && i < sn.classcap.getNumRows()) {
                for (int r = 0; r < K && r < sn.classcap.getNumCols(); r++) {
                    double ccap = sn.classcap.get(i, r);
                    if (!Double.isFinite(ccap) || ccap < 0 || ccap >= Integer.MAX_VALUE) {
                        continue;
                    }
                    double[] row = new double[nstate];
                    boolean anyc = false;
                    for (int s = 0; s < nstate; s++) {
                        if (coordStation[s] == i && coordClass[s] == r) {
                            row[s] = 1.0;
                            anyc = true;
                        }
                    }
                    if (anyc) {
                        addRow(rows, bs, rgs, sts, kls, sdg, labels, row, ccap, -1, i, r, false,
                                "station " + (i + 1) + " class " + (r + 1) + " buffer");
                    }
                }
            }
        }

        int ncand = rows.size();
        boolean[] keep = new boolean[ncand];
        Arrays.fill(keep, true);

        // A CAP THE POPULATION CANNOT REACH IS NOT A CAP, and dropping it here keeps
        // it out of the active-set loop, where it would be tested on every pass and
        // never bind while its multiplier stayed an unknown with no equation to pin
        // it. THE BOUND IS PER CLASS, and it has to be: the heaviest weight times
        // the whole population never prunes a PER-CLASS cap set to its own class
        // population -- exactly the row refreshCapacity derives at every station of
        // every closed model. An open class makes the reachable total unbounded.
        for (int c = 0; c < ncand; c++) {
            if (reach(rows.get(c), coordClass, coordStation, sn.njobs, K, nstate)
                    <= bs.get(c) + ZERO) {
                keep[c] = false;
            }
        }

        // A CAP ON COORDINATES THE IMMEDIATE REDUCTION FOLDED AWAY IS VACUOUS, NOT
        // MALFORMED. hide_immediate stochastic-complements an Immediate-rate
        // coordinate out of the event set -- the MMT transform's zero-service Join is
        // one, until the fork-join fixed point gives it a synchronisation delay -- and
        // the reduced drift then holds no mass there and has no event landing on it.
        // Left in, such a row reaches capacityGates with no gating event and is
        // reported as a limit the model cannot approach, which is the right message
        // for a station that really does hold jobs and the wrong one here: this
        // station holds none, so its buffer is satisfied identically.
        boolean[] gone = eliminatedCoords(terms, nstate);
        for (int c = 0; c < ncand; c++) {
            if (!keep[c]) {
                continue;
            }
            boolean live = false;
            double[] row = rows.get(c);
            for (int s = 0; s < nstate && !live; s++) {
                if (row[s] != 0.0 && !gone[s]) {
                    live = true;
                }
            }
            if (!live) {
                keep[c] = false;
            }
        }

        // THE SAME ROW TWICE IS A SINGULAR NEWTON SYSTEM, not a redundancy the least
        // squares absorbs: two identical rows both bind, each takes a multiplier,
        // and nothing distinguishes them. refreshCapacity derives classcap from cap,
        // so a single-class model declares the station total and the class buffer as
        // the same row. The TIGHTER bound survives; on a tie the region row does.
        for (int c = 0; c < ncand; c++) {
            if (!keep[c]) {
                continue;
            }
            for (int d = c + 1; d < ncand; d++) {
                if (!keep[d] || !sameRow(rows.get(c), rows.get(d), nstate)) {
                    continue;
                }
                boolean takeover = bs.get(d) < bs.get(c) - ZERO
                        || (FastMath.abs(bs.get(d) - bs.get(c)) <= ZERO && sdg.get(d) && !sdg.get(c));
                if (takeover) {
                    bs.set(c, bs.get(d));
                    rgs.set(c, rgs.get(d));
                    sts.set(c, sts.get(d));
                    kls.set(c, kls.get(d));
                    sdg.set(c, sdg.get(d));
                    labels.set(c, labels.get(d));
                }
                keep[d] = false;
            }
        }

        // A ROW ITS OWN PER-CLASS ROWS ALREADY IMPLY IS RANK, NOT INFORMATION. LINE
        // derives the station total from the per-class buffers, so a two-class
        // station capped 3 and 3 declares a total of 6 as well -- exactly the sum of
        // the two class rows. All three then bind together with rank 2, and the
        // multipliers are one arbitrary point of a line of solutions. Only an exact
        // implication is dropped, so a total TIGHTER than the sum of its parts
        // survives.
        for (int c = 0; c < ncand; c++) {
            if (!keep[c] || kls.get(c) != -1) {
                continue;
            }
            List<Integer> classes = new ArrayList<Integer>();
            for (int r = 0; r < K; r++) {
                for (int s = 0; s < nstate; s++) {
                    if (coordClass[s] == r && coordStation[s] >= 0 && rows.get(c)[s] > 0) {
                        classes.add(r);
                        break;
                    }
                }
            }
            if (classes.isEmpty()) {
                continue;
            }
            boolean implied = true;
            double budget = 0;
            for (int ci = 0; ci < classes.size() && implied; ci++) {
                int r = classes.get(ci);
                int part = -1;
                for (int d = 0; d < ncand; d++) {
                    if (!keep[d] || d == c || kls.get(d) != r) {
                        continue;
                    }
                    if (!rgs.get(d).equals(rgs.get(c)) || !sts.get(d).equals(sts.get(c))) {
                        continue;
                    }
                    boolean covers = true;
                    for (int s = 0; s < nstate; s++) {
                        if (coordClass[s] == r && coordStation[s] >= 0 && rows.get(c)[s] > 0
                                && rows.get(d)[s] < rows.get(c)[s] - ZERO) {
                            covers = false;
                            break;
                        }
                    }
                    if (covers) {
                        part = d;
                        break;
                    }
                }
                if (part < 0) {
                    implied = false;
                } else {
                    budget += bs.get(part);
                }
            }
            if (implied && budget <= bs.get(c) + ZERO) {
                keep[c] = false;
            }
        }

        List<Integer> kept = new ArrayList<Integer>();
        for (int c = 0; c < ncand; c++) {
            if (keep[c]) {
                kept.add(c);
            }
        }
        con.A = new Matrix(kept.size(), nstate);
        con.As = new Matrix(kept.size(), 0);
        con.b = new double[kept.size()];
        con.region = new int[kept.size()];
        con.station = new int[kept.size()];
        con.klassRow = new int[kept.size()];
        con.staged = new boolean[kept.size()];
        con.label = new String[kept.size()];
        for (int c = 0; c < kept.size(); c++) {
            int src = kept.get(c);
            for (int s = 0; s < nstate; s++) {
                con.A.set(c, s, rows.get(src)[s]);
            }
            con.b[c] = bs.get(src);
            con.region[c] = rgs.get(src);
            con.station[c] = sts.get(src);
            con.klassRow[c] = kls.get(src);
            con.staged[c] = sdg.get(src);
            con.label[c] = labels.get(src);
        }

        // WHICH STATION RULES ARE A CONSTRAINT ON THIS DRIFT, and which are a
        // different event set. THE RULE IS NOT WHAT DECIDES THE SEMANTICS -- the
        // class type is, exactly as State.arrivalIsLost decides it for every other
        // solver. So a waiting queue and a drop are BOTH constraints here (the
        // struct refresh declares DROP by default at a capped station reached by an
        // open class, so refusing it would refuse every open loss model); what is
        // refused is the rules that add STATE.
        //
        // A rule is only a contradiction where the cap can BIND, so this runs on the
        // surviving rows.
        if (sn.droprule != null && sn.stations != null && sn.jobclasses != null) {
            for (int c = 0; c < con.b.length; c++) {
                if (con.station[c] < 0 || con.station[c] >= sn.stations.size()) {
                    continue;
                }
                Map<JobClass, DropStrategy> perClass = sn.droprule.get(sn.stations.get(con.station[c]));
                if (perClass == null) {
                    continue;
                }
                for (int r = 0; r < K && r < sn.jobclasses.size(); r++) {
                    boolean weighs = false;
                    for (int s = 0; s < nstate; s++) {
                        if (coordStation[s] == con.station[c] && coordClass[s] == r
                                && con.A.get(c, s) > 0) {
                            weighs = true;
                            break;
                        }
                    }
                    if (!weighs) {
                        continue;
                    }
                    DropStrategy rule = perClass.get(sn.jobclasses.get(r));
                    if (rule != null && rule != DropStrategy.WaitingQueue && rule != DropStrategy.Drop) {
                        line_error(mfilename(new Object() {
                        }), String.format("Station %d applies %s to class %d, and its buffer binds. Only "
                                + "a waiting queue or a drop is a constraint on this drift: the first "
                                + "conserves the population and throttles the admission flow, the second "
                                + "discards the flow the cap will not take. BAS/BBS/RSRD add a "
                                + "blocked-server state to the upstream station and the retrial rules add "
                                + "an orbit, so each needs a different drift rather than an algebraic "
                                + "equation on this one. Use SolverCTMC, SolverJMT, SolverSSA or "
                                + "SolverLDES.", con.station[c] + 1, rule.toString(), r + 1));
                    }
                }
            }
        }
        return con;
    }

    private static void addRow(List<double[]> rows, List<Double> bs, List<Integer> rgs,
                               List<Integer> sts, List<Integer> kls, List<Boolean> sdg,
                               List<String> labels, double[] row, double b, int region,
                               int station, int klass, boolean staged, String label) {
        rows.add(row);
        bs.add(b);
        rgs.add(region);
        sts.add(station);
        kls.add(klass);
        sdg.add(staged);
        labels.add(label);
    }

    private static boolean sameRow(double[] a, double[] b, int nstate) {
        for (int s = 0; s < nstate; s++) {
            if (FastMath.abs(a[s] - b[s]) > ZERO) {
                return false;
            }
        }
        return true;
    }

    /**
     * The largest {@code row x} the population can produce, ignoring the coupling.
     * An upper bound is what is wanted: too loose only costs a constraint that stays
     * inactive, while too tight would discard a cap that does bind.
     */
    /**
     * The coordinates the immediate reduction folded away, where the reduced drift
     * holds no mass and no event lands.
     *
     * {@code terms.immediateAbsorb} is the projector the reduction returns: the
     * identity on a surviving coordinate and the absorption distribution on an
     * eliminated one, so a zero diagonal is exactly the eliminated case. It is null
     * when nothing was eliminated, where every coordinate survives.
     */
    private static boolean[] eliminatedCoords(FluidMomentTerms terms, int nstate) {
        boolean[] gone = new boolean[nstate];
        Matrix absorb = terms.immediateAbsorb;
        if (absorb == null || absorb.getNumRows() == 0) {
            return gone;
        }
        int n = Math.min(nstate, Math.min(absorb.getNumRows(), absorb.getNumCols()));
        for (int s = 0; s < n; s++) {
            gone[s] = absorb.get(s, s) == 0.0;
        }
        return gone;
    }

    private static double reach(double[] row, int[] coordClass, int[] coordStation,
                                Matrix njobs, int K, int nstate) {
        double total = 0;
        for (int r = 0; r < K; r++) {
            double w = 0;
            boolean any = false;
            for (int s = 0; s < nstate; s++) {
                if (coordClass[s] == r && coordStation[s] >= 0) {
                    w = FastMath.max(w, row[s]);
                    any = true;
                }
            }
            if (!any || w <= 0) {
                continue;
            }
            double nr = (njobs != null && r < njobs.length()) ? njobs.get(r) : Double.POSITIVE_INFINITY;
            if (!Double.isFinite(nr)) {
                return Double.POSITIVE_INFINITY;
            }
            total += w * nr;
        }
        return total;
    }

    private static boolean isFiniteCap(double v) {
        // FiniteCapacityRegion.UNBOUNDED is -1, and a negative cap is not a cap.
        return Double.isFinite(v) && v != -1.0 && v >= 0;
    }

    /**
     * Per cap and per event: is this event an admission the cap throttles, and where
     * does the stopped mass go. An admission is an event that pushes the constrained
     * quantity UP, read off {@code A D} rather than off the topology so that a
     * per-class or memory-weighted row picks out its own admissions with no extra
     * code.
     */
    static CapacityGates capacityGates(NetworkStruct sn, FluidMomentTerms terms,
                                       CapacityConstraints con) {
        int nevents = terms.eventIdx.length;
        int ncon = con.b.length;
        int nstate = terms.nstate;
        CapacityGates gates = new CapacityGates();
        gates.gate = new boolean[ncon][nevents];
        gates.held = new boolean[ncon][nevents];
        gates.loss = new boolean[ncon][nevents];
        if (ncon == 0) {
            return gates;
        }
        double tol = FastMath.sqrt(ZERO);
        gates.Dn = new Matrix(nstate, nevents);
        gates.DnExt = new Matrix(nstate, nevents);
        gates.Dp = new Matrix(nstate, nevents);
        boolean[] extCoord = new boolean[nstate];
        for (int i = 0; i < terms.M; i++) {
            if (terms.isExt[i]) {
                for (int s = 0; s < terms.stationBlock[i].length; s++) {
                    extCoord[terms.stationBlock[i][s]] = true;
                }
            }
        }
        for (int s = 0; s < nstate; s++) {
            for (int e = 0; e < nevents; e++) {
                double d = terms.D.get(s, e);
                if (d < 0) {
                    // A LOST ARRIVAL IS RETURNED TO THE SOURCE POOL: the EXT
                    // coordinate is a normalisation and not a population, so scaling
                    // only the arrival leg of a lost event would unbalance its row by
                    // exactly the loss.
                    if (extCoord[s]) {
                        gates.DnExt.set(s, e, d);
                    } else {
                        gates.Dn.set(s, e, d);
                    }
                }
                if (d > 0) {
                    gates.Dp.set(s, e, d);
                }
            }
        }
        for (int c = 0; c < ncon; c++) {
            boolean anyGate = false;
            for (int e = 0; e < nevents; e++) {
                double delta = 0;
                for (int s = 0; s < nstate; s++) {
                    delta += con.A.get(c, s) * terms.D.get(s, e);
                }
                if (delta <= tol) {
                    continue;
                }
                gates.gate[c][e] = true;
                anyGate = true;
                if (con.staged[c]) {
                    continue;
                }
                // the class is read off the coordinate the mass LANDS on, inside the
                // capped station: a class switch on entry would otherwise ask the
                // class the job is leaving behind whether it may be lost
                int k = -1;
                for (int s = 0; s < nstate; s++) {
                    if (gates.Dp.get(s, e) > tol && con.A.get(c, s) > 0) {
                        k = con.klass[s];
                        break;
                    }
                }
                boolean isopen = k >= 0 && k < sn.njobs.length() && !Double.isFinite(sn.njobs.get(k));
                gates.loss[c][e] = isopen;
                gates.held[c][e] = !isopen;
            }
            if (!anyGate) {
                final String lbl = con.label[c];
                line_error(mfilename(new Object() {
                }), String.format("No event increases %s, so the cap can never be approached and there "
                        + "is no admission flow for the constraint to throttle. This is a malformed "
                        + "limit rather than a solvable one.", lbl));
            }
        }
        return gates;
    }

    /**
     * The waiting room outside a capped region.
     *
     * <p>WHY THE CONSTRAINT ALONE IS NOT ENOUGH, FOR A REGION. Throttling a
     * region's admission events does hold its population at the cap, but it holds it
     * by slowing the UPSTREAM STATION'S COMPLETIONS -- an admission event IS that
     * station finishing a job -- so blocked mass piles up at a station it has
     * already finished being served by. A waiting queue means the job COMPLETES
     * upstream service and then waits; it is somewhere else, and the model needs
     * somewhere else to put it.</p>
     *
     * <p>A STATION BUFFER GETS NO ROOM, and that is not an omission: LINE disables
     * the upstream departure rather than moving the job out, so the blocked mass is
     * still at the upstream station and still counted there.</p>
     */
    static CapacityStaging capacityStaging(FluidMomentTerms terms, CapacityConstraints con) {
        int nevents = terms.eventIdx.length;
        int K = terms.K;
        int ncon = con.b.length;
        CapacityStaging stg = new CapacityStaging();
        stg.n = 0;
        stg.region = new int[0];
        stg.klass = new int[0];
        stg.adm = new boolean[nevents];
        stg.admRegion = new int[nevents];
        stg.admStage = new int[nevents];
        stg.gatedBy = new boolean[ncon][0];
        boolean anyStaged = false;
        for (int c = 0; c < ncon; c++) {
            anyStaged |= con.staged[c];
        }
        if (con.nregions == 0 || con.isEmpty() || !anyStaged) {
            return stg;
        }
        int nstate = terms.nstate;
        stg.Dn = new Matrix(nstate, nevents);
        stg.Dp = new Matrix(nstate, nevents);
        for (int s = 0; s < nstate; s++) {
            for (int e = 0; e < nevents; e++) {
                double d = terms.D.get(s, e);
                if (d < 0) {
                    stg.Dn.set(s, e, d);
                }
                if (d > 0) {
                    stg.Dp.set(s, e, d);
                }
            }
        }
        boolean[] stagedRegion = new boolean[con.nregions];
        for (int c = 0; c < ncon; c++) {
            if (con.staged[c] && con.region[c] >= 0) {
                stagedRegion[con.region[c]] = true;
            }
        }
        double tol = FastMath.sqrt(ZERO);
        int[][] idx = new int[con.nregions][K];
        List<Integer> region = new ArrayList<Integer>();
        List<Integer> klass = new ArrayList<Integer>();
        for (int f = 0; f < con.nregions; f++) {
            if (!stagedRegion[f]) {
                continue;
            }
            boolean anyMember = false;
            for (int s = 0; s < nstate; s++) {
                anyMember |= con.member[f][s];
            }
            if (!anyMember) {
                continue;
            }
            for (int e = 0; e < nevents; e++) {
                double delta = 0;
                for (int s = 0; s < nstate; s++) {
                    if (con.member[f][s]) {
                        delta += terms.D.get(s, e);
                    }
                }
                if (delta <= tol) {
                    continue;
                }
                int c = -1;
                for (int s = 0; s < nstate; s++) {
                    if (con.member[f][s] && stg.Dp.get(s, e) > tol) {
                        c = con.klass[s];
                        break;
                    }
                }
                if (c < 0 || c >= K) {
                    continue;
                }
                if (idx[f][c] == 0) {
                    region.add(f);
                    klass.add(c);
                    idx[f][c] = region.size();   // 1-based; 0 means none
                }
                stg.adm[e] = true;
                stg.admRegion[e] = f;
                stg.admStage[e] = idx[f][c] - 1;
            }
        }
        stg.n = region.size();
        stg.region = new int[stg.n];
        stg.klass = new int[stg.n];
        for (int j = 0; j < stg.n; j++) {
            stg.region[j] = region.get(j);
            stg.klass[j] = klass.get(j);
        }

        // WHICH ROWS GATE WHICH ROOM. A region-global cap gates every room of its
        // region, a per-class cap only the room of its class. A room gated by
        // several ACTIVE rows drains at the harmonic composition of their rates,
        // which is what lets two caps of one region bind at once.
        stg.gatedBy = new boolean[ncon][stg.n];
        for (int c = 0; c < ncon; c++) {
            if (!con.staged[c] || con.region[c] < 0) {
                continue;
            }
            for (int j = 0; j < stg.n; j++) {
                if (stg.region[j] != con.region[c]) {
                    continue;
                }
                double w = 0;
                for (int s = 0; s < nstate; s++) {
                    if (con.member[con.region[c]][s] && con.klass[s] == stg.klass[j]) {
                        w = FastMath.max(w, con.A.get(c, s));
                    }
                }
                stg.gatedBy[c][j] = w > 0;
            }
        }
        return stg;
    }

    /**
     * Extend every cap to the staging coordinates that hold mass INSIDE it.
     *
     * <p>A waiting room is outside the region it feeds, which is the whole point of
     * it -- but it is not outside every OTHER limit. Where two regions overlap, an
     * admission into the inner one is an INTERNAL move of the outer one: the job
     * leaves a station of the outer region, waits, and re-enters a station of the
     * same outer region, never having left it. Counting only the state coordinates
     * would take that mass out of the outer cap for as long as it waits, and the
     * Newton system that results is inconsistent rather than merely inexact.</p>
     */
    static void capacityExtend(CapacityConstraints con, CapacityStaging stg,
                               FluidMomentTerms terms) {
        int ncon = con.b.length;
        con.As = new Matrix(ncon, stg.n);
        if (ncon == 0 || stg.n == 0) {
            return;
        }
        int nstate = terms.nstate;
        int nevents = terms.eventIdx.length;
        double tol = FastMath.sqrt(ZERO);
        List<Set<Integer>> feeds = new ArrayList<Set<Integer>>();
        int[] dest = new int[stg.n];
        for (int j = 0; j < stg.n; j++) {
            feeds.add(new HashSet<Integer>());
            dest[j] = -1;
        }
        for (int e = 0; e < nevents; e++) {
            if (!stg.adm[e]) {
                continue;
            }
            int j = stg.admStage[e];
            for (int s = 0; s < nstate; s++) {
                if (terms.D.get(s, e) < -tol) {
                    feeds.get(j).add(s);
                }
                if (dest[j] < 0 && terms.D.get(s, e) > tol && con.member[stg.region[j]][s]) {
                    dest[j] = s;
                }
            }
        }
        for (int c = 0; c < ncon; c++) {
            for (int j = 0; j < stg.n; j++) {
                if (dest[j] < 0 || feeds.get(j).isEmpty()) {
                    continue;
                }
                double w = con.A.get(c, dest[j]);
                if (w <= 0) {
                    continue;
                }
                boolean allInside = true;
                for (Integer s : feeds.get(j)) {
                    if (con.A.get(c, s) <= 0) {
                        allInside = false;
                        break;
                    }
                }
                if (allInside) {
                    con.As.set(c, j, w);
                }
            }
        }
    }

    /**
     * The two legs of every event under the active caps, and the waiting rooms.
     *
     * <p>Shared by the steady-state residual and the transient right-hand side so
     * that the two solve the SAME model and not two spellings of it. What differs is
     * only what a STAGED cap's multiplier means: a drain RATE for the steady state,
     * where the room mass is pinned by the cap; the admitted FLOW for the transient,
     * because at the instant a region fills its room is EMPTY and a rate times zero
     * mass cannot hold the cap. The two rules agree at a fixed point, where the room
     * mass is proportional to its inflow.</p>
     *
     * <p>A held or lost cap composes as a PRODUCT of fractions either way, which is
     * what independent blocking gives and what keeps every active cap present in the
     * Jacobian. Several staged caps gating one room compose HARMONICALLY, because
     * the waits a job serves in turn add.</p>
     */
    static Legs legs(FluidMomentTerms terms, CapacityGates gates, CapacityStaging stg,
                     CapacityConstraints con, int[] active, double[] x, double[] sg,
                     Matrix r, double[] mult, boolean stagedFlow) {
        int nevents = terms.eventIdx.length;
        int ns = stg.n;
        int nact = active.length;
        double[] whole = new double[nevents];
        double[] entry = new double[nevents];
        Arrays.fill(whole, 1.0);
        Arrays.fill(entry, 1.0);
        double[] invTheta = new double[ns];
        double[] flowOf = new double[ns];
        boolean[] throttled = new boolean[ns];
        for (int k = 0; k < nact; k++) {
            int c = active[k];
            double m = mult[k];
            if (con.staged[c]) {
                for (int j = 0; j < ns; j++) {
                    if (stg.gatedBy[c][j]) {
                        throttled[j] = true;
                        if (!stagedFlow) {
                            invTheta[j] += m > ZERO ? 1.0 / m : Double.POSITIVE_INFINITY;
                        }
                    }
                }
            } else {
                for (int e = 0; e < nevents; e++) {
                    if (gates.held[c][e]) {
                        whole[e] *= m;
                    }
                    if (gates.loss[c][e]) {
                        entry[e] *= m;
                    }
                }
            }
        }
        boolean[] split = new boolean[nevents];
        for (int e = 0; e < nevents; e++) {
            if (ns > 0 && stg.adm[e] && throttled[stg.admStage[e]]) {
                split[e] = true;
            }
        }
        Legs out = new Legs();
        out.rup = new double[nevents];
        out.rin = new double[nevents];
        for (int e = 0; e < nevents; e++) {
            // A staged event is NOT suppressed upstream even when a held cap gates
            // it: the job completes upstream service into the waiting room, so the
            // held fraction moves to the room's exit leg below.
            out.rup[e] = r.get(e, 0) * (split[e] ? 1.0 : whole[e]);
            out.rin[e] = out.rup[e] * entry[e];
        }
        double[] inflow = new double[ns];
        double[] R = new double[ns];
        for (int e = 0; e < nevents; e++) {
            if (split[e]) {
                inflow[stg.admStage[e]] += out.rup[e];
                R[stg.admStage[e]] += out.rup[e];
            }
        }
        if (stagedFlow) {
            for (int k = 0; k < nact; k++) {
                int c = active[k];
                if (!con.staged[c]) {
                    continue;
                }
                double mass = 0;
                double tot = 0;
                int nrooms = 0;
                for (int j = 0; j < ns; j++) {
                    if (stg.gatedBy[c][j]) {
                        mass += sg[j];
                        tot += inflow[j];
                        nrooms++;
                    }
                }
                if (nrooms == 0) {
                    continue;
                }
                for (int j = 0; j < ns; j++) {
                    if (!stg.gatedBy[c][j]) {
                        continue;
                    }
                    double w;
                    if (mass > FINE) {
                        w = sg[j] / mass;
                    } else if (tot > FINE) {
                        w = inflow[j] / tot;
                    } else {
                        w = 1.0 / nrooms;
                    }
                    flowOf[j] += mult[k] * w;
                }
            }
        } else {
            for (int j = 0; j < ns; j++) {
                double theta = invTheta[j] > 0 ? 1.0 / invTheta[j] : 0.0;
                flowOf[j] = theta * sg[j];
            }
        }
        out.drain = flowOf.clone();
        double[] left = new double[ns];
        for (int e = 0; e < nevents; e++) {
            if (!split[e]) {
                continue;
            }
            int j = stg.admStage[e];
            double q = R[j] > FINE ? flowOf[j] * out.rup[e] / R[j] : 0.0;
            // the held fraction gates the room's EXIT: what it stops stays in the
            // room. The lost fraction gates the ARRIVAL: that mass leaves the room
            // and is destroyed.
            out.rin[e] = q * whole[e] * entry[e];
            left[j] += q * whole[e];
        }
        for (int j = 0; j < ns; j++) {
            if (R[j] > FINE) {
                out.drain[j] = left[j];
            }
        }
        out.ds = new double[ns];
        for (int j = 0; j < ns; j++) {
            // a waiting room with no cap above it must be EMPTY, not merely balanced
            out.ds[j] = throttled[j] ? inflow[j] - out.drain[j] : sg[j];
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // the residual system
    // -----------------------------------------------------------------------

    /** Everything the residual needs, gathered so the Newton can stay generic. */
    private static final class Ctx {
        FluidMomentTerms terms;
        Matrix C;
        double[] Nvec;
        int[] cidx;
        Matrix[] covblk;
        int nstate;
        CapacityConstraints con;
        CapacityStaging stg;
        CapacityGates gates;
        int[] active;
        boolean[] driftRows;
        Matrix clampT;      // the tangent space of the caps that clamp, or null
        int M;
    }

    /** The rate vector the metrics are read at: the flow that actually crosses. */
    private Matrix lastRates;

    /** The rate vector the events FIRE at, which is what the diffusion counts. */
    private Matrix lastFireRates;

    /**
     * Orthogonal projector onto the subspace the CLAMPING caps leave free.
     *
     * <p>A cap that holds the job upstream or loses it fixes its own combination of
     * the state for as long as it binds -- blocking answers the state instantaneously
     * -- so that combination does not fluctuate and the noise lives on the subspace
     * orthogonal to A. WITHOUT THIS THE LYAPUNOV SOLVE HAS NO SOLUTION AT ALL: with
     * the multiplier held constant the drift along the constrained direction is
     * neutral (an overloaded M/M/1/K sits on a whole line of equilibria), so the
     * Jacobian carries a zero eigenvalue there. A STAGED cap is not clamped and must
     * not be projected: its room is a state and supplies the restoring force.</p>
     */
    private static Matrix clampTangent(CapacityConstraints con, int[] active,
                                       FluidMomentTerms terms) {
        List<Integer> rows = new ArrayList<Integer>();
        for (int k = 0; k < active.length; k++) {
            if (!con.staged[active[k]]) {
                rows.add(active[k]);
            }
        }
        if (rows.isEmpty()) {
            return null;
        }
        int nc = terms.covIdx.length;
        Matrix R = new Matrix(rows.size(), nc);
        boolean any = false;
        for (int i = 0; i < rows.size(); i++) {
            for (int j = 0; j < nc; j++) {
                double v = con.A.get(rows.get(i), terms.covIdx[j]);
                R.set(i, j, v);
                any |= FastMath.abs(v) > ZERO;
            }
        }
        if (!any) {
            return null;
        }
        Matrix Rt = R.transpose();
        Matrix T = Matrix.eye(nc).sub(1, Rt.mult(R.mult(Rt).pinv()).mult(R));
        return T;
    }

    /**
     * The coupled algebraic system, stacked. Returns null when the closure cannot be
     * evaluated at this iterate, so the line search can back off; the caller
     * evaluates the seed with {@code quiet=false}, where a genuine failure must
     * surface.
     *
     * <p>THE UNKNOWNS ARE {@code [x; staging; sigma2; mult]}. Staging carries the
     * blocked mass a capped region will not admit yet; MULT is one multiplier per
     * ACTIVE cap, and which multiplier it is depends on where that cap sends the
     * mass it stops -- a FRACTION of the admissions allowed for a cap that holds or
     * loses, the RATE its waiting room drains at for one that stages.</p>
     */
    private double[] residual(double[] u, Ctx ctx, boolean quiet) {
        FluidMomentTerms terms = ctx.terms;
        int n = ctx.nstate;
        int ns = ctx.stg.n;
        int ncl = ctx.cidx.length;
        int nact = ctx.active.length;
        int nevents = terms.eventIdx.length;

        double[] x = Arrays.copyOfRange(u, 0, n);
        double[] sg = Arrays.copyOfRange(u, n, n + ns);
        double[] s2 = new double[ctx.M];
        // The Newton unknowns are unconstrained but a variance, a population and a
        // multiplier are not. The Newton projects the iterate into the feasible box,
        // so these clamps are guards rather than the mechanism.
        for (int j = 0; j < ncl; j++) {
            s2[ctx.cidx[j]] = FastMath.max(0, u[n + ns + j]);
        }
        double[] mult = new double[nact];
        for (int k = 0; k < nact; k++) {
            mult[k] = FastMath.max(0, u[n + ns + ncl + k]);
        }

        Matrix r;
        Legs lg;
        double[] s2new;
        try {
            r = terms.rates(x, s2, ctx.covblk);
            lg = legs(terms, ctx.gates, ctx.stg, ctx.con, ctx.active, x, sg, r, mult, false);
            Matrix A = terms.jacobian(x, s2, ctx.covblk);
            Matrix rfire = new Matrix(nevents, 1);
            for (int e = 0; e < nevents; e++) {
                rfire.set(e, 0, lg.rup[e]);
            }
            Matrix sigma = localLyapunov(A, rfire, terms, ctx.clampT);
            s2new = sigma2From(sigma, terms, ctx.cidx, ctx.M);
            this.lastFireRates = rfire;
        } catch (RuntimeException e) {
            if (quiet) {
                return null;
            }
            throw e;
        }
        Matrix rmet = new Matrix(nevents, 1);
        for (int e = 0; e < nevents; e++) {
            rmet.set(e, 0, lg.rin[e]);
        }
        this.lastRates = rmet;

        // THE THREE LEGS. The removal leaves at the rate the event FIRES and the
        // arrival lands at the rate mass actually ARRIVES; Dn + DnExt + Dp = D, so
        // this collapses to D rup wherever the two agree. A LOST arrival is returned
        // to the source pool rather than destroyed -- the EXT coordinate is a
        // normalisation, and its row is a real equation -- while a job lost leaving a
        // REAL station is destroyed, so that station's row keeps the firing rate.
        double[] drift = new double[n];
        for (int s = 0; s < n; s++) {
            double acc = 0;
            for (int e = 0; e < nevents; e++) {
                if (ctx.gates.Dn == null) {
                    acc += terms.D.get(s, e) * lg.rup[e];
                } else {
                    acc += ctx.gates.Dn.get(s, e) * lg.rup[e]
                            + ctx.gates.DnExt.get(s, e) * lg.rin[e]
                            + ctx.gates.Dp.get(s, e) * lg.rin[e];
                }
            }
            drift[s] = acc;
        }

        int ndrift = 0;
        for (int s = 0; s < n; s++) {
            if (ctx.driftRows[s]) {
                ndrift++;
            }
        }
        int nrows = ndrift + ns + ctx.C.getNumRows() + ncl + nact;
        double[] G = new double[nrows];
        int p = 0;
        for (int s = 0; s < n; s++) {
            if (ctx.driftRows[s]) {
                G[p++] = drift[s];
            }
        }
        for (int j = 0; j < ns; j++) {
            G[p++] = lg.ds[j];
        }
        for (int c = 0; c < ctx.C.getNumRows(); c++) {
            double acc = 0;
            for (int s = 0; s < n; s++) {
                acc += ctx.C.get(c, s) * x[s];
            }
            for (int j = 0; j < ns; j++) {
                acc += ctx.C.get(c, n + j) * sg[j];
            }
            G[p++] = acc - ctx.Nvec[c];
        }
        for (int j = 0; j < ncl; j++) {
            G[p++] = s2[ctx.cidx[j]] - s2new[ctx.cidx[j]];
        }
        for (int k = 0; k < nact; k++) {
            // THE UNDIFFERENTIATED CONSTRAINT, on purpose. At a fixed point every
            // flow already balances, so the differentiated form A D r = 0 is
            // satisfied by anything and pins no multiplier. A x = b does pin it,
            // through the dependence of the fixed point on the multiplier.
            double acc = 0;
            for (int s = 0; s < n; s++) {
                acc += ctx.con.A.get(ctx.active[k], s) * x[s];
            }
            for (int j = 0; j < ns; j++) {
                acc += ctx.con.As.get(ctx.active[k], j) * sg[j];
            }
            G[p++] = acc - ctx.con.b[ctx.active[k]];
        }
        return G;
    }

    /**
     * The multipliers that hold the active caps at this state, by small Newton.
     *
     * <p>Each active cap contributes one equation, {@code d/dt (A x + As s) = 0}, and
     * one unknown -- a fraction for a cap that holds or loses, an admitted flow for
     * one that stages -- so the system is square and small. This is what makes an
     * event RESTART consistent: RODAS needs the algebraic unknowns to satisfy their
     * equations at the initial point of an index-1 DAE. It is also the FEASIBILITY
     * test at an activation: a fraction above one means the cap would have to admit
     * more than arrives, so it is not binding after all.</p>
     */
    private double[] holdMultipliers(FluidMomentTerms terms, CapacityGates gates,
                                     CapacityStaging stg, CapacityConstraints con,
                                     int[] active, double[] x, double[] sg, double[] s2,
                                     Matrix[] covblk, double[] m0, boolean[] okOut) {
        double[] mult = m0.clone();
        if (active.length == 0) {
            okOut[0] = true;
            return mult;
        }
        double[] F = holdResidual(terms, gates, stg, con, active, x, sg, s2, covblk, mult);
        for (int it = 0; it < 40; it++) {
            if (infNorm(F) < FastMath.max(1e-12, 1e-10 * (infNorm(mult) + 1))) {
                break;
            }
            Matrix J = new Matrix(F.length, mult.length);
            for (int j = 0; j < mult.length; j++) {
                double h = FastMath.max(1e-7 * FastMath.abs(mult[j]), 1e-9);
                double[] mp = mult.clone();
                mp[j] += h;
                double[] Fp = holdResidual(terms, gates, stg, con, active, x, sg, s2, covblk, mp);
                for (int i = 0; i < F.length; i++) {
                    J.set(i, j, (Fp[i] - F[i]) / h);
                }
            }
            Matrix Fm = new Matrix(F.length, 1);
            for (int i = 0; i < F.length; i++) {
                Fm.set(i, 0, F[i]);
            }
            Matrix step;
            try {
                step = J.pinv().mult(Fm);
            } catch (RuntimeException e) {
                break;
            }
            double lam = 1.0;
            boolean stepped = false;
            for (int ls = 0; ls < 20; ls++) {
                double[] mn = new double[mult.length];
                for (int j = 0; j < mult.length; j++) {
                    mn[j] = FastMath.max(0, mult[j] - lam * step.get(j, 0));
                }
                double[] Fn = holdResidual(terms, gates, stg, con, active, x, sg, s2, covblk, mn);
                if (infNorm(Fn) < infNorm(F)) {
                    mult = mn;
                    F = Fn;
                    stepped = true;
                    break;
                }
                lam *= 0.5;
            }
            if (!stepped) {
                break;
            }
        }
        okOut[0] = infNorm(F) < 1e-6;
        return mult;
    }

    private double[] holdResidual(FluidMomentTerms terms, CapacityGates gates,
                                  CapacityStaging stg, CapacityConstraints con, int[] active,
                                  double[] x, double[] sg, double[] s2, Matrix[] covblk,
                                  double[] mult) {
        int n = terms.nstate;
        int nevents = terms.eventIdx.length;
        Matrix r = terms.rates(x, s2, covblk);
        Legs lg = legs(terms, gates, stg, con, active, x, sg, r, mult, true);
        double[] dx = new double[n];
        for (int s = 0; s < n; s++) {
            double acc = 0;
            for (int e = 0; e < nevents; e++) {
                if (gates.Dn == null) {
                    acc += terms.D.get(s, e) * lg.rup[e];
                } else {
                    acc += gates.Dn.get(s, e) * lg.rup[e] + gates.DnExt.get(s, e) * lg.rin[e]
                            + gates.Dp.get(s, e) * lg.rin[e];
                }
            }
            dx[s] = acc;
        }
        double[] F = new double[active.length];
        for (int k = 0; k < active.length; k++) {
            double acc = 0;
            for (int s = 0; s < n; s++) {
                acc += con.A.get(active[k], s) * dx[s];
            }
            for (int j = 0; j < stg.n; j++) {
                acc += con.As.get(active[k], j) * lg.ds[j];
            }
            F[k] = acc;
        }
        return F;
    }

    /**
     * Damped PROJECTED Newton with a finite-difference Jacobian and an Armijo
     * backtrack on the residual norm. The step is solved in LEAST SQUARES: the
     * drift block is rank deficient by exactly the number of conserved chains,
     * and the constraint rows restore that rank, so the stacked system is
     * consistent and overdetermined rather than square.
     *
     * <p>WHY PROJECTED, AND NOT MERELY CLAMPED. The unknowns past {@code nfree}
     * are a variance and an admission throttle, neither of which may go
     * negative, and the residual reads them through max(0,.). Clamping inside the
     * residual ALONE is a trap: once an iterate goes negative the residual stops
     * depending on it, so the finite-difference column is exactly zero, there is
     * no derivative to climb back on, and the unknown is pinned at the boundary
     * for good. Projecting the ITERATE keeps every evaluation inside the box,
     * where the forward difference across max(0,.) is live even at zero.</p>
     */
    private double[] newton(Ctx ctx, double[] u0, double tol, int maxit, int nfree,
                            int[] itOut, boolean[] convOut, double[] resOut) {
        double[] u = project(u0.clone(), nfree);
        double[] G = residual(u, ctx, false);
        double resnorm = infNorm(G);
        int it = 0;
        while (resnorm >= tol && it < maxit) {
            it++;
            int m = G.length;
            int nun = u.length;
            Matrix J = new Matrix(m, nun);
            for (int j = 0; j < nun; j++) {
                double h = FastMath.max(1e-7 * FastMath.abs(u[j]), 1e-9);
                double[] up = u.clone();
                up[j] += h;
                double[] Gp = residual(up, ctx, true);
                if (Gp == null || !allFinite(Gp)) {
                    double[] um = u.clone();
                    um[j] -= h;
                    double[] Gm = residual(um, ctx, true);
                    if (Gm == null || !allFinite(Gm)) {
                        continue;   // column left at zero; least squares absorbs it
                    }
                    for (int i = 0; i < m; i++) {
                        J.set(i, j, (G[i] - Gm[i]) / h);
                    }
                } else {
                    for (int i = 0; i < m; i++) {
                        J.set(i, j, (Gp[i] - G[i]) / h);
                    }
                }
            }
            // The minimum-norm least-squares step, which is what MATLAB's
            // backslash and numpy's lstsq both give on a full-column-rank
            // overdetermined system -- the case here, since the constraint rows
            // restore the rank the drift block is missing.
            Matrix Gm1 = new Matrix(m, 1);
            for (int i = 0; i < m; i++) {
                Gm1.set(i, 0, G[i]);
            }
            Matrix du;
            try {
                du = J.pinv().mult(Gm1);
            } catch (RuntimeException e) {
                break;
            }
            double lam = 1.0;
            boolean stepped = false;
            for (int ls = 0; ls < 25; ls++) {
                double[] un = new double[nun];
                for (int j = 0; j < nun; j++) {
                    un[j] = u[j] - lam * du.get(j, 0);
                }
                // project BEFORE evaluating, so the accepted point and the
                // residual that measured it are the same feasible point
                un = project(un, nfree);
                double[] Gn = residual(un, ctx, true);
                if (Gn != null && allFinite(Gn) && infNorm(Gn) < resnorm * (1 - 1e-4 * lam)) {
                    u = un;
                    G = Gn;
                    resnorm = infNorm(Gn);
                    stepped = true;
                    break;
                }
                lam *= 0.5;
            }
            if (!stepped) {
                break;   // no descent along this direction; report the iterate
            }
        }
        itOut[0] = it;
        convOut[0] = resnorm < tol;
        resOut[0] = resnorm;
        return u;
    }

    /** Onto the feasible box: the state is free, the variances are not. */
    private static double[] project(double[] u, int nfree) {
        for (int i = nfree; i < u.length; i++) {
            u[i] = FastMath.max(0, u[i]);
        }
        return u;
    }

    private static boolean allFinite(double[] v) {
        for (int i = 0; i < v.length; i++) {
            if (!Double.isFinite(v[i])) {
                return false;
            }
        }
        return true;
    }

    private static double infNorm(double[] v) {
        double m = 0;
        for (int i = 0; i < v.length; i++) {
            m = FastMath.max(m, FastMath.abs(v[i]));
        }
        return m;
    }

    /**
     * One Lyapunov solve on the coordinates that carry a real population.
     *
     * <p>Sigma is LINEAR in itself for a held x, so this is a SOLVE and not an
     * iteration -- which is why Sigma stays out of the Newton unknowns and only
     * the M numbers it projects onto go in. R IS THE FIRING RATE VECTOR when a cap
     * is active, so the diffusion matrix counts the events that actually happen.</p>
     *
     * <p>CLAMPT, when present, is the tangent space of the caps that CLAMP: a cap
     * that holds the job upstream or loses it fixes its own combination of the state
     * while it binds, so that combination does not fluctuate. Projecting the jump
     * directions is enough to state the reduced problem, because FluidLyapunov
     * restricts everything to range(D) already. See {@link #clampTangent}.</p>
     */
    private static Matrix localLyapunov(Matrix A, Matrix r, FluidMomentTerms terms,
                                        Matrix clampT) {
        int[] idx = terms.covIdx;
        int nc = idx.length;
        int nevents = terms.eventIdx.length;
        Matrix Dc = new Matrix(nc, nevents);
        for (int i = 0; i < nc; i++) {
            for (int e = 0; e < nevents; e++) {
                Dc.set(i, e, terms.D.get(idx[i], e));
            }
        }
        if (clampT != null) {
            Dc = clampT.mult(Dc);
        }
        Matrix Qc = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                double acc = 0;
                for (int e = 0; e < nevents; e++) {
                    acc += Dc.get(i, e) * r.get(e, 0) * Dc.get(j, e);
                }
                Qc.set(i, j, acc);
            }
        }
        Matrix Ac = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Ac.set(i, j, A.get(idx[i], idx[j]));
            }
        }
        Matrix Sc = FluidLyapunov.solve(Ac, Qc, Dc, 0).sigma;
        Matrix sigma = new Matrix(terms.nstate, terms.nstate);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                sigma.set(idx[i], idx[j], Sc.get(i, j));
            }
        }
        return sigma;
    }

    private static double[] sigma2From(Matrix sigma, FluidMomentTerms terms, int[] cidx, int M) {
        double[] s2 = new double[M];
        for (int j = 0; j < cidx.length; j++) {
            int[] blk = terms.stationBlock[cidx[j]];
            double acc = 0;
            for (int a = 0; a < blk.length; a++) {
                for (int b = 0; b < blk.length; b++) {
                    acc += sigma.get(blk[a], blk[b]);
                }
            }
            s2[cidx[j]] = FastMath.max(0, acc);
        }
        return s2;
    }

    private static double blockSum(Matrix sigma, int[] blk) {
        double acc = 0;
        for (int a = 0; a < blk.length; a++) {
            for (int b = 0; b < blk.length; b++) {
                acc += sigma.get(blk[a], blk[b]);
            }
        }
        return acc;
    }

    // -----------------------------------------------------------------------
    // the entry point
    // -----------------------------------------------------------------------

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        long t0 = System.nanoTime();
        int M = sn.nstations;
        int K = sn.nclasses;

        FluidMomentTerms terms = new FluidMomentTerms(sn, options);
        terms.checkSupported();

        // The simultaneous solve carries one Lyapunov solve per residual
        // evaluation and takes a finite-difference Jacobian over
        // nstate + nclosable unknowns, so its cost is cubic per evaluation and
        // quartic overall. Affordable at the scale the moment methods run at,
        // but the crossover is lower than the 200 minnormal permits, so it gets
        // its own limit.
        int maxstate = (options.config == null) ? 0 : options.config.dae_maxstate;
        if (maxstate <= 0) {
            maxstate = 100;
        }
        if (terms.nstate > maxstate) {
            line_error(mfilename(new Object() {
            }), String.format("The dae method solves a %d-unknown algebraic system with a "
                    + "finite-difference Jacobian, above the limit of %d set by "
                    + "options.config.dae_maxstate. Raise that limit, or use "
                    + "options.method=\"minnormal\" for the same closure by successive "
                    + "substitution.", terms.nstate, maxstate));
        }
        // DPS and GPS close on the covariance BETWEEN a station's class
        // coordinates, not on the station total, so their closure state is a
        // matrix block rather than the scalar this solves for. Carrying those
        // blocks as Newton unknowns puts the quartic term back; carrying them as
        // a chord reintroduces the alternation this method exists to remove.
        for (int i = 0; i < M; i++) {
            if (terms.sched[i] == SchedStrategy.DPS || terms.sched[i] == SchedStrategy.GPS) {
                line_error(mfilename(new Object() {
                }), "The dae method closes on the per-station variance only, but this model has a "
                        + "DPS or GPS station whose share closes on the covariance BETWEEN its class "
                        + "coordinates. Use options.method=\"minnormal\", which carries those blocks "
                        + "through its outer iteration.");
            }
        }

        // Caps first: staging depends on them, and conservation must count the
        // blocked mass staging holds. GATES decides, per cap and per event, where
        // the mass a cap stops actually goes -- held upstream, lost, or staged.
        CapacityConstraints con = capacityConstraints(sn, terms);
        int ncon = con.b.length;
        CapacityGates gates = capacityGates(sn, terms, con);
        CapacityStaging stg = capacityStaging(terms, con);
        capacityExtend(con, stg, terms);

        // -- population conservation, as equations ---------------------------
        Matrix C = conservationRows(sn, terms, stg);
        double[] Nvec = this.lastNvec;
        if (C.getNumRows() > 0) {
            double leak = 0;
            for (int c = 0; c < C.getNumRows(); c++) {
                for (int e = 0; e < terms.eventIdx.length; e++) {
                    double acc = 0;
                    for (int s = 0; s < terms.nstate; s++) {
                        acc += C.get(c, s) * terms.D.get(s, e);
                    }
                    leak = FastMath.max(leak, FastMath.abs(acc));
                }
            }
            if (leak > FastMath.sqrt(ZERO)) {
                line_error(mfilename(new Object() {
                }), String.format("The event set does not conserve a closed chain: the largest "
                        + "population leak per unit rate is %g, where it must be zero. The "
                        + "conservation constraint would contradict the drift rather than complete "
                        + "it.", leak));
            }
        }

        // -- which stations have a variance that enters the drift ------------
        // A delay or a source has no min() to close. A station that cannot fill
        // its servers has min(n,c) = n on its whole support, so closing it is
        // not an improvement but an error. Those are held at zero and are not
        // unknowns, which keeps the Newton system as small as the closure is.
        List<Integer> closable = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (terms.isExt[i] || terms.sched[i] == SchedStrategy.INF) {
                continue;
            }
            if (!Double.isFinite(terms.S.get(i, 0)) || terms.minExact[i]
                    || terms.stationBlock[i].length == 0) {
                continue;
            }
            closable.add(i);
        }
        int[] cidx = new int[closable.size()];
        for (int j = 0; j < cidx.length; j++) {
            cidx[j] = closable.get(j);
        }
        Matrix[] covblk = new Matrix[M];   // always null: DPS/GPS refused above

        // A SOURCE COORDINATE HAS NO FIXED POINT, so its drift row is not an
        // equation and must not be asked to vanish. In the closing representation
        // an open class's departures are routed back onto the EXT pseudo-station,
        // so that coordinate ACCUMULATES the mass that left the system: its drift
        // is the throughput, permanently, and 0 = D r is unsatisfiable there by
        // construction. Dropping those rows is exact rather than a concession,
        // because the EXT rate factor of a single-phase source is 1 identically
        // -- FluidMomentTerms refuses a multi-phase source outright -- so no
        // other row depends on the value of that coordinate.
        boolean[] driftRows = new boolean[terms.nstate];
        Arrays.fill(driftRows, true);
        for (int i = 0; i < M; i++) {
            if (terms.isExt[i]) {
                int[] blk = terms.stationBlock[i];
                for (int j = 0; j < blk.length; j++) {
                    driftRows[blk[j]] = false;
                }
            }
        }

        // -- seed ------------------------------------------------------------
        // Newton needs a point in the basin, not an answer. One first-order
        // solve supplies it, at the cost of the single integration this method
        // exists to avoid repeating twenty times.
        SolverOptions meanopt = options.copy();
        meanopt.method = "closing";
        meanopt.config.moment_sigma2 = new double[M];
        meanopt.config.moment_cov = new Matrix[M];
        solver_fluid(sn, meanopt, result);
        double[] x = new double[terms.nstate];
        for (int j = 0; j < terms.nstate; j++) {
            x[j] = this.xvec_it.get(0, j);
        }
        double[] x0 = new double[terms.nstate];
        for (int j = 0; j < terms.nstate; j++) {
            x0[j] = this.xvec_t.get(0, j);
        }
        // and reset the EXT coordinates to the unit normalisation mass the layout
        // intends: the seed integration ran to a horizon, so what it left there
        // is however much work happened to cross the system, not a state.
        for (int j = 0; j < terms.nstate; j++) {
            if (!driftRows[j]) {
                x[j] = x0[j];
            }
        }

        // THE VARIANCE IS SEEDED POSITIVE, which is why this route has no kink
        // probe. sigma2 = 0 is where min(n,c) has no derivative and a saturated
        // model's first-order fixed point sits exactly there. The station mean is
        // an O(N) starting value in the right units and costs no Lyapunov solve.
        double[] sigma2cur = new double[M];
        for (int j = 0; j < cidx.length; j++) {
            int[] blk = terms.stationBlock[cidx[j]];
            double acc = 0;
            for (int a = 0; a < blk.length; a++) {
                acc += x[blk[a]];
            }
            sigma2cur[cidx[j]] = FastMath.max(FINE, acc);
        }

        double tol = options.tol;
        if (!(tol > 0) || !Double.isFinite(tol)) {
            tol = FINE;
        }
        int newtonMax = FastMath.max(50, options.iter_max > 0 ? options.iter_max : 0);

        // -- the active-set loop ---------------------------------------------
        // A capacity constraint is an inequality, and an inequality has no
        // residual to hand a Newton solver -- only the caps that actually bind
        // become equations. Each pass is a complete simultaneous solve, so the
        // loop iterates over WHICH caps bind, not over the closure.
        //
        // ONE MULTIPLIER PER ACTIVE ROW, not per region. A region's waiting room
        // used to carry a single drain rate, so two caps of one region were two
        // equalities against one control and the case was refused by name. A room
        // gated by several active rows now drains at the harmonic composition of
        // their rates and a held cap composes as a product of fractions.
        //
        // THE SET IS WARM-STARTED FROM THE SEED, not empty. An empty start asks for
        // the UNCONSTRAINED fixed point first, and there are models where that point
        // does not exist: an overloaded M/M/1/K has no equilibrium without its cap.
        // THE SET STARTS EMPTY, and the first pass therefore asks for the
        // UNCONSTRAINED fixed point -- which is what makes the answer of a model
        // whose caps bind independent of how far outside the seed happened to land.
        // It is also not always solvable: an overloaded M/M/1/K has NO equilibrium
        // without its cap, so that pass fails rather than converging and the caps
        // the iterate violates are seeded into the set instead (below). Doing that
        // up front for every model looked cheaper and changed answers: a region of 8
        // over two identical queues settled 7.43/0.57 from the clipped seed and 4/4
        // from the unconstrained one, and only the second is the reference's.
        int[] active = new int[0];
        boolean seededFromFailure = false;
        int[] activeSolved = new int[0];
        int iters = 0;
        double[] mult = new double[0];
        double[] sg = new double[stg.n];
        int asetMax = FastMath.max(4, 2 * ncon + 2);
        boolean conv = false;
        boolean clamped = false;
        double resn = Double.POSITIVE_INFINITY;
        double[] u = null;
        Ctx ctx = null;
        // the last iterate that WAS a fixed point, and the answer read off it: a pass
        // that fails to converge leaves an iterate that is not a fixed point of
        // anything, and reading the next active set off it is how one bad pass turns
        // into a walk through unrelated capacity combinations
        double[] xOk = x.clone();
        double[] sigma2Ok = sigma2cur.clone();
        double[] bestX = null;
        double[] bestSg = null;
        double[] bestSigma2 = null;
        double[] bestMult = null;
        int[] bestActive = null;
        double bestRes = Double.POSITIVE_INFINITY;
        boolean bestClamped = false;
        boolean bestFeasible = false;
        int aset = 0;
        for (aset = 0; aset < asetMax; aset++) {
            boolean hasClamp = false;
            for (int a = 0; a < active.length; a++) {
                hasClamp |= !con.staged[active[a]];
            }
            // Seed each waiting room with the mass that does not fit, and every
            // multiplier at unity -- an unthrottled fraction for a held or lost cap,
            // and the drain rate the region route has always started from.
            double[] sg0 = new double[stg.n];
            for (int k = 0; k < active.length; k++) {
                int c = active[k];
                if (!con.staged[c]) {
                    continue;
                }
                double cur = 0;
                for (int s = 0; s < terms.nstate; s++) {
                    cur += con.A.get(c, s) * x[s];
                }
                double excess = FastMath.max(0, cur - con.b[c]);
                int cnt = 0;
                for (int j = 0; j < stg.n; j++) {
                    if (stg.gatedBy[c][j]) {
                        cnt++;
                    }
                }
                if (cnt > 0) {
                    for (int j = 0; j < stg.n; j++) {
                        if (stg.gatedBy[c][j]) {
                            sg0[j] = excess / cnt;
                        }
                    }
                }
            }
            double[] u0 = new double[terms.nstate + stg.n + cidx.length + active.length];
            System.arraycopy(x, 0, u0, 0, terms.nstate);
            System.arraycopy(sg0, 0, u0, terms.nstate, stg.n);
            for (int j = 0; j < cidx.length; j++) {
                u0[terms.nstate + stg.n + j] = sigma2cur[cidx[j]];
            }
            for (int k = 0; k < active.length; k++) {
                u0[terms.nstate + stg.n + cidx.length + k] = 1.0;
            }

            // THE CLAMPED COVARIANCE IS A FALLBACK, NOT THE DEFAULT. A cap that holds
            // or loses fixes its own combination of the state, so the honest LNA puts
            // no fluctuation there -- but the truth is neither that nor the
            // unprojected variance: the population under a cap follows a TRUNCATED
            // distribution. The unprojected solve is what every other fluid method
            // computes, so it runs first; the projection is tried only when the
            // unprojected system has no stationary covariance at all, which is the
            // neutral case an overloaded loss station produces. Deciding once per
            // PASS matters: a projection switching between iterates would give Newton
            // a discontinuous system.
            boolean[] attempts = hasClamp ? new boolean[]{false, true} : new boolean[]{false};
            u = u0;
            RuntimeException nohyp = null;
            for (int ia = 0; ia < attempts.length; ia++) {
                ctx = new Ctx();
                ctx.terms = terms;
                ctx.C = C;
                ctx.Nvec = Nvec;
                ctx.cidx = cidx;
                ctx.covblk = covblk;
                ctx.nstate = terms.nstate;
                ctx.con = con;
                ctx.stg = stg;
                ctx.gates = gates;
                ctx.active = active;
                ctx.driftRows = driftRows;
                ctx.M = M;
                ctx.clampT = attempts[ia] ? clampTangent(con, active, terms) : null;
                int[] itOut = new int[1];
                boolean[] convOut = new boolean[1];
                double[] resOut = new double[1];
                try {
                    u = newton(ctx, u0, tol, newtonMax, terms.nstate, itOut, convOut, resOut);
                } catch (RuntimeException e) {
                    if (attempts[ia] || !hasClamp) {
                        nohyp = e;
                        break;
                    }
                    continue;
                }
                iters += itOut[0];
                conv = convOut[0];
                resn = resOut[0];
                clamped = attempts[ia];
                if (conv) {
                    break;
                }
            }
            if (nohyp != null) {
                // THE UNCONSTRAINED FIXED POINT NEED NOT EXIST. An overloaded open
                // station has no equilibrium until its buffer bounds it, so the pass
                // that asks for one fails and the caps the iterate violates -- or
                // left non-finite, which no comparison catches -- are seeded into the
                // active set instead. Once: a second failure with caps already bound
                // is the model's answer, not a starting point to improve.
                List<Integer> cand = new ArrayList<Integer>();
                for (int c = 0; c < ncon; c++) {
                    boolean inActive = false;
                    for (int a = 0; a < active.length; a++) {
                        inActive |= active[a] == c;
                    }
                    if (inActive) {
                        continue;
                    }
                    double acc = 0;
                    boolean bad = false;
                    for (int s = 0; s < terms.nstate; s++) {
                        if (con.A.get(c, s) > 0 && !Double.isFinite(x[s])) {
                            bad = true;
                        }
                        acc += con.A.get(c, s) * x[s];
                    }
                    if (bad || acc > con.b[c] + FastMath.max(1e-9, tol)) {
                        cand.add(c);
                    }
                }
                if (seededFromFailure || cand.isEmpty()) {
                    throw nohyp;
                }
                for (int ci = 0; ci < cand.size(); ci++) {
                    int c = cand.get(ci);
                    double val = 0;
                    int cnt = 0;
                    for (int s = 0; s < terms.nstate; s++) {
                        val += con.A.get(c, s) * x[s];
                        if (con.A.get(c, s) > 0) {
                            cnt++;
                        }
                    }
                    if (cnt == 0 || con.b[c] <= 0) {
                        continue;
                    }
                    if (Double.isFinite(val) && val > con.b[c]) {
                        double f = con.b[c] / val;
                        for (int s = 0; s < terms.nstate; s++) {
                            if (con.A.get(c, s) > 0) {
                                x[s] *= f;
                            }
                        }
                    } else if (!Double.isFinite(val)) {
                        for (int s = 0; s < terms.nstate; s++) {
                            if (con.A.get(c, s) > 0) {
                                x[s] = con.b[c] / cnt;
                            }
                        }
                    }
                }
                for (int s = 0; s < terms.nstate; s++) {
                    if (!Double.isFinite(x[s])) {
                        x[s] = 0;
                    }
                }
                List<Integer> merged = new ArrayList<Integer>();
                for (int a = 0; a < active.length; a++) {
                    merged.add(active[a]);
                }
                merged.addAll(cand);
                java.util.Collections.sort(merged);
                active = new int[merged.size()];
                for (int a = 0; a < merged.size(); a++) {
                    active[a] = merged.get(a);
                }
                seededFromFailure = true;
                continue;
            }

            System.arraycopy(u, 0, x, 0, terms.nstate);
            sg = Arrays.copyOfRange(u, terms.nstate, terms.nstate + stg.n);
            sigma2cur = new double[M];
            for (int j = 0; j < cidx.length; j++) {
                sigma2cur[cidx[j]] = FastMath.max(0, u[terms.nstate + stg.n + j]);
            }
            mult = new double[active.length];
            for (int k = 0; k < active.length; k++) {
                mult[k] = FastMath.max(0, u[terms.nstate + stg.n + cidx.length + k]);
            }
            // the set that produced THIS x and these multipliers, which is what the
            // metrics are read at; `active` below is the set to try NEXT
            activeSolved = active;
            if (ncon == 0) {
                break;
            }

            double[] slack = new double[ncon];
            boolean feasible = true;
            for (int c = 0; c < ncon; c++) {
                double acc = 0;
                for (int s = 0; s < terms.nstate; s++) {
                    acc += con.A.get(c, s) * x[s];
                }
                for (int j = 0; j < stg.n; j++) {
                    acc += con.As.get(c, j) * sg[j];
                }
                slack[c] = con.b[c] - acc;
                feasible &= slack[c] >= -FastMath.max(1e-9, tol);
            }
            if (conv) {
                xOk = x.clone();
                sigma2Ok = sigma2cur.clone();
                bestX = x.clone();
                bestSg = sg.clone();
                bestSigma2 = sigma2cur.clone();
                bestMult = mult.clone();
                bestActive = active.clone();
                bestRes = resn;
                bestClamped = clamped;
                bestFeasible = feasible;
            }
            List<Integer> violated = new ArrayList<Integer>();
            for (int c = 0; c < ncon; c++) {
                boolean inActive = false;
                for (int a = 0; a < active.length; a++) {
                    inActive |= active[a] == c;
                }
                if (!inActive && slack[c] < -FastMath.max(1e-9, tol)) {
                    violated.add(c);
                }
            }
            // THE RELEASE SIGNAL IS THE MULTIPLIER'S OWN UNITS, and the two kinds do
            // not share them. A held or lost cap throttles by a FRACTION, so one that
            // came back above one was holding the flow down for no reason. A staged
            // cap throttles by a RATE, which has no such scale -- there the signal is
            // a waiting room with no blocked mass at all.
            List<Integer> released = new ArrayList<Integer>();
            for (int k = 0; k < active.length; k++) {
                int c = active[k];
                if (con.staged[c]) {
                    double held = 0;
                    int rooms = 0;
                    for (int j = 0; j < stg.n; j++) {
                        if (stg.gatedBy[c][j]) {
                            held += sg[j];
                            rooms++;
                        }
                    }
                    if (rooms == 0 || held < 1e-9) {
                        released.add(c);
                    }
                } else if (mult[k] > 1.0 + FastMath.max(1e-9, tol)) {
                    released.add(c);
                }
            }
            if (!conv) {
                // A FAILED PASS SAYS NOTHING ABOUT WHICH CAPS BIND. Its release
                // signal is still information, but its state is not, so no row is
                // ADDED from it and the next pass restarts from the last point that
                // was a fixed point.
                violated.clear();
                x = xOk.clone();
                sigma2cur = sigma2Ok.clone();
                if (released.isEmpty()) {
                    break;
                }
            }
            if (violated.isEmpty() && released.isEmpty()) {
                break;
            }
            List<Integer> next = new ArrayList<Integer>();
            for (int a = 0; a < active.length; a++) {
                if (!released.contains(active[a])) {
                    next.add(active[a]);
                }
            }
            for (int v = 0; v < violated.size(); v++) {
                if (!next.contains(violated.get(v))) {
                    next.add(violated.get(v));
                }
            }
            java.util.Collections.sort(next);
            active = new int[next.size()];
            for (int a = 0; a < next.size(); a++) {
                active[a] = next.get(a);
            }
        }
        if (ncon > 0 && aset == asetMax) {
            line_warning(mfilename(new Object() {
            }), "The active set did not settle: the same capacity constraints kept binding and "
                    + "releasing. The reported point satisfies the last set tried.\n");
        }
        // A CONVERGED POINT BEATS THE LAST ITERATE. The loop can end on a pass that
        // did not converge -- a cap the closure cannot hold at any multiplier is
        // added, fails and is released for as long as the loop runs -- and the last
        // iterate of such a pass is not a fixed point of anything.
        if (bestX != null && !conv) {
            x = bestX.clone();
            sg = bestSg.clone();
            sigma2cur = bestSigma2.clone();
            mult = bestMult.clone();
            activeSolved = bestActive.clone();
            resn = bestRes;
            clamped = bestClamped;
            conv = true;
            if (!bestFeasible) {
                StringBuilder names = new StringBuilder();
                for (int c = 0; c < ncon; c++) {
                    double acc = 0;
                    for (int s = 0; s < terms.nstate; s++) {
                        acc += con.A.get(c, s) * x[s];
                    }
                    for (int j = 0; j < stg.n; j++) {
                        acc += con.As.get(c, j) * sg[j];
                    }
                    if (acc > con.b[c] + FastMath.max(1e-9, tol)) {
                        if (names.length() > 0) {
                            names.append(", ");
                        }
                        names.append(con.label[c]);
                    }
                }
                final String over = names.toString();
                line_warning(mfilename(new Object() {
                }), String.format("No fixed point of the closure satisfies %s. The closure wants more "
                        + "jobs there than the cap allows and no admission multiplier holds it: the "
                        + "reported point is the converged fixed point of the model without that cap, "
                        + "and it exceeds the cap. Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES "
                        + "for this model.\n", over));
            }
        }
        if (!conv) {
            // Report rather than return a point that is not a fixed point. The
            // substitution route would simply have stopped at outer_max with no
            // indication, which is the failure mode this replaces.
            line_warning(mfilename(new Object() {
            }), String.format("The simultaneous closure solve stopped at residual %.3e after %d "
                    + "Newton steps without reaching %.3e. The reported point is the last "
                    + "iterate.\n", resn, iters, tol));
        }

        // the variance as it entered the DRIFT
        double[] sigma2drift = sigma2cur.clone();
        Matrix[] covdrift = covblk.clone();

        // RMET is the flow that actually CROSSES each event, which is what the
        // throughput table must report, and the FIRING rate is what the diffusion
        // counts: a lost arrival fires and lands nowhere, so the two differ by
        // exactly the loss.
        ctx = new Ctx();
        ctx.terms = terms;
        ctx.C = C;
        ctx.Nvec = Nvec;
        ctx.cidx = cidx;
        ctx.covblk = covblk;
        ctx.nstate = terms.nstate;
        ctx.con = con;
        ctx.stg = stg;
        ctx.gates = gates;
        ctx.active = activeSolved;
        ctx.driftRows = driftRows;
        ctx.M = M;
        ctx.clampT = clamped ? clampTangent(con, activeSolved, terms) : null;
        double[] uFinal = new double[terms.nstate + stg.n + cidx.length + activeSolved.length];
        System.arraycopy(x, 0, uFinal, 0, terms.nstate);
        System.arraycopy(sg, 0, uFinal, terms.nstate, stg.n);
        for (int j = 0; j < cidx.length; j++) {
            uFinal[terms.nstate + stg.n + j] = sigma2cur[cidx[j]];
        }
        for (int k = 0; k < activeSolved.length; k++) {
            uFinal[terms.nstate + stg.n + cidx.length + k] = mult[k];
        }
        residual(uFinal, ctx, false);
        Matrix r = this.lastRates;
        Matrix sigma = localLyapunov(terms.jacobian(x, sigma2cur, covblk), this.lastFireRates,
                terms, ctx.clampT);

        // -- transient -------------------------------------------------------
        this.sigmat = null;
        this.qVart = null;
        this.tvar = null;
        boolean wantTran = options.timespan != null && options.timespan.length > 1
                && Double.isFinite(options.timespan[1]);
        if (wantTran) {
            // UNDER A CAP THIS IS A HYBRID DAE, integrated segment by segment with
            // the binding set updated at each located crossing: what the steady
            // state settles once with an active-set loop, the trajectory settles
            // again at every fill and every drain.
            transient_(terms, C, Nvec, cidx, covblk, sigma2drift, options, tol, x0, M, K,
                    con, gates, stg);
        }

        // -- performance measures --------------------------------------------
        // Read back from the same event representation that defines the drift,
        // so throughput balances flow at the fixed point. The THROTTLED rates, so
        // the throughput reported at a blocked region is the flow that actually
        // crosses it rather than the nominal one.
        Matrix gfac = terms.factors(x, sigma2drift, covdrift);
        Matrix Seff = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double v = terms.S.get(i, 0);
            if (terms.lldscaling != null && terms.lldscaling.getNumRows() > i) {
                for (int j = 0; j < terms.lldscaling.getNumCols(); j++) {
                    v = FastMath.max(v, terms.lldscaling.get(i, j));
                }
            }
            Seff.set(i, 0, v);
        }
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length == 0) {
                    continue;
                }
                double q = 0;
                double g = 0;
                for (int j = 0; j < blk.length; j++) {
                    q += x[blk[j]];
                    g += gfac.get(blk[j], 0);
                }
                // A Source holds no jobs: its coordinates are the arrival
                // process's phase indicator, not a population.
                result.QN.set(i, k, terms.isExt[i] ? 0.0 : q);
                result.UN.set(i, k, terms.sched[i] == SchedStrategy.INF ? q : g / Seff.get(i, 0));
                double tn = 0;
                for (int e = 0; e < terms.eventIdx.length; e++) {
                    if (terms.evIsDeparture[e] && terms.evStation[e] == i && terms.evClass[e] == k) {
                        tn += r.get(e, 0);
                    }
                }
                result.TN.set(i, k, tn);
                // See SolverFluid: TN is zero only to the integrator's accuracy.
                if (tn > GlobalConstants.Zero) {
                    result.RN.set(i, k, result.QN.get(i, k) / tn);
                }
            }
        }

        // UTILIZATION MUST BE READ FROM THE FLOW THAT ACTUALLY CROSSES once a cap
        // binds. Away from a constraint the in-service fluid sum(gfac)/s and the
        // carried utilization T/(mu s) are the same number, because the drift
        // balances; under an ACTIVE cap they are not -- the multiplier throttles
        // the departures (TN) and leaves the in-service fluid (gfac) alone, so the
        // two columns disagreed: a closed tandem capped at 1 reported Util 0.688 at
        // a station whose own Tput/mu was 0.550, and the exact answer is neither.
        // Every other solver reports the CARRIED utilization (Util = X*D, the
        // utilization law), and SolverCTMC gives 0.444 for the same station, so that
        // is the convention the throttled point has to keep. Unconstrained runs are
        // bit-identical: the branch is entered only when the active set is non-empty.
        // see _kb/06-solver-catalog.md
        if (activeSolved.length > 0) {
            for (int i = 0; i < M; i++) {
                if (terms.sched[i] == SchedStrategy.INF || terms.isExt[i]) {
                    continue;
                }
                for (int k = 0; k < K; k++) {
                    if (terms.classBlock[i * K + k].length == 0) {
                        continue;
                    }
                    double mu = sn.rates.get(i, k);
                    if (Double.isFinite(mu) && mu > 0) {
                        result.UN.set(i, k, result.TN.get(i, k) / (mu * Seff.get(i, 0)));
                    }
                }
            }
        }

        // -- transient measures ----------------------------------------------
        // On the DAE route sigma2 VARIES along the trajectory, so the rate
        // factors are read at the variance the trajectory HAD at each instant.
        // minnormal reads its whole transient at the single stationary variance,
        // which is the quasi-steady-state approximation of this.
        int nt = this.xvec_t.getNumRows();
        result.QNt = new Matrix[M][K];
        result.UNt = new Matrix[M][K];
        result.TNt = new Matrix[M][K];
        Matrix Gt = new Matrix(nt, terms.nstate);
        Matrix Rt = new Matrix(nt, terms.eventIdx.length);
        double[] xs = new double[terms.nstate];
        for (int s = 0; s < nt; s++) {
            for (int j = 0; j < terms.nstate; j++) {
                xs[j] = this.xvec_t.get(s, j);
            }
            double[] s2s = sigma2drift;
            if (this.sigmat != null && s < this.sigmat.length) {
                s2s = sigma2From(this.sigmat[s], terms, cidx, M);
            }
            Matrix gs = terms.factors(xs, s2s, covdrift);
            Matrix rs = terms.rates(xs, s2s, covdrift);
            for (int j = 0; j < terms.nstate; j++) {
                Gt.set(s, j, gs.get(j, 0));
            }
            for (int e = 0; e < terms.eventIdx.length; e++) {
                Rt.set(s, e, rs.get(e, 0));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                result.QNt[i][k] = new Matrix(nt, 1);
                result.UNt[i][k] = new Matrix(nt, 1);
                result.TNt[i][k] = new Matrix(nt, 1);
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length == 0) {
                    continue;
                }
                for (int s = 0; s < nt; s++) {
                    double q = 0;
                    double g = 0;
                    for (int j = 0; j < blk.length; j++) {
                        q += this.xvec_t.get(s, blk[j]);
                        g += Gt.get(s, blk[j]);
                    }
                    result.QNt[i][k].set(s, 0, terms.isExt[i] ? 0.0 : q);
                    result.UNt[i][k].set(s, 0,
                            terms.sched[i] == SchedStrategy.INF ? q : g / Seff.get(i, 0));
                    double tn = 0;
                    for (int e = 0; e < terms.eventIdx.length; e++) {
                        if (terms.evIsDeparture[e] && terms.evStation[e] == i
                                && terms.evClass[e] == k) {
                            tn += Rt.get(s, e);
                        }
                    }
                    result.TNt[i][k].set(s, 0, tn);
                }
            }
        }
        if (this.tvar != null) {
            result.t = this.tvar;
        }

        this.sigmaMatrix = sigma;
        this.sigma2 = sigma2cur;
        this.sigma2Drift = sigma2drift;
        this.classBlock = terms.classBlock;
        this.stationBlock = terms.stationBlock;
        this.outerIters = iters;
        this.residual = resn;
        this.converged = conv;
        this.qVar = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length > 0) {
                    this.qVar.set(i, k, FastMath.max(0, blockSum(sigma, blk)));
                }
            }
        }
        double cerr = 0;
        for (int c = 0; c < C.getNumRows(); c++) {
            double acc = 0;
            for (int s = 0; s < terms.nstate; s++) {
                acc += C.get(c, s) * x[s];
            }
            for (int j = 0; j < stg.n; j++) {
                acc += C.get(c, terms.nstate + j) * sg[j];
            }
            cerr = FastMath.max(cerr, FastMath.abs(acc - Nvec[c]));
        }
        this.conservation = cerr;

        // What the capacity limits did: which bound, how much mass each holds
        // outside its region, and the multiplier each settled at -- a fraction of
        // the admissions allowed for a cap that holds the job upstream or loses it,
        // a drain rate for one that stages it. STAGING IS REPORTED HERE AND NOT
        // FOLDED INTO QN, because a STAGED job is at no station -- the same choice
        // LDES makes. A HELD job is not staged and IS folded in, at the upstream
        // station where the reference counts it.
        this.capacityLabel = con.label;
        this.capacityB = new Matrix(con.b.length, 1);
        this.capacityValue = new Matrix(con.b.length, 1);
        for (int c = 0; c < con.b.length; c++) {
            this.capacityB.set(c, 0, con.b[c]);
            double acc = 0;
            for (int s = 0; s < terms.nstate; s++) {
                acc += con.A.get(c, s) * x[s];
            }
            for (int j = 0; j < stg.n; j++) {
                acc += con.As.get(c, j) * sg[j];
            }
            this.capacityValue.set(c, 0, acc);
        }
        this.capacityActive = activeSolved;
        this.capacityStaged = con.staged;
        this.capacityRegion = con.region;
        this.capacityStation = con.station;
        this.staging = new Matrix(stg.n, 1);
        double blk2 = 0;
        for (int j = 0; j < stg.n; j++) {
            this.staging.set(j, 0, sg[j]);
            blk2 += sg[j];
        }
        this.blocked = blk2;
        this.stagingRegion = stg.region;
        this.stagingClass = stg.klass;
        this.drain = new Matrix(mult.length, 1);
        for (int k = 0; k < mult.length; k++) {
            this.drain.set(k, 0, mult[k]);
        }

        result.method = "dae";
        result.runtime = (System.nanoTime() - t0) / 1.0e9;
    }

    /** Carries Nvec out of {@link #conservationRows}, which returns the matrix. */
    private double[] lastNvec = new double[0];

    /**
     * Population conservation, one row per CLOSED chain.
     *
     * <p>An open chain has no conserved population and contributes nothing. The
     * EXT coordinates are excluded because the closing representation holds unit
     * mass there as a normalisation constant, not as a job count.</p>
     */
    private Matrix conservationRows(NetworkStruct sn, FluidMomentTerms terms, CapacityStaging stg) {
        int M = terms.M;
        int K = terms.K;
        int nstate = terms.nstate;
        int[] coordClass = new int[nstate];
        int[] coordStation = new int[nstate];
        Arrays.fill(coordStation, -1);
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                int kic = (int) terms.Kic.get(i, c);
                if (kic > 0) {
                    int lo = (int) terms.qIndices.get(i, c);
                    for (int s = lo; s < lo + kic; s++) {
                        coordClass[s] = c;
                        coordStation[s] = i;
                    }
                }
            }
        }
        List<double[]> rows = new ArrayList<double[]>();
        List<Double> nlist = new ArrayList<Double>();
        if (sn.chains != null) {
            for (int ch = 0; ch < sn.chains.getNumRows(); ch++) {
                List<Integer> inch = new ArrayList<Integer>();
                double nch = 0;
                boolean finite = true;
                for (int k = 0; k < K; k++) {
                    if (sn.chains.get(ch, k) == 0) {
                        continue;
                    }
                    inch.add(k);
                    double nj = sn.njobs.get(k);
                    if (!Double.isFinite(nj)) {
                        finite = false;
                    }
                    nch += nj;
                }
                if (inch.isEmpty() || !finite || nch <= 0) {
                    continue;   // open chain: nothing is conserved
                }
                double[] row = new double[nstate + stg.n];
                boolean any = false;
                for (int s = 0; s < nstate; s++) {
                    if (coordStation[s] < 0 || terms.isExt[coordStation[s]]) {
                        continue;
                    }
                    if (!inch.contains(coordClass[s])) {
                        continue;
                    }
                    row[s] = 1.0;
                    any = true;
                }
                // A JOB IN THE WAITING QUEUE IS STILL IN THE CHAIN. Leaving the
                // staging coordinates out of this row would let the constraint
                // balance while the blocked mass quietly left the model.
                for (int j = 0; j < stg.n; j++) {
                    if (inch.contains(stg.klass[j])) {
                        row[nstate + j] = 1.0;
                    }
                }
                if (!any) {
                    continue;
                }
                rows.add(row);
                nlist.add(nch);
            }
        }
        Matrix C = new Matrix(rows.size(), nstate + stg.n);
        this.lastNvec = new double[rows.size()];
        for (int c = 0; c < rows.size(); c++) {
            for (int s = 0; s < nstate + stg.n; s++) {
                C.set(c, s, rows.get(c)[s]);
            }
            this.lastNvec[c] = nlist.get(c);
        }
        return C;
    }

    /**
     * The event functions, all of them, as one vector: a cap that is NOT active is
     * watched for REACHING its bound, one that IS active for stopping to bind -- a
     * fraction above one, or a room that has emptied.
     *
     * <p>ARMED IS A LATCH, and it has to be: a staged cap activates with an EMPTY
     * room, so "the room emptied" is true at the instant it starts binding and the
     * release would fire immediately. The latch only ever goes false to true.</p>
     */
    private static double[] capEvents(CapacityConstraints con, CapacityStaging stg, int[] active,
                                      Set<Integer> armed, double[] z, int nstate, int oSg, int oM,
                                      int ncon, int nstg) {
        double[] g = new double[ncon];
        for (int c = 0; c < ncon; c++) {
            boolean isActive = false;
            for (int a = 0; a < active.length; a++) {
                isActive |= active[a] == c;
            }
            if (isActive) {
                if (con.staged[c]) {
                    double mass = 0;
                    for (int j = 0; j < nstg; j++) {
                        if (stg.gatedBy[c][j]) {
                            mass += z[oSg + j];
                        }
                    }
                    g[c] = armed.contains(c) ? mass : 1.0;
                } else {
                    g[c] = 1.0 - z[oM + c];
                }
            } else {
                double val = 0;
                for (int s = 0; s < nstate; s++) {
                    val += con.A.get(c, s) * z[s];
                }
                for (int j = 0; j < nstg; j++) {
                    val += con.As.get(c, j) * z[oSg + j];
                }
                g[c] = con.b[c] - val;
            }
        }
        return g;
    }

    /**
     * The transient closure as an index-1 DAE, integrated with a SINGULAR mass
     * matrix by the vendored RODAS.
     *
     * <p>UNDER A CAP THIS IS A HYBRID DAE, and it is integrated as one rather than
     * refused. The system SWITCHES every time a cap starts or stops binding, so the
     * horizon is covered by SEGMENTS: each segment is one index-1 DAE with a fixed
     * set of binding caps, the segment ends at a located crossing, and the next one
     * starts from that state with the set updated. Integrating with the set frozen
     * instead would silently report the unconstrained trajectory through a cap the
     * model declares.</p>
     *
     * <p>THE MULTIPLIER IS AN ALGEBRAIC UNKNOWN HERE, one per cap, carried in the
     * state vector with a zero mass row, and its equation is the DIFFERENTIATED
     * constraint {@code d/dt (A x + As s) = 0}. The undifferentiated form
     * {@code A x = b} would be an index-2 DAE -- the multiplier appears only after
     * one differentiation -- and RODAS solves index 1. A cap that is not binding
     * keeps its unknown pinned at the inert value (one for a fraction, zero for a
     * flow) by an equation of its own, so the state layout is the same on every
     * segment.</p>
     *
     * <p>A STAGED cap's unknown is the admitted FLOW and not the drain rate the
     * steady state solves for: at the instant a region fills its room is EMPTY, so a
     * rate times zero mass cannot hold the cap.</p>
     */
    private void transient_(final FluidMomentTerms terms, final Matrix C, final double[] Nvec,
                            final int[] cidx, final Matrix[] covblk, final double[] sigma2ss,
                            SolverOptions options, double tol, double[] x0In, int M, int K,
                            final CapacityConstraints con, final CapacityGates gates,
                            final CapacityStaging stg) {
        final int nstate = terms.nstate;
        final int[] idx = terms.covIdx;
        final int nc = idx.length;
        final int ncon = con.b.length;
        final int nstg = stg.n;
        int maxcov = (options.config == null) ? 0 : options.config.dae_maxcov;
        if (maxcov <= 0) {
            maxcov = 25;
        }
        // The covariance adds nc^2 differential states and the Jacobian is formed
        // by finite differences over all of them, so this grows as nc^4. Above the
        // cap the mean is still integrated as a DAE -- conservation stays an
        // equation -- but the variance is held at its stationary value, which is
        // what minnormal does for the whole transient anyway.
        final boolean withcov = nc > 0 && nc <= maxcov;
        if (nc > maxcov) {
            line_warning(mfilename(new Object() {
            }), String.format("The transient covariance would add %d differential states, above "
                    + "the limit of %d set by options.config.dae_maxcov. Integrating the mean as a "
                    + "DAE with the variance held at its stationary value; raise that limit for a "
                    + "time-varying second moment.\n", nc * nc, maxcov * maxcov));
        }

        double[] x0 = x0In.clone();
        final int oSg = nstate;
        final int oM = nstate + nstg;
        final int oCov = nstate + nstg + ncon;
        final int nz = oCov + (withcov ? nc * nc : 0);
        final double[] inert = new double[ncon];
        for (int c = 0; c < ncon; c++) {
            inert[c] = con.staged[c] ? 0.0 : 1.0;   // a flow is inert at zero
        }

        // A STATE ABOVE A CAP IS NOT A STATE THE MODEL CAN BE IN. The initial point
        // comes from the caller, so it may sit outside a cap; holding the cap from
        // there would freeze the violation for the whole horizon, since the equation
        // says the constrained quantity does not MOVE. Start on the cap, and move the
        // excess rather than dropping it: conservation is an algebraic row here, so
        // an initial state that does not satisfy it is an INCONSISTENT
        // initialisation and no index-1 solver may be handed one.
        double[] sg0 = new double[nstg];
        double[] excess = new double[FastMath.max(ncon, 1)];
        List<Integer> over = new ArrayList<Integer>();
        for (int c = 0; c < ncon; c++) {
            double val = 0;
            for (int s = 0; s < nstate; s++) {
                val += con.A.get(c, s) * x0[s];
            }
            if (val > con.b[c] + FastMath.max(1e-9, tol) && con.b[c] > 0) {
                double tot = 0;
                for (int s = 0; s < nstate; s++) {
                    if (con.A.get(c, s) > 0) {
                        tot += x0[s];
                    }
                }
                excess[c] = tot * (1.0 - con.b[c] / val);
                for (int s = 0; s < nstate; s++) {
                    if (con.A.get(c, s) > 0) {
                        x0[s] *= con.b[c] / val;
                    }
                }
                over.add(c);
            }
        }

        // THE INITIAL MODE. A cap the initial state sits ON is already binding, so
        // the first segment must carry it. Feasibility decides: a multiplier above
        // one means the drift is pulling the constrained quantity DOWN.
        double[] mInit = inert.clone();
        Set<Integer> active = new java.util.TreeSet<Integer>();
        if (!over.isEmpty()) {
            int[] trial = new int[over.size()];
            double[] m0 = new double[over.size()];
            for (int a = 0; a < over.size(); a++) {
                trial[a] = over.get(a);
                m0[a] = con.staged[trial[a]] ? 0.0 : 1.0;
            }
            boolean[] okOut = new boolean[1];
            double[] mm = holdMultipliers(terms, gates, stg, con, trial, x0, sg0, sigma2ss,
                    covblk, m0, okOut);
            if (okOut[0]) {
                boolean allFeasible = true;
                for (int a = 0; a < trial.length; a++) {
                    if (!con.staged[trial[a]] && mm[a] > 1.0 + 1e-9) {
                        allFeasible = false;
                    }
                }
                if (allFeasible) {
                    for (int a = 0; a < trial.length; a++) {
                        active.add(trial[a]);
                        mInit[trial[a]] = mm[a];
                    }
                }
            }
        }
        for (int oi = 0; oi < over.size(); oi++) {
            int c = over.get(oi);
            if (excess[c] <= 0) {
                continue;
            }
            int rooms = 0;
            if (nstg > 0 && con.staged[c] && active.contains(c)) {
                for (int j = 0; j < nstg; j++) {
                    if (stg.gatedBy[c][j]) {
                        rooms++;
                    }
                }
            }
            if (rooms > 0) {
                for (int j = 0; j < nstg; j++) {
                    if (stg.gatedBy[c][j]) {
                        sg0[j] += excess[c] / rooms;
                    }
                }
                continue;
            }
            // otherwise onto the coordinates that FEED the capped stations, which is
            // where a held job waits
            boolean[] pool = new boolean[nstate];
            boolean anyFeeder = false;
            for (int s = 0; s < nstate; s++) {
                pool[s] = con.A.get(c, s) <= 0;
            }
            if (gates.Dn != null) {
                boolean[] feeders = new boolean[nstate];
                for (int e = 0; e < terms.eventIdx.length; e++) {
                    if (!gates.gate[c][e]) {
                        continue;
                    }
                    for (int s = 0; s < nstate; s++) {
                        if (gates.Dn.get(s, e) < 0 || gates.DnExt.get(s, e) < 0) {
                            feeders[s] = true;
                        }
                    }
                }
                for (int s = 0; s < nstate; s++) {
                    anyFeeder |= feeders[s] && pool[s];
                }
                if (anyFeeder) {
                    for (int s = 0; s < nstate; s++) {
                        pool[s] = pool[s] && feeders[s];
                    }
                }
            }
            double w = 0;
            int cnt = 0;
            for (int s = 0; s < nstate; s++) {
                if (pool[s]) {
                    w += x0[s];
                    cnt++;
                }
            }
            if (w > FINE) {
                for (int s = 0; s < nstate; s++) {
                    if (pool[s]) {
                        x0[s] += excess[c] * x0[s] / w;
                    }
                }
            } else if (cnt > 0) {
                for (int s = 0; s < nstate; s++) {
                    if (pool[s]) {
                        x0[s] += excess[c] / cnt;
                    }
                }
            }
        }
        if (!over.isEmpty()) {
            StringBuilder names = new StringBuilder();
            for (int oi = 0; oi < over.size(); oi++) {
                if (names.length() > 0) {
                    names.append(", ");
                }
                names.append(con.label[over.get(oi)]);
            }
            final String lbls = names.toString();
            line_warning(mfilename(new Object() {
            }), String.format("The initial state holds more jobs than %s allows; the transient starts "
                    + "on the cap, with the excess where the model would hold it -- in the waiting "
                    + "room, or at the stations feeding the cap.\n", lbls));
        }

        // One differential equation per closed chain is redundant -- the rows of D
        // sum to zero there -- so one is REPLACED by the constraint rather than added
        // to it. The row dropped is the one carrying the most mass at t=0.
        final int[] algrow = new int[C.getNumRows()];
        for (int c = 0; c < C.getNumRows(); c++) {
            int best = -1;
            double bestv = -1;
            for (int s = 0; s < nstate; s++) {
                if (C.get(c, s) == 0) {
                    continue;
                }
                boolean taken = false;
                for (int d = 0; d < c; d++) {
                    taken |= algrow[d] == s;
                }
                if (taken) {
                    continue;
                }
                if (x0[s] > bestv) {
                    bestv = x0[s];
                    best = s;
                }
            }
            if (best < 0) {
                line_error(mfilename(new Object() {
                }), "A conserved chain has no free coordinate to carry its algebraic row.");
            }
            algrow[c] = best;
        }

        final Matrix Dc = new Matrix(nc, terms.eventIdx.length);
        for (int i = 0; i < nc; i++) {
            for (int e = 0; e < terms.eventIdx.length; e++) {
                Dc.set(i, e, terms.D.get(idx[i], e));
            }
        }
        // Sigma(0) = 0 is the consistent initialisation, and the physically right
        // one: the population at t=0 is a known deterministic state.
        double[] z = new double[nz];
        System.arraycopy(x0, 0, z, 0, nstate);
        System.arraycopy(sg0, 0, z, oSg, nstg);
        System.arraycopy(mInit, 0, z, oM, ncon);

        final double[] sigma2ssF = sigma2ss;
        final List<Double> ts = new ArrayList<Double>();
        final List<double[]> zs = new ArrayList<double[]>();
        final Set<Integer> armed = new java.util.HashSet<Integer>();
        for (Integer c : active) {
            if (con.staged[c]) {
                double mass = 0;
                for (int j = 0; j < nstg; j++) {
                    if (stg.gatedBy[c][j]) {
                        mass += sg0[j];
                    }
                }
                if (mass > 1e-8) {
                    armed.add(c);
                }
            }
        }
        this.capacitySwitches = new ArrayList<double[]>();

        double rtol = FastMath.max(tol, 1e-10);
        double atol = FastMath.max(tol, 1e-12);
        double tcur = options.timespan[0];
        final double t1v = options.timespan[1];
        int segMax = 4 * ncon + 8;
        boolean frozen = false;
        int seg = 0;
        for (seg = 0; seg < segMax; seg++) {
            final int[] act = new int[active.size()];
            int ai = 0;
            for (Integer c : active) {
                act[ai++] = c;
            }
            final boolean[] gatedRooms = new boolean[nstg];
            for (int a = 0; a < act.length; a++) {
                if (con.staged[act[a]]) {
                    for (int j = 0; j < nstg; j++) {
                        gatedRooms[j] |= stg.gatedBy[act[a]][j];
                    }
                }
            }
            Rodas.Fcn fcn = new Rodas.Fcn() {
                public void eval(double t, double[] zz, double[] out) {
                    double[] xx = Arrays.copyOfRange(zz, 0, nstate);
                    double[] ss = Arrays.copyOfRange(zz, oSg, oSg + nstg);
                    double[] mm = Arrays.copyOfRange(zz, oM, oM + ncon);
                    double[] s2;
                    Matrix Sc = null;
                    if (withcov) {
                        Sc = new Matrix(nc, nc);
                        for (int i = 0; i < nc; i++) {
                            for (int j = 0; j < nc; j++) {
                                // the equation preserves symmetry; rounding does not
                                Sc.set(i, j, 0.5 * (zz[oCov + i * nc + j] + zz[oCov + j * nc + i]));
                            }
                        }
                        Matrix sig = new Matrix(nstate, nstate);
                        for (int i = 0; i < nc; i++) {
                            for (int j = 0; j < nc; j++) {
                                sig.set(idx[i], idx[j], Sc.get(i, j));
                            }
                        }
                        s2 = sigma2From(sig, terms, cidx, terms.M);
                    } else {
                        s2 = sigma2ssF;
                    }
                    Matrix rr = terms.rates(xx, s2, covblk);
                    double[] mAct = new double[act.length];
                    for (int a = 0; a < act.length; a++) {
                        mAct[a] = mm[act[a]];
                    }
                    Legs lg = legs(terms, gates, stg, con, act, xx, ss, rr, mAct, true);
                    double[] dx = new double[nstate];
                    for (int s = 0; s < nstate; s++) {
                        double acc = 0;
                        for (int e = 0; e < terms.eventIdx.length; e++) {
                            if (gates.Dn == null) {
                                acc += terms.D.get(s, e) * lg.rup[e];
                            } else {
                                acc += gates.Dn.get(s, e) * lg.rup[e]
                                        + gates.DnExt.get(s, e) * lg.rin[e]
                                        + gates.Dp.get(s, e) * lg.rin[e];
                            }
                        }
                        dx[s] = acc;
                        out[s] = acc;
                    }
                    // the algebraic rows: the mass matrix has zeroed these, so what
                    // is written here is the CONSTRAINT RESIDUAL, not a derivative
                    for (int c = 0; c < algrow.length; c++) {
                        double acc = 0;
                        for (int s = 0; s < nstate; s++) {
                            acc += C.get(c, s) * xx[s];
                        }
                        for (int j = 0; j < nstg; j++) {
                            acc += C.get(c, nstate + j) * ss[j];
                        }
                        out[algrow[c]] = acc - Nvec[c];
                    }
                    for (int j = 0; j < nstg; j++) {
                        out[oSg + j] = gatedRooms[j] ? lg.ds[j] : ss[j];
                    }
                    for (int c = 0; c < ncon; c++) {
                        boolean isActive = false;
                        for (int a = 0; a < act.length; a++) {
                            isActive |= act[a] == c;
                        }
                        if (isActive) {
                            double acc = 0;
                            for (int s = 0; s < nstate; s++) {
                                acc += con.A.get(c, s) * dx[s];
                            }
                            for (int j = 0; j < nstg; j++) {
                                acc += con.As.get(c, j) * lg.ds[j];
                            }
                            out[oM + c] = acc;
                        } else {
                            out[oM + c] = mm[c] - inert[c];
                        }
                    }
                    if (withcov) {
                        Matrix A = terms.jacobian(xx, s2, covblk);
                        Matrix Ac = new Matrix(nc, nc);
                        for (int i = 0; i < nc; i++) {
                            for (int j = 0; j < nc; j++) {
                                Ac.set(i, j, A.get(idx[i], idx[j]));
                            }
                        }
                        for (int i = 0; i < nc; i++) {
                            for (int j = 0; j < nc; j++) {
                                double acc = 0;
                                for (int l = 0; l < nc; l++) {
                                    acc += Ac.get(i, l) * Sc.get(l, j) + Sc.get(i, l) * Ac.get(j, l);
                                }
                                for (int e = 0; e < terms.eventIdx.length; e++) {
                                    acc += Dc.get(i, e) * lg.rup[e] * Dc.get(j, e);
                                }
                                out[oCov + i * nc + j] = acc;
                            }
                        }
                    }
                }
            };

            // The mass matrix: the identity with one row zeroed per closed chain,
            // per un-gated room and per multiplier. Diagonal, hence MLMAS = MUMAS = 0
            // and RODAS's banded path (IJOB 3).
            final double[] diag = new double[nz];
            Arrays.fill(diag, 1.0);
            for (int c = 0; c < algrow.length; c++) {
                diag[algrow[c]] = 0.0;
            }
            for (int j = 0; j < nstg; j++) {
                diag[oSg + j] = gatedRooms[j] ? 1.0 : 0.0;
            }
            for (int c = 0; c < ncon; c++) {
                diag[oM + c] = 0.0;
            }
            Rodas.Mas mas = new Rodas.Mas() {
                public void eval(double[][] am) {
                    for (int j = 0; j < nz; j++) {
                        am[0][j] = diag[j];
                    }
                }
            };

            final double[] gprev = capEvents(con, stg, act, armed, z, nstate, oSg, oM, ncon, nstg);
            final int[] hitIdx = new int[]{-1};
            final double[] hitT = new double[]{0};
            final double[][] hitZ = new double[1][];
            Rodas.Solout solout = new Rodas.Solout() {
                public int eval(int nr, double xold, double xcur, double[] y, Rodas.Dense dns) {
                    if (ncon == 0) {
                        ts.add(xcur);
                        zs.add(y.clone());
                        return 0;
                    }
                    // arm a staged release only once its room has really filled, or
                    // the release would fire at the activation instant, where the room
                    // is empty by construction. The latch only ever goes false -> true.
                    for (int a = 0; a < act.length; a++) {
                        int c = act[a];
                        if (con.staged[c] && !armed.contains(c)) {
                            double mass = 0;
                            for (int j = 0; j < nstg; j++) {
                                if (stg.gatedBy[c][j]) {
                                    mass += y[oSg + j];
                                }
                            }
                            if (mass > 1e-8) {
                                armed.add(c);
                            }
                        }
                    }
                    double[] gcur = capEvents(con, stg, act, armed, y, nstate, oSg, oM, ncon, nstg);
                    int cross = -1;
                    for (int c = 0; c < ncon; c++) {
                        if (gprev[c] > 0 && gcur[c] <= 0) {
                            cross = c;
                            break;
                        }
                    }
                    System.arraycopy(gcur, 0, gprev, 0, ncon);
                    if (cross < 0) {
                        ts.add(xcur);
                        zs.add(y.clone());
                        return 0;
                    }
                    // BISECT ON THE INTERPOLANT rather than on the integration: RODAS
                    // carries a third-order interpolant over the step just accepted,
                    // so locating the crossing costs no extra step.
                    double lo = xold;
                    double hi = xcur;
                    double[] ymid = new double[nz];
                    for (int bit = 0; bit < 60; bit++) {
                        double mid = 0.5 * (lo + hi);
                        for (int i = 0; i < nz; i++) {
                            ymid[i] = dns.value(i, mid);
                        }
                        double[] gmid = capEvents(con, stg, act, armed, ymid, nstate, oSg, oM,
                                ncon, nstg);
                        if (gmid[cross] > 0) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                        if (hi - lo <= 1e-12 * FastMath.max(1.0, FastMath.abs(hi))) {
                            break;
                        }
                    }
                    double[] yhit = new double[nz];
                    for (int i = 0; i < nz; i++) {
                        yhit[i] = dns.value(i, hi);
                    }
                    // THE OVERSHOOTING STEP IS NOT REPORTED: the path stops AT the cap
                    // and the next segment starts there.
                    hitIdx[0] = cross;
                    hitT[0] = hi;
                    hitZ[0] = yhit;
                    ts.add(hi);
                    zs.add(yhit.clone());
                    return -1;
                }
            };

            Rodas.Options ro = new Rodas.Options();
            ro.ijac = 0;
            ro.mljac = nz;
            ro.mujac = nz;
            ro.ifcn = 1;
            ro.imas = 1;
            ro.mas = mas;
            ro.mlmas = 0;
            ro.mumas = 0;
            ro.iout = 1;
            ro.solout = solout;
            try {
                Rodas.Result res = Rodas.integrate(nz, fcn, tcur, z, t1v, 1e-6,
                        new double[]{rtol}, new double[]{atol}, 0, ro);
                if (res.idid != Rodas.IDID_SUCCESS && hitIdx[0] < 0) {
                    throw new RuntimeException("RODAS returned idid=" + res.idid);
                }
            } catch (RuntimeException e) {
                line_warning(mfilename(new Object() {
                }), String.format("The transient DAE failed (%s); reporting the seed trajectory and "
                        + "the stationary variance instead.\n", e.getMessage()));
                return;
            }

            if (hitIdx[0] < 0) {
                break;   // the horizon was reached with this set of caps binding
            }
            tcur = hitT[0];
            z = hitZ[0].clone();
            int c = hitIdx[0];
            double[] xh = Arrays.copyOfRange(z, 0, nstate);
            double[] sgh = Arrays.copyOfRange(z, oSg, oSg + nstg);
            double[] s2h = sigma2ss;
            if (withcov) {
                Matrix sig = new Matrix(nstate, nstate);
                for (int i = 0; i < nc; i++) {
                    for (int j = 0; j < nc; j++) {
                        sig.set(idx[i], idx[j],
                                0.5 * (z[oCov + i * nc + j] + z[oCov + j * nc + i]));
                    }
                }
                s2h = sigma2From(sig, terms, cidx, terms.M);
            }
            if (active.contains(c)) {
                active.remove(c);
                armed.remove(c);
                z[oM + c] = inert[c];
                this.capacitySwitches.add(new double[]{tcur, c, 0});   // 0: release
            } else {
                // ACTIVATING NEEDS THE MULTIPLIER THAT HOLDS THE CAP, both because
                // the restart must satisfy the algebraic rows and because that
                // multiplier IS the feasibility test.
                Set<Integer> trialSet = new java.util.TreeSet<Integer>(active);
                trialSet.add(c);
                int[] trial = new int[trialSet.size()];
                double[] m0 = new double[trial.length];
                int ti = 0;
                for (Integer k : trialSet) {
                    trial[ti] = k;
                    m0[ti] = active.contains(k) ? z[oM + k] : (con.staged[k] ? 0.0 : 1.0);
                    ti++;
                }
                boolean[] okOut = new boolean[1];
                double[] mm = holdMultipliers(terms, gates, stg, con, trial, xh, sgh, s2h,
                        covblk, m0, okOut);
                boolean feasible = okOut[0];
                for (int a = 0; a < trial.length; a++) {
                    if (!con.staged[trial[a]] && mm[a] > 1.0 + 1e-9) {
                        feasible = false;
                    }
                }
                if (feasible) {
                    active = trialSet;
                    for (int a = 0; a < trial.length; a++) {
                        z[oM + trial[a]] = mm[a];
                    }
                    this.capacitySwitches.add(new double[]{tcur, c, 1});   // 1: activate
                } else {
                    // a crossing whose cap cannot be held: stop rather than freeze a
                    // violation into the algebraic rows
                    this.capacitySwitches.add(new double[]{tcur, c, 2});   // 2: unheld
                    frozen = true;
                    break;
                }
            }
            if (tcur >= t1v - 1e-12) {
                break;
            }
        }
        if (seg >= segMax) {
            frozen = true;
        }
        if (frozen) {
            final int segF = seg;
            final double tcurF = tcur;
            line_warning(mfilename(new Object() {
            }), String.format("The transient stopped switching after %d segments at t=%g of %g: the "
                    + "reported trajectory ends there rather than continuing with a capacity "
                    + "constraint that does not hold. Shorten options.timespan, or use "
                    + "SolverCTMC/SolverJMT/SolverSSA/SolverLDES.\n", segF + 1, tcurF, t1v));
        }

        int nt = ts.size();
        this.xvec_t = new Matrix(nt, nstate);
        for (int s = 0; s < nt; s++) {
            for (int j = 0; j < nstate; j++) {
                this.xvec_t.set(s, j, zs.get(s)[j]);
            }
        }
        this.tvar = new Matrix(nt, 1);
        for (int s = 0; s < nt; s++) {
            this.tvar.set(s, 0, ts.get(s));
        }
        if (withcov) {
            this.sigmat = new Matrix[nt];
            for (int s = 0; s < nt; s++) {
                Matrix sig = new Matrix(nstate, nstate);
                for (int i = 0; i < nc; i++) {
                    for (int j = 0; j < nc; j++) {
                        sig.set(idx[i], idx[j],
                                0.5 * (zs.get(s)[nstate + i * nc + j] + zs.get(s)[nstate + j * nc + i]));
                    }
                }
                this.sigmat[s] = sig;
            }
            this.qVart = new Matrix[M * K];
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    Matrix v = new Matrix(nt, 1);
                    int[] blk = terms.classBlock[i * K + k];
                    for (int s = 0; s < nt; s++) {
                        v.set(s, 0, blk.length == 0 ? 0.0
                                : FastMath.max(0, blockSum(this.sigmat[s], blk)));
                    }
                    this.qVart[i * K + k] = v;
                }
            }
        }
    }
}
