/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;


import jline.lang.Element;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ServerType;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Element of a LayeredNetwork model
 */
public class LayeredNetworkElement extends Element {
    public static final int ACTIVITY = 3;
    public static final int CALL = 4;
    public static final int ENTRY = 2;
    public static final int HOST = 0;
    public static final int PROCESSOR = 0;
    public static final int TASK = 1;
    public LayeredNetwork model;

    /** One admission constraint row declared by operand name. */
    public static class LinConRow {
        public final List<String> names;
        public final double[] coeffs;
        public final double cap;

        public LinConRow(List<String> names, double[] coeffs, double cap) {
            this.names = names;
            this.coeffs = coeffs;
            this.cap = cap;
        }
    }

    /** Matrix(C,K): admission constraint matrix on this server's layer station. */
    public Matrix linConA;
    /** Matrix(C,1): admission constraint capacities. */
    public Matrix linConB;
    /** Rows declared by operand name, resolved in LayeredNetwork.getStruct. */
    public List<LinConRow> linConRows = new ArrayList<LinConRow>();

    /** Vector alpha(n): rate scaling of this server's layer station when it holds n jobs. */
    public Matrix lldScaling;
    /** Handle beta(n): product-form class-dependent scaling, n counted over this server's operands. */
    public SerializableFunction<Matrix, Matrix> lcdScaling;
    /** Peak rate scaling per operand of lcdScaling, normalizes Util = T*S/peak. */
    public Matrix lcdScalingPeak;
    /** Handle eta(n): non-product-form joint-dependent scaling, n counted over this server's operands. */
    public SerializableFunction<Matrix, Matrix> ljdScaling;
    /** Peak rate scaling per operand of ljdScaling, normalizes Util = T*S/peak. */
    public Matrix ljdScalingPeak;

    /**
     * Heterogeneous server pools with a compatibility graph over this server's
     * operands. Each entry is one declared pool; SolverLN lowers the whole list
     * to the activated-server rate of {@code SnCompatRate}, carried onto the
     * layer station as a joint dependence.
     */
    public List<ServerPool> serverPools = new ArrayList<ServerPool>();

    /** One declared pool, its operands held by name until getStruct resolves them. */
    public static class ServerPool implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public String name;
        public double count;
        public double rate;
        public List<String> compatible;

        public ServerPool(String name, double count, double rate, List<String> compatible) {
            this.name = name;
            this.count = count;
            this.rate = rate;
            this.compatible = compatible;
        }
    }

    public LayeredNetworkElement(String name) {
        super(name);
    }

    /**
     * <p>Appends one admission constraint row naming its operands, so the meaning
     * does not depend on declaration order:</p>
     *
     * <pre>
     *   t2.addConstraint(Arrays.asList(e2, e3), new double[]{1, 1}, 2); // n(E2) + n(E3) &lt;= 2
     * </pre>
     *
     * <p>Operands are the entries of a Task, or the tasks of a Host. Names are
     * resolved against the model in LayeredNetwork.getStruct, where an operand
     * that does not belong to this server is an error rather than a silent
     * mis-mapping.</p>
     *
     * @param operands entries of this task, or tasks on this host
     * @param coeffs   one coefficient per operand, or null for all ones
     * @param cap      right-hand side of the row
     */
    public void addConstraint(List<? extends LayeredNetworkElement> operands, double[] coeffs, double cap) {
        if (operands == null || operands.isEmpty()) {
            throw new IllegalArgumentException("Admission constraint requires at least one operand.");
        }
        List<String> names = new ArrayList<String>();
        for (int k = 0; k < operands.size(); k++) {
            LayeredNetworkElement op = operands.get(k);
            if (op == null) {
                throw new IllegalArgumentException("Admission constraint operand " + (k + 1) + " is null.");
            }
            names.add(op.getName());
        }
        addConstraintByName(names, coeffs, cap);
    }

    /**
     * Name-based form of {@link #addConstraint(List, double[], double)}, for
     * operands not yet held as handles.
     *
     * @param names operand names, entries of this task or tasks on this host
     * @param coeffs one coefficient per operand, or null for all ones
     * @param cap    right-hand side of the row
     */
    public void addConstraintByName(List<String> names, double[] coeffs, double cap) {
        if (!(this instanceof Task) && !(this instanceof Host)) {
            throw new IllegalArgumentException("Admission constraints can only be set on a Task or a Host, which are the only elements that become server stations in a layer.");
        }
        if (names == null || names.isEmpty()) {
            throw new IllegalArgumentException("Admission constraint requires at least one operand.");
        }
        int n = names.size();
        double[] c;
        if (coeffs == null || coeffs.length == 0) {
            c = new double[n];
            for (int k = 0; k < n; k++) {
                c[k] = 1.0;
            }
        } else if (coeffs.length == 1 && n > 1) {
            c = new double[n];
            for (int k = 0; k < n; k++) {
                c[k] = coeffs[0];
            }
        } else {
            c = coeffs.clone();
        }
        if (c.length != n) {
            throw new IllegalArgumentException("Admission constraint has " + n + " operands but " + c.length + " coefficients.");
        }
        boolean allZero = true;
        for (int k = 0; k < n; k++) {
            if (!Double.isFinite(c[k]) || c[k] < 0) {
                throw new IllegalArgumentException("Admission constraint coefficients must be finite and non-negative.");
            }
            if (c[k] != 0) {
                allZero = false;
            }
        }
        if (allZero) {
            throw new IllegalArgumentException("Admission constraint has all-zero coefficients, which constrains nothing.");
        }
        Set<String> distinct = new HashSet<String>(names);
        if (distinct.size() != n) {
            throw new IllegalArgumentException("Admission constraint names the same operand more than once; give it a single combined coefficient instead.");
        }
        if (!Double.isFinite(cap) || cap < 1) {
            throw new IllegalArgumentException("Admission constraint capacity must be a finite scalar of at least 1.");
        }
        this.linConRows.add(new LinConRow(new ArrayList<String>(names), c, cap));
    }

    /**
     * <p>Raw form of addConstraint, for programmatic construction. Declares
     * A*n &lt;= b on the station that represents this server in its layer, where n
     * counts the jobs in service or queueing at that station.</p>
     *
     * <p>Columns of A are indexed positionally by the entries of a Task, or by
     * the tasks of a Host, in declaration order, so the mapping shifts if an
     * entry is added later; prefer addConstraint, which names its operands. Only
     * the column count is checked, in LayeredNetwork.getStruct, since entries may
     * be added after this call. Rows from both forms are concatenated.</p>
     *
     * @param A constraint matrix, one column per entry or task
     * @param b capacity vector, one entry per row of A
     */
    public void setConstraint(Matrix A, Matrix b) {
        if (!(this instanceof Task) && !(this instanceof Host)) {
            throw new IllegalArgumentException("Admission constraints can only be set on a Task or a Host, which are the only elements that become server stations in a layer.");
        }
        if (A == null || b == null || A.isEmpty() || b.isEmpty()) {
            throw new IllegalArgumentException("Constraint matrix A and capacity vector b must be non-empty.");
        }
        Matrix bcol = b;
        if (b.getNumRows() == 1 && b.getNumCols() > 1) {
            bcol = b.transpose();
        }
        if (A.getNumRows() != bcol.getNumRows()) {
            throw new IllegalArgumentException("A and b must have matching number of rows.");
        }
        for (int i = 0; i < A.getNumRows(); i++) {
            boolean allZero = true;
            for (int j = 0; j < A.getNumCols(); j++) {
                double v = A.get(i, j);
                if (!Double.isFinite(v) || v < 0) {
                    throw new IllegalArgumentException("Constraint matrix A must be finite and non-negative.");
                }
                if (v != 0) {
                    allZero = false;
                }
            }
            if (allZero) {
                throw new IllegalArgumentException("Constraint matrix A has an all-zero row, which constrains nothing.");
            }
            double bv = bcol.get(i, 0);
            if (!Double.isFinite(bv) || bv < 1) {
                throw new IllegalArgumentException("Capacity vector b must be finite and at least 1.");
            }
        }
        this.linConA = A;
        this.linConB = bcol;
    }

    /**
     * Positional constraint pair declared on this element, before name resolution.
     *
     * @return a two-element array holding A and b, either of which may be null
     */
    public Matrix[] getLinearConstraints() {
        return new Matrix[]{this.linConA, this.linConB};
    }

    /**
     * Whether this element declares any admission constraint, in either form.
     *
     * @return true if a positional pair or at least one named row is present
     */
    public boolean hasLinearConstraints() {
        return (this.linConA != null && this.linConB != null) || !this.linConRows.isEmpty();
    }

    /**
     * <p>Sets the service-rate scaling of the station that represents this server
     * in its layer: alpha[n] applies when that station holds n jobs in total, as
     * in Queue.setLoadDependence. The scaling multiplies the station rate on top
     * of its multiplicity, so a multi-server host applies min(n,m)*alpha[n].</p>
     *
     * @param alpha row vector of positive scalings, indexed by station population
     */
    public void setLoadDependence(Matrix alpha) {
        assertRateDependent("Load");
        if (alpha == null || alpha.isEmpty()) {
            throw new IllegalArgumentException("Load-dependence scalings must be non-empty.");
        }
        if (alpha.getNumCols() == 1 && alpha.getNumRows() > 1) {
            alpha = alpha.transpose(); // the station reads alpha as a row, as Queue.setLoadDependence does
        }
        for (int i = 0; i < alpha.length(); i++) {
            double v = alpha.get(i);
            if (!Double.isFinite(v) || v <= 0) {
                throw new IllegalArgumentException("Load-dependence scalings must be finite and positive.");
            }
        }
        this.lldScaling = alpha;
    }

    /**
     * <p>Sets a class-dependent service-rate scaling on this server's layer
     * station. The handle takes the per-operand population vector of this server:
     * entry j counts the jobs held on behalf of operand j, which is task j of a
     * Host or entry j of a Task, in the same tasksof/entriesof order as the
     * columns of setConstraint. It returns a 1x1 scaling shared by every operand,
     * or a per-operand row vector.</p>
     *
     * <p>Product form holds only where an operand occupies the layer station
     * through a single job class; otherwise SolverLN emits the equivalent joint
     * dependence, which is numerically identical but carries no exactness
     * guarantee.</p>
     *
     * @param beta             maps the per-operand population vector to a rate scaling
     * @param peakRatePerOperand 1x1 (broadcast) or per-operand peak rate scaling
     */
    public void setClassDependence(SerializableFunction<Matrix, Matrix> beta, Matrix peakRatePerOperand) {
        assertRateDependent("Class");
        assertDependenceHandle(beta, peakRatePerOperand, "Class");
        this.lcdScaling = beta;
        this.lcdScalingPeak = peakRatePerOperand;
    }

    /**
     * <p>Sets a joint-dependent (non-product-form) service-rate scaling on this
     * server's layer station. The handle reads the per-operand population vector
     * arbitrarily (e.g. min(n[0],c)), so solvers treat it as an approximation.
     * Operand order and the required peak rate are as in
     * {@link #setClassDependence}.</p>
     *
     * @param eta                maps the per-operand population vector to a rate scaling
     * @param peakRatePerOperand 1x1 (broadcast) or per-operand peak rate scaling
     */
    public void setJointDependence(SerializableFunction<Matrix, Matrix> eta, Matrix peakRatePerOperand) {
        assertRateDependent("Joint");
        assertDependenceHandle(eta, peakRatePerOperand, "Joint");
        if (!this.serverPools.isEmpty()) {
            throw new IllegalArgumentException(this.getName() + " already declares server pools, "
                    + "which are themselves a rate law, so it cannot also take a joint dependence.");
        }
        this.ljdScaling = eta;
        this.ljdScalingPeak = peakRatePerOperand;
    }

    /**
     * Declares one pool of {@code serverType.getNumOfServers()} identical
     * servers, each running at {@code serverType.getRate()}, eligible only for
     * the operands the ServerType names.
     *
     * <p>The operands are the tasks of a Processor, or the entries of a Task.
     * They are resolved against the model in {@code LayeredNetwork.getStruct},
     * where an operand that does not belong to this server is an error rather
     * than a silent mis-mapping, exactly as for an admission constraint.
     *
     * <p>SolverLN lowers the whole declaration to the activated-server rate of
     * {@link jline.api.sn.SnCompatRate}, carried onto the layer station as a
     * joint dependence, so the pools are an APPROXIMATION in a layer for the
     * same reason {@code setJointDependence} is.
     *
     * @param serverType the pool to declare
     */
    public void addServerType(ServerType serverType) {
        assertRateDependent("Compatibility");
        if (this.ljdScaling != null) {
            throw new IllegalArgumentException(this.getName() + " already declares a joint "
                    + "dependence, so it cannot also declare server pools, which are a rate law "
                    + "of their own.");
        }
        if (serverType.getNumOfServers() < 1) {
            throw new IllegalArgumentException("Server pool '" + serverType.getName()
                    + "' must hold at least one server");
        }
        List<String> compat = serverType.getCompatibleOperands();
        if (compat == null || compat.isEmpty()) {
            throw new IllegalArgumentException("Server pool '" + serverType.getName()
                    + "' is compatible with no operand, so it can never serve");
        }
        for (int k = 0; k < this.serverPools.size(); k++) {
            if (this.serverPools.get(k).name.equals(serverType.getName())) {
                throw new IllegalArgumentException("Server pool '" + serverType.getName()
                        + "' is already declared on " + this.getName());
            }
        }
        Set<String> seen = new HashSet<String>(compat);
        if (seen.size() != compat.size()) {
            throw new IllegalArgumentException("Server pool '" + serverType.getName()
                    + "' names the same operand more than once");
        }
        serverType.setId(this.serverPools.size());
        this.serverPools.add(new ServerPool(serverType.getName(),
                serverType.getNumOfServers(), serverType.getRate(),
                new ArrayList<String>(compat)));
    }

    /** Declared compatibility pools. */
    public List<ServerPool> getServerTypes() {
        return this.serverPools;
    }

    /** Whether this element declares compatibility pools. */
    public boolean hasServerPools() {
        return !this.serverPools.isEmpty();
    }

    /**
     * Whether this element declares any service-rate dependence.
     *
     * @return true if a load, class or joint dependence is present
     */
    public boolean hasRateDependence() {
        return this.lldScaling != null || this.lcdScaling != null || this.ljdScaling != null
                || !this.serverPools.isEmpty();
    }

    /** Rejects elements that do not become a rate-scalable layer station. */
    private void assertRateDependent(String what) {
        SchedStrategy sched;
        if (this instanceof Task) {
            sched = ((Task) this).getScheduling();
        } else if (this instanceof Host) {
            sched = ((Host) this).getScheduling();
        } else {
            throw new IllegalArgumentException(what + "-dependence can only be set on a Task or a Host, which are the only elements that become server stations in a layer.");
        }
        if (sched != SchedStrategy.PS && sched != SchedStrategy.FCFS) {
            throw new IllegalArgumentException(what + "-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) servers, but " + this.getName() + " is scheduled " + SchedStrategy.toText(sched) + ".");
        }
    }

    /** Common validation of a class- or joint-dependence declaration. */
    private static void assertDependenceHandle(SerializableFunction<Matrix, Matrix> f, Matrix peak, String what) {
        if (f == null) {
            throw new IllegalArgumentException(what + " dependence must be specified through a function handle.");
        }
        if (peak == null || peak.isEmpty()) {
            throw new IllegalArgumentException(what + " dependence requires an explicit peak rate: pass a 1x1 matrix (identical peak for every operand) or a per-operand vector.");
        }
        for (int i = 0; i < peak.length(); i++) {
            double v = peak.get(i);
            if (!Double.isFinite(v) || v <= 0) {
                throw new IllegalArgumentException("peakRatePerOperand must be finite and positive.");
            }
        }
    }

}
