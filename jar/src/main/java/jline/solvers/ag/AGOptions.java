/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.util.ArrayList;
import java.util.List;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

/**
 * Options of the agent-based (RCAT) solver.
 *
 * <p>Beside the truncation level of an open agent, these carry the EXECUTION
 * BACKEND of the reversed-rate fixed point. The backend decides who evaluates an
 * agent, never what the agent evaluates to: agent k reads the rest of the model
 * only through the scalar reversed rates x, and the sweep is Jacobi, so the
 * agent order is immaterial and every backend walks the same iterates.</p>
 */
public class AGOptions extends SolverOptions {

    /** One agent per (station, class); solve them one after another. */
    public static final String EXEC_SERIAL = "serial";
    /**
     * Fan the agents of a sweep out over a local thread pool.
     *
     * Named "parallel", with "para" accepted as an alias -- the same pair
     * SolverSSA's replica analyzer answers to, so one spelling convention covers
     * both solvers. It was called "threads" until 2026-08-19; that name is no
     * longer accepted, and {@link AgExec#isParallel(String)} is the only place
     * either spelling is recognised.
     */
    public static final String EXEC_PARALLEL = "parallel";
    /** The accepted alias of {@link #EXEC_PARALLEL}, as SolverSSA spells it. */
    public static final String EXEC_PARA = "para";
    /** Partition the agents over remote ag-worker processes. */
    public static final String EXEC_CLUSTER = "cluster";

    /**
     * Truncation level of an OPEN agent's queue-length dimension. A closed class
     * uses its own population instead, so this bounds only the open agents;
     * 'inapinf' ignores it and solves them on the infinite state space.
     */
    public int maxStates = 100;

    /** One of {@link #EXEC_SERIAL}, {@link #EXEC_PARALLEL} (or {@link #EXEC_PARA}), {@link #EXEC_CLUSTER}. */
    public String exec = EXEC_SERIAL;

    /**
     * Thread-pool size for {@link #EXEC_PARALLEL}; 0 means one thread per
     * available processor. Pinned rather than derived per sweep so that a run is
     * reproducible on a machine whose load changes under it.
     */
    public int nworkers = 0;

    /** Worker addresses ("host:port") for {@link #EXEC_CLUSTER}. */
    public List<String> endpoints = new ArrayList<String>();

    /**
     * Seconds to wait on a worker before solving its agents locally instead. A
     * lost worker is never fatal: any agent can be solved anywhere given x, so a
     * cluster run degrades to a slower run and never to a wrong one.
     */
    public double workerTimeout = 30.0;

    public AGOptions() {
        super(SolverType.AG);
        this.method = "default";
        this.iter_max = 100;
        // MATLAB's SolverOptions('AG') sets iter_max and lets iter_tol fall
        // through to the global 1e-4. The reversed-rate fixed point converges
        // only linearly, so a tighter tolerance stops ~1e-4 further along and
        // the codebases disagree on a goldened example; keep them equal.
        this.iter_tol = 1e-4;
    }

    /**
     * Preserve the AG fields across the defensive copy the Solver constructor
     * takes.
     *
     * <p>WITHOUT THIS OVERRIDE THE WHOLE EXECUTION AXIS IS SILENTLY INERT.
     * {@code Solver} stores {@code options.copy()}, and the base copy returns a
     * plain {@link SolverOptions}, so {@link #exec}, {@link #endpoints},
     * {@link #maxStates}, {@link #nworkers} and {@link #workerTimeout} would all
     * be dropped between construction and the solve -- a run asked for on the
     * cluster would quietly execute serially and report nothing. Worse, a test
     * comparing two backends would then compare two SERIAL runs and pass
     * vacuously, which is how this was found.</p>
     */
    @Override
    public SolverOptions copy() {
        AGOptions out = new AGOptions();
        SolverOptions base = super.copy();
        out.verbose = base.verbose;
        out.samples = base.samples;
        out.seed = base.seed;
        out.timespan = base.timespan;
        out.iter_max = base.iter_max;
        out.iter_tol = base.iter_tol;
        out.tol = base.tol;
        out.cutoff = base.cutoff;
        out.init_sol = base.init_sol;
        out.remote = base.remote;
        out.remote_endpoint = base.remote_endpoint;
        out.cache = base.cache;
        out.keep = base.keep;
        out.force = base.force;
        out.hide_immediate = base.hide_immediate;
        out.lang = base.lang;
        out.method = base.method;
        out.stiff = base.stiff;
        out.config = base.config;

        out.maxStates = this.maxStates;
        out.exec = this.exec;
        out.nworkers = this.nworkers;
        out.endpoints = new ArrayList<String>(this.endpoints);
        out.workerTimeout = this.workerTimeout;
        return out;
    }

    public AGOptions maxStates(int maxStates) {
        this.maxStates = maxStates;
        return this;
    }

    public AGOptions exec(String exec) {
        this.exec = exec;
        return this;
    }

    public AGOptions nworkers(int nworkers) {
        this.nworkers = nworkers;
        return this;
    }

    public AGOptions endpoint(String hostPort) {
        this.endpoints.add(hostPort);
        return this;
    }

    public AGOptions workerTimeout(double seconds) {
        this.workerTimeout = seconds;
        return this;
    }
}
