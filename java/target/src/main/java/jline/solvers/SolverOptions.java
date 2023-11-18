// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.lang.constant.VerboseLevel;
import jline.solvers.ssa.strategies.TauLeapingConfig;
import jline.util.Matrix;
import jline.lang.constant.SolverType;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Double.POSITIVE_INFINITY;

import odesolver.LSODA;

public class SolverOptions {

  public static class Config {

    public String highvar; // TODO: enum?
    public String multiserver; // TODO: enum?
    public String np_priority; // TODO: enum?
    public List<Double> pstar; // For p-norm smoothing in SolverFluid
    public String fork_join;
    public TauLeapingConfig tau_leaping;
    public String merge;
    public String compress;
    public int space_max;
    public boolean interlocking;
  }

  public static class ODESolvers {
    public double odeminstep;
    public double odemaxstep;

    public FirstOrderIntegrator fastODESolver;
    public FirstOrderIntegrator accurateODESolver;
    public LSODA fastStiffODESolver;
    public LSODA accurateStiffODESolver;
  }

  public boolean cache;
  public double cutoff;
  public Config config;
  public boolean force;
  public boolean hide_immediate;
  public Matrix init_sol;
  public int iter_max;
  public double iter_tol;
  public double tol;
  public boolean keep;
  public String lang;
  public String method;
  public boolean remote;
  public String remote_endpoint;
  public ODESolvers odesolvers;
  public int samples;
  public int seed;
  public boolean stiff;
  public double[] timespan;
  public VerboseLevel verbose;

  public SolverOptions() {

    // Solver Default Options
    this.cache = true;
    this.cutoff = POSITIVE_INFINITY;
    this.config = new Config();
    this.config.tau_leaping = null;
    this.config.pstar = new ArrayList<>();
    this.config.fork_join = "default";
    this.force = false;
    this.hide_immediate = true; // Hide immediate transitions if possible
    this.init_sol = new Matrix(0, 0);
    this.iter_max = 10;
    this.iter_tol = 0.0001; // Convergence tolerance to stop iterations
    this.tol = 0.0001; // Tolerance for all other uses
    this.keep = false;
    this.lang = "java";
    this.method = "default";
    this.remote = false;
    this.remote_endpoint = "127.0.0.1";

    this.odesolvers = new ODESolvers();
    this.odesolvers.odeminstep = 0.001;
    //this.odeMinStep = 0.00000001;
    this.odesolvers.odemaxstep = POSITIVE_INFINITY;
    this.odesolvers.fastODESolver = null; // TODO
    this.odesolvers.accurateODESolver =
            new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
    this.odesolvers.fastStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3); // TODO
    this.odesolvers.accurateStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);

    this.samples = 10000;
    this.seed = Math.toIntExact(Math.round((Math.random() * (1e6 - 1)) + 1));
    this.stiff = true;
    this.timespan = new double[2];
    this.timespan[0] = POSITIVE_INFINITY;
    this.timespan[1] = POSITIVE_INFINITY;
    this.verbose = VerboseLevel.STD;
  }

  public SolverOptions(SolverType solverType) {
    this();
    if (solverType==null) {
      return;
    }

    // Solver-specific Defaults
    switch (solverType) {
      case Env:
        this.iter_max = 100;
        this.verbose = VerboseLevel.SILENT;
        break;
      case Fluid:
        this.config.highvar = "none";
        this.iter_max = 5;
        this.timespan[0] = 0;
        break;
      case LN:
        this.config.interlocking = true;
        this.config.multiserver = "default";
        this.iter_tol = 0.05;
        this.iter_max = 100;
        break;
      case LQNS:
        this.keep = true;
        break;
      case MAM:
        this.iter_max = 100;
        break;
      case MVA:
        this.iter_max = 1000;
        this.config.highvar = "none";
        this.config.multiserver = "default";
        this.config.np_priority = "default";
        this.config.fork_join = "default";
        break;
      case NC:
        this.samples = 100000;
        this.config.highvar = "interp";
        break;
      case SSA:
        this.timespan[0] = 0;
        this.timespan[1] = POSITIVE_INFINITY;
        this.config.tau_leaping = null; // null stands for disabled
        break;
      case JMT:
        //TODO add configs
        break;
      default: // Global options unless overridden by a solver
    }
  }

  public void setODEMinStep(double odeMinStep) {
    this.odesolvers.odeminstep = odeMinStep;
    this.odesolvers.fastODESolver = null; // TODO
    this.odesolvers.accurateODESolver =
            new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
    this.odesolvers.fastStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3); // TODO
    this.odesolvers.accurateStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);
  }

  public void setODEMaxStep(double odeMaxStep) {
    this.odesolvers.odemaxstep = odeMaxStep;
    this.odesolvers.fastODESolver = null; // TODO
    this.odesolvers.accurateODESolver =
            new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
    this.odesolvers.fastStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3); // TODO
    this.odesolvers.accurateStiffODESolver =
            new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);
  }
}
