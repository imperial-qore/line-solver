// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.lang.JobClass;
import jline.lang.nodes.Station;

import java.util.Map;

// Class for handles for the mean performance metrics
public class SolverHandles {

  public SolverHandles(Map<Station, Map<JobClass, Metric>> Q, Map<Station, Map<JobClass, Metric>> U, Map<Station, Map<JobClass, Metric>> R, Map<Station, Map<JobClass, Metric>> T, Map<Station, Map<JobClass, Metric>> A) {
    this.Q = Q;
    this.U = U;
    this.R = R;
    this.T = T;
    this.A = A;
  }

  public Map<Station, Map<JobClass, Metric>> Q;
  public Map<Station, Map<JobClass, Metric>> U;
  public Map<Station, Map<JobClass, Metric>> R;
  public Map<Station, Map<JobClass, Metric>> T;
  public Map<Station, Map<JobClass, Metric>> A;
  public Map<Station, Map<JobClass, Metric>> Qt;
  public Map<Station, Map<JobClass, Metric>> Ut;
  public Map<Station, Map<JobClass, Metric>> Tt;

  public static class Metric {

    public String type;
    public Station station;
    public JobClass jobClass;
    public boolean isDisabled;
    public boolean isTransient;
  }
}
