// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.odes;

import jline.util.Matrix;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.Map;
import java.util.Objects;

import static java.lang.Math.min;
import static jline.lang.constant.SchedStrategy.DPS;
import static jline.lib.KPCToolbox.*;

public class ClosingAndStateDepMethodsODE implements FirstOrderDifferentialEquations {
  private final NetworkStruct sn;
  private final Map<Station, Map<JobClass, Matrix>> mu;
  private final Map<Station, Map<JobClass, Matrix>> phi;
  private final Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc;
  private final Matrix rt;
  private final Matrix nservers;
  private final SolverOptions options;
  private final int numDimensions;

  public ClosingAndStateDepMethodsODE(
      NetworkStruct sn,
      Map<Station, Map<JobClass, Matrix>> mu,
      Map<Station, Map<JobClass, Matrix>> phi,
      Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc,
      Matrix rt,
      Matrix S,
      SolverOptions options,
      int numDimensions) {
    this.sn = sn;
    this.mu = mu;
    this.phi = phi;
    this.proc = proc;
    this.rt = rt;
    this.nservers = S;
    this.options = options;
    this.numDimensions = numDimensions;
  }

  public ClosingAndStateDepMethodsODE(
      NetworkStruct sn,
      Map<Station, Map<JobClass, Matrix>> mu,
      Map<Station, Map<JobClass, Matrix>> phi,
      Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc,
      Matrix rt,
      Matrix S,
      SolverOptions options) {
    this(sn, mu, phi, proc, rt, S, options, options.init_sol.length());
  }

  private void setNextJump(Matrix jumps, int completionIdx, int startIdx) {

    int jumpsCols = jumps.getNumCols();
    jumps.expandMatrix(jumps.getNumRows(), jumpsCols + 1, jumps.getNumElements() + 2);
    jumps.set(completionIdx, jumpsCols, -1); // type c in stat i completes service
    jumps.set(startIdx, jumpsCols, 1); // type c job starts in stat j
  }

  private Matrix calculateJumps(boolean[][] enabled, Matrix qIndices, Matrix Kic) {

    int M = sn.nstations; // Number of stations
    int K = mu.get(sn.stations.get(0)).size(); // Number of classes
    int jumpsRows = (int) Kic.elementSum();
    Matrix jumps =
        new Matrix(jumpsRows, 0); // Returns state changes triggered by all the events

    for (int i = 0; i < M; i++) { // state changes from departures in service phases 2
      for (int c = 0; c < K; c++) {
        if (enabled[i][c]) {
          int xic = (int) qIndices.get(i, c); //  index of x_ic
          for (int j = 0; j < M; j++) {
            for (int l = 0; l < K; l++) {
              if (rt.get(i * K + c, j * K + l) > 0) {
                int xjl = (int) qIndices.get(j, l); // index of x_jl
                for (int ki = 0; ki < Kic.get(i, c); ki++) { // job can leave from any phase in i
                  for (int kj = 0; kj < Kic.get(j, l); kj++) { // job can start from any phase in j
                    setNextJump(jumps, xic + ki, xjl + kj);
                  }
                }
              }
            }
          }
        }
      }
    }

    for (int i = 0; i < M; i++) { // state changes: "next service phase" transition
      for (int c = 0; c < K; c++) {
        if (enabled[i][c]) {
          int xic = (int) qIndices.get(i, c);
          for (int ki = 0; ki < Kic.get(i, c) - 1; ki++) {
            for (int kip = 0; kip < Kic.get(i, c); kip++) {
              if (ki != kip) {
                setNextJump(jumps, xic + ki, xic + kip);
              }
            }
          }
        }
      }
    }

    return jumps;
  }

  private void calculateRateBaseAndEventIdxs(
      boolean[][] enabled,
      Matrix qIndices,
      Matrix Kic,
      Matrix rateBase,
      Matrix eventIdx) {

    int M = sn.nstations; // Number of stations
    int K = mu.get(sn.stations.get(0)).size(); // Number of classes
    int rateIdx = 0;

    // State changes from departures in service phases 2...
    for (int i = 0; i < M; i++) {
      for (int c = 0; c < K; c++) {
        if (enabled[i][c]) {
          for (int j = 0; j < M; j++) {
            for (int l = 0; l < K; l++) {
              Matrix pie;
              if (proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).isEmpty()) {
                pie = new Matrix(1, 1, 1);
                pie.set(0, 0, 1);
              } else {
                Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);

                pie = map_pie(D0,D1);
              }
              if (rt.get(i * K + c, j * K + l) > 0) {
                for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                  for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                    rateBase.set(
                        rateIdx,
                        0,
                        phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                            * mu.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                            * rt.get(i * K + c, j * K + l)
                            * pie.get(0, kjl));
                    eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                    rateIdx++;
                  }
                }
              }
            }
          }
        }
      }
    }

    // State changes from "next service phase" transition in phases 2...
    for (int i = 0; i < M; i++) {
      for (int c = 0; c < K; c++) {
        if (enabled[i][c]) {
          for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
            for (int kicp = 0; kicp < Kic.get(i, c); kicp++) {
              if (kicp != kicIdx) {
                rateBase.set(
                    rateIdx,
                    0,
                    proc.get(sn.stations.get(i))
                        .get(sn.jobclasses.get(c))
                        .get(0)
                        .get(kicIdx, kicp));
                eventIdx.set(rateIdx, 0, qIndices.get(i, c) + kicIdx);
                rateIdx++;
              }
            }
          }
        }
      }
    }
  }

  private Matrix calculatedxdtClosingMethod(
      double[] x,
      Matrix w,
      boolean[][] enabled,
      Matrix qIndices,
      Matrix Kic,
      Matrix allJumps,
      Matrix rateBase,
      Matrix eventIdx) {

    int M = sn.nstations; // Number of stations
    int K = mu.get(sn.stations.get(0)).size(); // Number of classes

    // Basic vector valid for INF and PS case min(ni, nservers(i)) = ni
    Matrix rates = new Matrix(x.length, 1);
    for (int i = 0; i < x.length; i++) {
      rates.set(i, 0, x[i]);
    }

    for (int i = 0; i < M; i++) {
      switch (sn.sched.get(sn.stations.get(i))) {
        case INF:
          break;
        case EXT:
          // This is treated by a delay except that we require mass conservation in the local
          // population
          for (int k = 0; k < K; k++) {
            int idxIni = (int) qIndices.get(i, k);
            int idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
            if (enabled[i][k]) {
              // Keep total mass 1 into the source for all classes at all
              // times, not needed for idxIni+1:idxEnd as rates is initialized equal to x
              double tmpSum = 0;
              for (int idx = idxIni + 1; idx < idxEnd; idx++) {
                tmpSum += x[idx];
              }
              rates.set(idxIni, 0, 1 - tmpSum);
            }
          }
          break;
        case PS:
        case FCFS:
          int idxIni = (int) qIndices.get(i, 0);
          int idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
          double ni = 0;
          for (int idx = idxIni; idx < idxEnd; idx++) {
            ni += x[idx];
          }
          if (ni > nservers.get(i, 0)) { // case min = ni handled by rates = x
            for (int idx = idxIni; idx < idxEnd; idx++) {
              rates.set(idx, 0, x[idx] / ni * nservers.get(i, 0));
            }
          }
          break;
        case DPS:
          double sumWI = w.sumRows(i);
          for (int col = 0; col < K; col++) {
            w.set(i, col, w.get(i, col) / sumWI);
          }
          sumWI = w.sumRows(i);
          ni = sumWI / K;

          for (int k = 0; k < K; k++) {
            idxIni = (int) qIndices.get(i, k);
            idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
            if (enabled[i][k]) {
              double tmpSum = 0;
              for (int idx = idxIni; idx < idxEnd; idx++) {
                tmpSum += x[idx] * w.get(i, k);
              }
              ni += tmpSum;
            }
          }

          for (int k = 0; k < K; k++) {
            idxIni = (int) qIndices.get(i, k);
            idxEnd = (int) qIndices.get(i, k) + (int) Kic.get(i, k);
            if (enabled[i][k]) {
              for (int idx = idxIni; idx < idxEnd; idx++) {
                // Not needed for idxIni+1:idxEnd as rates is initialised equal to x
                rates.set(idx, 0, w.get(i, k) * x[idx] / ni * nservers.get(i, 0));
              }
            }
          }
      }
    }

    int numEventIndices = eventIdx.length();
    Matrix newRates = new Matrix(numEventIndices, 1);
    for (int i = 0; i < numEventIndices; i++) {
      newRates.set(i, 0, rates.get((int) eventIdx.get(i, 0), 0));
    }
    newRates.elementMult(rateBase, newRates);
    return allJumps.mult(newRates, null);
  }

  private Matrix calculatedxdtStateDepMethod(
          double[] x, boolean[][] enabled, Matrix qIndices, Matrix Kic, Matrix w) {

    int M = sn.nstations; // Number of stations
    int K = mu.get(sn.stations.get(0)).size(); // Number of classes
    Matrix dxdt = new Matrix(x.length, 1);

    for (int i = 0; i < M; i++) {
      switch (sn.sched.get(sn.stations.get(i))) {
        case INF:
          // Phase changes
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                  if (kicIdx != kic_p) {
                    double rate =
                        proc.get(sn.stations.get(i))
                            .get(sn.jobclasses.get(c))
                            .get(0)
                            .get(kicIdx, kic_p);
                    dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                    dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                  }
                }
              }
            }
          }
          // Service completions
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int j = 0; j < M; j++) {
                for (int l = 0; l < K; l++) {
                  int xjl = (int) qIndices.get(j, l);
                  if (enabled[j][l]) {
                    Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                    Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                    Matrix pie =
                        map_pie(D0,D1);
                    if (rt.get(i * K + c, j * K + l) > 0) {
                      for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                        for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                          if (j != i) {
                            double rate =
                                phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                    * mu.get(sn.stations.get(i))
                                        .get(sn.jobclasses.get(c))
                                        .get(kicIdx, 0)
                                    * rt.get(i * K + c, j * K + l)
                                    * pie.get(0, kjl);
                            dxdt.set(
                                xic + kicIdx,
                                0,
                                dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                            dxdt.set(
                                xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          break;

        case EXT:
          // TODO
          System.err.println(
              "State dependent ODE method does not support open models. Try with default method.");
          break;

        case PS:
          int idxIni = (int) qIndices.get(i, 0);
          int idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
          double ni = 0;
          for (int idx = idxIni; idx < idxEnd; idx++) {
            ni += x[idx];
          }
          // Phase changes
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                  if (kicIdx != kic_p) {
                    double rate =
                        proc.get(sn.stations.get(i))
                            .get(sn.jobclasses.get(c))
                            .get(0)
                            .get(kicIdx, kic_p);
                    if (ni > sn.nservers.get(i, 0)) {
                      dxdt.set(
                          xic + kicIdx,
                          0,
                          dxdt.get(xic + kicIdx, 0)
                              - (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                      dxdt.set(
                          xic + kic_p,
                          0,
                          dxdt.get(xic + kic_p, 0)
                              + (x[xic + kicIdx] * rate * sn.nservers.get(i, 0) / ni));
                    } else {
                      dxdt.set(
                          xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                      dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                    }
                  }
                }
              }
            }
          }
          // Service completions
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int j = 0; j < M; j++) {
                for (int l = 0; l < K; l++) {
                  int xjl = (int) qIndices.get(j, l);
                  if (enabled[j][l]) {
                    Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                    Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                    Matrix pie = map_pie(D0,D1);

                    if (rt.get(i * K + c, j * K + l) > 0) {
                      for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                        for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                          double rate =
                              phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                  * mu.get(sn.stations.get(i))
                                      .get(sn.jobclasses.get(c))
                                      .get(kicIdx, 0)
                                  * rt.get(i * K + c, j * K + l)
                                  * pie.get(0, kjl);
                          if (ni > sn.nservers.get(i, 0)) {
                            rate *= sn.nservers.get(i, 0) / ni;
                          }
                          dxdt.set(
                              xic + kicIdx,
                              0,
                              dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                          dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          break;

        case FCFS:
          idxIni = (int) qIndices.get(i, 0);
          idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
          ni = 0;
          for (int idx = idxIni; idx < idxEnd; idx++) {
            ni += x[idx];
          }
          double wni = 0.001;
          for (int c = 0; c < K; c++) {
            for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
              if (enabled[i][c]) {
                int xic = (int) qIndices.get(i, c);
                w.set(
                    c,
                    kicIdx,
                    -1
                        / proc.get(sn.stations.get(i))
                            .get(sn.jobclasses.get(c))
                            .get(0)
                            .get(kicIdx, kicIdx));
                wni += w.get(c, kicIdx) * x[xic + kicIdx];
              }
            }
          }
          // Phase changes
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                  if (kicIdx != kic_p) {
                    double rate =
                        proc.get(sn.stations.get(i))
                                .get(sn.jobclasses.get(c))
                                .get(0)
                                .get(kicIdx, kic_p)
                            * min(ni, sn.nservers.get(i, 0))
                            * w.get(c, kicIdx)
                            / wni;
                    dxdt.set(xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                    dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                  }
                }
              }
            }
          }
          // Service completions
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int j = 0; j < M; j++) {
                for (int l = 0; l < K; l++) {
                  int xjl = (int) qIndices.get(j, l);
                  if (enabled[j][l]) {
                    Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                    Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                    Matrix pie = map_pie(D0,D1);

                    if (rt.get(i * K + c, j * K + l) > 0) {
                      for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                        for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                          double rate =
                              phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                  * mu.get(sn.stations.get(i))
                                      .get(sn.jobclasses.get(c))
                                      .get(kicIdx, 0)
                                  * rt.get(i * K + c, j * K + l)
                                  * pie.get(0, kjl)
                                  * min(ni, sn.nservers.get(i, 0))
                                  * w.get(c, kicIdx)
                                  / wni;
                          dxdt.set(
                              xic + kicIdx,
                              0,
                              dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                          dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          break;

        case DPS:
          double sumWI = w.sumRows(i);
          for (int col = 0; col < w.getNumCols(); col++) {
            w.set(i, col, w.get(i, col) / sumWI);
          }
          idxIni = (int) qIndices.get(i, 0);
          idxEnd = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
          wni = 0;
          for (int idx = idxIni; idx < idxEnd; idx++) {
            wni += x[idx];
          }
          // Phase changes
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int kicIdx = 0; kicIdx < Kic.get(i, c) - 1; kicIdx++) {
                for (int kic_p = 0; kic_p < Kic.get(i, c); kic_p++) {
                  if (kicIdx != kic_p) {
                    double rate =
                        proc.get(sn.stations.get(i))
                            .get(sn.jobclasses.get(c))
                            .get(0)
                            .get(kicIdx, kic_p);
                    if (wni > sn.nservers.get(i, 0)) {
                      dxdt.set(
                          xic + kicIdx,
                          0,
                          dxdt.get(xic + kicIdx, 0)
                              - (x[xic + kicIdx]
                                  * rate
                                  * sn.nservers.get(i, 0)
                                  * w.get(c, kicIdx)
                                  / wni));
                      dxdt.set(
                          xic + kic_p,
                          0,
                          dxdt.get(xic + kic_p, 0)
                              + (x[xic + kicIdx]
                                  * rate
                                  * sn.nservers.get(i, 0)
                                  * w.get(c, kicIdx)
                                  / wni));
                    } else {
                      dxdt.set(
                          xic + kicIdx, 0, dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                      dxdt.set(xic + kic_p, 0, dxdt.get(xic + kic_p, 0) + (x[xic + kicIdx] * rate));
                    }
                  }
                }
              }
            }
          }
          // Service completions
          for (int c = 0; c < K; c++) {
            if (enabled[i][c]) {
              int xic = (int) qIndices.get(i, c);
              for (int j = 0; j < M; j++) {
                for (int l = 0; l < K; l++) {
                  int xjl = (int) qIndices.get(j, l);
                  if (enabled[j][l]) {
                    Matrix D0 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(0);
                    Matrix D1 = proc.get(sn.stations.get(j)).get(sn.jobclasses.get(l)).get(1);
                    Matrix pie = map_pie(D0,D1);

                    if (rt.get(i * K + c, j * K + l) > 0) {
                      for (int kicIdx = 0; kicIdx < Kic.get(i, c); kicIdx++) {
                        for (int kjl = 0; kjl < Kic.get(j, l); kjl++) {
                          double rate =
                              phi.get(sn.stations.get(i)).get(sn.jobclasses.get(c)).get(kicIdx, 0)
                                  * mu.get(sn.stations.get(i))
                                      .get(sn.jobclasses.get(c))
                                      .get(kicIdx, 0)
                                  * rt.get(i * K + c, j * K + l)
                                  * pie.get(0, kjl);
                          if (wni > sn.nservers.get(i, 0)) {
                            rate *= sn.nservers.get(i, 0) * w.get(c, kicIdx) / wni;
                          }
                          dxdt.set(
                              xic + kicIdx,
                              0,
                              dxdt.get(xic + kicIdx, 0) - (x[xic + kicIdx] * rate));
                          dxdt.set(xjl + kjl, 0, dxdt.get(xjl + kjl, 0) + (x[xic + kicIdx] * rate));
                        }
                      }
                    }
                  }
                }
              }
            }
          }
      }
    }

    return dxdt;
  }

  @Override
  public int getDimension() {
    return numDimensions;
  }

  @Override
  public void computeDerivatives(double t, double[] x, double[] dxdt)
      throws MaxCountExceededException, DimensionMismatchException {

    int M = sn.nservers.length(); // Number of stations
    int K = mu.get(sn.stations.get(0)).size(); // Number of classes

    Matrix w = new Matrix(M, K);
    boolean[][] enabled = new boolean[M][K]; // Indicates whether a class is served at a station
    for (int i = 0; i < M; i++) {
      for (int k = 0; k < K; k++) {
        w.set(i, k, 1);
        enabled[i][k] = false;
      }
    }
    Matrix qIndices = new Matrix(M, K);
    Matrix Kic = new Matrix(M, K);
    int cumSum = 0;

    for (int i = 0; i < M; i++) {
      Station station = sn.stations.get(i);
      for (int c = 0; c < K; c++) {
        JobClass jobClass = sn.jobclasses.get(c);
        int numPhases = 0;
        int numNans = 0;
        int rows = mu.get(station).get(jobClass).getNumRows();
        int cols = mu.get(station).get(jobClass).getNumCols();
        for (int row = 0; row < rows; row++) {
          for (int col = 0; col < cols; col++) {
            if (Double.isNaN(mu.get(station).get(jobClass).get(row, col))) {
              numNans++;
            }
          }
        }
        if ((numNans != mu.get(station).get(jobClass).getNumElements())
            && !mu.get(station).get(jobClass).isEmpty()) {
          numPhases = mu.get(station).get(jobClass).length();
          enabled[i][c] = true;
        }
        qIndices.set(i, c, cumSum);
        Kic.set(i, c, numPhases);
        cumSum += numPhases;
      }

      if (sn.sched.get(station) == DPS) {
        for (int k = 0; k < K; k++) {
          w.set(i, k, sn.schedparam.get(i, k));
        }
      }
    }

    Matrix dxdtTmp;
    if (Objects.equals(options.method, "statedep")) {
      dxdtTmp = calculatedxdtStateDepMethod(x, enabled, qIndices, Kic, w);
    } else {
      // Determine all the jumps and save them for later use
      Matrix allJumps = calculateJumps(enabled, qIndices, Kic);
      // Determines a vector with the fixed part of the rates and defines the indexes that
      // correspond to the events that occur
      Matrix rateBase = new Matrix(allJumps.getNumCols(), 1);
      Matrix eventIdx = new Matrix(allJumps.getNumCols(), 1);
      calculateRateBaseAndEventIdxs(enabled, qIndices, Kic, rateBase, eventIdx);
      dxdtTmp =
          calculatedxdtClosingMethod(x, w, enabled, qIndices, Kic, allJumps, rateBase, eventIdx);
    }

    for (int i = 0; i < dxdt.length; i++) {
      dxdt[i] = dxdtTmp.get(i);
    }
  }
}
