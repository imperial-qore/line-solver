// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.lang;
import jline.api.CTMC;
import jline.lang.distributions.MarkovianDistribution;
import jline.util.Matrix;


import java.util.*;
import static jline.lib.KPCToolbox.*;
import static jline.lib.M3A.*;


/**
 * An environment model defined by a collection of network sub-models coupled with an environment transition rule
 * that selects the active sub-model.
 */
public class Env extends Ensemble {

  public interface ResetQueueLengthsFunction {
    Matrix reset(Matrix input);
  }

  public interface ResetEnvRatesFunction {
    MarkovianDistribution reset(
        MarkovianDistribution originalDist,
        Matrix QExit,
        Matrix UExit,
        Matrix TExit);
  }

  public final MarkovianDistribution[][] env;
  private final String[] names;
  private final String[] types;
  private final Network[] models;

  // Markovian representation of each stage transition
  public Map<Integer, Matrix>[][] proc;
  public Map<Integer, Matrix>[] holdTime; // Holding times
  public Matrix probEnv; // Steady-stage probability of the environment
  public Matrix probOrig; // Probability that a request originated from phase
  public ResetQueueLengthsFunction[][] resetQLFun; // Function implementing the reset policy
  public ResetEnvRatesFunction[][] resetEnvRatesFun;

  @SuppressWarnings("unchecked")
  public Env(String name, int numStages) {
    super(new ArrayList<>());
    this.setName(name);
    this.env = new MarkovianDistribution[numStages][numStages];
    this.names = new String[numStages];
    this.types = new String[numStages];
    this.models = new Network[numStages];
    this.proc = new HashMap[numStages][numStages];
    this.holdTime = new HashMap[numStages];
    this.probEnv = new Matrix(0, 0);
    this.probOrig = new Matrix(0, 0);
    this.resetQLFun = new ResetQueueLengthsFunction[numStages][numStages];
    this.resetEnvRatesFun = new ResetEnvRatesFunction[numStages][numStages];
  }

  public void addStage(int stageIdx, String name, String type, Network model) {
    this.names[stageIdx] = name;
    this.types[stageIdx] = type;
    this.models[stageIdx] = model;
    if (stageIdx > 0 && model.getNumberOfStatefulNodes() != models[0].getNumberOfStatefulNodes()) {
      throw new RuntimeException(
          "Unsupported feature. Random environment stages must map to networks with identical number of stateful nodes.");
    }
    this.ensemble.add(stageIdx, model);
  }

  public void addTransition(int fromStageIdx, int toStageIdx, MarkovianDistribution distrib) {
    this.addTransition(fromStageIdx, toStageIdx, distrib, input -> input);
  }

  public void addTransition(
      int fromStageIdx,
      int toStageIdx,
      MarkovianDistribution distrib,
      ResetQueueLengthsFunction resetFun) {
    this.env[fromStageIdx][toStageIdx] = distrib;
    this.resetQLFun[fromStageIdx][toStageIdx] = resetFun;
  }
  @SuppressWarnings("unchecked")
  public void init() {
    int E = this.models.length;
    Matrix Pemb = new Matrix(E, E);

    // Analyse holding times
    Map<Integer, Matrix>[][] emmap = new HashMap[E][E];
    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        // Multiclass MMAP representation
        if (this.env[e][h] == null) {
          Matrix zero = new Matrix(1, 1);
          Map<Integer, Matrix> zeroMap = new HashMap<>();
          zeroMap.put(0, zero);
          zeroMap.put(1, zero);
          emmap[e][h] = zeroMap;
        } else {
          emmap[e][h] = this.env[e][h].getRepres();
        }
        for (int j = 0; j < E; j++) {
          emmap[e][h].put(j + 2, emmap[e][h].get(1).clone());
          if (j != h) {
            emmap[e][h].get(j + 2).zero();
          }
        }
      }

      holdTime[e] = new HashMap<>();
      for (int i = 0; i < emmap[e][e].size(); i++) {
        holdTime[e].put(i, emmap[e][e].get(i).clone());
      }

      for (int h = 0; h < E; h++) {
        if (h != e) {
          this.holdTime[e].put(0, this.holdTime[e].get(0).krons(emmap[e][h].get(0)));
          for (int j = 1; j < E + 2; j++) {
            this.holdTime[e].put(j, this.holdTime[e].get(j).krons(emmap[e][h].get(j)));
            Matrix ones = new Matrix(holdTime[e].get(j).length(), 1);
            ones.ones();
            Matrix completionRates = holdTime[e].get(j).mult(ones, new Matrix(0, 0));
            holdTime[e].get(j).zero();
            for (int row = 0; row < completionRates.getNumRows(); row++) {
              holdTime[e].get(j).set(row, 0, completionRates.get(row, 0));
            }
          }
          holdTime[e] = mmap_normalize(holdTime[e]);
        }
      }
      // Completion rates for the different transitions
      Matrix countLambda = mmap_count_lambda(holdTime[e]);
      double sumCountLambda = countLambda.sumRows(0);
      for (int col = 0; col < Pemb.getNumCols(); col++) {
        Pemb.set(e, col, countLambda.get(0, col) / sumCountLambda);
      }
    }
    this.proc = emmap;

    Matrix lambda = new Matrix(1, E);
    Matrix A = new Matrix(E, E);
    Matrix I = Matrix.eye(E);
    for (int e = 0; e < E; e++) {
      lambda.set(
          0,
          e,
          1
              / map_mean( holdTime[e].get(0), holdTime[e].get(1)));
      for (int h = 0; h < E; h++) {
        A.set(e, h, -lambda.get(0, e) * (I.get(e, h) - Pemb.get(e, h)));
      }
    }

    int countLambdaValuesLEQZero = 0;
    for (int col = 0; col < lambda.length(); col++) {
      if (lambda.get(0, col) <= 0) {
        countLambdaValuesLEQZero++;
      }
    }
    if (countLambdaValuesLEQZero == 0) {
      this.probEnv = CTMC.ctmc_solve(A);
      this.probOrig = new Matrix(E, E);
      for (int e = 0; e < E; e++) {
        for (int h = 0; h < E; h++) {
          probOrig.set(h, e, probEnv.get(0, h) * lambda.get(0, h) * Pemb.get(h, e));
        }
        if (probEnv.get(0, e) > 0) {
          double probOrigSumCol = probOrig.sumCols(e);
          for (int row = 0; row < E; row++) {
            probOrig.set(row, e, probOrig.get(row, e) / probOrigSumCol);
          }
        }
      }
    }
  }

  public void printStageTable() {
    // TODO: implementation - note return type should likely not be void
    throw new RuntimeException("printStageTable() has not yet been implemented in JLINE.");
  }

  // NOTE: the following LINE methods have not been migrated to JLINE
  // a) getEnv() - env has been made public instead, therefore no need for getter
  // b) setEnv() - env has been made public instead, therefore no need for setter
  // c) setStageName() - appears to be legacy code
  // d) setStageType() - appears to be legacy code
  // e) copyElement - overrides an unimplemented method in "Copyable", and is unused
}
