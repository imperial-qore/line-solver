package jline.lang.state;

import java.io.Serializable;
import java.util.*;

import jline.examples.ClosedModel;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.*;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.Matrix;
import jline.lang.NetworkStruct;
import jline.util.SerializableFunction;

import static jline.lang.constant.SchedStrategy.*;

/**
 * Class modeling the state of Stateful nodes
 */
public class State implements Serializable {

  /*
     The state of the network is described by:
         1. For each node.. there is a :
             a. Initial State
             b. Prior State

  */

  public final Map<StatefulNode, Matrix> initialState;
  public final Map<StatefulNode, Matrix> priorInitialState;

  public State(
          Map<StatefulNode, Matrix> initialState, Map<StatefulNode, Matrix> priorInitialState) {
    this.initialState = initialState;
    this.priorInitialState = priorInitialState;
  }



  public static class StateMarginalStatistics {

    public Matrix ni;
    public Matrix nir;
    public Matrix sir;
    public List<Matrix> kir;

    public StateMarginalStatistics(
            Matrix ni, Matrix nir, Matrix sir, List<Matrix> kir) {
      this.ni = ni;
      this.nir = nir;
      this.sir = sir;
      this.kir = kir;
    }
  }

  public static StateMarginalStatistics toMarginal(
      NetworkStruct sn,
      int ind,
      Matrix state_i,
      Matrix phasesz,
      Matrix phaseshift,
      Matrix space_buf,
      Matrix space_srv,
      Matrix space_var) {
    // Cached node, currently not support in JLINE
    if (sn.isstation.get(ind, 0) == 0 && sn.isstateful.get(ind, 0) > 0) {
      throw new RuntimeException("Not implemented"); // TODO: not implemented
    }

    int R = sn.nclasses;
    int ist = (int) sn.nodeToStation.get(0, ind);

    if (phasesz == null) {
      phasesz = new Matrix(1, sn.phasessz.getNumCols());
      Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), phasesz, 0, 0);
    }
    if (phaseshift == null) {
      phaseshift = new Matrix(1, sn.phaseshift.getNumCols());
      Matrix.extract(
          sn.phaseshift, ist, ist + 1, 0, sn.phaseshift.getNumCols(), phaseshift, 0, 0);
    }

    boolean isExponential = (phasesz.elementMax() == 1);

    if (space_var == null) {
      int col = (int) sn.nvars.sumRows(ind);
      space_var = new Matrix(state_i.getNumRows(), col);
      Matrix.extract(
          state_i,
          0,
          state_i.getNumRows(),
          state_i.getNumCols() - col,
          state_i.getNumCols(),
          space_var,
          0,
          0);
    }
    if (space_srv == null) {
      int sumPhasesz = (int) phasesz.elementSum();
      int sumNvars = (int) sn.nvars.sumRows(ind);
      space_srv = new Matrix(state_i.getNumRows(), sumPhasesz);
      Matrix.extract(
          state_i,
          0,
          state_i.getNumRows(),
          state_i.getNumCols() - sumPhasesz - sumNvars,
          state_i.getNumCols() - sumNvars,
          space_srv,
          0,
          0);
    }
    if (space_buf == null) {
      int col =
          state_i.getNumCols() - (int) (phasesz.elementSum() + sn.nvars.sumRows(ind));
      space_buf = new Matrix(state_i.getNumRows(), col);
      Matrix.extract(state_i, 0, state_i.getNumRows(), 0, col, space_buf, 0, 0);
    }

    Matrix nir = new Matrix(state_i.getNumRows(), R);
    Matrix sir = new Matrix(state_i.getNumRows(), R);
    List<Matrix> kir = new ArrayList<>();
    if (isExponential) {
      sir = space_srv;
      kir.add(space_srv);
    } else {
      // Initialize kir
      for (int i = 0; i < phasesz.elementMax(); i++)
        kir.add(new Matrix(state_i.getNumRows(), R));

      for (int r = 0; r < R; r++) {
        for (int k = 0; k < phasesz.get(r); k++) {
          // kir(:,r,k) = space_srv(:,phaseshift(r)+k);
          Matrix tmp_kir = kir.get(k);
          Matrix.extract(
              space_srv,
              0,
              space_srv.getNumRows(),
              (int) phaseshift.get(r) + k,
              (int) phaseshift.get(r) + k + 1,
              tmp_kir,
              0,
              r);

          // sir(:,r) = sir(:,r) + kir(:,r,k);
          for (int i = 0; i < sir.getNumRows(); i++)
            sir.set(i, r, sir.get(i, r) + tmp_kir.get(i, r));
        }
      }
    }

    switch (sn.sched.get(sn.stations.get(ist))) {
      case INF:
      case PS:
      case DPS:
      case GPS:
        nir = sir.clone();
        break;
      case EXT:
        nir.fill(Double.POSITIVE_INFINITY);
        break;
      case FCFS:
      case HOL:
      case LCFS:
        for (int r = 0; r < R; r++) {
          Matrix sumval;
          if (space_buf.getNumRows() == 1
              && space_buf.getNumCols() == 1
              && space_buf.get(0, 0) == 0) {
            sumval = new Matrix(1, 1);
            sumval.set(0, 0, 0);
          } else {
            // +1 since loop starts from 0 but classes start from 1
            sumval = space_buf.countEachRow(r+1);
          }
          for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r) + sumval.get(i));
        }
        break;
      case LCFSPR:
        if (space_buf.length() > 1) {
          // space_buf = space_buf(1:2:end);
          Matrix sub_space_buf =
              new Matrix(1, (space_buf.getNumCols() * space_buf.getNumRows() + 1) / 2);
          for (int i = 0; i < sub_space_buf.getNumCols(); i++)
            sub_space_buf.set(0, i, space_buf.get(i * 2));

          for (int r = 0; r < R; r++) {
            Matrix sumval = sub_space_buf.countEachRow(r);
            for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r) + sumval.get(i));
          }
        } else {
          nir = sir.clone();
        }
        break;
      case SIRO:
      case SEPT:
      case LEPT:
        for (int r = 0; r < R; r++) {
          for (int i = 0; i < sir.getNumRows(); i++)
            nir.set(i, r, sir.get(i, r) + space_buf.get(i, r));
        }
        break;
      default:
        for (int r = 0; r < R; r++) {
          for (int i = 0; i < sir.getNumRows(); i++) nir.set(i, r, sir.get(i, r));
        }
    }

    if (sn.nodetypes.get(ind) != NodeType.Place) {
      for (int r = 0; r < R; r++) {
        if (Double.isNaN(sn.rates.get(ist, r))) {
          for (int i = 0; i < nir.getNumRows(); i++) {
            nir.remove(i, r);
            sir.remove(i, r);
          }
          for (int k = 0; k < phasesz.get(r); k++) {
            for (int j = 0; j < kir.size(); j++) {
              kir.get(j).remove(r, k);
            }
          }
        }
      }
    }

    Matrix ni = nir.sumRows();

    return new StateMarginalStatistics(ni, nir, sir, kir);
  }

public static StateMarginalStatistics toMarginalAggr(NetworkStruct sn,
                                              int ind,
                                              Matrix state_i,
                                              Matrix K,
                                              Matrix Ks,
                                              Matrix space_buf,
                                              Matrix space_srv,
                                              Matrix space_var) {

  int ist = (int) sn.nodeToStation.get(ind);
  int R = sn.nclasses;

  if (sn.isstation.get(ind, 0) == 0 && sn.isstateful.get(ind, 0) > 0) {
    throw new RuntimeException("Not implemented"); // TODO: not implemented
  }


  if (K == null) {
    K = new Matrix(1, sn.phasessz.getNumCols());
    Matrix.extract(sn.phasessz, ist, ist + 1, 0, sn.phasessz.getNumCols(), K, 0, 0);
  }
  if (Ks == null) {
    Ks = new Matrix(1, sn.phaseshift.getNumCols());
    Matrix.extract(
            sn.phaseshift, ist, ist + 1, 0, sn.phaseshift.getNumCols(), Ks, 0, 0);
  }

  if (space_var == null) {
    int col = (int) sn.nvars.sumRows(ind);
    space_var = new Matrix(state_i.getNumRows(), col);
    Matrix.extract(
            state_i,
            0,
            state_i.getNumRows(),
            state_i.getNumCols() - col,
            state_i.getNumCols(),
            space_var,
            0,
            0);
  }

  if (space_srv == null) {
    int sumK = (int) K.elementSum();
    int sumNvars = (int) sn.nvars.sumRows(ind);
    int numCols = sumNvars - sumK;
    space_srv = new Matrix(state_i.getNumRows(), numCols);
    Matrix.extract(
            state_i,
            0,
            state_i.getNumRows(),
            state_i.getNumCols() - sumK,
            state_i.getNumCols() - sumNvars,
            space_srv,
            0,
            0);
  }

  if (space_buf == null) {
    int col =
            state_i.getNumCols() - (int) (K.elementSum());
    space_buf = new Matrix(state_i.getNumRows(), col);
    Matrix.extract(state_i, 0, state_i.getNumRows(), 0, col, space_buf, 0, 0);
  }
  Matrix nir = new Matrix(state_i.getNumRows(), R);
  nir.zero();
  for (int r = 0; r < R; r++) {
    for (int k = 0; k < K.get(r); k++) {
      Matrix tmp_kir = new Matrix(state_i.getNumRows(), 1);
      Matrix.extract(
              space_srv,
              0,
              space_srv.getNumRows(),
              (int) Ks.get(r) + k,
              (int) Ks.get(r) + k + 1,
              tmp_kir,
              0,
              0);
      for (int i = 0; i < nir.getNumRows(); i++) {
        nir.set(i, r, nir.get(i, r) + tmp_kir.get(i));
      }
    }
  }

  // MATLAB LINE does not handle LCFSPR, INF, PS, DPS, GPS cases
  // nir: class-r jobs at the station
  switch (sn.sched.get(sn.stations.get(ist))) {
    case EXT:
      nir.fill(Double.POSITIVE_INFINITY);
      break;
    case FCFS:
    case HOL:
    case LCFS:
      for (int r = 0; r < R; r++) {
        Matrix sumval;
        if (space_buf.getNumRows() == 1
                && space_buf.getNumCols() == 1
                && space_buf.get(0, 0) == 0) {
          sumval = new Matrix(1, 1);
          sumval.set(0, 0, 0);
        } else {
          // +1 since loop starts from 0 but classes start from 1
          sumval = space_buf.countEachRow(r + 1);
        }
        for (int i = 0; i < nir.getNumRows(); i++) nir.set(i, r, nir.get(i, r) + sumval.get(i));
      }
      break;
    case SIRO:
    case SEPT:
    case LEPT:
      for (int r = 0; r < R; r++) {
        for (int i = 0; i < nir.getNumRows(); i++)
          nir.set(i, r, nir.get(i, r) + space_buf.get(i, r));
      }
      break;
    default:
      // do nothing, no-op case
  }

  for (int r = 0; r < R; r++) {
    if (Double.isNaN(sn.rates.get(ist, r))) { // if station disabled
      for (int i = 0; i < nir.getNumRows(); i++) {
        nir.remove(i, r);
      }
    }
  }

  Matrix ni = nir.sumRows();
  return new StateMarginalStatistics(ni, nir, null, null);

}


  public static Matrix fromMarginalAndRunning(NetworkStruct sn, int ind, Matrix n, Matrix s) {
    return fromMarginalAndRunning(sn, ind, n, s, true);
  }

  public static Matrix fromMarginalAndRunning(Network sn, int ind, Matrix n, Matrix s) {
    return fromMarginalAndRunning(sn.getStruct(true), ind, n, s, true);
  }

  public static Matrix fromMarginalAndRunning(NetworkStruct sn, int ind, Matrix n, Matrix s, boolean optionsForce) {
    int ist = (int) sn.nodeToStation.get(ind);
    int isf = (int) sn.nodeToStateful.get(ind);

    // generate one initial state such that the marginal queue-lengths are as in vector n
    // n(r): number of jobs at the station in class r
    // s(r): jobs of class r that are running
    int R = sn.nclasses;
    Matrix S = sn.nservers;
    Matrix K = new Matrix(1, R);

    for (int r = 0; r < R; r++) {
      if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
        K.set(0, r, 0);
      } else {
        K.set(0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
      }
    }
    Matrix state = new Matrix(0, 0);
    Matrix space = new Matrix(0, 0);
    LinkedList<Integer> exceeded = new LinkedList<>();
    for (int i = 0; i < sn.classcap.getNumCols(); i++) {
      if (n.get(0, i) > sn.classcap.get(ist, i)) {
        exceeded.add(i);
      }
    }
    if (!exceeded.isEmpty()) {
      for (Integer r : exceeded) {
        if (!sn.proc.isEmpty()
                && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()
                && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
          System.err.format(
                  "State vector at station %d exceeds the class capacity. Some service classes are disabled.\n",
                  ist);
        } else {
          System.err.format("State vector at station %d exceeds the class capacity.\n", ist);
        }
      }
      return space;
    }

    if (sn.nservers.get(ist, 0) > 0
            && s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols()) > sn.nservers.get(ist, 0)) {
      return space;
    }

    if ((sn.nodetypes.get(ind) == NodeType.Queue)
            || (sn.nodetypes.get(ind) == NodeType.Delay)
            || (sn.nodetypes.get(ind) == NodeType.Source)) {
      switch (sn.sched.get(sn.stations.get(ist))) {
        case EXT:
            // TODO: complete
          break;

        case INF:
        case PS:
        case DPS:
        case GPS:
          // in these policies we only track the jobs in the servers
          for (int r = 0; r < R; r++) {
            Matrix init = spaceClosedSingle(K.get(r), n.get(r));
            state = Matrix.decorate(state, init);
          }
          space = Matrix.decorate(space, state);
          break;
        case SIRO:
        case LEPT:
        case SEPT:
          // in these policies we track an un-ordered buffer and the jobs in the servers
          // we build a list of job classes in the node with repetition
          if (n.elementSum() <= S.get(ist)) {
            for (int r = 0; r < R; r++) {
              Matrix init = spaceClosedSingle(K.get(r), n.get(r));
              state = Matrix.decorate(state, init);
            }
            Matrix newStates = new Matrix(state.getNumRows(), R);
            newStates.zero();
            newStates = newStates.concatCols(state);
            space = Matrix.decorate(space, newStates);
          }  else {
            Matrix si = s.clone();
            Matrix mi_buf = n.repmat(si.getNumRows(), 1).sub(1, si); // jobs of class r in buffer
            for (int k = 0; k < si.getNumRows(); k++) {
              Matrix kstate = new Matrix(0, 0);
              for (int r = 0; r < R; r++) {
                Matrix init = spaceClosedSingle(K.get(r), si.get(k, r));
                kstate = Matrix.decorate(kstate, init);
              }
              state = Matrix.extractRows(mi_buf, k, k + 1, null).repmat(kstate.getNumRows(), 1).concatCols(kstate);
              if (space.isEmpty()) {
                space = state.clone();
              } else {
                space = Matrix.concatRows(space, state, null);
              }
            }
          }
          break;
        case FCFS:
        case HOL:
        case LCFS:
          double sizeEstimator = Maths.multinomialln(n.sub(1, s));
          sizeEstimator = Math.round(sizeEstimator/Math.log(10));
          if (sizeEstimator > 2) {
            if (!optionsForce) {
              System.err.format("State space size is very large: 1e%d states. Stopping execution. " +
                      "Set options.force=true to bypass this control.\n",Math.round(sizeEstimator/Math.log(10)));
            }
          }
          if (n.elementSum() == 0) {
            space = new Matrix(1, (int) (1 + K.elementSum()));
            space.zero();
            return space;
          }
          // in these policies we track an ordered buffer and
          // the jobs in the servers

          // build list of job classes in the buffer, with repetition
          Matrix inbuf = new Matrix(0,0);
          for (int r = 0; r < R; r++) {
            if (n.get(0, r) > s.get(r)) {
              int numNewCols = (int) (n.get(0,r) - s.get(0,r));
              Matrix newInBuf = new Matrix(1, inbuf.getNumCols() + numNewCols);
              for (int i = 0; i < inbuf.getNumCols(); i++) {
                newInBuf.set(0,i,inbuf.get(0,i));
              }
              for (int i = inbuf.getNumCols(); i < newInBuf.getNumCols(); i++) {
                newInBuf.set(0,i,r+1);
              }
              inbuf = newInBuf.clone();
            }
          }

          // gen permutation of their positions in the fcfs buffer
          Matrix mi = Maths.uniquePerms(inbuf);
          if (mi.isEmpty()) {
            Matrix mi_buf = new Matrix(1, (int) Maths.max(0, n.elementSum() - S.get(ist)));
            state = new Matrix(1, R);
            state.zero();
            state = Matrix.decorate(state, mi_buf.concatCols(state));
          } else {
            // mi_buf: class of job in buffer position i (0=empty)
            Matrix mi_buf = new Matrix(0,0);
            double sumN = n.elementSum();
            double sums = s.elementSum();
            if (sumN > sums) {
              mi_buf = new Matrix(mi.getNumRows(), (int) sumN - (int) sums);
              for (int row = 0; row < mi_buf.getNumRows(); row++) {
                for (int col = 0; col < sumN - sums; col++) {
                  mi_buf.set(row, col, mi.get(row, col));
                }
              }
            } else {
              mi_buf = new Matrix(1,1);
              mi_buf.set(0,0,0);
            }

            // si: number of class r jobs that are running
            Matrix si = s.clone();
            for (int b = 0; b < mi_buf.getNumRows(); b++) {
              for (int k = 0; k < si.getNumRows(); k++) {
                Matrix kstate = new Matrix(0,0);
                for (int r = 0; r < R; r++) {
                  Matrix init = spaceClosedSingle(K.get(r), si.get(k,r));
                  kstate = Matrix.decorate(kstate, init);
                }
                Matrix miBufRep = Matrix.extractRows(mi_buf,b,b+1,null).repmat(kstate.getNumRows(),1);
                miBufRep = miBufRep.concatCols(kstate);
                if (state.isEmpty()) {
                  state = miBufRep;
                } else {
                  state = Matrix.concatRows(state, miBufRep, null);
                }
              }
            }
          }
          space = state;
          break;

        case LCFSPR:
          double sizeEstimatorLPR = Maths.multinomialln(n.sub(1,s));
          sizeEstimatorLPR = Math.round(sizeEstimatorLPR/Math.log(10));
          if (sizeEstimatorLPR > 2) {
            if (!optionsForce) {
              System.err.format("State space size is very large: 1e%d states. Stopping execution. Set options = true," +
                      "to bypass this control.\n", Math.round(sizeEstimatorLPR/Math.log(10)));
            }
          }
          if (n.elementSum() == 0) {
            space = new Matrix(1, (int) (1+K.elementSum()));
            space.zero();
            return space;
          }
          // in these policies we track an ordered buffer and the jobs in the servers

          // build list of job classes in the buffer with repetition
          Matrix inbufLpr = new Matrix(0,0);
          for (int r = 0; r < R; r++) {
            if (n.get(r) > s.get(r)) {
              int numNewCols = (int) (n.get(r) - s.get(r));
              Matrix newInBuf = new Matrix(1, inbufLpr.getNumCols() + numNewCols);
              for (int i = 0; i < inbufLpr.getNumCols(); i++) {
                newInBuf.set(0,i,inbufLpr.get(0,i));
              }
              for (int i = inbufLpr.getNumCols(); i < newInBuf.getNumCols(); i++) {
                newInBuf.set(0,i,r+1);
              }
              inbufLpr = newInBuf.clone();
            }
          }

          // gen permutation of their positions in the FCFS buffer
          Matrix miLpr = Maths.uniquePerms(inbufLpr);
          if (miLpr.isEmpty()) {
            Matrix mi_buf = new Matrix(1, (int) Math.max(0, n.elementSum() - S.get(ist)));
            mi_buf.zero();
            state = new Matrix(1, R);
            state.zero();
            state = Matrix.decorate(state, mi_buf.concatCols(state));
          } else {
            // mi_buf: class of job in buffer position i (0=empty)
            Matrix mi_buf = new Matrix(0,0);
            if (n.elementSum() > s.elementSum()) {
              double sumN = n.elementSum();
              double sums = s.elementSum();
              mi_buf = new Matrix(miLpr.getNumRows(), (int) sumN - (int) sums);
              for (int row = 0; row < mi_buf.getNumRows(); row++) {
                for (int col = 0; col < sumN - sums; col++) {
                  mi_buf.set(row, col, miLpr.get(row, col));
                }
              }
            } else {
              mi_buf = new Matrix(1,1);
              mi_buf.set(0,0,0);
            }

            // si: number of class r jobs that are running
            Matrix si = s.clone();
            for (int b = 0; b < mi_buf.getNumRows(); b++) {
              for (int k = 0; k < si.getNumRows(); k++) {
                Matrix kState = new Matrix(0,0);
                for (int r = 0; r < R; r++) {
                  Matrix init = spaceClosedSingle(K.get(r), si.get(k,r));
                  kState = Matrix.decorate(kState, init);
                }
                Matrix bkState = new Matrix(0,0);
                Matrix jobsInBuffer = Matrix.extractRows(mi_buf, b, b+1, null);
                for (int j = 0; j < jobsInBuffer.length(); j++) {
                  double job = jobsInBuffer.get(j);
                  if (job > 0) {
                    List<Double> phasesJRange = new ArrayList<>();
                    for (double i = 1; i <= K.get((int) job - 1); i++) {
                      phasesJRange.add(i);
                    }
                    // no transpose as constructor makes column vector
                    bkState = Matrix.decorate(bkState, new Matrix(phasesJRange));
                  } else {
                    bkState = new Matrix(1,1);
                    bkState.zero();
                  }
                }
                Matrix bufStateTmp = Matrix.decorate(Matrix.extractRows(mi_buf, b, b + 1, null), bkState);
                // here we interleave positions of class and phases in buffer
                Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                bufState.zero();

                // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                int colForBufStateTmp = 0;
                for (int row = 0; row < bufState.getNumRows(); row++) {
                  for (int col = 0; col < bufState.getNumCols(); col += 2) {
                    if (colForBufStateTmp < mi_buf.getNumCols()) {
                      bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                      colForBufStateTmp++;
                    }
                  }
                  colForBufStateTmp = 0;
                }
                colForBufStateTmp = mi_buf.getNumCols();
                for (int row = 0; row < bufState.getNumRows(); row++) {
                  for (int col = 1; col < bufState.getNumCols(); col += 2) {
                    if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                      bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                      colForBufStateTmp++;
                    }
                  }
                  colForBufStateTmp = mi_buf.getNumCols();
                }
                if (state.isEmpty()) {
                  state = Matrix.decorate(bufState, kState);
                } else {
                  state = Matrix.concatRows(state, Matrix.decorate(bufState, kState), null);
                }
              }
            }
          }
          space = state;
          break;

        case SJF:
        case LJF:
          // in these policies the state space includes continuous random variables
          // for the service times
          System.err.format("The scheduling policy does not admit a discrete state space");

      }




    } else if (sn.nodetypes.get(ind) == NodeType.Cache) {
      // TODO: handle cache node
    }

     //Required to sort empty state as first
    List<Matrix> uniqueRows = new ArrayList<>();
    for (int i = 0; i < space.getNumRows(); i++) {
      Matrix tmp = new Matrix(1, space.getNumCols());
      Matrix tmp2 = new Matrix(1, space.getNumCols());
      Matrix.extractRows(space, i, i + 1, tmp);
      boolean unique = true;
      for (int j = i + 1; j < space.getNumRows(); j++) {
        Matrix.extractRows(space, j, j + 1, tmp2);
        if (tmp.isEqualTo(tmp2)) {
          unique = false;
        }
      }
      if (unique) {
        uniqueRows.add(tmp);
      }
    }

    Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
    // this ensures that states where jobs start in phase 1 are first, which is used eg
    // in SSA
    int row = 0;
    Comparator<Matrix> lexico = (mat1, mat2) -> {
      for (int col = 0; col < mat1.getNumCols(); col++) {
        int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
        if (result != 0) {
          return result;
        }
      }
      return 0;
    };

    uniqueRows.sort(lexico);

    for (int i = uniqueRows.size() - 1; i >= 0; i--) {
      for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
        newSpace.set(row, j, uniqueRows.get(i).get(0, j));
      }
      row++;
    }

    return newSpace;
  }









  public static Matrix fromMarginalAndStarted(
          NetworkStruct sn, int ind, Matrix n, Matrix s) {
    return fromMarginalAndStarted(sn, ind, n, s, true);
  }

  public static Matrix fromMarginalAndStarted(
          Network network, int ind, Matrix n, Matrix s) {
    return fromMarginalAndStarted(network.getStruct(true), ind, n, s, true);
  }

  public static Matrix fromMarginalAndStarted(
          NetworkStruct sn, int ind, Matrix n, Matrix s, Boolean optionsForce) {

    // generate one initial state such that the marginal queue-lengths are as in vector n
    // n(r): number of jobs at the station in class r
    // s(r): jobs of class r that are running
    int R = sn.nclasses;
    Matrix S = sn.nservers;
    int ist = (int) sn.nodeToStation.get(0, ind);

    Matrix K = new Matrix(1, R);
    for (int r = 0; r < R; r++) {
      if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
        K.set(0, r, 0);
      } else {
        K.set(0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
      }
    }

    Matrix state = new Matrix(0, 0);
    Matrix space = new Matrix(0, 0);
    LinkedList<Integer> exceeded = new LinkedList<>();
    for (int i = 0; i < sn.classcap.getNumCols(); i++) {
      if (n.get(0, i) > sn.classcap.get(ist, i)) {
        exceeded.add(i);
      }
    }
    if (!exceeded.isEmpty()) {
      for (Integer r : exceeded) {
        if (!sn.proc.isEmpty()
            && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()
            && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
          System.err.format(
              "State vector at station %d exceeds the class capacity. Some service classes are disabled.\n",
              ist);
        } else {
          System.err.format("State vector at station %d exceeds the class capacity.\n", ist);
        }
      }
      return space;
    }

    if (sn.nservers.get(ist, 0) > 0
        && s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols()) > sn.nservers.get(ist, 0)) {
      return space;
    }

    // Generate local-state space
    if ((sn.nodetypes.get(ind) == NodeType.Queue)
        || (sn.nodetypes.get(ind) == NodeType.Delay)
        || (sn.nodetypes.get(ind) == NodeType.Source)) {
      switch (sn.sched.get(sn.stations.get(ist))) {
        case EXT:
          for (int r = 0; r < R; r++) {
            Matrix init = State.spaceClosedSingle(K.get(0, r), 0);
            if (Double.isInfinite(sn.njobs.get(0, r))) {
              if ((!sn.proc.isEmpty())
                  && (!sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty())
                  && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) {
                init.set(0, 0, 0); // class is not processed at this source
              } else {
                // init the job generation
                init.set(0, 0, 1);
              }
            }
            state = Matrix.decorate(state, init);
          }
          space = Matrix.decorate(space, state);
          Matrix ones = new Matrix(space.getNumRows(), 1);
          ones.ones();
          Matrix infBuffer =
                  ones.mult(new Matrix(1,1).fromArray2D(new double[][]{{Double.POSITIVE_INFINITY}}));
          space = infBuffer.concatCols(space);
          break;
        case INF:
        case DPS:
        case GPS:
        case PS:
          // In these policies we only track the jobs in the servers
          for (int r = 0; r < R; r++) {
            Matrix init = State.spaceClosedSingle(K.get(0, r), 0);
            init.set(0, 0, n.get(0, r));
            state = Matrix.decorate(state, init);
          }
          space = Matrix.decorate(space, state);
          break;

        case SIRO:
        case LEPT:
        case SEPT:
          // In these policies we track an un-ordered buffer and the jobs in the servers build list
          // of job classes in the node, with repetition
          if (n.elementSum() <= S.get(ist)) {
            for (int r = 0; r < R; r++) {
              Matrix init = spaceClosedSingle(K.get(0,r), 0);
              init.set(0,0,n.get(0,r));
              state = Matrix.decorate(state, init);
            }
            Matrix newStates = new Matrix(state.getNumRows(), R);
            newStates.zero();
            newStates = newStates.concatCols(state);
            space = Matrix.decorate(space, newStates);
          } else {
            Matrix si = s.clone();
            Matrix mi_buf = n.repmat(si.getNumRows(),1).sub(1, si); // jobs of class r in buffer
            for (int k = 0; k < si.getNumRows(); k++) {
              Matrix kstate = new Matrix(0,0);
              for (int r = 0; r < R; r++) {
                Matrix init = spaceClosedSingle(K.get(0,r), 0);
                init.set(0,0,si.get(k,r));
                kstate = Matrix.decorate(kstate, init);
              }
              state = Matrix.extractRows(mi_buf, k, k+1, null).repmat(kstate.getNumRows(),1).concatCols(kstate);
              if (space.isEmpty()) {
                space = state.clone();
              } else {
                space = Matrix.concatRows(space, state, null);
              }
            }
          }
          break;
        case FCFS:
        case HOL:
        case LCFS:
          Matrix inbuf = new Matrix(0,0);
          double sizeEstimator = 0;
          Matrix mi = new Matrix(0,0);
          Matrix mi_buf = new Matrix(0,0);
          Matrix mi_srv = new Matrix(0,0);
          if (n.elementSum() == 0) {
            space = new Matrix(1, (int) (1+Maths.max(R, K.elementSum())));
            return space;
          }

          // In these policies we track an ordered buffer and the jobs in the servers
          // build list of job classes in the buffer, with repetition
          inbuf = new Matrix(0, 0);
          for (int r = 0; r < R; r++) {
            if (n.get(0, r) > 0) {
              int numNewCols = (int) (n.get(0,r) - s.get(0,r));
              Matrix newInBuf = new Matrix(1, inbuf.getNumCols() + numNewCols);
              for (int i = 0; i < inbuf.getNumCols(); i++) {
                newInBuf.set(0,i,inbuf.get(0,i));
              }
              for (int i = inbuf.getNumCols(); i < newInBuf.getNumCols(); i++) {
                newInBuf.set(0,i,r+1);
              }
              inbuf = newInBuf.clone();
            }
          }

          sizeEstimator = Maths.multinomialln(n);
          sizeEstimator = Math.round(sizeEstimator / Math.log(10));
          if (sizeEstimator > 2) {
            if (!optionsForce) {
              System.err.format(
                      "State space size is very large: 1e%f states. Cannot generate valid state space. Initializing station %d from a default state.\n",
                      sizeEstimator, ind);
              state = inbuf.clone();
              return state;
            }
          }

          // Gen permutation of their positions in the FCFS buffer
          mi = Maths.uniquePerms(inbuf);
          double sumN = n.sumSubMatrix(0, n.getNumRows(), 0, n.getNumCols());
          double sumS = s.sumSubMatrix(0, s.getNumRows(), 0, s.getNumCols());
          if (mi.isEmpty()) {
            mi_buf = new Matrix(1, (int) Math.max(1, sumN - S.get(ist, 0)));
            mi_buf.zero();
            state = new Matrix(1, (int) K.sumSubMatrix(0, K.getNumRows(), 0, K.getNumCols()));
            state.zero();
            Matrix newState = new Matrix(1, mi_buf.getNumCols() + state.getNumCols());
            for (int i = 0; i < mi_buf.getNumCols(); i++) {
              newState.set(0, i, mi_buf.get(0, i));
            }
            for (int i = mi_buf.getNumCols(); i < newState.getNumCols(); i++) {
              newState.set(0, i, state.get(0, i - mi_buf.getNumCols()));
            }
            state = newState.clone();
          } else {
            // mi_buf: class of job in buffer position i (0 = empty)
            if (sumN > sumS) {
              mi_buf = new Matrix(mi.getNumRows(), (int) sumN - (int) sumS);
              for (int row = 0; row < mi_buf.getNumRows(); row++) {
                for (int col = 0; col < sumN - sumS; col++) {
                  mi_buf.set(row, col, mi.get(row, col));
                }
              }
            } else {
              mi_buf = new Matrix(1, 1);
              mi_buf.set(0, 0, 0);
            }
          }

          // mi_srv: class of jobs running in the server of i
          mi_srv = new Matrix(0, 0);
          for (int r = 0; r < R; r++) {
            Matrix new_mi_srv = new Matrix(1, mi_srv.getNumCols() + (int) s.get(0, r));
            for (int i = 0; i < mi_srv.getNumCols(); i++) {
              new_mi_srv.set(0, i, mi_srv.get(0, i));
            }
            for (int i = mi_srv.getNumCols(); i < new_mi_srv.getNumCols(); i++) {
              new_mi_srv.set(0, i, r);
            }
            mi_srv = new_mi_srv.clone();
          }

          // si: number of class r jobs that are running
          Matrix si = s.clone();
          for (int b = 0; b < mi_buf.getNumRows(); b++) {
            for (int k = 0; k < si.getNumRows(); k++) {
              Matrix kState = new Matrix(0, 0);
              for (int r = 0; r < R; r++) {
                Matrix init = State.spaceClosedSingle(K.get(0, r), 0);
                init.set(0, 0, si.get(k, r));
                kState = Matrix.decorate(kState, init);
              }
              Matrix miBufRep = Matrix.extractRows(mi_buf,b,b+1,null).repmat(kState.getNumRows(),1);
              miBufRep = miBufRep.concatCols(kState);
              if (state.isEmpty()) {
                state = miBufRep;
              } else {
                state = Matrix.concatRows(state, miBufRep, null);
              }

            }
          }
          space = state;
          break;
        case LCFSPR:
          Matrix inbuf_lpr = new Matrix(0,0);
          double sizeEstimator_lpr = 0;
          Matrix mi_lpr = new Matrix(0,0);
          Matrix mi_buf_lpr = new Matrix(0,0);
          Matrix mi_srv_lpr = new Matrix(0,0);
          if (n.elementSum() == 0) {
            space = new Matrix(1, (int) (1 + K.elementSum()));
            return space;
          }
            // in this policy we track an ordered buffer and the jobs in the servers
            // build list of job classes in the buffer, with repetition

            inbuf_lpr = new Matrix(0, 0);
            for (int r = 0; r < R; r++) {
              if (n.get(0, r) > 0) {
                int numNewCols = (int) (n.get(0,r) - s.get(0,r));
                Matrix newInBuf = new Matrix(1, inbuf_lpr.getNumCols() + numNewCols);
                for (int i = 0; i < inbuf_lpr.getNumCols(); i++) {
                  newInBuf.set(0,i,inbuf_lpr.get(0,i));
                }
                for (int i = inbuf_lpr.getNumCols(); i < newInBuf.getNumCols(); i++) {
                  newInBuf.set(0,i,r+1);
                }
                inbuf_lpr = newInBuf.clone();
              }
            }
            sizeEstimator_lpr = Maths.multinomialln(n);
            sizeEstimator_lpr = Math.round(sizeEstimator_lpr/Math.log(10));
            if (sizeEstimator_lpr > 2) {
              if (!optionsForce) {
                System.err.format(
                        "State space size is very large: 1e%f states. Cannot generate valid state space. Initializing station %d from a default state.\n",
                        sizeEstimator_lpr, ind);
                state = inbuf_lpr.clone();
                return state;
              }
            }

            // gen permutation of their positions in the buffer
            mi_lpr = Maths.uniquePerms(inbuf_lpr);
            mi_buf_lpr = new Matrix(0,0);
            if (mi_lpr.isEmpty()) {
              mi_buf_lpr = new Matrix (1, (int) Maths.max(1, n.elementSum() - S.get(ist)));
              state = new Matrix(1, (int) K.elementSum());
              state.zero();
              state = mi_buf_lpr.concatCols(state);
            } else {
              // mi_buf: class of job in buffer position i (0 = empty)
              if (n.elementSum() > s.elementSum()) {
                mi_buf_lpr = new Matrix(mi_lpr.getNumRows(), (int) (n.elementSum() - s.elementSum()));
                for (int row = 0; row < mi_buf_lpr.getNumRows(); row++) {
                  for (int col = 0; col < mi_buf_lpr.getNumCols(); col++) {
                    mi_buf_lpr.set(row, col, mi_lpr.get(row, col));
                  }
                }
              } else {
                mi_buf_lpr = new Matrix(1,1);
                mi_buf_lpr.zero();
              }
            }

            // mi_srv: class of jobs running in the server of i
            mi_srv_lpr = new Matrix(0,0);
            for (int r = 0; r < R; r++) {
              if (n.get(0, r) > 0) {
                Matrix newMiSrv =
                        new Matrix(1, mi_srv_lpr.getNumCols() + (int) s.get(0, r));
                for (int i = 0; i < mi_srv_lpr.getNumCols(); i++) {
                  newMiSrv.set(0, i, mi_srv_lpr.get(0, i));
                }
                for (int i = mi_srv_lpr.getNumCols(); i < newMiSrv.getNumCols(); i++) {
                  newMiSrv.set(0, i, r+1);
                }
                mi_srv_lpr = newMiSrv.clone();
              }
            }


            // si: number of class r jobs that are running
            Matrix si_lpr = s.clone();

            for (int b = 0; b < mi_buf_lpr.getNumRows(); b++) {
              for (int k = 0; k < si_lpr.getNumRows(); k++) {
                Matrix kState = new Matrix(0,0);
                for (int r = 0; r < R; r++) {
                  Matrix init = spaceClosedSingle(K.get(r), 0);
                  init.set(0,0, si_lpr.get(k, r));
                  kState = Matrix.decorate(kState, init);
                }
                Matrix bkState = new Matrix(0,0);
                Matrix jobsInBuffer = Matrix.extractRows(mi_buf_lpr, b, b+1, null);
                for (int j = 0; j < jobsInBuffer.length(); j++) {
                  double job = jobsInBuffer.get(j);
                  if (job > 0) {
                    List<Double> phasesJRange = new ArrayList<>();
                    for (double i = 1; i <= K.get((int) job - 1); i++) {
                      phasesJRange.add(i);
                    }
                    // no transpose as constructor makes column vector
                    bkState = Matrix.decorate(bkState, new Matrix(phasesJRange));
                  } else {
                    bkState = new Matrix(1,1);
                    bkState.zero();
                  }
                }
                Matrix bufStateTmp = Matrix.decorate(Matrix.extractRows(mi_buf_lpr, b, b + 1, null), bkState);
                // here we interleave positions of class and phases in buffer
                Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                bufState.zero();

                // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                int colForBufStateTmp = 0;
                for (int row = 0; row < bufState.getNumRows(); row++) {
                  for (int col = 0; col < bufState.getNumCols(); col += 2) {
                    if (colForBufStateTmp < mi_buf_lpr.getNumCols()) {
                      bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                      colForBufStateTmp++;
                    }
                  }
                  colForBufStateTmp = 0;
                }
                colForBufStateTmp = mi_buf_lpr.getNumCols();
                for (int row = 0; row < bufState.getNumRows(); row++) {
                  for (int col = 1; col < bufState.getNumCols(); col += 2) {
                    if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                      bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                      colForBufStateTmp++;
                    }
                  }
                  colForBufStateTmp = mi_buf_lpr.getNumCols();
                }
                if (state.isEmpty()) {
                  state = Matrix.decorate(bufState, kState);
                } else {
                  state = Matrix.concatRows(state, Matrix.decorate(bufState, kState), null);
                }
              }
            }
            space = state;
            break;
        case SJF:
        case LJF:
          // In these policies the state space includes continuous random variables for the service
          // times in these policies we only track the jobs in the servers
          for (int r = 0; r < R; r++) {
            Matrix init = spaceClosedSingle(K.get(0,r), 0);
            init.set(0,0,n.get(0, r));
            state = Matrix.decorate(state, init);
          }
          space = Matrix.decorate(space, state);
          System.err.format("The schedyling policy does not admit a discrete state space");
          break;
      }

    } else if (sn.nodetypes.get(ind) == NodeType.Cache) {
      // TODO: finish remainder of implementation, lines 264 to 272
      System.out.println(
              "Warning: unimplemented code reached in fromMarginalAndStarted for Cache Nodes");
    } else if (sn.nodetypes.get(ind) == NodeType.Join ||sn.nodetypes.get(ind) == NodeType.Transition
            || sn.nodetypes.get(ind) == NodeType.Place) {
      space = new Matrix(1, 1);
      space.zero();
    }

    // Required to sort empty state as first
    List<Matrix> uniqueRows = new ArrayList<>();
    for (int i = 0; i < space.getNumRows(); i++) {
      Matrix tmp = new Matrix(1, space.getNumCols());
      Matrix tmp2 = new Matrix(1, space.getNumCols());
      Matrix.extractRows(space, i, i + 1, tmp);
      boolean unique = true;
      for (int j = i + 1; j < space.getNumRows(); j++) {
        Matrix.extractRows(space, j, j + 1, tmp2);
        if (tmp.isEqualTo(tmp2)) {
          unique = false;
        }
      }
      if (unique) {
        uniqueRows.add(tmp);
      }
    }

    Comparator<Matrix> lexico = (mat1, mat2) -> {
      for (int col = 0; col < mat1.getNumCols(); col++) {
        int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
        if (result != 0) {
          return result;
        }
      }
      return 0;
    };

    uniqueRows.sort(lexico);
    Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
    // this ensures that states where jobs start in phase 1 are first, which is used eg
    // in SSA
    int row = 0;
    for (int i = uniqueRows.size() - 1; i >= 0; i--) {
      for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
        newSpace.set(row, j, uniqueRows.get(i).get(0, j));
      }
      row++;
    }

    return newSpace;
  }

  private static Matrix spaceClosedSingle(double M, double N) {

    if (M != 0) {
      return Maths.multiChoose(M, N);
    }
    return new Matrix(0, 0);
  }

  public static boolean isValid(Network sn, Matrix n, Matrix s) {
    return isValid(sn.getStruct(true), n, s);
  }

  // TODO: may need to create another polymorphic version where n is a cell array/list of JLMs
  // Currently no need at this stage as not using as part of SolverFluid implementation
  public static boolean isValid(NetworkStruct sn, Matrix n, Matrix s) {

    // n(r): number of jobs at the station in class r
    // s(r): jobs of class r that are running

    if (n.isEmpty() & !s.isEmpty()) {
      return false;
    }

    // Note that at this point in LINE there is a check whether n has more than one state possible
    // by checking whether n is a cell, not a single matrix. I believe this is not possible
    // currently as state has been implemented as a Hashmap, so no duplicate keys, so only one state
    // possible. If this changes i.e. possible to have more than one possible starting state, then
    // would need to implement this part of the method (lines 21 through 29).

    int R = sn.nclasses;
    Matrix K = new Matrix(1, R);
    K.zero();

    for (int ist = 0; ist < sn.nstations; ist++) {
      for (int r = 0; r < R; r++) {
        K.set(0, r, sn.phases.get(ist, r));
        if (sn.nodetypes.get((int) sn.stationToNode.get(0, ist)) != NodeType.Place) {
          if (!sn.proc.isEmpty()
              && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()
              && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()
              && n.get(ist, r) > 0) { // if disabled
            return false;
          }
        }
      }

      for (int j = 0; j < n.getNumCols(); j++) {
        if (n.get(ist, j) > sn.classcap.get(ist, j)) {
          System.out.println("n: " + n.get(ist, j) + " classcap: " + sn.classcap.get(ist, j));
          System.err.format(
              "Station %d is in a state with more jobs than its allowed capacity. ", ist);
          return false;
        }
      }
    }

    if (!s.isEmpty()) {
      for (int ist = 0; ist < sn.nstations; ist++) {
        if (sn.nservers.get(ist, 0) > 0) {
          // If more running jobs than servers
          if (s.sumRows(ist) > sn.nservers.get(ist, 0)) {
            // Don't flag invalid if PS
            SchedStrategy schedStrat = sn.sched.get(sn.stations.get(ist));
            if (schedStrat == SchedStrategy.FCFS
                || schedStrat == SIRO
                || schedStrat == LCFS
                || schedStrat == HOL) {
              return false;
            }
          }
          // if more running jobs than jobs at the node
          for (int row = 0; row < n.getNumRows(); row++) {
            for (int col = 0; col < n.getNumCols(); col++) {
              if (n.get(row, col) < s.get(row, col)) {
                return false;
              }
            }
          }
        }
      }
    }

    for (int nc = 0; nc < sn.nchains; nc++) {
      double njobs_chain = 0;
      LinkedList<Integer> chainsIdx = new LinkedList<>();
      for (int i = 0; i < sn.chains.getNumCols(); i++) {
        if (sn.chains.get(nc, i) > 0) {
          chainsIdx.add(i);
          njobs_chain += sn.njobs.get(0, i);
        }
      }
      double statejobs_chain = 0;
      if (!Double.isInfinite(njobs_chain)) {
        for (int i = 0; i < n.getNumRows(); i++) {
          for (Integer idx : chainsIdx) {
            statejobs_chain += n.get(i, idx);
          }
        }
        if (Math.abs(1 - (njobs_chain / statejobs_chain)) > 0.0001) {
          System.err.format(
              "Chain %d is initialized with an incorrect number of jobs: %f instead of %f.",
              nc, statejobs_chain, njobs_chain);
          return false;
        }
      }
    }

    return true;
  }


  public static Matrix fromMarginal(NetworkStruct sn, int ind, Matrix n) {

    // Generate states such that the marginal queue-lengths are as in vector n
    // n(r): number of jobs at the station in class r
    int R = sn.nclasses;
    Matrix S = sn.nservers;
    Matrix state = new Matrix(0, 0);
    Matrix space = new Matrix(0, 0);

    if (sn.isstation.get(ind) == 1) {
      boolean mapOrMMPP2ProcTypes = false;
      for (int i = 0; i < sn.proctype.get(sn.stations.get(ind)).size(); i++) {
        if (sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(i)) == ProcessType.MAP
            || sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(i))
                == ProcessType.MMPP2) {
          mapOrMMPP2ProcTypes = true;
        }
      }
      if (mapOrMMPP2ProcTypes) {
        if (sn.nservers.get(ind) > 1) {
          throw new RuntimeException("Multiserver MAP stations are not supported.");
        }
        if (sn.sched.get(sn.stations.get(ind)) == FCFS && !sn.nodetypes.contains(NodeType.Source)) {
          throw new RuntimeException("Non-FCFS MAP stations are not supported.");
        }
      }
    }

    int ist = (int) sn.nodeToStation.get(ind);
    int isf = (int) sn.nodeToStateful.get(ind);

    if (sn.isstateful.get(ind, 0) == 1 && sn.isstation.get(ind, 0) == 0) {
      for (int r = 0; r < R; r++) {
        Matrix init_r = spaceClosedSingle(1, n.get(r));
        state = Matrix.decorate(state, init_r);
      }
      return Matrix.decorate(space, state);
    }

    Matrix phases = new Matrix(1, R);
    for (int r = 0; r < R; r++) {
      if (sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()) {
        phases.set(0, r, 0);
      } else {
        phases.set(
            0, r, sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).length());
      }
    }

    boolean anyNGreaterThanClassCap = false;
    for (int i = 0; i < sn.classcap.getNumCols(); i++) {
      if (n.get(i) > sn.classcap.get(ist, i)) {
        anyNGreaterThanClassCap = true;
      }
    }
    if (sn.sched.get(sn.stations.get(ind)) != EXT && anyNGreaterThanClassCap) {
      return space;
    }

    // Generate local-state space
    switch (sn.nodetypes.get(ind)) {
      case Queue:
      case Delay:
      case Source:
        switch (sn.sched.get(sn.stations.get(ist))) {
          case EXT:
            // source node case, treated as an infinite pool of jobs in buffer and a server for each class
            for (int r = 0; r < R; r++) {
              Matrix init_r = new Matrix(0,0);
              if (!sn.proc.isEmpty() && !sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).isEmpty()
                && sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).get(0).hasNaN()) { // if service disabled
                init_r = new Matrix(1, (int) phases.get(r));
                init_r.zero();
              } else {
                init_r = State.spaceClosedSingle(phases.get(r), 1);
              }
              state = Matrix.decorate(state, init_r);
            }
            space = Matrix.decorate(space, state); // server part

            Matrix ones = new Matrix(space.getNumRows(),1);
            ones.ones();
            Matrix infBuffer = ones.mult(new Matrix(1,1).fromArray2D(new double[][]{{Double.POSITIVE_INFINITY}}));
            // Attach infinite buffer to all state spaces containing job distribution across servers
           space = infBuffer.concatCols(space);
           break;
          case INF:
          case PS:
          case DPS:
          case GPS:
            // In these policies we only track the jobs in the servers
            for (int r = 0; r < R; r++) {
              Matrix init_r = spaceClosedSingle(phases.get(r), n.get(r));
              state = Matrix.decorate(state, init_r);
            }
            space = Matrix.decorate(space, state);
            break;
          case SIRO:
          case LEPT:
          case SEPT:
            // In these policies we track an un-ordered buffer and the jobs in the servers.
            // We build list of job classes in the node, with repetition
            if (n.elementSum() <= S.get(ist)) {
              // buffer will be empty as we have enough servers to handle all tasks in parallel
              for (int r = 0; r < R; r++) {
                Matrix init_r = State.spaceClosedSingle(phases.get(r), n.get(r));
                state = Matrix.decorate(state, init_r);
              }
              Matrix newStates = new Matrix(state.getNumRows(), R);
              newStates.zero();
              newStates = newStates.concatCols(state);
              space = Matrix.decorate(space, newStates);
            } else {
              Matrix si = Maths.multiChooseCon(n, S.get(ist)); // jobs of class r that are running
              Matrix mi_buf = n.repmat(si.getNumRows(),1).sub(1, si);
              for (int k = 0; k < si.getNumRows(); k++) {
                // determine number of class r jobs running in phase j
                Matrix kstate = new Matrix(0, 0);
                for (int r = 0; r < R; r++) {
                  Matrix init_r = State.spaceClosedSingle(phases.get(r), si.get(k,r));
                  kstate = Matrix.decorate(kstate, init_r);
                }
                state = Matrix.extractRows(mi_buf, k, k+1, null).repmat(kstate.getNumRows(),1).concatCols(kstate);
                if (space.isEmpty()) {
                  space = state.clone();
                } else {
                  space = Matrix.concatRows(space, state, null);
                }
              }
            }
            break;
          case FCFS:
          case HOL:
          case LCFS:
            Matrix vi = new Matrix(0,0);
            Matrix mi = new Matrix(0,0);
            double sizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) +
                    Maths.factln(sn.cap.get(ist));
            sizeEstimator = Math.round(sizeEstimator/Math.log(10));
            if (sizeEstimator > 3) {
              // TODO: Options force and line warning. Line warning commented in MATLAB
            }
            if (n.elementSum() == 0) {
              space = new Matrix(1, (int) (1+phases.elementSum()));
              if (!sn.nodetypes.get(ind).equals(NodeType.Source)) {
                for (int r = 0; r < R; r++) {
                  switch (sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(ind))) {
                    case MAP:
                    case MMPP2:
                      List<Double> phasesRange = new ArrayList<>();
                      for (double i = 1; i <= sn.phases.get(ind, r); i++) {
                        phasesRange.add(i);
                      }
                      // no transpose as constructor creates a column vector
                      space = Matrix.decorate(space, new Matrix(phasesRange));
                  }
                }
              }
            }
            vi = new Matrix(0, 0);
            for (int r = 0; r < R; r++) {
              if (n.get(0, r) > 0) {
                Matrix newVi =
                        new Matrix(1, vi.getNumCols() + (int) n.get(0, r));
                for (int i = 0; i < vi.getNumCols(); i++) {
                  newVi.set(0, i, vi.get(0, i));
                }
                for (int i = vi.getNumCols(); i < newVi.getNumCols(); i++) {
                  newVi.set(0, i, r+1);
                }
                vi = newVi.clone();
              }
            }
            // gen permutation of their positions in the waiting buffer
            mi = Maths.uniquePerms(vi);
            Matrix mi_buf = new Matrix(0,0);
            // now generate server states
            if (mi.isEmpty()) {
              mi_buf = new Matrix(1, (int) Math.max(0, n.elementSum() - S.get(ist)));
              mi_buf.zero();
              state = new Matrix(1, R);
              state.zero();
              state = Matrix.decorate(state, mi_buf.concatCols(state));
            } else {
              int numCols = (int) Maths.min(n.elementSum(), sn.cap.get(ist));
              Matrix miClone = mi.clone();
              mi = new Matrix(mi.getNumRows(), numCols);

              for (int row = 0; row < miClone.getNumRows(); row++) {
                for (int col = miClone.getNumCols() - numCols; col < miClone.getNumCols(); col++) {
                  mi.set(row, col, miClone.get(row, col));
                }
              }
              mi = uniqueAndSort(mi);

              // mi_buf: class of job in buffer position i (0=empty)
              int numColumnsRight = (int) Maths.max((mi.getNumCols()-S.get(ist)), 0);
              Matrix right = new Matrix(mi.getNumRows(), numColumnsRight);
              Matrix.extract(mi, 0, mi.getNumRows(), 0,
                      numColumnsRight, right, 0, 0);
              double x= sn.cap.get(ist);
              int numColumnsLeft = (int) Maths.max(0, (Maths.min(n.elementSum(), sn.cap.get(ist)) - S.get(ist)
                      - right.getNumCols()));
              Matrix left = new Matrix(mi.getNumRows(), numColumnsLeft);
              mi_buf =  left.concatCols(right);
              if (mi_buf.isEmpty()) {
                mi_buf = new Matrix(mi.getNumRows(), 1);
                mi_buf.zero();
              }
              // mi_srv: class of job running in server i
              int numColsSrv = (int) Maths.max(S.get(ist), 1);
              Matrix miSrv = new Matrix(mi.getNumRows(), numColsSrv);

              int colForMiSrv = 0;
              for (int row = 0; row < miSrv.getNumRows(); row++) {
                colForMiSrv = 0;
                for (int col = mi.getNumCols() - numColsSrv; col < mi.getNumCols(); col++) {
                  miSrv.set(row, colForMiSrv, mi.get(row, col));
                  colForMiSrv++;
                }
              }

              // si: number of class r jobs that are running
              Matrix si = new Matrix(miSrv.getNumRows(), R);
              for (int k = 0; k < mi.getNumRows(); k++) {
                Matrix miSrvKRow = Matrix.extractRows(miSrv,k,k+1,null);
                Matrix histRow = Maths.hist(miSrvKRow, 1, R);
                for (int j = 0; j < R; j++) {
                  si.set(k, j, histRow.get(j));
                }
              }

              for (int k = 0; k < si.getNumRows(); k++) {
                // determine number of class r jobs running in phase
                // j in server state mi_srv(k,:) and build state
                Matrix kState = new Matrix(0,0);
                kState.zero();
                Matrix map_cols = new Matrix(0, 0);
                for (int r = 0; r < R; r++) {
                  Matrix init_r = spaceClosedSingle(phases.get(r), si.get(k, r));
                  // TODO: MAP, MMP2 case, lines 248-259
                  kState = Matrix.decorate(kState, init_r).clone();
                  // TODO: MAP, MMP2 case, lines 261-263
                }
                // TODO: modify kState wrt map_cols, line 265
                Matrix miBufreplicated = Matrix.extractRows(mi_buf,k,k+1,null).repmat(kState.getNumRows(),1);
                miBufreplicated = miBufreplicated.concatCols(kState).clone();

                if (state.isEmpty()) {
                  state = miBufreplicated;
                } else {
                  state = Matrix.concatRows(state, miBufreplicated, null);
                }
              }
            }
            space = state;
            break;
          case LCFSPR:
            Matrix vi_lpr = new Matrix(0,0);
            Matrix mi_lpr = new Matrix(0,0);
              // sum(n) - 1 due to Maths.factln including + 1
              double lcfsprSizeEstimator = Maths.multinomialln(n) - Maths.factln(n.elementSum() - 1) +
                      Maths.factln(sn.cap.get(ist));
              lcfsprSizeEstimator = Math.round(lcfsprSizeEstimator/Math.log(10));
              if (lcfsprSizeEstimator > 3) {
                // TODO: Options force and line warning. Line warning commented in MATLAB
              }

              if (n.elementSum() == 0) {
                Matrix newSpace = new Matrix(1, (int) (1+phases.elementSum()));
                newSpace.zero();
                space = newSpace;
                return space;
              }
              // Similar to FCFS/HOL/LCFS case we track an ordered buffer and the jobs in the servers
              // but in this case due to pre-emption jobs in buffer can be in not initial phase

              // build list of job classes in the node, with repetition

              vi_lpr = new Matrix(0, 0);
              for (int r = 0; r < R; r++) {
                if (n.get(0, r) > 0) {
                  Matrix newVi =
                          new Matrix(1, vi_lpr.getNumCols() + (int) n.get(0, r));
                  for (int i = 0; i < vi_lpr.getNumCols(); i++) {
                    newVi.set(0, i, vi_lpr.get(0, i));
                  }
                  for (int i = vi_lpr.getNumCols(); i < newVi.getNumCols(); i++) {
                    newVi.set(0, i, r+1);
                  }
                  vi_lpr = newVi.clone();
                }
              }
              // gen permutation of their positions in the waiting buffer
              mi_lpr = Maths.uniquePerms(vi_lpr);
              // now generate server states
              if (mi_lpr.isEmpty()) {
                Matrix mi_buf_lpr =new Matrix(1, (int) Maths.max(0, n.elementSum() - S.get(ist)));
                mi_buf_lpr.zero();
                state = new Matrix(1, R);
                state.zero();
                state = Matrix.decorate(state, mi_buf_lpr.concatCols(state));
              } else {
                int numCols = (int) Maths.min(n.elementSum(), sn.cap.get(ist));
                Matrix miClone = mi_lpr.clone();
                mi_lpr = new Matrix(mi_lpr.getNumRows(), numCols);

                for (int row = 0; row < miClone.getNumRows(); row++) {
                  for (int col = miClone.getNumCols() - numCols; col < miClone.getNumCols(); col++) {
                    mi_lpr.set(row, col, miClone.get(row, col));
                  }
                }
                // mi_buf: class of job in buffer position i (0 = empty)
                int numColumnsRight = (int) Maths.max((mi_lpr.getNumCols() - S.get(ist)), 0);
                Matrix right = new Matrix(mi_lpr.getNumRows(), numColumnsRight);
                Matrix.extract(mi_lpr, 0, mi_lpr.getNumRows(), 0,
                        numColumnsRight, right, 0, 0);
                int numColumnsLeft = (int) Maths.max(0, (Maths.min(n.elementSum(), sn.cap.get(ist)) - S.get(ist)
                        - right.getNumCols()));
                Matrix left = new Matrix(mi_lpr.getNumRows(), numColumnsLeft);
                Matrix mi_buf_lpr = left.concatCols(right);
                if (mi_buf_lpr.isEmpty()) {
                  mi_buf_lpr = new Matrix(mi_lpr.getNumRows(), 1);
                  mi_buf_lpr.zero();
                }
                // miSrv: class of job running in server i
                int numColsSrv = (int) Maths.max(S.get(ist), 1);
                Matrix miSrv = new Matrix(mi_lpr.getNumRows(), numColsSrv);

                int colForMiSrv = 0;
                for (int row = 0; row < miSrv.getNumRows(); row++) {
                  colForMiSrv = 0;
                  for (int col = mi_lpr.getNumCols() - numColsSrv; col < mi_lpr.getNumCols(); col++) {
                    miSrv.set(row, colForMiSrv, mi_lpr.get(row, col));
                    colForMiSrv++;
                  }
                }

                // si: number of class r jobs that are running
                Matrix si = new Matrix(miSrv.getNumRows(), R);
                for (int k = 0; k < mi_lpr.getNumRows(); k++) {
                  Matrix miSrvKRow = Matrix.extractRows(miSrv, k, k + 1, null);
                  Matrix histRow = Maths.hist(miSrvKRow, 1, R);
                  for (int j = 0; j < R; j++) {
                    si.set(k, j, histRow.get(j));
                  }
                }
                for (int k = 0; k < si.getNumRows(); k++) {
                  // determine number of class r jobs running in phase j
                  // in server state miSrv(k, :) and build state
                  Matrix kState = new Matrix(0, 0);
                  for (int r = 0; r < R; r++) {
                    kState = Matrix.decorate(kState, State.spaceClosedSingle(phases.get(r), si.get(k, r)));
                  }
                  // generate job phases for all buffer states since we have pre-emption
                  Matrix bkState = new Matrix(0, 0);

                  Matrix jobsInBuffer = Matrix.extractRows(mi_buf_lpr, k, k + 1, null);
                  for (int j = 0; j < jobsInBuffer.length(); j++) {
                    double job = jobsInBuffer.get(j);
                    if (job > 0) {
                      List<Double> phasesJRange = new ArrayList<>();
                      for (double i = 1; i <= phases.get((int) job - 1); i++) {
                        phasesJRange.add(i);
                      }
                      // no transpose as constructor makes column vector
                      bkState = Matrix.decorate(bkState, new Matrix(phasesJRange));
                    } else {
                      bkState = new Matrix(1, 1);
                      bkState.zero();
                    }
                  }
                  Matrix bufStateTmp = Matrix.decorate(Matrix.extractRows(mi_buf_lpr, k, k + 1, null), bkState);
                  // here we interleave positions of class and phases in buffer
                  Matrix bufState = new Matrix(bufStateTmp.getNumRows(), bufStateTmp.getNumCols());
                  bufState.zero();

                  // bufstateTmp has classses followrd by phases. here we interleave the classes and phases
                  int colForBufStateTmp = 0;
                  for (int row = 0; row < bufState.getNumRows(); row++) {
                    for (int col = 0; col < bufState.getNumCols(); col += 2) {
                      if (colForBufStateTmp < mi_buf_lpr.getNumCols()) {
                        bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                        colForBufStateTmp++;
                      }
                    }
                    colForBufStateTmp = 0;
                  }
                  colForBufStateTmp = mi_buf_lpr.getNumCols();
                  for (int row = 0; row < bufState.getNumRows(); row++) {
                    for (int col = 1; col < bufState.getNumCols(); col += 2) {
                      if (colForBufStateTmp < bufStateTmp.getNumCols()) {
                        bufState.set(row, col, bufStateTmp.get(row, colForBufStateTmp));
                        colForBufStateTmp++;
                      }
                    }
                    colForBufStateTmp = mi_buf_lpr.getNumCols();
                  }
                  if (state.isEmpty()) {
                    state = Matrix.decorate(bufState, kState);
                  } else {
                    state = Matrix.concatRows(state, Matrix.decorate(bufState, kState), null);
                  }
                }
              }
            space = state;
            break;
          case SJF:
          case LJF:
            // in these policies the state space includes continuous
            // random variables for the service times
            throw new RuntimeException("The scheduling policy does not admit a discrete state space.");
        }

        for (int r = 0; r < R; r++) {
          if (sn.routing.get(sn.stations.get(ind)).get(sn.jobclasses.get(r))
                  == RoutingStrategy.RROBIN
              || sn.routing.get(sn.stations.get(ind)).get(sn.jobclasses.get(r))
                  == RoutingStrategy.WRROBIN) {
            System.out.println("Unimplemented code reached in NetworkState fromMarginal 8");
          }
        }
        break;
      case Cache:
        System.out.println("Unimplemented code reached in NetworkState fromMarginal 9");
        break;
    }

    // Required to sort empty state as first
    List<Matrix> uniqueRows = new ArrayList<>();
    for (int i = 0; i < space.getNumRows(); i++) {
      Matrix tmp = new Matrix(1, space.getNumCols());
      Matrix tmp2 = new Matrix(1, space.getNumCols());
      Matrix.extractRows(space, i, i + 1, tmp);
      boolean unique = true;
      for (int j = i + 1; j < space.getNumRows(); j++) {
        Matrix.extractRows(space, j, j + 1, tmp2);
        if (tmp.isEqualTo(tmp2)) {
          unique = false;
        }
      }
      if (unique) {
        uniqueRows.add(tmp);
      }
    }

    Comparator<Matrix> lexico = (mat1, mat2) -> {
      for (int col = 0; col < mat1.getNumCols(); col++) {
        int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
        if (result != 0) {
          return result;
        }
      }
      return 0;
    };

    uniqueRows.sort(lexico);
    Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
    // So that states with jobs in phase 1 comes earlier
    int row = 0;
    for (int i = uniqueRows.size() - 1; i >= 0; i--) {
      for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
        newSpace.set(row, j, uniqueRows.get(i).get(0, j));
      }
      row++;
    }
//    return space;
    return newSpace;
  }


  // "unique" function in matlab sorts in lexicographic, replicated here
  private static Matrix uniqueAndSort(Matrix space) {
    List<Matrix> uniqueRows = new ArrayList<>();
    for (int i = 0; i < space.getNumRows(); i++) {
      Matrix tmp = new Matrix(1, space.getNumCols());
      Matrix tmp2 = new Matrix(1, space.getNumCols());
      Matrix.extractRows(space, i, i + 1, tmp);
      boolean unique = true;
      for (int j = i + 1; j < space.getNumRows(); j++) {
        Matrix.extractRows(space, j, j + 1, tmp2);
        if (tmp.isEqualTo(tmp2)) {
          unique = false;
        }
      }
      if (unique) {
        uniqueRows.add(tmp);
      }
    }

    Comparator<Matrix> lexico = (mat1, mat2) -> {
      for (int col = 0; col < mat1.getNumCols(); col++) {
        int result = Integer.compare((int) mat1.get(0, col), (int) mat2.get(0, col));
        if (result != 0) {
          return result;
        }
      }
      return 0;
    };

    uniqueRows.sort(lexico);
    Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
    // So that states with jobs in phase 1 comes earlier
    for (int i = 0; i <uniqueRows.size(); i++) {
      for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
        newSpace.set(i, j, uniqueRows.get(i).get(0, j));
      }
    }

    return newSpace;
  }

  public static void main(String[] args) {

//    ARRIVALS:
//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, );
//    inspace.fromArray2D(new int[][]{{1,2,1,1}});
//    EventResult result = afterEvent(sn, 0, inspace, EventType.ARV, 0, false);
//
//
//    System.out.println(result.outspace);
//    System.out.println(result.outrate);
//    System.out.println(result.outprob);

//
//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{0,0,1,1,0,1,1}});
//    EventResult result = afterEvent(sn, 2, inspace, EventType.ARV, 3, false);
//
//

    // DEPARTURES:
//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 4);
//    inspace.fromArray2D(new int[][]{{4,8,3,2}});
//    EventResult result = afterEvent(sn, 0, inspace, EventType.DEP, 3, false);


//    NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
//    Matrix inspace = new Matrix(1, 3);
//    inspace.fromArray2D(new int[][]{{1,0,0}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.DEP, 0, false);

//    NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
//    Matrix inspace = new Matrix(1, 3);
//    inspace.fromArray2D(new int[][]{{1,0,0}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.DEP, 0, false);

    // FCFS
//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{4,2,1,1,1,1,1}});
//    EventResult result = afterEvent(sn, 2, inspace, EventType.DEP, 2, false);

//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{4,2,1,1,1,1,1}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.DEP, 0, false);

//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{0,2,1,1,0,1,1}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.DEP, 1, false);
//
//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{0,0,0,1,0,0,1}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.DEP, 1, false);

//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{0,0,0,1,0,0,1}});
//    EventResult result = afterEvent(sn, 2, inspace, EventType.DEP, 3, false);

//    NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//    Matrix inspace = new Matrix(1, 7);
//    inspace.fromArray2D(new int[][]{{0,0,1,1,0,1,1}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.PHASE, 1, false);

//    NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
//    Matrix inspace = new Matrix(1, 4);
//    inspace.fromArray2D(new int[][]{{1,0,1,0}});
//    EventResult result = afterEvent(sn, 0, inspace, EventType.PHASE, 0, false);

//    NetworkStruct sn = ClosedModel.ex2_line().getStruct(false);
//    Matrix inspace = new Matrix(1, 3);
//    inspace.fromArray2D(new int[][]{{1,1,1}});
//    EventResult result = afterEvent(sn, 1, inspace, EventType.PHASE, 0, true);


    NetworkStruct sn = ClosedModel.ex4_line().getStruct(false);
    Matrix inspace = new Matrix(1, 7);
    inspace.fromArray2D(new int[][]{{0,2,0,0,0,0,0}});

    // time this call
    long startTime = System.nanoTime();
    EventResult result = afterEvent(sn, 1, inspace, EventType.ARV, 0, true);
    long endTime = System.nanoTime();
    long elapsedTimeSec = (endTime - startTime) / 1_000_000; // Convert nanoseconds to milli

//
//    System.out.println("Outspace: ");
//    System.out.println(result.outspace);
//    System.out.println("Outrate: ");
//    System.out.println(result.outrate);
//    System.out.println("Outprob");
//    System.out.println(result.outprob);
    System.out.println(result.outspace);
    System.out.println(elapsedTimeSec);


  }


  public static EventResult afterEvent(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation) {


    int M = sn.nstations;
    int R = sn.nclasses;
    Matrix S = sn.nservers;
    Matrix phasessz = sn.phasessz;
    Matrix phaseshift = sn.phaseshift;
    Map<Station, Map<JobClass, Matrix>> pie = sn.pie;
    Matrix outspace = new Matrix(0,0);
    Matrix outrate = new Matrix(0,0);
    Matrix outprob = new Matrix(1,1);
    outprob.fill(1);


//    double isf = sn.nodeToStateful.get(ind);

      boolean ismkvmod = false;
      Matrix ismkvmodclass = new Matrix(0,0);
      if (sn.isstation.get(ind) == 1) {
          for (int i = 0; i < sn.proctype.get(sn.stations.get(ind)).size(); i++) {
              if (sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(i)) == ProcessType.MAP
                      || sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(i))
                      == ProcessType.MMPP2) {
                  ismkvmod = true;
              }
          }
          ismkvmodclass = new Matrix(R,1);
          ismkvmodclass.zero();
          for (int r = 0; r < R; r++) {
              if (sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(r)) == ProcessType.MAP
                      || sn.proctype.get(sn.stations.get(ind)).get(sn.jobclasses.get(r))
                      == ProcessType.MMPP2) {
                  ismkvmodclass.set(r, 0, 1);
              }
          }
      }


    Matrix lldscaling = sn.lldscaling;
    int lldlimit = 0;
    if (lldscaling.isEmpty()) {
      lldlimit = (int) Maths.max(sn.nclosedjobs, 1);
      lldscaling = new Matrix(M, lldlimit);
      lldscaling.ones();
    } else {
      lldlimit = lldscaling.getNumCols();
    }

    Map<Station, SerializableFunction<Matrix, Double>> cdscaling = sn.cdscaling;
    if (cdscaling == null || cdscaling.isEmpty()) {
      cdscaling = new HashMap<>();
      for (Station s : sn.stations) {
        // make a serializabl function that takes in a matrix and returns 1
        cdscaling.put(s, (Matrix x) -> 1.0);
      }
    }

    boolean hasOnlyExp = false; // true if all service processes are exponential
    int ist = (int) sn.nodeToStation.get(ind);
    Matrix K = Matrix.extractRows(phasessz, ist, ist + 1, null);
    Matrix Ks = Matrix.extractRows(phaseshift, ist, ist + 1, null);
    if (K.elementMax() == 1) { // ie no multi phase service, all are exponential
      hasOnlyExp = true;
    }
    Map<Station, Map<JobClass, Matrix>> mu = sn.mu;
    Map<Station, Map<JobClass, Matrix>> phi = sn.phi;

    Map<Station, Map<JobClass, Map<Integer, Matrix>>> proc = sn.proc;
    Matrix capacity = sn.cap;
    Matrix classcap = sn.classcap;


    double V = 0;

    // for a stateless node:
    Matrix spaceVar = new Matrix(0,0);
    Matrix spaceSrv = new Matrix(0,0);
    Matrix spaceBuf = new Matrix(0,0);


    if (sn.isstation.get(ind) == 1) {
      if (K.get(jobClass) == 0) {
        return new EventResult(outspace, outrate, outprob);
      }
      V = Matrix.extractRows(sn.nvars, ind, ind + 1, null).elementSum();
      int inspaceRows = inspace.getNumRows();

      // local state variables
      spaceVar = Matrix.extract(inspace, 0, inspaceRows, (int) (inspace.getNumCols() - V), inspace.getNumCols());

      spaceSrv = Matrix.extract(inspace, 0, inspaceRows, (int) (inspace.getNumCols() - K.elementSum() - V),
              (int) (inspace.getNumCols() - V)); // server state

      int spaceBufCols = (int) (inspace.getNumCols() - K.elementSum() - V);
      spaceBuf = Matrix.extract(inspace, 0, inspaceRows, 0, spaceBufCols); // buffer state


    } else if (sn.isstateful.get(ind, 0) == 1) {
      V = Matrix.extractRows(sn.nvars, ind, ind + 1, null).elementSum();
      int inspaceRows = inspace.getNumRows();

      int spaceVarCols = (int) V;
      spaceVar = new Matrix(inspaceRows, spaceVarCols);
      Matrix.extract(inspace, (int) (inspace.getNumCols() - V + 1), inspace.getNumCols(), 0,
              inspaceRows, spaceVar, 0,0); // local state variables

      int spaceSrvCols = R;
      spaceSrv = new Matrix(inspaceRows, spaceSrvCols);
      Matrix.extract(inspace, (int) (inspace.getNumCols() - R - V + 1), (int) (inspace.getNumCols() - V), 0,
              inspaceRows, spaceSrv, 0,0 ); // server state

    }
    if (sn.isstation.get(ind) == 1) {
      switch(event) {
        case ARV:
          // return if no space to accept the arrival, otherwise check scheduling strategy
          StateMarginalStatistics stats = toMarginalAggr(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
          Matrix ni = stats.ni;
          Matrix nir = stats.nir;
          Matrix pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
          outprob = new Matrix(0,0);
          Matrix outprobK = new Matrix(0, 0);
          for (int kentry = 0; kentry < K.get(jobClass); kentry++) {
            Matrix spaceVarK = spaceVar.clone();
            Matrix spaceSrvK = spaceSrv.clone();
            Matrix spaceBufK = spaceBuf.clone();
            switch (sn.sched.get(sn.stations.get(ist))) {
              case EXT: // source, can receive any "virtual" arrival from the sink as long as it is from an open class
                if (Double.isInfinite(sn.njobs.get(jobClass))) {
                  outspace = inspace.clone();
                  outrate = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                  outrate.zero();
                  outprob = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                  outprob.ones();
                }
                break;
              case PS:
              case INF:
              case DPS:
              case GPS:
                // due to nature of these policies, a new job enters service immediately.
                boolean exceedsClassCap = false;
                double istCap = classcap.get(ist, jobClass);
                int col = (int) (Ks.get(jobClass) + kentry);
                for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                  if (spaceSrvK.get(row, col) >= istCap) {
                    exceedsClassCap = true;
                  }
                }
                if (!exceedsClassCap) {
                  // increment spacesrvk by one
                  for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                    spaceSrvK.set(row, col, spaceSrvK.get(row, col) + 1);
                  }
                  outprobK = new Matrix(spaceSrvK.getNumRows(), spaceSrvK.getNumRows());
                  outprobK.fill(pentry.get(kentry));
                } else {
                  outprobK = new Matrix(spaceSrvK.getNumRows(), spaceSrvK.getNumRows());
                  outprobK.zero();
                }
                break;
              case SIRO:
              case SEPT:
              case LEPT:
                throw new RuntimeException("SIRO/SEPT/LEPT scheduling not supported");

              case FCFS:
              case HOL:
              case LCFS:
                // find states with all servers busy
                // if MAP service, when empty restart from the phase stored in spaceVar for this class
                // sn.nvars(ind,1:class)):

                Matrix nVarsAtInd = Matrix.extract(sn.nvars,  ind, ind + 1, 0, jobClass + 1);
                int sum = (int) nVarsAtInd.elementSum();
                if (ismkvmodclass.get(jobClass) == 0 || (ismkvmodclass.get(jobClass) == 1 && kentry == spaceVar.get(sum))) {
                  if (ismkvmodclass.get(jobClass) == 1) {
                    pentry.zero();
                    pentry.set(kentry, 1);
                  }
                  // construct all_busy_srv, a matrix with 0s where sum of that row in space_srv_k is >= S.get(ist) and 1s where its <
                  Matrix all_busy_srv = new Matrix(spaceSrvK.getNumRows(), 1);
                  for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                    Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                    int rowSum = (int) row.elementSum();
                    if (rowSum >= S.get(ist)) {
                      all_busy_srv.set(i, 0, 1);
                    } else {
                      all_busy_srv.set(i, 0, 0);
                    }
                  }

                  // find and modify states with an idle server
                  Matrix idle_srv = new Matrix(spaceSrvK.getNumRows(), 1);
                  for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                    Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                    int rowSum = (int) row.elementSum();
                    if (rowSum < S.get(ist)) {
                      idle_srv.set(i, 0, 1);
                    } else {
                      idle_srv.set(i, 0, 0);
                    }
                  }
                  // job enters service, increments idle server terms in space_srv_k
                  int colSrvK = (int) (spaceSrvK.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                  for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                    if (idle_srv.get(row, 0) == 1) {
                      spaceSrvK.set(row, colSrvK, spaceSrvK.get(row, colSrvK) + 1);
                    }
                  }
                  // this section dynamically grows the number of elements in the buffer
                  //  ni is an Mx1
                  // check if there is space in buffer
                  boolean spaceInBuffer = false;
                  for (int i = 0; i < ni.getNumRows(); i++) {
                    if (ni.get(i, 0) < capacity.get(ist)) {
                      spaceInBuffer = true;
                    }
                  }
                  if (spaceInBuffer) {
                    boolean spaceForClass = false;
                    for (int i = 0; i < nir.getNumRows(); i++) {
                      if (nir.get(i, jobClass) < classcap.get(ist, jobClass)) {
                        spaceForClass = true;
                      }
                    }
                    if (spaceForClass) {
                      // there is room for the job, check if buffer has empty slots. If not, append job slot
                      boolean emptySlots = false;
                      for (int i = 0; i < spaceBufK.getNumRows(); i++) {
                        if (spaceBufK.get(i, 0) == 0) {
                          emptySlots = true;
                        }
                      }
                      if (!emptySlots) {
                        // append job slot
                        Matrix left = new Matrix(spaceBufK.getNumRows(), 1);
                        left.zero();
                        spaceBufK = Matrix.concatColumns(left, spaceBufK, null);
                      }
                    }
                  }
                  // get position of first empty slot in buffer
                  Matrix empty_slots = new Matrix(all_busy_srv.getNumRows(), 1);
                  empty_slots.fill(-1);
                  int spaceBufKCols = spaceBufK.getNumCols();
                  if (spaceBufKCols == 0) {
                    //set empty_slots(all_busy_srv) = false;
                    for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                      if (all_busy_srv.get(i, 0) == 1) {
                        empty_slots.set(i, 0, 0);
                      }
                    }
                  } else if (spaceBufKCols == 1) {
                    // set empty_slots(all_busy_srv) = space_buf_k(all_busy_srv,:)==0;
                    for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                      if (all_busy_srv.get(i, 0) == 1) {
                        if (spaceBufK.get(i, 0) == 0) {
                          empty_slots.set(i, 0, 1);
                        } else {
                          empty_slots.set(i, 0, 0);
                        }
                      }
                    }

                  } else {
                    Matrix spaceBufKAllBusySrv = new Matrix(0, 0);
                    for (int i = 0; i < spaceBufK.getNumRows(); i++) {
                      if (all_busy_srv.get(i, 0) == 1) {
                        if (spaceBufKAllBusySrv.isEmpty()) {
                          spaceBufKAllBusySrv = Matrix.extractRows(spaceBufK, i, i + 1, null);
                        } else {
                          spaceBufKAllBusySrv = Matrix.concatRows(spaceBufKAllBusySrv, Matrix.extractRows(spaceBufK, i, i + 1, null), null);
                        }
                      }
                    }

                    for (int row = 0; row < spaceBufKAllBusySrv.getNumRows(); row++) {
                      for (int c = 0; c < spaceBufKAllBusySrv.getNumCols(); c++) {
                        if (spaceBufKAllBusySrv.get(row, c) == 0) {
                          spaceBufKAllBusySrv.set(row, c, 1);
                        } else {
                          spaceBufKAllBusySrv.set(row, c, 0);
                        }
                      }
                    }
                    Matrix sizeSpaceBufK = new Matrix(1, spaceBufK.getNumCols());
                    // TODO: check if this should be i instead of i+1
                    for (int i = 0; i < sizeSpaceBufK.getNumCols(); i++) {
                      sizeSpaceBufK.set(0, i, i+1);
                    }
                    if (!spaceBufKAllBusySrv.isEmpty() && !sizeSpaceBufK.isEmpty()) {
                      Matrix elementMult = spaceBufKAllBusySrv.elementMult(sizeSpaceBufK, null);
                      Matrix max = new Matrix(elementMult.getNumRows(), 1);
                      for (int i = 0; i < elementMult.getNumRows(); i++) {
                        max.set(i, 0, Matrix.extractRows(elementMult, i, i+1, null).elementMax());
                      }
                      int max_ind = 0;
                      for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                        if (all_busy_srv.get(i, 0) == 1) {
                          empty_slots.set(i, 0, max.get(max_ind));
                          max_ind++;
                        }
                      }
                    }
                  }
                  // ignore states where buffer has no empty slots
                  // set wbuf_empty = empty_slots >0
                  Matrix wbuf_empty = new Matrix(empty_slots.getNumRows(), 1);
                  boolean space_available = false;
                  for (int i = 0; i < empty_slots.getNumRows(); i++) {
                    if (empty_slots.get(i, 0) > 0) {
                      wbuf_empty.set(i, 0, 1);
                      space_available = true;
                    } else {
                      wbuf_empty.set(i, 0, 0);
                    }
                  }
                  if (space_available) {
                    // space_srv_k set to the rows of s[ace_srv_k where wbuf_empty = 1
                    Matrix spaceSrvKTmp = new Matrix(0, 0);
                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                      if (wbuf_empty.get(i, 0) == 1) {
                        if (spaceSrvKTmp.isEmpty()) {
                          spaceSrvKTmp = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                        } else {
                          spaceSrvKTmp = Matrix.concatRows(spaceSrvKTmp, Matrix.extractRows(spaceSrvK, i, i + 1, null), null);
                        }
                      }
                    }
                    spaceSrvK = spaceSrvKTmp;
                    Matrix spaceBufKTmp = new Matrix(0, 0);
                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                      if (wbuf_empty.get(i, 0) == 1) {
                        if (spaceBufKTmp.isEmpty()) {
                          spaceBufKTmp = Matrix.extractRows(spaceBufK, i, i + 1, null);
                        } else {
                          spaceBufKTmp = Matrix.concatRows(spaceBufKTmp, Matrix.extractRows(spaceBufK, i, i + 1, null), null);
                        }
                      }
                    }
                    spaceBufK = spaceBufKTmp;
                    Matrix spaceVarKTmp = new Matrix(0, 0);
                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                      if (wbuf_empty.get(i, 0) == 1) {
                        if (spaceVarKTmp.isEmpty()) {
                          spaceVarKTmp = Matrix.extractRows(spaceVarK, i, i + 1, null);
                        } else {
                          spaceVarKTmp = Matrix.concatRows(spaceVarKTmp, Matrix.extractRows(spaceVarK, i, i + 1, null), null);
                        }
                      }
                    }
                    spaceVarK = spaceVarKTmp;
                    Matrix emptySlotsTmp = new Matrix(0, 0);
                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                      if (wbuf_empty.get(i, 0) == 1) {
                        if (emptySlotsTmp.isEmpty()) {
                          emptySlotsTmp = Matrix.extractRows(empty_slots, i, i + 1, null);
                        } else {
                          emptySlotsTmp = Matrix.concatRows(emptySlotsTmp, Matrix.extractRows(empty_slots, i, i + 1, null), null);
                        }
                      }
                    }
                    empty_slots = emptySlotsTmp;
                    Matrix dims = new Matrix(1,2);
                    dims.set(0,0, spaceBufK.getNumRows());
                    dims.set(0,1, spaceBufK.getNumCols());

                    Matrix row_indices = new Matrix(1, spaceBufK.getNumRows());
                    for (int i = 0; i < row_indices.getNumCols(); i++) {
                      row_indices.set(0, i, i);
                    }

                    // need col indices to be = empty slots, but decrement as matrix 0 indexed
                    Matrix col_indices = empty_slots.clone();
                    for (int r = 0; r < col_indices.getNumRows(); r++) {
                      for (int c = 0; c < col_indices.getNumCols(); c++) {
                        col_indices.set(r, c, col_indices.get(r, c) - 1);
                      }
                    }
                    col_indices = col_indices.transpose();
                    List<Integer> indices = Maths.sub2ind(dims, row_indices, col_indices);
                    for (Integer n : indices) {
                      // jobClass + 1 since final form needs jobs to be 1 indexed
                      spaceBufK.set(n,  jobClass + 1);
                    }

                  }
                  outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                  outprobK.fill(pentry.get(kentry));
                } else {
                  outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                  outprobK.zero(); // zero probability event
                }
                break;
              case LCFSPR:
                // TODO: implement
                throw new RuntimeException("LCFSPR scheduling not supported");
            }
            // form the new state
            Matrix outspaceKTmp = Matrix.concatColumns(spaceBufK, spaceSrvK, null);
            Matrix outspaceK = Matrix.concatColumns(outspaceKTmp, spaceVarK, null);
            // remove states where new arrival violates capacity or cutoff constraints
            StateMarginalStatistics oi_oir = toMarginalAggr(sn, ind, outspaceK, K, Ks, spaceBufK, spaceSrvK, spaceVarK);
            Matrix oi = oi_oir.ni;
            Matrix oir = oi_oir.nir;

            Matrix en_o = new Matrix(oi.getNumRows(), 1);
            for (int row = 0; row < oi.getNumRows(); row++) {
              Matrix m = new Matrix(oi.getNumRows(), 1);
              m.fill(capacity.get(ist));
              boolean violates = false;
              for (int col = 0; col < oi.getNumCols(); col++) {
                if (m.get(row, 0) < oi.get(row, col)) {
                  violates = true;
                }
              }

              if (classcap.get(ist, jobClass) >= oir.get(row, jobClass) || !violates) {
                en_o.set(row, 0, 1);
              }
            }

            // need to extract all rows of outspace_k where en_o is true
            Matrix outspace_k_en_o = new Matrix(0, 0);
            for (int row = 0; row < en_o.getNumRows(); row++) {
              if (en_o.get(row, 0) == 1) {
                if (outspace_k_en_o.isEmpty()) {
                  outspace_k_en_o = Matrix.extractRows(outspaceK, row, row + 1, null);
                } else {
                  outspace_k_en_o = Matrix.concatRows(outspace_k_en_o, Matrix.extractRows(outspaceK, row, row + 1, null), null);
                }
              }
            }

            Matrix outprob_k_en_o = new Matrix(0, 0);
            for (int row = 0; row < en_o.getNumRows(); row++) {
              if (en_o.get(row, 0) == 1) {
                if (outprob_k_en_o.isEmpty()) {
                  outprob_k_en_o = Matrix.extractRows(outprobK, row, row + 1, null);
                } else {
                  outprob_k_en_o = Matrix.concatRows(outprob_k_en_o, Matrix.extractRows(outprobK, row, row + 1, null), null);
                }
              }
            }

            if (outspace.getNumCols() > outspace_k_en_o.getNumCols()) {
              Matrix zeros = new Matrix(1, outspace.getNumCols() - outspace_k_en_o.getNumCols());
              zeros.zero();
              Matrix bottom = Matrix.concatColumns(zeros, outspace_k_en_o, null);
              outspace = Matrix.concatRows(outspace, bottom, null);
            } else if (outspace.getNumCols() < outspace_k_en_o.getNumCols()) {
              Matrix zeros = new Matrix(outspace.getNumRows(), outspace_k_en_o.getNumCols() - outspace.getNumCols());
              zeros.zero();
              Matrix top = Matrix.concatColumns(zeros, outspace, null);
              outspace = Matrix.concatRows(top, outspace_k_en_o, null);
            } else {
              outspace = Matrix.concatRows(outspace, outspace_k_en_o, null);
            }
            Matrix newRates = new Matrix(outspace_k_en_o.getNumRows(), 1);
            newRates.fill(-1);
            outrate = Matrix.concatRows(outrate, newRates, null);
            outprob = Matrix.concatRows(outprob, outprob_k_en_o, null);
          }
          if (isSimulation) {
            if (outprob.getNumRows() > 1) {
              Matrix cum_sum = outprob.cumsumViaCol();
              Matrix sum_by_col = outprob.sumCols();
              Matrix cum_prob = Matrix.scale_mult(cum_sum, 1.0 / sum_by_col.get(0, 0));


              int firing_ctr = 0;
              double rand = Math.random();
              // we need the indicies where rand is bigger than cum_prob
              for (int row = 0; row < cum_prob.getNumRows(); row++) {
                if (rand > cum_prob.get(row, 0)) {
                  firing_ctr = row;
                }
              }
              outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
              outrate = new Matrix(1, 1);
              outrate.set(0, 0, -1);
              outprob = new Matrix(1, 1);
              outprob.set(0, 0, 1);
            }
          }
          break;
        case DEP:
          boolean busy = false;
          for (int row = 0; row < spaceSrv.getNumRows(); row++) {
            for (int col = (int) Ks.get(jobClass); col < (Ks.get(jobClass) + K.get(jobClass)); col++) {
              if (spaceSrv.get(row, col) > 0) {
                busy = true;
              }
            }
          }
          if (busy) {
            SchedStrategy strategy = sn.sched.get(sn.stations.get(ist));
            Matrix sir = new Matrix(0, 0);
            List<Matrix> kir = new ArrayList<>();
            if (hasOnlyExp && (strategy == PS || strategy == INF || strategy == DPS || strategy == GPS)) {
              nir = spaceSrv.clone();
              // set ni to sum of nir row-wise
              ni = nir.sumRows();
              sir = nir.clone();
              kir.add(sir.clone());
            } else {
              StateMarginalStatistics dep_stats = toMarginal(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
              ni = dep_stats.ni;
              nir = dep_stats.nir;
              sir = dep_stats.sir;
              kir = dep_stats.kir;
            }

            if (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass)) == RoutingStrategy.RROBIN) {
              // TODO: implement
            }

            if (sir.get(jobClass) > 0) {
              outprob = new Matrix(0, 0);
              for (int k = 0; k < K.get(jobClass); k++) {
                spaceSrv = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V)); // server state
                int spaceBufCols = (int) (inspace.getNumCols() - K.elementSum() - V);
                spaceBuf = Matrix.extract(inspace, 0, inspace.getNumRows(), 0, spaceBufCols); // buffer state
                Matrix rate = new Matrix(spaceSrv.getNumRows(), 1);
                rate.zero();
                Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                boolean en_set = false;
                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                  if (spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) > 0) {
                    en.set(row, 0, 1);
                    en_set = true;
                  } else {
                    en.set(row, 0, 0);
                  }
                }
                if (en_set) {
                  switch (sn.sched.get(sn.stations.get(ist))) {
                    case EXT:
                      // source, can produce an arrival from phase-k as long as it is from an open class
                      if (Double.isInfinite(sn.njobs.get(jobClass))) {
                        pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
                        for (int kentry = 0; kentry < K.get(jobClass); kentry++) {
                          Matrix.extract(inspace, (int) (inspace.getNumCols() - K.elementSum() - V + 1), (int) (inspace.getNumCols() - V), 0,
                                  inspace.getNumRows(), spaceSrv, 0, 0); // server state

                          // record a departure
                          for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                            }
                          }
                          // record a new job arriving
                          for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              spaceSrv.set(row, (int) (Ks.get(jobClass) + kentry), spaceSrv.get(row, (int) (Ks.get(jobClass) + kentry)) + 1);
                            }
                          }
                          // extract all rows of spaceBuf where en==1 rto space_buf_en
                          Matrix inspaceEn = new Matrix(0, 0);
                          Matrix spaceBufEn = new Matrix(0, 0);
                          Matrix spaceSrvEn = new Matrix(0, 0);
                          Matrix spaceVarEn = new Matrix(0, 0);
                          for (int row = 0; row < en.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              if (inspaceEn.isEmpty()) {
                                inspaceEn = Matrix.extractRows(inspace, row, row + 1, null);
                              } else {
                                inspaceEn = Matrix.concatRows(inspaceEn, Matrix.extractRows(inspace, row, row + 1, null), null);
                              }
                            }
                          }
                          for (int row = 0; row < en.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              if (spaceBufEn.isEmpty()) {
                                spaceBufEn = Matrix.extractRows(spaceBuf, row, row + 1, null);
                              } else {
                                spaceBufEn = Matrix.concatRows(spaceBufEn, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                              }

                            }
                          }
                          for (int row = 0; row < en.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              if (spaceSrvEn.isEmpty()) {
                                spaceSrvEn = Matrix.extractRows(spaceSrv, row, row + 1, null);
                              } else {
                                spaceSrvEn = Matrix.concatRows(spaceSrvEn, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                              }
                            }
                          }
                          for (int row = 0; row < en.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              if (spaceVarEn.isEmpty()) {
                                spaceVarEn = Matrix.extractRows(spaceVar, row, row + 1, null);
                              } else {
                                spaceVarEn = Matrix.concatRows(spaceVarEn, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                              }
                            }
                          }


                          Matrix left_bottom = Matrix.concatColumns(spaceBufEn, spaceSrvEn, null);
                          Matrix bottom = Matrix.concatColumns(left_bottom, spaceVarEn, null);
                          outspace = Matrix.concatRows(outspace, bottom, null);
                          if (isinf(ni)) {
                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                            Matrix outrate_bottom = new Matrix(inspace.getNumRows(), 1);
                            outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, lldlimit - 1) *
                                    pentry.get(kentry) * mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k)
                                    * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k));
                            Matrix ones = new Matrix(inspaceEn.getNumRows(), 1);
                            ones.ones();
                            outrate_bottom = outrate_bottom.mult(ones);
                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                          } else {
                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                            Matrix outrate_bottom = new Matrix(inspace.getNumRows(), 1);
                            // TODO: check usage of ni
                            outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit))
                                    * pentry.get(kentry) * mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k)
                                    * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k));
                            Matrix ones = new Matrix(inspaceEn.getNumRows(), 1);
                            ones.ones();
                            outrate_bottom = outrate_bottom.mult(ones);
                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                          }
                          Matrix outprob_bottom = new Matrix(spaceBufEn.getNumRows(), 1);
                          outprob_bottom.ones();
                          outprob = Matrix.concatRows(outprob, outprob_bottom, null);

                        }
                      }
                      break;
                    case INF:
                      // record a departure
                      for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        if (en.get(row) == 1) {
                          spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                        }
                      }

                      // TODO: check these 2 loops
                      // get kir(en, class, k)
                      Matrix kirEnClassK = new Matrix(0, 0);
                      for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                        if (en.get(l_ind, 0) == 1) {
                          if (kirEnClassK.isEmpty()) {
                              kirEnClassK = new Matrix(1,1);
                              kirEnClassK.set(0,0, kir.get(k).get(l_ind, jobClass));
                          } else {
                            Matrix new_elem = new Matrix(1,1);
                            new_elem.set(0,0, kir.get(k).get(l_ind, jobClass));
                            kirEnClassK = Matrix.concatRows(kirEnClassK, new_elem, null);
                          }
                        }
                      }
                      for (int l_ind = 0; l_ind < kirEnClassK.getNumRows(); l_ind++) {
                        if (rate.isEmpty()) {
                          rate = new Matrix(1, 1);
                          rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *
                                  phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *
                                  kirEnClassK.get(l_ind, 0));
                        } else {
                          Matrix new_elem = new Matrix(1, 1);
                          new_elem.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *
                                  phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *
                                  kirEnClassK.get(l_ind));
                          if (l_ind < rate.getNumElements()) {
                            // replacing an existing element in rate
                            rate.set(l_ind, new_elem.get(0,0));
                          } else {
                            // expand rate accordingly
                            rate = Matrix.concatRows(rate, new_elem, null);
                          }
                        }
                      }
                      // if state unchanged, add with rate 0
                      Matrix spaceBufEn = new Matrix(0, 0);
                      Matrix spaceSrvEn = new Matrix(0, 0);
                      for (int row = 0; row < en.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          if (spaceBufEn.isEmpty()) {
                            spaceBufEn = Matrix.extractRows(spaceBuf, row, row + 1, null);
                          } else {
                            spaceBufEn = Matrix.concatRows(spaceBufEn, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                          }
                          if (spaceSrvEn.isEmpty()) {
                            spaceSrvEn = Matrix.extractRows(spaceSrv, row, row + 1, null);
                          } else {
                            spaceSrvEn = Matrix.concatRows(spaceSrvEn, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                          }
                        }
                      }
                      Matrix space_var_last = Matrix.extractRows(spaceVar, spaceVar.getNumRows() - 1, spaceVar.getNumRows(), null);
                      Matrix left_bottom = Matrix.concatColumns(spaceBufEn, spaceSrvEn, null);
                      Matrix bottom = Matrix.concatColumns(left_bottom, space_var_last, null);
                      outspace = Matrix.concatRows(outspace, bottom, null);
                      Matrix rateEn = new Matrix(0, 0);
                      for (int row = 0; row < en.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          if (rateEn.isEmpty()) {
                            rateEn = Matrix.extractRows(rate, row, row + 1, null);
                          } else {
                            rateEn = Matrix.concatRows(rateEn, Matrix.extractRows(rate, row, row + 1, null), null);
                          }
                        }
                      }
                      if (isinf(ni)) {
                        // hit limited load-dependence
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                        Matrix outrate_bottom = Matrix.scale_mult(rateEn, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      } else {
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                        Matrix outrate_bottom = Matrix.scale_mult(rateEn, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      }
                      Matrix outprob_bottom = new Matrix(rateEn.getNumRows(), 1);
                      outprob_bottom.ones();
                      outprob = Matrix.concatRows(outprob, outprob_bottom, null);
                      break;

                    case PS:
                      // record departure
                      for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                        }
                      }

                      Matrix kirEnClassKPs = new Matrix(0, 0);
                      for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                        if (en.get(l_ind) == 1) {
                          if (kirEnClassKPs.isEmpty()) {
                            kirEnClassKPs = new Matrix(1,1);
                            kirEnClassKPs.set(0,0, kir.get(k).get(l_ind, jobClass));
                          } else {
                            Matrix new_elem = new Matrix(1,1);
                            new_elem.set(0,0, kir.get(k).get(l_ind, jobClass));
                            kirEnClassKPs = Matrix.concatRows(kirEnClassKPs, new_elem, null);
                          }
                        }
                      }

                      // assume active event
                      for (int l_ind = 0; l_ind < kirEnClassKPs.getNumRows(); l_ind++) {
                        if (rate.isEmpty()) {
                          rate = new Matrix(1, 1);
                          rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *
                                  phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) *(
                                  kirEnClassKPs.get(l_ind) / ni.get(l_ind)) * Maths.min(ni.get(l_ind), S.get(ist)));
                        } else {
                          Matrix new_elem = new Matrix(1, 1);
                          double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                          double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                          new_elem.set(0, 0, mu_value *
                                  phi_value *(
                                  kirEnClassKPs.get(l_ind) / ni.get(l_ind)) * Maths.min(ni.get(l_ind), S.get(ist)));
                          if (l_ind < rate.getNumElements()) {
                            // replacing an existing element in rate
                            rate.set(l_ind, new_elem.get(0,0));
                          } else {
                            // expand rate accordingly
                            rate = Matrix.concatRows(rate, new_elem, null);
                          }
                        }
                      }

                      Matrix spaceBufEnPs = new Matrix(0, 0);
                      Matrix spaceSrvEnPs = new Matrix(0, 0);
                      Matrix spaceVarEn = new Matrix(0, 0);
                      for (int row = 0; row < en.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          if (spaceBufEnPs.isEmpty()) {
                            spaceBufEnPs = Matrix.extractRows(spaceBuf, row, row + 1, null);
                          } else {
                            spaceBufEnPs = Matrix.concatRows(spaceBufEnPs, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                          }
                          if (spaceSrvEnPs.isEmpty()) {
                            spaceSrvEnPs = Matrix.extractRows(spaceSrv, row, row + 1, null);
                          } else {
                            spaceSrvEnPs = Matrix.concatRows(spaceSrvEnPs, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                          }
                          if (spaceVarEn.isEmpty()) {
                            spaceVarEn = Matrix.extractRows(spaceVar, row, row + 1, null);
                          } else {
                            spaceVarEn = Matrix.concatRows(spaceVarEn, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                          }
                        }
                      }

                      Matrix left_bottom_ps = Matrix.concatColumns(spaceBufEnPs, spaceSrvEnPs, null);
                      Matrix bottom_ps = Matrix.concatColumns(left_bottom_ps, spaceVarEn, null);
                      outspace = Matrix.concatRows(outspace, bottom_ps, null);
                      Matrix rateEnPs = new Matrix(0, 0);
                      for (int row = 0; row < en.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          if (rateEnPs.isEmpty()) {
                            rateEnPs = Matrix.extractRows(rate, row, row + 1, null);
                          } else {
                            rateEnPs = Matrix.concatRows(rateEnPs, Matrix.extractRows(rate, row, row + 1, null), null);
                          }
                        }
                      }

                      if (isinf(ni)) {
                        // hit limited load-dependence
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                        Matrix outrate_bottom = Matrix.scale_mult(rateEnPs, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      } else {
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                        Matrix outrate_bottom = Matrix.scale_mult(rateEnPs, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      }
                      Matrix outprob_bottom_ps = new Matrix(rateEnPs.getNumRows(), 1);
                      outprob_bottom_ps.ones();
                      outprob = Matrix.concatRows(outprob, outprob_bottom_ps, null);
                      break;

                    case FCFS:
                      for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        if (en.get(row, 0) == 1) {
                          spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                        }
                      }
                      // set en_wbuf to states with jobs in buffer
                      Matrix enWbuf = new Matrix(en.getNumRows(), 1);
                      for (int row = 0; row < en.getNumRows(); row++) {
                        if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                          enWbuf.set(row, 0, 1);
                        } else {
                          enWbuf.set(row, 0, 0);
                        }
                      }
                      for (int kdest = 0; kdest < K.get(jobClass); kdest++) {
                        Matrix space_buf_kd = spaceBuf.clone();
                        Matrix space_var_kd = spaceVar.clone();
                        if (ismkvmodclass.get(jobClass) == 1) {
                          // set space_var_kd(en, sum(sn.nvars(ind,1:class)) = kdest
                          Matrix nvar_ind = new Matrix(1, jobClass);
                          Matrix.extract(sn.nvars, ind, ind + 1, 0, jobClass, nvar_ind, 0, 0);
                          int nvar_sum = (int) nvar_ind.elementSum();
                          for (int row = 0; row < space_var_kd.getNumRows(); row++) {
                            if (en.get(row, 0) == 1) {
                              space_var_kd.set(row, (int) nvar_sum, kdest);
                            }
                          }
                        }
                        Matrix rate_kd = rate.clone();

                        // set rate_kd(en) = proc{ist}{class}{2}(k,kdest).*kir(en,class,k);


                        Matrix kirEnClassKFcfs = new Matrix(0, 0);
                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                          if (en.get(l_ind) == 1) {
                            if (kirEnClassKFcfs.isEmpty()) {
                              kirEnClassKFcfs = new Matrix(1,1);
                              kirEnClassKFcfs.set(0,0, kir.get(k).get(l_ind, jobClass));
                            } else {
                              Matrix new_elem = new Matrix(1,1);
                              new_elem.set(0,0, kir.get(k).get(l_ind, jobClass));
                              kirEnClassKFcfs = Matrix.concatRows(kirEnClassKFcfs, new_elem, null);
                            }
                          }
                        }


                        for (int row = 0; row < rate_kd.getNumRows(); row++) {
                          if (en.get(row, 0) == 1) {
                            double v = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(1).get(k, kdest) *
                                    kirEnClassKFcfs.get(row);
                            rate_kd.set(row, 0, v);
                          }
                        }
                        Matrix en_wobuf = new Matrix(enWbuf.getNumRows(), enWbuf.getNumCols());
                        // set all elems in en_wobuf to 1 where the same elem in enWBuf is 0 and vice versa
                        boolean anyStateNoJobs = false;
                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                          if (enWbuf.get(row, 0) == 0) {
                            en_wobuf.set(row, 0, 1);
                            anyStateNoJobs = true;
                          } else {
                            en_wobuf.set(row, 0, 0);
                          }
                        }
                        Matrix rate_kd_no_jobs = new Matrix(0, 0);
                        for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                          if (en_wobuf.get(row, 0) == 1) {
                            if (rate_kd_no_jobs.isEmpty()) {
                              rate_kd_no_jobs = Matrix.extractRows(rate_kd, row, row + 1, null);
                            } else {
                              rate_kd_no_jobs = Matrix.concatRows(rate_kd_no_jobs, Matrix.extractRows(rate_kd, row, row + 1, null), null);
                            }
                          }
                        }
                        if (anyStateNoJobs) {
                          // set outspace = [outspace; space_buf_kd(en_wobuf,:), space_srv(en_wobuf,:), space_var_kd(en_wobuf,:)];
                          Matrix space_buf_kd_no_jobs = new Matrix(0, 0);
                          for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                            if (en_wobuf.get(row, 0) == 1) {
                              if (space_buf_kd_no_jobs.isEmpty()) {
                                space_buf_kd_no_jobs = Matrix.extractRows(space_buf_kd, row, row + 1, null);
                              } else {
                                space_buf_kd_no_jobs = Matrix.concatRows(space_buf_kd_no_jobs, Matrix.extractRows(space_buf_kd, row, row + 1, null), null);
                              }
                            }
                          }
                          Matrix space_srv_no_jobs = new Matrix(0, 0);
                          for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                            if (en_wobuf.get(row, 0) == 1) {
                              if (space_srv_no_jobs.isEmpty()) {
                                space_srv_no_jobs = Matrix.extractRows(spaceSrv, row, row + 1, null);
                              } else {
                                space_srv_no_jobs = Matrix.concatRows(space_srv_no_jobs, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                              }
                            }
                          }
                          Matrix space_var_kd_no_jobs = new Matrix(0, 0);
                          for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                            if (en_wobuf.get(row, 0) == 1) {
                              if (space_var_kd_no_jobs.isEmpty()) {
                                space_var_kd_no_jobs = Matrix.extractRows(space_var_kd, row, row + 1, null);
                              } else {
                                space_var_kd_no_jobs = Matrix.concatRows(space_var_kd_no_jobs, Matrix.extractRows(space_var_kd, row, row + 1, null), null);
                              }
                            }
                          }

                          Matrix left_bottom_fcfs = Matrix.concatColumns(space_buf_kd_no_jobs, space_srv_no_jobs, null);
                          Matrix bottom_fcfs = Matrix.concatColumns(left_bottom_fcfs, space_var_kd_no_jobs, null);
                          outspace = Matrix.concatRows(outspace, bottom_fcfs, null);
                          // if all jobs (ni) are Infinite: hit limited load-dependence
                          if (isinf(ni)) {
                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                            // must multiply w rate_kd(en_wobuf, :)

                            Matrix outrate_bottom_fcfs = Matrix.scale_mult(rate_kd_no_jobs, cdscalingIst * lld);
                            outrate = Matrix.concatRows(outrate, outrate_bottom_fcfs, null);
                          } else {
                            // set outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_kd(en_wobuf,:)];
                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);

//                            Matrix lldscaling_ist = new Matrix(0, 0);
//                            for (int col = 0; col < lldscaling.getNumCols(); col++) {
//                              if (lldscaling.get(ist, col) == 1) {
//                                if (lldscaling_ist.isEmpty()) {
//                                  lldscaling_ist = Matrix.extractColumn(lldscaling, col, null);
//                                } else {
//                                  lldscaling_ist = Matrix.concatColumns(lldscaling_ist, Matrix.extractColumn(lldscaling, col, null), null);
//                                }
//                              }
//                            }
                            double lldscaling_ist = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit));

                            //perform element wise multiplication between cdscalingIst, lldscaling_ist and rate_kd(en_wobuf,:)

                            Matrix outrate_bottom_fcfs = Matrix.scale_mult(rate_kd_no_jobs, cdscalingIst * lldscaling_ist);
                            outrate = Matrix.concatRows(outrate, outrate_bottom_fcfs, null);
                          }
                        }
                          // now process states with jobs in buffer
                          Matrix outprob_bottom_fcfs = new Matrix(rate_kd_no_jobs.getNumRows(), 1);
                          outprob_bottom_fcfs.ones();
                          outprob = Matrix.concatRows(outprob, outprob_bottom_fcfs, null);
                          boolean any_jobs_in_buffer = false;
                          for (int row = 0; row < enWbuf.getNumRows(); row++) {
                            if (enWbuf.get(row, 0) == 1) {
                              any_jobs_in_buffer = true;
                            }
                          }
                          if (any_jobs_in_buffer) { // if there is any state with jobs in the buffer
                            // get class of job at head
                            Matrix space_buf_kd_last = Matrix.extractColumn(space_buf_kd, space_buf_kd.getNumCols() - 1, null);
                            // from space_buf_kd_last, extract all rows where en_wbuf = 1 into a matrix start_svc_class
                            Matrix start_svc_class = new Matrix(0, 0);
                            for (int row = 0; row < enWbuf.getNumRows(); row++) {
                              if (enWbuf.get(row, 0) == 1) {
                                if (start_svc_class.isEmpty()) {
                                  start_svc_class = Matrix.extractRows(space_buf_kd_last, row, row + 1, null);
                                } else {
                                  start_svc_class = Matrix.concatRows(start_svc_class, Matrix.extractRows(space_buf_kd_last, row, row + 1, null), null);
                                }
                              }
                            }
                            // if all elements of start_svc_lass are bigger than 0 set boolean x to true
                            boolean all_elems_bigger_than_zero = true;
                            for (int row = 0; row < start_svc_class.getNumRows(); row++) {
                              if (start_svc_class.get(row, 0) <= 0) {
                                all_elems_bigger_than_zero = false;
                              }
                            }
                            if (all_elems_bigger_than_zero) {
                              // update input buffer
                              // set space_buf_kd(en_wbuf,:) = [zeros(sum(en_wbuf),1),space_buf_kd(en_wbuf,1:end-1)];
                              Matrix left = new Matrix((int) enWbuf.elementSum(), 1);
                              left.zero();

                              Matrix space_buf_kd_end_removed = Matrix.extract(space_buf_kd, 0, space_buf_kd.getNumRows(), 0, space_buf_kd.getNumCols() - 1);

                              // extract into a matrix "right" all rows of space_buf_kd_end_removed where enWbuf = 1
                              Matrix right = new Matrix(0, 0);
                              for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                if (enWbuf.get(row, 0) == 1) {
                                  if (right.isEmpty()) {
                                    right = Matrix.extractRows(space_buf_kd_end_removed, row, row + 1, null);
                                  } else {
                                    right = Matrix.concatRows(right, Matrix.extractRows(space_buf_kd_end_removed, row, row + 1, null), null);
                                  }
                                }
                              }

                              Matrix new_space_buf_kd = Matrix.concatColumns(left, right, null);
                              // where en_wbuf = 1, take row of that index in "new_space_buf_kd" and assign it to that row in space_buf_kd
                              for (int row = 0; row < space_buf_kd.getNumRows(); row++) {
                                if (enWbuf.get(row) == 1) {
                                  // write the row-th row in new_space_buf_kd into this row
                                  for (int col = 0; col < space_buf_kd.getNumCols(); col++) {
                                    space_buf_kd.set(row, col, new_space_buf_kd.get(row, col));
                                  }
                                }
                              }


                              // TODO: use values in start_svc_class as indices into ismkvmodclass and see if all idxed values are true
                              boolean all_ismkvmodclass_true = false;
//                              for (int row = 0; row < start_svc_class.getNumRows(); row++) {
//                                if (ismkvmodclass.get((int) start_svc_class.get(row)) == 0) {
//                                  all_ismkvmodclass_true = false;
//                                }
//                              }
                              double kentry_range = 0;
                              Matrix pentry_svc_class = new Matrix(0,0);
                              if (all_ismkvmodclass_true) {
                                // TODO: markov-modulated case
                              } else {
                                pentry_svc_class = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) start_svc_class.get(0,0) - 1));
                                kentry_range = K.get((int) (start_svc_class.get(0,0) - 1)); // -1 since classes 1 indexed, matrix 0-indexed
                              }
                              for (int kentry = 0; kentry < kentry_range; kentry++) {
                                // increment all values in space_srv at rows where enWbuf is 1 at the column = Ks(start_svc_class.get(0,0)+kentry)
                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                  if (enWbuf.get(row, 0) == 1) {
                                    spaceSrv.set(row, (int) (Ks.get((int) start_svc_class.get(0,0) - 1) + kentry), spaceSrv.get(row, (int) (Ks.get((int) start_svc_class.get(0,0) - 1) + kentry)) + 1);
                                  }
                                }
                                // extract 3 matrices: space_buf_kd with only the rows where en istrue, space_srv with only the rows where en is true, and space_var_kd with only the rows where en is true
                                Matrix space_buf_kd_en = new Matrix(0, 0);
                                for (int row = 0; row < en.getNumRows(); row++) {
                                  if (en.get(row, 0) == 1) {
                                    if (space_buf_kd_en.isEmpty()) {
                                      space_buf_kd_en = Matrix.extractRows(space_buf_kd, row, row + 1, null);
                                    } else {
                                      space_buf_kd_en = Matrix.concatRows(space_buf_kd_en, Matrix.extractRows(space_buf_kd, row, row + 1, null), null);
                                    }
                                  }
                                }
                                Matrix space_srv_en = new Matrix(0, 0);
                                for (int row = 0; row < en.getNumRows(); row++) {
                                  if (en.get(row, 0) == 1) {
                                    if (space_srv_en.isEmpty()) {
                                      space_srv_en = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                    } else {
                                      space_srv_en = Matrix.concatRows(space_srv_en, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                    }
                                  }
                                }
                                Matrix space_var_kd_en = new Matrix(0, 0);
                                for (int row = 0; row < en.getNumRows(); row++) {
                                  if (en.get(row, 0) == 1) {
                                    if (space_var_kd_en.isEmpty()) {
                                      space_var_kd_en = Matrix.extractRows(space_var_kd, row, row + 1, null);
                                    } else {
                                      space_var_kd_en = Matrix.concatRows(space_var_kd_en, Matrix.extractRows(space_var_kd, row, row + 1, null), null);
                                    }
                                  }
                                }
                                Matrix left_bottom_outspace = Matrix.concatColumns(space_buf_kd_en, space_srv_en, null);
                                Matrix bottom_outspace = Matrix.concatColumns(left_bottom_outspace, space_var_kd_en, null);
                                outspace = Matrix.concatRows(outspace, bottom_outspace, null);

                                Matrix rate_k = rate_kd.clone();
                                // multiply each element in rate_k in rows (across all columns) where enWbuf is one by pentry_svc_class(kentry)
                                for (int row = 0; row < rate_k.getNumRows(); row++) {
                                  if (enWbuf.get(row, 0) == 1) {
                                    for (int col = 0; col < rate_k.getNumCols(); col++) {
                                      rate_k.set(row, col, rate_k.get(row, col) * pentry_svc_class.get(kentry));
                                    }
                                  }
                                }


                                Matrix rate_k_en = new Matrix(0, 0);
                                for (int row = 0; row < en.getNumRows(); row++) {
                                  if (en.get(row, 0) == 1) {
                                    if (rate_k_en.isEmpty()) {
                                      rate_k_en = Matrix.extractRows(rate_k, row, row + 1, null);
                                    } else {
                                      rate_k_en = Matrix.concatRows(rate_k_en, Matrix.extractRows(rate_k, row, row + 1, null), null);
                                    }
                                  }
                                }

                                if (isinf(ni)) {
                                  // hit limited load-dependence
                                  double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                  double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                  Matrix outrate_bottom = Matrix.scale_mult(rate_k_en, cdscalingIst * lld);
                                  outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                } else {
                                  double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                  double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                                  Matrix outrate_bottom = Matrix.scale_mult(rate_k_en, cdscalingIst * lld);
                                  outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                }
                                // extract rows of rate_kd using en as indices
                                Matrix rate_kd_en = new Matrix(0, 0);
                                for (int row = 0; row < en.getNumRows(); row++) {
                                  if (en.get(row, 0) == 1) {
                                    if (rate_kd_en.isEmpty()) {
                                      rate_kd_en = Matrix.extractRows(rate_kd, row, row + 1, null);
                                    } else {
                                      rate_kd_en = Matrix.concatRows(rate_kd_en, Matrix.extractRows(rate_kd, row, row + 1, null), null);
                                    }
                                  }
                                }
                                Matrix outprob_cur = new Matrix(rate_kd_en.getNumRows(), 1);
                                outprob_cur.ones();

                                // at the indices where outrate == 0, set outprob_cur to 0. make outprob_cur bigger if needed
                                for (int row = 0; row < outrate.getNumRows(); row++) {
                                  if (outrate.get(row) == 0) {
                                    // may need to expand outprob_cur. add a new row to outprob_cur with one column containing value 0
                                    if (row < outprob_cur.getNumElements()) {
                                      outprob_cur.set(row, 0);
                                    } else {
                                      // extending outprob_cur
                                      Matrix new_elem = new Matrix(1, 1);
                                      new_elem.set(0, 0, 0);
                                      outprob_cur = Matrix.concatColumns(outprob_cur, new_elem, null);
                                    }

                                  }
                                }



                                outprob = Matrix.concatRows(outprob, outprob_cur.transpose(), null);
                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                  if (enWbuf.get(row, 0) == 1) {
                                    spaceSrv.set(row, (int) (Ks.get((int) start_svc_class.get(0,0) - 1) + kentry), spaceSrv.get(row, (int) (Ks.get((int) start_svc_class.get(0,0) - 1) + kentry)) - 1);
                                  }
                                }
                              }
                            }
                          }
                        }
                      // if state unchanged still add with rate 0
                      break;

                    case DPS:
                      throw new RuntimeException("Scheduling strategy (DPS) not implemented");


                    case GPS:
                      throw new RuntimeException("Scheduling strategy (GPS) not implemented");
                    case HOL:
                      throw new RuntimeException("Scheduling strategy (HOL) not implemented");

                    case LCFS:
                      throw new RuntimeException("Scheduling strategy (LCFSP) not implemented");

                    case LCFSPR:
                      throw new RuntimeException("Scheduling strategy (LCFSPR) not implemented");

                    case SIRO:
                      throw new RuntimeException("Scheduling strategy (SIRO) not implemented");

                    case SEPT:
                    case LEPT:
                      throw new RuntimeException("Scheduling strategy (SEPT/LEPT) not implemented");
                    default:
                      throw new RuntimeException(String.format("Scheduling strategy %s is not supported",
                              sn.sched.get(sn.stations.get(ind)).toString()));


                  }


                }
              }
              if (isSimulation) {
                if (outspace.getNumRows() > 1) {
                  Matrix tot_rate = outrate.sumCols();
                  Matrix cum_sum = outrate.cumsumViaCol();
                  Matrix cum_rate = Matrix.scale_mult(cum_sum, (double) 1 / tot_rate.get(0,0));
                  int firing_ctr = 0;
                  double rand = Math.random();
                  // we need the indicies where rand is bigger than cum_prob
                  for (int row = 0; row < cum_rate.getNumRows(); row++) {
                    if (rand > cum_rate.get(row)) {
                      firing_ctr = row;
                    }
                  }
                  outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                  double outrate_val = outrate.elementSum();
                  outrate = new Matrix(1,1);
                  outrate.set(0, 0, outrate_val);
                  outprob = Matrix.extractRows(outprob, firing_ctr, firing_ctr + 1, null);
                }

              }
            }
          }
          break;
          case PHASE:
            outspace = new Matrix(0,0);
            outrate = new Matrix(0,0);
            outprob = new Matrix(0,0);
            StateMarginalStatistics stateMarginalStatistics = toMarginal(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
            ni = stateMarginalStatistics.ni;
            nir = stateMarginalStatistics.nir;
            List<Matrix> kir = stateMarginalStatistics.kir;

            if (nir.get(jobClass) > 0) {
              for (int k = 0; k < K.get(jobClass); k++) {
                // set en = space_srv(:,Ks(class)+k) > 0;
                // set en to a matrix which has a 1 if that row in space_srv in column Ks(class) + k is bigger than 0, and a 0 if not
                Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                en.zero();
                boolean any_en = false;
                for (int row = 0; row < en.getNumRows(); row++) {
                  if (spaceSrv.get(row, (int) Ks.get(jobClass) + k) > 0) {
                    en.set(row, 0, 1);
                    any_en = true;
                  }
                }

                if (any_en) {
                  for (int kdest = 0; kdest < K.get(jobClass); kdest++) {
                    if (kdest != k) {
                      Matrix rate = new Matrix(1,1);

                      Matrix spaceSrvK = new Matrix(0, 0);
                      for (int i = 0; i < en.getNumRows(); i++) {
                        if (en.get(i, 0) == 1) {
                          if (spaceSrvK.isEmpty()) {
                            spaceSrvK = Matrix.extractRows(spaceSrv, i, i + 1, null);
                          } else {
                            spaceSrvK = Matrix.concatRows(spaceSrvK, Matrix.extractRows(spaceSrv, i, i + 1, null), null);
                          }
                        }
                      }

                      Matrix spaceBufK = new Matrix(0, 0);
                      for (int i = 0; i < en.getNumRows(); i++) {
                        if (en.get(i, 0) == 1) {
                          if (spaceBufK.isEmpty()) {
                            spaceBufK = Matrix.extractRows(spaceBuf, i, i + 1, null);
                          } else {
                            spaceBufK = Matrix.concatRows(spaceBufK, Matrix.extractRows(spaceBuf, i, i + 1, null), null);
                          }
                        }
                      }

                      Matrix spaceVarK = new Matrix(0, 0);
                      for (int i = 0; i < en.getNumRows(); i++) {
                        if (en.get(i, 0) == 1) {
                          if (spaceVarK.isEmpty()) {
                            spaceVarK = Matrix.extractRows(spaceVar, i, i + 1, null);
                          } else {
                            spaceVarK = Matrix.concatRows(spaceVarK, Matrix.extractRows(spaceVar, i, i + 1, null), null);
                          }
                        }
                      }

                      // TODO: markov-modulated case: lines 712-714

                      for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                        spaceSrvK.set(row, (int) (Ks.get(jobClass) + k), spaceSrvK.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                      }
                      for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                        spaceSrvK.set(row, (int) (Ks.get(jobClass) + kdest), spaceSrvK.get(row, (int) (Ks.get(jobClass) + kdest)) + 1);
                      }

                      switch (sn.sched.get(sn.stations.get(ist))) {
                        case EXT:
                          rate.set(0,0,(int) proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest));
                          break;
                        case INF:
                          double proc_value_inf = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                          double kir_value_inf = kir.get(k).get(jobClass);
                          rate.set(0,0, proc_value_inf * kir_value_inf);
                          break;
                        case PS:
                          double proc_value_ps = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                          double kir_value_ps = kir.get(k).get(jobClass);
                          Matrix numerator = new Matrix(1,1);
                          numerator.set(0,0,proc_value_ps * kir_value_ps);
                          Matrix denom = new Matrix(1,1);
                          double ni_value = ni.get(0);
                          denom.set(0,0, ni_value * Maths.min(ni_value, S.get(ist)));
                          rate = numerator.element_divide(denom);
                          break;
                        case DPS:
                          throw new RuntimeException("Scheduling strategy (DPS) not implemented");
                        case GPS:
                          throw new RuntimeException("Scheduling strategy (GPS) not implemented");

                        case FCFS:
                        case HOL:
                        case LCFS:
                        case LCFSPR:
                        case SIRO:
                        case SEPT:
                        case LEPT:
                             double proc_value = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                             double kir_value = kir.get(k).get(jobClass);
                             rate.set(0,0, proc_value * kir_value);
                             break;
                      }

                      // if class cannot be served locally, rate = NaN since mu{i, class} = NaN
                      // TODO: lines 744-750
                      if (isinf(ni)) {
                        // hit limited load-dependence
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                        Matrix outrate_bottom = Matrix.scale_mult(rate, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      } else {
                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                        Matrix outrate_bottom = Matrix.scale_mult(rate, cdscalingIst * lld);
                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                      }
                      Matrix outspace_bottom_left = Matrix.concatColumns(spaceBufK, spaceSrvK, null);
                      Matrix bottom = Matrix.concatColumns(outspace_bottom_left, spaceVarK, null);
                      outspace = Matrix.concatRows(outspace, bottom, null);
                      Matrix outprob_bottom = new Matrix(rate.getNumRows(), 1);
                      outprob_bottom.ones();
                      outprob = Matrix.concatRows(outprob, outprob_bottom, null);
                    }
                  }
                }
              }
              if (isSimulation) {
                if (outspace.getNumRows() > 1) {
                  Matrix tot_rate = outrate.sumCols();
                  Matrix cum_sum = outrate.cumsumViaCol();
                  Matrix cum_rate = Matrix.scale_mult(cum_sum, (double) 1 / tot_rate.get(0,0));
                  int firing_ctr = 0;
                  double rand = Math.random();
                  // we need the indicies where rand is bigger than cum_prob
                  for (int row = 0; row < cum_rate.getNumRows(); row++) {
                    if (rand > cum_rate.get(row)) {
                      firing_ctr = row;
                    }
                  }
                  outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                  double outrate_val = outrate.elementSum();
                  outrate.set(0, 0, outrate_val);
                  outprob = Matrix.extractRows(outprob, firing_ctr, firing_ctr + 1, null);
                }

              }

            }
            break;
          }



    } else if (sn.isstateful.get(ind) == 1) {

    }

    return new EventResult(outspace, outrate, outprob);

  }


  public static boolean isinf(Matrix matrix) {
    for (int row = 0; row < matrix.getNumRows(); row++) {
      for (int col = 0; col < matrix.getNumCols(); col++) {
        if (!Double.isInfinite(matrix.get(row, col))) {
          return false;
        }
      }
    }
    return true;
  }














}


