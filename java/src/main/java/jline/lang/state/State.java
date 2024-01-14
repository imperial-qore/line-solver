package jline.lang.state;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import jline.examples.ClosedModel;
import jline.examples.OpenModel;
import jline.examples.RandomEnvironment;
import jline.lang.Network;
import jline.lang.nodes.StatefulNode;
import jline.util.Maths;
import jline.util.Matrix;
import jline.lang.constant.ProcessType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;

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
            sumval = space_buf.countEachRow(r);
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


  public static Matrix fromMarginalAndRunning(NetworkStruct sn, int ind, Matrix n, Matrix s) {
    return fromMarginalAndStarted(sn, ind, n, s, true);
  }

  public static Matrix fromMarginalAndRunning(Network sn, int ind, Matrix n, Matrix s) {
    return fromMarginalAndStarted(sn.getStruct(true), ind, n, s, true);
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
          break;

        case LCFSPR:
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
    // TODO: check from this point (line 33)
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
    // TODO: down to this point, line 46

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
          // TODO: check from here
          if (n.sumSubMatrix(0, n.getNumRows(), 0, n.getNumCols()) == 0) {
            space =
                new Matrix(
                    1, 1 + Math.max((int) K.sumSubMatrix(0, K.getNumRows(), 0, K.getNumCols()), R));
            space.zero();
            return space;
          }

          // In these policies we track an ordered buffer and the jobs in the servers
          // build list of job classes in the buffer, with repetition
          Matrix inbuf = new Matrix(0, 0);
          for (int r = 0; r < R; r++) {
            if (n.get(0, r) > 0) {
              Matrix newInbuf =
                  new Matrix(1, inbuf.getNumCols() + (int) n.get(0, r) - (int) s.get(0, r));
              for (int i = 0; i < inbuf.getNumCols(); i++) {
                newInbuf.set(0, i, inbuf.get(0, i));
              }
              for (int i = inbuf.getNumCols(); i < newInbuf.getNumCols(); i++) {
                newInbuf.set(0, i, r);
              }
              inbuf = newInbuf.clone();
            }
          }

          double sizeEstimator = 1; // TODO: sizeEstimator = multinomialln(n);
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
          Matrix mi = Maths.uniquePerms(inbuf);
          Matrix mi_buf;
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
              mi_buf = new Matrix(1, 1, 1);
              mi_buf.set(0, 0, 0);
            }
          }

          // mi_srv: class of jobs running in the server of i
          Matrix mi_srv = new Matrix(0, 0);
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
              Matrix newState =
                  new Matrix(
                      1, state.getNumCols() + mi_buf.getNumCols() + kState.getNumCols());
              for (int row = 0; row < mi_buf.getNumRows(); row++) {
                for (int col = 0; col < state.getNumCols(); col++) {
                  newState.set(row, col, state.get(row, col));
                }
                for (int col = state.getNumCols();
                    col < state.getNumCols() + mi_buf.getNumCols();
                    col++) {
                  newState.set(row, col, mi_buf.get(b, col - state.getNumCols()));
                }
                for (int col = state.getNumCols() + mi_buf.getNumCols();
                    col < newState.getNumCols();
                    col++) {
                  newState.set(
                      row, col, kState.get(row, col - state.getNumCols() - mi_buf.getNumCols()));
                }
              }
              state = newState.clone();
            }
          }

          space = state.clone();
          // TODO: check down to here
          break;

        case LCFSPR:
          // TODO: finish remainder of implementation, lines 167 to 245
          System.out.println(
              "Warning: unimplemented code reached in fromMarginalAndStarted  for LCFSPR strategies");
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
            || sn.nodetypes.get(ind) == NodeType.Place  ) {
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

    public static void main(String[] args) {
//      NetworkStruct sn = ClosedModel.ex3().getStruct(true);
//      Matrix m = new Matrix(1, 3);
//      m.fromArray2D(new int[][]{{1,1,1}});
//
//      Matrix s = new Matrix(1, 3);
//      s.fromArray2D(new int[][]{{1,1,0}});
//      Matrix res = fromMarginalAndStarted(sn, 1, m, s);
//      System.out.println(res.getNumRows());
//      System.out.println(res);

//      NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//      Matrix m = new Matrix(1, 4);
//      m.fromArray2D(new int[][]{{2,1,1,1}});
//
//      Matrix s = new Matrix(1,4);
//      s.fromArray2D(new int[][]{{2,0,1,0}});
//      Matrix res = fromMarginal(sn, 1, m);
//      System.out.println(res.getNumRows());
//      System.out.println(res);

//
//      NetworkStruct sn = ClosedModel.ex6().getStruct(false);
//      Matrix m = new Matrix(1, 2);
//      m.fromArray2D(new int[][]{{5,6}});
//
//      Matrix s = new Matrix(1, 2);
//      s.fromArray2D(new int[][]{{1, 0}});
//
//      Matrix res = fromMarginalAndStarted(sn, 1, m, s);
//      System.out.println(res.getNumRows());
//      System.out.println(res);

//        NetworkStruct sn = ClosedModel.ex7_fcfs().getStruct(false);
//        Matrix m = new Matrix(1, 2);
//        m.fromArray2D(new int[][]{{0,2}});
//        Matrix res = fromMarginal(sn, 1, m);
//        System.out.println(res.getNumRows());
//        System.out.println(res);

//
//        NetworkStruct sn = ClosedModel.ex9().getStruct(true);
//        Matrix m = new Matrix(1, 1);
//        m.fromArray2D(new int[][]{{1}});
//        Matrix res = fromMarginal(sn, 1, m);
//        System.out.println(res.getNumRows());
//        System.out.println(res);

//      Matrix m = new Matrix(3,1);
//
//      m.fromArray2D(new double[][]{{Double.POSITIVE_INFINITY}, {3},{3}});
//      System.out.println(m);
//
//      NetworkStruct sn = ClosedModel.ex4().getStruct(false);
//      System.out.println(sn.nservers);
//
//      Matrix ones = new Matrix(3, 1);
//      ones.ones();
//      Matrix infBuffer = ones.mult(new Matrix(1,1).fromArray2D(new double[][]{{Double.POSITIVE_INFINITY}}));
//      System.out.println(infBuffer);

//      NetworkStruct sn = OpenModel.ex4().getStruct(true);
//      Matrix m = new Matrix(1, 1);
//      m.fromArray2D(new int[][]{{1}});
//      Matrix res = fromMarginal(sn, 1, m);
//      System.out.println(res.getNumRows());
//      System.out.println(res);

//        NetworkStruct sn = ClosedModel.ex7_lcfspr().getStruct(false);
//        Matrix m = new Matrix(1, 2);
//        m.fromArray2D(new int[][]{{3,3}});
//        Matrix res = fromMarginal(sn, 1, m);
//        System.out.println(res.getNumRows());
//        System.out.println(res);


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
          case LCFSPR:
            Matrix vi = new Matrix(0,0);
            Matrix mi = new Matrix(0,0);
            if (sn.sched.get(sn.stations.get(ist)) == LCFSPR) {

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
              // now generate server states
              if (mi.isEmpty()) {
                Matrix mi_buf =new Matrix(1, (int) Maths.max(0, n.elementSum() - S.get(ist)));
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
                // mi_buf: class of job in buffer position i (0 = empty)
                int numColumnsRight = (int) Maths.max((mi.getNumCols() - S.get(ist)), 0);
                Matrix right = new Matrix(mi.getNumRows(), numColumnsRight);
                Matrix.extract(mi, 0, mi.getNumRows(), 0,
                        numColumnsRight, right, 0, 0);
                int numColumnsLeft = (int) Maths.max(0, (Maths.min(n.elementSum(), sn.cap.get(ist)) - S.get(ist)
                        - right.getNumCols()));
                Matrix left = new Matrix(mi.getNumRows(), numColumnsLeft);
                Matrix mi_buf = left.concatCols(right);
                if (mi_buf.isEmpty()) {
                  mi_buf = new Matrix(mi.getNumRows(), 1);
                  mi_buf.zero();
                }
                Matrix miBufKState = new Matrix(0, 0);
                // miSrv: class of job running in server i
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

                  Matrix jobsInBuffer = Matrix.extractRows(mi_buf, k, k + 1, null);
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
                  Matrix bufStateTmp = Matrix.decorate(Matrix.extractRows(mi_buf, k, k + 1, null), bkState);
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
                  }
                  if (state.isEmpty()) {
                    state = Matrix.decorate(bufState, kState);
                  } else {
                    state = Matrix.concatRows(state, Matrix.decorate(bufState, kState), null);
                  }
                }
              }
            } else {
              // FCFS, HOL, LCFS case
                // sum(n) - 1 due to Maths.factln including + 1
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

    Matrix newSpace = new Matrix(uniqueRows.size(), uniqueRows.get(0).getNumCols());
    // So that states with jobs in phase 1 comes earlier
    int row = 0;
    for (int i = uniqueRows.size() - 1; i >= 0; i--) {
      for (int j = 0; j < uniqueRows.get(0).getNumCols(); j++) {
        newSpace.set(row, j, uniqueRows.get(i).get(0, j));
      }
      row++;
    }

    return newSpace;
  }
}
