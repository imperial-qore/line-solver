package jline.api;

import java.util.*;


import jline.lang.constant.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.util.Matrix;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.util.Maths;

/**
 * APIs for macros to process NetworkStruct objects
 */
public class SN {

    /**
     * Reculate the visits to each node
     *
     * @param sn      - NetworkStruct object for the queueing network model. The object will be updated by the method.
     * @param chains  - chain membership for each class
     * @param rt      - routing table among stations
     * @param rtnodes - routing table among nodes
     * @return updated sn structure
     */
    // TODO: clone sn at start and do not change
    public static NetworkStruct snRefreshVisits(NetworkStruct sn, Matrix chains, Matrix rt, Matrix rtnodes) {
        int I = sn.nnodes;
        int M = sn.nstateful;
        int K = sn.nclasses;
        Matrix refstat = sn.refstat;
        int nchains = sn.nchains;

        /* Obtain chain characteristics */
        Map<Integer, Matrix> inchain = sn.inchain;
        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            double val = refstat.get((int) inchain_c.get(0, 0), 0);
            for (int col = 1; col < inchain_c.getNumCols(); col++) {
                int row = (int) inchain_c.get(0, col);
                if (val != refstat.get(row, 0))
                    refstat.set(row, 0, val);
                //throw new RuntimeException("Classes within chain have different reference station");
            }
        }

        /* Transfer inchain to List<Integer> in order to reduce the time of type conversion (double -> int) which is time consuming) */
        Map<Integer, List<Integer>> new_inchain = new HashMap<Integer, List<Integer>>();
        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            List<Integer> inchain_c_list = new ArrayList<Integer>();
            for (int i = 0; i < inchain_c.getNumCols(); i++)
                inchain_c_list.add((int) inchain_c.get(i));
            new_inchain.put(c, inchain_c_list);
        }

        /* Generate visits */
        Map<Integer, Matrix> visits = new HashMap<Integer, Matrix>();
        for (int c = 0; c < nchains; c++) {
            List<Integer> inchain_c = new_inchain.get(c);
            List<Integer> cols = new ArrayList<Integer>();    //If use JLineMatrix, there would be more data type transfer in Pchain creation
            for (int i = 0; i < M; i++) {
                for (int ik = 0; ik < inchain_c.size(); ik++) {
                    cols.add(i * K + inchain_c.get(ik));
                }
            }

            //Pchain = rt(cols,cols);
            Matrix Pchain = new Matrix(cols.size(), cols.size());
            for (int row = 0; row < cols.size(); row++) {
                for (int col = 0; col < cols.size(); col++) {
                    Pchain.set(row, col, rt.get(cols.get(row), cols.get(col)));
                }
            }


            //visited = sum(Pchain,2) > 0;
            Matrix visited = new Matrix(Pchain.getNumRows(), 1);
            int countTrue = 0;
            for (int row = 0; row < Pchain.getNumRows(); row++) {
                if (Pchain.sumRows(row) > 0) {
                    countTrue++;
                    visited.set(row, 0, 1.0);
                }
            }

            //alpha_visited = dtmc_solve(Pchain(visited,visited));
            Matrix input = new Matrix(countTrue, countTrue);
            int row_input = 0, col_input = 0;
            for (int row = 0; row < visited.getNumRows(); row++) {
                if (visited.get(row, 0) > 0) {
                    for (int col = 0; col < visited.getNumRows(); col++) {
                        if (row == col || visited.get(col, 0) > 0)
                            input.set(row_input, (col_input++) % countTrue, Pchain.get(row, col));
                    }
                    row_input++;
                }
            }

            Matrix alpha_visited = DTMC.dtmc_solve(input);

            //alpha = zeros(1,M*K); alpha(visited) = alpha_visited;
            Matrix alpha = new Matrix(1, M * K);
            int idx = 0;
            for (int row = 0; row < visited.getNumRows(); row++) {
                if (visited.get(row, 0) > 0)
                    alpha.set(0, row, alpha_visited.get(0, idx++));
            }

            Matrix visits_c = new Matrix(M, K);
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < inchain_c.size(); k++) {
                    visits_c.set(i, inchain_c.get(k), alpha.get(0, i * inchain_c.size() + k));
                }
            }

            //visits{c} = visits{c} / sum(visits{c}(refstat(inchain{c}(1)),inchain{c}));
            double sum = 0;
            int row = (int) sn.stationToStateful.get((int) refstat.get(inchain_c.get(0)));
            for (int i = 0; i < inchain_c.size(); i++) {
                sum += visits_c.get(row, inchain_c.get(i));
            }
            Matrix visits_c_divide = new Matrix(0, 0);
            visits_c.divide(sum, visits_c_divide, true);

            visits_c_divide.abs();
            visits.put(c, visits_c_divide);
        }

        /* Generate node visits */
        Map<Integer, Matrix> nodeVisits = new HashMap<Integer, Matrix>();
        for (int c = 0; c < nchains; c++) {
            List<Integer> inchain_c = new_inchain.get(c);
            List<Integer> nodes_cols = new ArrayList<Integer>();    //If use JLineMatrix, there would be more data type transfer in Pchain creation
            for (int i = 0; i < I; i++) {
                for (int ik = 0; ik < inchain_c.size(); ik++) {
                    nodes_cols.add(i * K + inchain_c.get(ik));
                }
            }

            Matrix nodes_Pchain = new Matrix(nodes_cols.size(), nodes_cols.size());
            for (int row = 0; row < nodes_cols.size(); row++) {
                for (int col = 0; col < nodes_cols.size(); col++) {
                    nodes_Pchain.set(row, col, rtnodes.get(nodes_cols.get(row), nodes_cols.get(col)));
                }
            }

            Matrix nodes_visited = new Matrix(nodes_Pchain.getNumRows(), 1);
            int countTrue = 0;
            for (int row = 0; row < nodes_Pchain.getNumRows(); row++) {
                if (nodes_Pchain.sumRows(row) > 0) {
                    countTrue++;
                    nodes_visited.set(row, 0, 1.0);
                }
            }

            Matrix input = new Matrix(countTrue, countTrue);
            int row_input = 0, col_input = 0;
            for (int row = 0; row < nodes_visited.getNumRows(); row++) {
                if (nodes_visited.get(row, 0) > 0) {
                    for (int col = 0; col < nodes_visited.getNumRows(); col++) {
                        if (row == col || nodes_visited.get(col, 0) > 0)
                            input.set(row_input, (col_input++) % countTrue, nodes_Pchain.get(row, col));
                    }
                    row_input++;
                }
            }
            Matrix nodes_alpha_visited = DTMC.dtmc_solve(input);

            Matrix nodes_alpha = new Matrix(1, I * K);
            int idx = 0;
            for (int row = 0; row < nodes_visited.getNumRows(); row++) {
                if (nodes_visited.get(row, 0) > 0)
                    nodes_alpha.set(0, row, nodes_alpha_visited.get(0, idx++));
            }

            Matrix node_visits_c = new Matrix(I, K);
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < inchain_c.size(); k++) {
                    node_visits_c.set(i, inchain_c.get(k), nodes_alpha.get(0, i * inchain_c.size() + k));
                }
            }

            double sum = 0;
            int row = (int) refstat.get(inchain_c.get(0));
            for (int i = 0; i < inchain_c.size(); i++) {
                sum += node_visits_c.get(row, inchain_c.get(i));
            }
            Matrix node_visits_c_divide = new Matrix(0, 0);
            node_visits_c.divide(sum, node_visits_c_divide, true);

            node_visits_c_divide.removeNegative();
            node_visits_c_divide.removeNaN();
            nodeVisits.put(c, node_visits_c_divide);
        }

        /* Save result in sn */
        sn.visits = visits;
        sn.nodevisits = nodeVisits;
        sn.isslc = new Matrix(sn.nclasses, 1);

        return sn;
    }

    /**
     * Calculate the parameters at class and chain level for a queueing network model
     *
     * @param sn - NetworkStruct object for the queueing network model.
     * @return queueing network parameters
     */
    // TODO: the method seems not to return the chain parameters and to use the wrong return class
    public static snGetProductFormChainParamsReturn snGetProductFormChainParams(NetworkStruct sn) {
        snGetProductFormChainParamsReturn ret1 = snGetProductFormParams(sn);
        Matrix lambda = ret1.lambda;
        Matrix mu = ret1.mu;
        ArrayList<Integer> queueIndex, delayIndex;
        HashSet<Integer> ignoreIndex;
        queueIndex = new ArrayList<>();
        delayIndex = new ArrayList<>();
        ignoreIndex = new HashSet<>();
        for (int i = 0; i < sn.nodetypes.size(); i++) {
            switch (sn.nodetypes.get(i)) {
                case Queue:
                    queueIndex.add(i);
                    break;
                case Delay:
                    delayIndex.add(i);
                    break;
                case Source:
                case Join:
                    ignoreIndex.add((int) sn.nodeToStation.get(i));
                    break;
            }
        }
        snGetDemandsChainReturn ret2 = snGetDemandsChain(sn);
        Matrix Dchain = ret2.Lchain;
        Matrix Vchain = ret2.Vchain;
        Matrix Nchain = ret2.Nchain;
        Matrix lambda_chains = new Matrix(1, sn.nchains);

        Matrix D_chains = new Matrix(queueIndex.size(), sn.nchains);
        Matrix Z_chains = new Matrix(delayIndex.size(), sn.nchains);

        for (int c = 0; c < sn.nchains; c++) {
            double lambdaSum = 0;
            Matrix inc = sn.inchain.get(c);
            for (int i = 0; i < inc.getNumCols(); i++) {
                if (!Double.isNaN(lambda.get((int) inc.get(i)))) {
                    lambdaSum += lambda.get((int) inc.get(i));
                }
            }
            lambda_chains.set(0, c, lambdaSum);
            for (int i = 0; i < queueIndex.size(); i++) {
                D_chains.set(i, c, Dchain.get((int) sn.nodeToStation.get(queueIndex.get(i)), c));
            }
            for (int i = 0; i < delayIndex.size(); i++) {
                Z_chains.set(i, c, Dchain.get((int) sn.nodeToStation.get(delayIndex.get(i)), c));
            }
        }
        Matrix S = new Matrix(queueIndex.size(), 1);
        for (int i = 0; i < queueIndex.size(); i++) {
            S.set(i, 0, sn.nservers.get((int) sn.nodeToStation.get(queueIndex.get(i))));
        }
        snGetProductFormChainParamsReturn ret = new snGetProductFormChainParamsReturn();
        ret.S = S;
        ret.lambda = lambda_chains;
        ret.N = Nchain;
        Vchain.removeRows(ignoreIndex);
        ret.D = D_chains;
        ret.Z = Z_chains;
        ret.mu = mu;
        ret.V = Vchain;
        if (ret.Z.isEmpty()) {
            ret.Z = new Matrix(ret.N.getNumRows(), ret.N.getNumCols());
        }
        return ret;
    }

    /**
     * Calculate the parameters at class level for a queueing network model
     *
     * @param sn - NetworkStruct object for the queueing network model.
     * @return queueing network parameters
     */
    // TODO: the method seems not to use the wrong name for return class
    public static snGetProductFormChainParamsReturn snGetProductFormParams(NetworkStruct sn) {
        int R = sn.nclasses;
        Matrix N = sn.njobs;
        ArrayList<Integer> queueIndices = new ArrayList<>();
        ArrayList<Integer> delayIndices = new ArrayList<>();
        int sourceIndex = -1;
        for (int i = 0; i < sn.nodetypes.size(); i++) {
            if (sn.nodetypes.get(i) == NodeType.Queue) {
                queueIndices.add(i);
            } else if (sn.nodetypes.get(i) == NodeType.Delay) {
                delayIndices.add(i);
            } else if (sn.nodetypes.get(i) == NodeType.Source) {
                sourceIndex = i;
            }
        }
        int Mq = queueIndices.size();
        int Mz = delayIndices.size();

        Matrix lambda = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(r))) {
                lambda.set(0, r, sn.rates.get((int) sn.nodeToStation.get(sourceIndex), r));
            }
        }
        Matrix S = new Matrix(queueIndices.size(), 1);
        double Smax = Double.MIN_VALUE;
        for (int i = 0; i < queueIndices.size(); i++) {
            S.set(i, 0, sn.nservers.get((int) sn.nodeToStation.get(queueIndices.get(i))));
            if (Double.isFinite(S.get(i, 0)) && S.get(i, 0) > Smax)
                Smax = S.get(i, 0);
        }
        Matrix D = new Matrix(Mq, R);
        double Nct = 0;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                if (Double.isFinite(N.get(i, j)))
                    Nct += N.get(i, j);
            }
        }
        Matrix mu = Matrix.ones(Mq, (int) (Math.ceil(Nct) + Smax));
        for (int i = 0; i < Mq; i++) {
            for (int r = 0; r < R; r++) {
                int c = 0;
                for (; c < sn.chains.getNumRows(); c++) {
                    if (sn.chains.get(c, r) != 0)
                        break;
                }
                if (sn.refclass.get(c) > 0) {
                    D.set(i, r, sn.visits.get(c).get((int) sn.nodeToStateful.get(queueIndices.get(i)), r) / sn.rates.get((int) sn.nodeToStation.get(queueIndices.get(i)), r) / sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), (int) sn.refclass.get(c)));
                } else {
                    D.set(i, r, sn.visits.get(c).get((int) sn.nodeToStateful.get(queueIndices.get(i)), r) / sn.rates.get((int) sn.nodeToStation.get(queueIndices.get(i)), r));
                }
            }
            for (int j = 0; j < mu.getNumCols(); j++) {
                mu.set(i, j, Maths.min(j + 1, sn.nservers.get((int) sn.nodeToStation.get(queueIndices.get(i)))));
            }
        }
        Matrix Z = new Matrix((int) Maths.max(1, Mz), R);
        for (int i = 0; i < Mz; i++) {
            for (int r = 0; r < R; r++) {
                int c = 0;
                for (; c < sn.chains.getNumRows(); c++) {
                    if (sn.chains.get(c, r) != 0)
                        break;
                }
                if (sn.refclass.get(c) > 0) {
                    Z.set(r, sn.visits.get(c).get((int) sn.nodeToStateful.get(delayIndices.get(i)), r) / sn.rates.get((int) sn.nodeToStation.get(delayIndices.get(i)), r) / sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), (int) sn.refclass.get(c)));
                } else {
                    Z.set(r, sn.visits.get(c).get((int) sn.nodeToStateful.get(delayIndices.get(i)), r) / sn.rates.get((int) sn.nodeToStation.get(delayIndices.get(i)), r));
                }
            }
        }
        Matrix V = null;
        for (Integer i : sn.visits.keySet()) {
            Matrix val = sn.visits.get(i);
            if (V == null) {
                V = val;
            } else {
                V = V.add(1, val);
            }
        }
        snGetProductFormChainParamsReturn ret = new snGetProductFormChainParamsReturn();
        ret.lambda = lambda;
        ret.D = D;
        ret.N = N;
        ret.Z = Z;
        ret.mu = mu;
        ret.S = S;
        ret.V = V;
        return ret;
    }

    /**
     * Calculate new queueing network parameters after aggregating classes into chains
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return chain parameters
     */
    public static snGetDemandsChainReturn snGetDemandsChain(NetworkStruct sn) {
        int M = sn.nstations;
        int K = sn.nclasses;
        int C = sn.nchains;
        Matrix N = sn.njobs;

        Matrix scv = sn.scv.clone();
        scv.apply(Double.NaN, 1, "equal");

        Matrix ST = new Matrix(0, 0);
        sn.rates.divide(1, ST, false);
        ST.removeNaN();

        Matrix alpha = new Matrix(M, K);
        Matrix Vchain = new Matrix(M, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            if (sn.refclass.get(0, c) > -1) {
                for (int i = 0; i < M; i++) {
                    //Vchain(i,c) = sum(sn.visits{c}(i,inchain)) / sum(sn.visits{c}(sn.refstat(inchain(1)),sn.refclass(c)));
                    Matrix visits = sn.visits.get(c);
                    double res = 0;
                    int iIdx = (int) sn.stationToStateful.get(i);
                    for (int col = 0; col < inchain.getNumCols(); col++)
                        res += visits.get(iIdx, (int) inchain.get(0, col));
                    Vchain.set(i, c, res / visits.get((int) sn.stationToStateful.get((int) sn.refstat.get((int) inchain.get(0, 0), 0)), (int) sn.refclass.get(0, c)));
                    //alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k) / sum(sn.visits{c}(i,inchain));
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int k = (int) inchain.get(0, col);
                        alpha.set(i, k, alpha.get(i, k) + visits.get(i, k) / res);
                    }
                }
            } else {
                for (int i = 0; i < M; i++) {
                    //Vchain(i,c) = sum(sn.visits{c}(i,inchain)) / sum(sn.visits{c}(sn.refstat(inchain(1)),inchain));
                    Matrix visits = sn.visits.get(c);
                    double res1 = 0, res2 = 0;
                    int refIdx = (int) sn.stationToStateful.get((int) sn.refstat.get((int) inchain.get(0, 0), 0));
                    int iIdx = (int) sn.stationToStateful.get(i);
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        res1 += visits.get(iIdx, idx);
                        res2 += visits.get(refIdx, idx);
                    }
                    Vchain.set(i, c, res1 / res2);
                    //alpha(i,k) = alpha(i,k) + sn.visits{c}(i,k) / sum(sn.visits{c}(i,inchain));
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int k = (int) inchain.get(0, col);
                        alpha.set(i, k, alpha.get(i, k) + visits.get(iIdx, k) / res1);
                    }
                }
            }
        }

        Vchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        Vchain.apply(Double.NaN, 0, "equal");
        for (int c = 0; c < C; c++) {
            double val = Vchain.get((int) sn.refstat.get((int) sn.inchain.get(c).get(0, 0), 0), c);
            for (int i = Vchain.getColIndexes()[c]; i < Vchain.getColIndexes()[c + 1]; i++)
                Vchain.getNonZeroValues()[i] /= val;
        }
        alpha.apply(Double.POSITIVE_INFINITY, 0, "equal");
        alpha.apply(Double.NaN, 0, "equal");
        alpha.apply(GlobalConstants.Zero, 0, "less");

        Matrix Lchain = new Matrix(M, C);
        Matrix STchain = new Matrix(M, C);
        Matrix SCVchain = new Matrix(M, C);
        Matrix Nchain = new Matrix(1, C);
        Matrix refstatchain = new Matrix(C, 1);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            //Nchain(c) = sum(N(inchain)); isOpenChain = any(isinf(N(inchain)));
            boolean isOpenChain = false;
            double sum = 0;
            for (int col = 0; col < inchain.getNumCols(); col++) {
                sum += N.get((int) inchain.get(0, col));
                if (Double.isInfinite(sum)) {
                    isOpenChain = true;
                    break;
                }
            }
            Nchain.set(0, c, sum);

            for (int i = 0; i < M; i++) {
                sum = 0;
                if (isOpenChain && i == sn.refstat.get((int) inchain.get(0, 0), 0)) {
                    //STchain(i,c) = 1 / sumfinite(sn.rates(i,inchain));
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        double val = sn.rates.get(i, (int) inchain.get(0, col));
                        if (Double.isFinite(val))
                            sum += val;
                    }
                    STchain.set(i, c, 1 / sum);
                } else {
                    //STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        sum += ST.get(i, idx) * alpha.get(i, idx);
                    }
                    STchain.set(i, c, sum);
                }
                Lchain.set(i, c, Vchain.get(i, c) * STchain.get(i, c));
                //alphachain = sum(alpha(i,inchain(isfinite(SCV(i,inchain))))');
                double alphachain = 0;
                for (int col = 0; col < inchain.getNumCols(); col++) {
                    int idx = (int) inchain.get(0, col);
                    double val = scv.get(i, idx);
                    if (Double.isFinite(val))
                        alphachain += alpha.get(i, idx);
                }
                if (alphachain > 0) {
                    sum = 0;
                    for (int col = 0; col < inchain.getNumCols(); col++) {
                        int idx = (int) inchain.get(0, col);
                        sum += scv.get(i, idx) * alpha.get(i, idx);
                    }
                    SCVchain.set(i, c, sum / alphachain);
                }
            }
            refstatchain.set(c, 0, sn.refstat.get((int) inchain.get(0, 0), 0));
            for (int col = 1; col < inchain.getNumCols(); col++) {
                int classIdx = (int) inchain.get(0, col);
                if (sn.refstat.get(classIdx, 0) != refstatchain.get(c, 0))
                    throw new RuntimeException("Class have different reference station");
            }
        }
        Lchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        Lchain.apply(Double.NaN, 0, "equal");
        STchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        STchain.apply(Double.NaN, 0, "equal");
        return new snGetDemandsChainReturn(Lchain, STchain, Vchain, alpha, Nchain, SCVchain, refstatchain);
    }

    /**
     * Calculate class-based performance metrics for a queueing network based on performance measures of its chains
     *
     * @param sn      - NetworkStruct object for the queueing network model
     * @param Lchain  - service demands per chain
     * @param ST      - mean service times per class
     * @param STchain - mean service times per chain
     * @param Vchain  - mean visits per chain
     * @param alpha   - class aggregation coefficients
     * @param Qchain  - mean queue-lengths per chain
     * @param Uchain  - mean utilization per chain
     * @param Rchain  - mean response time per chain
     * @param Tchain  - mean throughput per chain
     * @param Cchain  - mean system response time per chain
     * @param Xchain  - mean system throughput per chain
     * @return chain performance metrics
     */
    public static snDeaggregateChainResultsReturn snDeaggregateChainResults(NetworkStruct sn, Matrix Lchain, Matrix ST, Matrix STchain, Matrix Vchain,
                                                                            Matrix alpha, Matrix Qchain, Matrix Uchain, Matrix Rchain, Matrix Tchain, Matrix Cchain, Matrix Xchain) {

        if (ST == null || ST.isEmpty()) {
            ST = new Matrix(0, 0);
            sn.rates.divide(1.0, ST, false);
            ST.removeNaN();
        }

        if (Cchain != null && !Cchain.isEmpty())
            throw new RuntimeException("Cchain input to snDeaggregateChainResults not yet supported");

        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix X = new Matrix(1, K);
        Matrix U = new Matrix(M, K);
        Matrix Q = new Matrix(M, K);
        Matrix T = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        Matrix C = new Matrix(1, K);

        int idxSink = 0;
        for (Node nodeIter : sn.nodes) {
            if (nodeIter instanceof Sink)
                idxSink = nodeIter.getNodeIdx();
        }

        Matrix Vsinktmp = new Matrix(sn.nodevisits.get(0).getNumRows(), sn.nodevisits.get(0).getNumCols());
        for (int i = 0; i < sn.nodevisits.size(); i++) {
            Vsinktmp = Vsinktmp.add(1, sn.nodevisits.get(i));
        }
        Matrix Vsink = Matrix.extractRows(Vsinktmp, idxSink, idxSink + 1, null);
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            double sum = 0;
            for (int idx = 0; idx < inchain_c.getNumCols(); idx++)
                sum += sn.njobs.get((int) inchain_c.get(idx));
            for (int idx = 0; idx < inchain_c.getNumCols(); idx++) {
                int k = (int) inchain_c.get(0, idx);
                if (Double.isInfinite(sum))
                    X.set(0, k, Xchain.get(0, c) * Vsink.get(0, k));
                else
                    X.set(0, k, Xchain.get(0, c) * alpha.get((int) sn.refstat.get(k, 0), k));
                for (int i = 0; i < M; i++) {
                    if (Uchain == null || Uchain.isEmpty()) {
                        if (Double.isInfinite(sn.nservers.get(i, 0)))
                            U.set(i, k, ST.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k));
                        else
                            U.set(i, k, ST.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k) / sn.nservers.get(i, 0));
                    } else {
                        if (Double.isInfinite(sn.nservers.get(i, 0)))
                            U.set(i, k, ST.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k));
                        else
                            U.set(i, k, Uchain.get(i, c) * alpha.get(i, k));
                    }

                    if (Lchain.get(i, c) > 0) {
                        if (Qchain != null && !Qchain.isEmpty())
                            Q.set(i, k, Qchain.get(i, c) * alpha.get(i, k));
                        else
                            Q.set(i, k, Rchain.get(i, c) * ST.get(i, k) / STchain.get(i, c) * Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c) * alpha.get(i, k));
                        T.set(i, k, Tchain.get(i, c) * alpha.get(i, k));
                        R.set(i, k, Q.get(i, k) / T.get(i, k));
                    } else {
                        T.remove(i, k);
                        R.remove(i, k);
                        Q.remove(i, k);
                    }
                }
                C.set(0, k, sn.njobs.get(0, k) / X.get(0, k));
            }
        }

        Q.abs();
        R.abs();
        X.abs();
        U.abs();
        T.abs();
        C.abs();
        Q.removeNaN();
        Q.apply(Double.POSITIVE_INFINITY, 0, "equal");
        R.removeNaN();
        R.apply(Double.POSITIVE_INFINITY, 0, "equal");
        X.removeNaN();
        X.apply(Double.POSITIVE_INFINITY, 0, "equal");
        U.removeNaN();
        U.apply(Double.POSITIVE_INFINITY, 0, "equal");
        T.removeNaN();
        T.apply(Double.POSITIVE_INFINITY, 0, "equal");
        C.removeNaN();
        C.apply(Double.POSITIVE_INFINITY, 0, "equal");

        return new snDeaggregateChainResultsReturn(Q, U, R, T, C, X);
    }

    /**
     * Checks if the network uses an identical scheduling strategy at every station
     *
     * @param sn       - NetworkStruct object for the queueing network model
     * @param strategy - Scheduling strategy
     * @return boolean
     */
    public static boolean snHasHomogeneousScheduling(NetworkStruct sn, SchedStrategy strategy) {
        int stratCount = 0;
        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == strategy)
                stratCount++;
        }
        return stratCount == sn.sched.size();
    }

    /**
     * Checks if the network uses class priorities
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasPriorities(NetworkStruct sn) {
        for (int i = 0; i < sn.classprio.getNumRows(); i++) {
            for (int j = 0; j < sn.classprio.getNumCols(); j++) {
                if (sn.classprio.get(i, j) > 0)
                    return true;
            }
        }
        return false;
    }

    /**
     * Checks if the network uses fork and/or join nodes
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasForkJoin(NetworkStruct sn) {
        for (int i = 0; i < sn.fj.getNumRows(); i++) {
            for (int j = 0; j < sn.fj.getNumCols(); j++) {
                if (sn.fj.get(i, j) > 0)
                    return true;
            }
        }
        return false;
    }

    /**
     * Checks if the network satisfies product-form assumptions except multiclass heterogeneous FCFS
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasProductFormExceptMultiClassHeterExpFCFS(NetworkStruct sn) {
        boolean ret = true;
        for (int i = 0; i < sn.sched.size(); i++) {
            ret = ret && (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF || sn.sched.get(sn.stations.get(i)) == SchedStrategy.PS ||
                    sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS || sn.sched.get(sn.stations.get(i)) == SchedStrategy.LCFSPR ||
                    sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT);
        }
        ret = ret && !snHasPriorities(sn);
        ret = ret && !snHasForkJoin(sn);

        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                for (int j = 0; j < sn.scv.getNumCols(); j++) {
                    if (Double.isFinite(sn.scv.get(i, j)) && sn.scv.get(i, j) > 0) {
                        ret = ret && (sn.scv.get(i, j) > 1 - GlobalConstants.FineTol) && (sn.scv.get(i, j) < 1 + GlobalConstants.FineTol);
                    }
                }
            }
        }

        return ret;
    }

    /**
     * Checks if the network has a station with load-dependent service process
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasLoadDependence(NetworkStruct sn) {
        return sn.lldscaling.getNumCols() > 0;
    }

    /**
     * Checks if the network has one or more open classes
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasOpenClasses(NetworkStruct sn) {
        for (int i = 0; i < sn.njobs.getNumRows(); i++) {
            for (int j = 0; j < sn.njobs.getNumCols(); j++) {
                if (Double.isInfinite(sn.njobs.get(i, j)))
                    return true;
            }
        }
        return false;
    }

    /**
     * Checks if the network uses class-switching
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasClassSwitching(NetworkStruct sn) {
        return sn.nclasses != sn.nchains;
    }

    /**
     * Checks if the network has a known product-form solution
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasProductForm(NetworkStruct sn) {
        boolean ret = true;
        for (int i = 0; i < sn.sched.size(); i++) {
            ret = ret && (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF || sn.sched.get(sn.stations.get(i)) == SchedStrategy.PS ||
                    sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS || sn.sched.get(sn.stations.get(i)) == SchedStrategy.LCFSPR ||
                    sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT);
        }
        ret = ret && !snHasMultiClassHeterFCFS(sn);
        ret = ret && !snHasPriorities(sn);
        ret = ret && !snHasForkJoin(sn);
        return ret;
    }

    /**
     * Checks if the network has one or more stations with multiclass heterogeneous FCFS
     *
     * @param sn - NetworkStruct object for the queueing network model
     * @return boolean
     */
    public static boolean snHasMultiClassHeterFCFS(NetworkStruct sn) {
        for (int i = 0; i < sn.sched.size(); i++) {
            if (sn.sched.get(sn.stations.get(i)) != SchedStrategy.FCFS) {
                continue;
            }
            Matrix row = Matrix.extractRows(sn.rates, i, i + 1, null);
            if (row.elementMax() - row.elementMin() > 0) {
                return true;
            }
        }
        return false;
    }

    public static class snGetProductFormChainParamsReturn {
        public Matrix lambda;
        public Matrix D;
        public Matrix N;
        public Matrix Z;
        public Matrix mu;
        public Matrix S;
        public Matrix V;
    }

    public static class snGetDemandsChainReturn {
        public Matrix Lchain;
        public Matrix STchain;
        public Matrix Vchain;
        public Matrix alpha;
        public Matrix Nchain;
        public Matrix SCVchain;
        public Matrix refstatchain;

        public snGetDemandsChainReturn(Matrix lchain, Matrix sTchain, Matrix vchain, Matrix alpha,
                                       Matrix nchain, Matrix sCVchain, Matrix refstatchain) {
            this.Lchain = lchain;
            this.STchain = sTchain;
            this.Vchain = vchain;
            this.alpha = alpha;
            this.Nchain = nchain;
            this.SCVchain = sCVchain;
            this.refstatchain = refstatchain;
        }
    }

    public static class snDeaggregateChainResultsReturn {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix T;
        public Matrix C;
        public Matrix X;

        public snDeaggregateChainResultsReturn(Matrix q, Matrix u, Matrix r, Matrix t,
                                               Matrix c, Matrix x) {
            Q = q;
            U = u;
            R = r;
            T = t;
            C = c;
            X = x;
        }
    }

    public static boolean snHasFCFS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasDPS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasGPS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasINF(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasPS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasRAND(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasHOL(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasLCFS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasSEPT(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasLEPT(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasSJF(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasLJF(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasMultiClassFCFS(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasMultiServer(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasSingleChain(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasMultiChain(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasSingleClass(NetworkStruct sn) {
        // TODO
        return false;
    }

    public static boolean snHasMultiClass(NetworkStruct sn) {
        // TODO
        return false;
    }

}
