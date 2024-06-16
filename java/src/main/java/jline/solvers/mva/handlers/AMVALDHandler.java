package jline.solvers.mva.handlers;

import jline.api.PFQN;
import jline.api.SN;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.SerializableFunction;
import jline.util.Matrix;
import jline.util.Pair;

import java.util.*;

/**
 * Handler for the solver_amvald function
 */
public class AMVALDHandler implements MVASolverHandler {
    /*
    * Note: there are some discrepancies between this and solver_amvald in LINE. There is an error around 0.0001 for
    * example_loadDependent_1 and around 0.001 for example_loadDependent_2. TODO: fix this discrepancy
    * */
    public SolverMVAResult solve(NetworkStruct sn, SolverOptions options) {
        long startTime = System.currentTimeMillis();

        if (options == null)
            options = new SolverOptions(SolverType.MVA);

        SN.snGetDemandsChainReturn res = SN.snGetDemandsChain(sn);
        Matrix Lchain = res.Lchain, STchain = res.STchain, Vchain = res.Vchain, alpha = res.alpha, Nchain = res.Nchain, SCVchain = res.SCVchain, refstatchain = res.refstatchain;

        int M = sn.nstations;
        int K = sn.nchains;
        double Nt = 0;
        for(int col = 0; col < Nchain.getNumCols(); col++) {
            if (Double.isFinite(Nchain.get(0, col)))
                Nt += Nchain.get(0, col);
        }
        double tol = options.iter_tol;
        Matrix nservers = sn.nservers;

        Matrix Uchain = new Matrix(M, K);
        Matrix Tchain = new Matrix(M, K);
        Matrix Rchain = new Matrix(M, K);
        Matrix Cchain_s = new Matrix(1, K);

        Matrix Qchain = options.init_sol.clone();
        if ((Qchain == null) || (Qchain.isEmpty())) {
            Qchain = new Matrix(M, K);
            //Qchain = Qchain ./ repmat(sum(Qchain,1),size(Qchain,1),1) .* repmat(Nchain,size(Qchain,1),1);
            for(int row = 0; row < M; row++) {
                for(int col = 0; col < K; col++)
                    Qchain.set(row, col, (1.0/M)*Nchain.get(0, col));
            }
            Qchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
            //Qchain(refstatchain(isinf(Nchain)))=0;
        }

        List<Integer> nnzclasses = new ArrayList<Integer>();
        Matrix Xchain = new Matrix(1, STchain.getNumCols());
        for(int r = 0; r < Nchain.getNumCols(); r++) {
            if (Double.isInfinite(Nchain.get(0, r)))
                Xchain.set(0, r, 1.0/STchain.get((int) refstatchain.get(r,0), r));
            else
                Xchain.set(0, r, 1.0/STchain.sumCols(r));

            if (Nchain.get(0,r) > 0)
                nnzclasses.add(r);
        }

        for(int k = 0; k < M; k++) {
            for(Integer r : nnzclasses) {
                if (Double.isInfinite(nservers.get(k,0)))
                    Uchain.set(k, r, Vchain.get(k,r) * STchain.get(k,r) * Xchain.get(0,r));
                else
                    Uchain.set(k, r, (Vchain.get(k,r) * STchain.get(k,r) * Xchain.get(0,r))/nservers.get(k,0));
            }
        }

        if (options.config.np_priority == null)
            options.config.np_priority = "default";
        if (options.config.multiserver == null)
            options.config.multiserver = "default";
        if (options.config.highvar == null)
            options.config.highvar = "default";

        if (options.method.equals("default")) {
            if (Nt <= 2)
                options.method = "bs";
            else
                options.method = "lin";
        }

        //Use list JLineMatrix to represent 3-D Matrix
        List<Matrix> gamma = new ArrayList<Matrix>();
        Matrix tau = new Matrix(K, K);
        switch(options.method) {
            case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                for(int i = 0; i < K; i++)
                    gamma.add(new Matrix(K, M));
                break;
            default:
                gamma.add(new Matrix(K, M));
        }

        /* Main Loop */
        double omicron = 0.5; // under-relaxation parameter
        double outer_iter = 0;
        Matrix QchainOuter_1 = Qchain.clone();
        Matrix XchainOuter_1 = Xchain.clone();
        Matrix UchainOuter_1 = Uchain.clone();
        Matrix STeff = null;
        while ((outer_iter < 2 || Qchain.sub(1, QchainOuter_1).elementMaxAbs() > tol) && (outer_iter <= options.iter_max)) {
            outer_iter++;

            QchainOuter_1 = Qchain.clone();
            XchainOuter_1 = Xchain.clone();
            UchainOuter_1 = Uchain.clone();

            if (Double.isFinite(Nt) && Nt > 0) {
                if (options.method.equals("default") || options.method.equals("aql") ||
                        options.method.equals("qdaql") || options.method.equals("lin") || options.method.equals("qdlin")){
                    /* Iteration at population N-1_s */
                    for(int s = 0; s < K; s++){
                        if(Double.isFinite(Nchain.get(0,s))) {
                            double iter_s = 0;
                            Matrix Nchain_s = Matrix.oner(Nchain, new ArrayList<Integer>(Collections.singletonList(s)));
                            Matrix Qchain_s = Qchain.clone();
                            Qchain_s.scale((Nt-1)/Nt,  Qchain_s);
                            Matrix Xchain_s = Xchain.clone();
                            Xchain_s.scale((Nt-1)/Nt, Xchain_s);
                            Matrix Uchain_s = Uchain.clone();
                            Uchain_s.scale((Nt-1)/Nt, Uchain_s);
                            Matrix Qchain_s_1 = Qchain_s.clone();
                            Matrix Xchain_s_1 = Xchain_s.clone();
                            Matrix Uchain_s_1 = Uchain_s.clone();

                            int count_inf = 0;
                            double Nt_s = 0;
                            Matrix deltaclass = new Matrix(Nchain_s.getNumRows(), Nchain_s.getNumCols());
                            List<Integer> ocl = new ArrayList<Integer>();
                            List<Integer> ccl = new ArrayList<Integer>();
                            List<Integer> nnzclasses_s = new ArrayList<Integer>();
                            for(int col = 0; col < Nchain_s.getNumCols(); col++) {
                                double val = Nchain_s.get(0, col);
                                if (Double.isInfinite(val)) {
                                    count_inf++;
                                    deltaclass.set(0, col, 1.0);
                                    ocl.add(col);
                                } else {
                                    Nt_s += val;
                                    deltaclass.set(0, col, (val-1)/val);
                                    if (val > 0)
                                        ccl.add(col);
                                }

                                if (val > 0)
                                    nnzclasses_s.add(col);
                            }
                            double delta = count_inf == Nchain_s.getNumCols() ? 1 : (Nt_s - 1) / Nt_s;

                            //Use List<Integer> instead of JLineMatrix to store the index to avoid type transformation
                            Map<Integer, List<Integer>> nnzclasses_eprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            Map<Integer, List<Integer>> nnzclasses_hprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            Map<Integer, List<Integer>> nnzclasses_ehprio = new HashMap<Integer, List<Integer>>(nnzclasses_s.size());
                            for(Integer r : nnzclasses_s) {
                                double prio = sn.classprio.get(0, r);
                                List<Integer> eprio_list = new ArrayList<Integer>();
                                List<Integer> hprio_list = new ArrayList<Integer>();
                                for(int i = 0; i < sn.classprio.getNumCols(); i++) {
                                    if (Double.compare(prio, sn.classprio.get(0,i)) == 0)
                                        eprio_list.add(i);
                                    else if (Double.compare(prio, sn.classprio.get(0,i)) < 0)
                                        hprio_list.add(i);
                                }
                                List<Integer> eprio_common = new ArrayList<Integer>(nnzclasses_s);
                                List<Integer> hprio_common = new ArrayList<Integer>(nnzclasses_s);
                                eprio_common.retainAll(eprio_list);
                                hprio_common.retainAll(hprio_list);
                                nnzclasses_eprio.put(r, eprio_common);
                                nnzclasses_hprio.put(r, hprio_common);

                                Set<Integer> ehprio_common = new LinkedHashSet<Integer>(eprio_common);
                                ehprio_common.addAll(hprio_common);
                                nnzclasses_ehprio.put(r, new ArrayList<Integer>(ehprio_common));
                            }
                            while ((iter_s < 2 || Qchain_s.sub(1, Qchain_s_1).elementMaxAbs() > tol) && (iter_s <= options.iter_max)) {
                                iter_s++;

                                Qchain_s_1 = Qchain_s.clone();
                                Xchain_s_1 = Xchain_s.clone();
                                Uchain_s_1 = Uchain_s.clone();

                                Pair<Matrix, Matrix> ret = solver_amvald_forward(gamma, tau, Qchain_s_1, Xchain_s_1, Uchain_s_1, STchain, Vchain, Nchain_s, SCVchain,
                                        Nt_s, delta, deltaclass, ocl, ccl, nnzclasses_s, nnzclasses_eprio, nnzclasses_hprio, nnzclasses_ehprio, sn, options);
                                Matrix Wchain_s = ret.getLeft();
                                Matrix STeff_s = ret.getRight();

                                for(Integer r : nnzclasses) {
                                    if (Wchain_s.sumCols(r) == 0) {
                                        Xchain_s.remove(0, r);
                                    } else {
                                        if (Double.isInfinite(Nchain_s.get(0, r))){
                                            double val = 0;
                                            for(int row = 0; row < Vchain.getNumRows(); row++)
                                                val += Vchain.get(row, r) * Wchain_s.get(row, r);
                                            Cchain_s.set(0, r, val);
                                        } else if (Nchain.get(0, r) == 0) {
                                            Xchain_s.remove(0, r);
                                            Cchain_s.remove(0, r);
                                        } else {
                                            double val = 0;
                                            for(int row = 0; row < Vchain.getNumRows(); row++)
                                                val += Vchain.get(row, r) * Wchain_s.get(row, r);
                                            Cchain_s.set(0, r, val);
                                            Xchain_s.set(0, r, omicron * Nchain_s.get(0, r) / Cchain_s.get(0, r) + (1-omicron) * Xchain_s_1.get(0, r));
                                        }
                                    }
                                    for (int k = 0; k < M; k++) {
                                        //Rchain_s(k,r) = Vchain(k,r) * Wchain_s(k,r); NOT USED
                                        //Tchain_s(k,r) = Xchain_s(r) * Vchain(k,r); NOT USED
                                        Qchain_s.set(k, r, omicron * Xchain_s.get(0, r) * Vchain.get(k, r) * Wchain_s.get(k, r) + (1-omicron) * Qchain_s_1.get(k, r));
                                        Uchain_s.set(k, r, omicron * Vchain.get(k, r) * STeff_s.get(k, r) * Xchain_s.get(0, r) + (1-omicron) * Uchain_s_1.get(k, r));
                                    }
                                }
                            }

                            switch (options.method) {
                                case "default":
                                case "lin":
                                    for(int k = 0; k < M; k++) {
                                        for(Integer r : nnzclasses) {
                                            if (Double.isFinite(Nchain.get(0, r)) && Nchain_s.get(0, r) > 0)
                                                gamma.get(r).set(s, k, Qchain_s_1.get(k,r)/Nchain_s.get(0,r) - QchainOuter_1.get(k,r)/Nchain.get(0,r));
                                        }
                                    }
                                    break;
                                default:
                                    for(int k = 0; k < M; k++) {
                                        gamma.get(0).set(s, k, Qchain_s_1.sumRows(k)/(Nt-1) - QchainOuter_1.sumRows(k)/Nt);
                                    }
                            }

                            for(Integer r : nnzclasses) {
                                tau.set(s, r, Xchain_s_1.get(0,r) - XchainOuter_1.get(0,r));
                            }
                        }
                    }
                }
            }

            double inner_iter = 0;
            Matrix Qchain_inner = Qchain.clone();
            Matrix Xchain_inner = Xchain.clone();
            Matrix Uchain_inner = Uchain.clone();

            int count_inf = 0;
            double Nt_inner = 0;
            Matrix deltaclass = new Matrix(Nchain.getNumRows(), Nchain.getNumCols());
            List<Integer> ocl = new ArrayList<Integer>();
            List<Integer> ccl = new ArrayList<Integer>();
            List<Integer> nnzclasses_inner = new ArrayList<Integer>();
            for(int col = 0; col < Nchain.getNumCols(); col++) {
                double val = Nchain.get(0, col);
                if (Double.isInfinite(val)) {
                    count_inf++;
                    deltaclass.set(0, col, 1.0);
                    ocl.add(col);
                } else {
                    Nt_inner += val;
                    deltaclass.set(0, col, (val-1)/val);
                    if (val > 0)
                        ccl.add(col);
                }

                if (val > 0)
                    nnzclasses_inner.add(col);
            }
            double delta = count_inf == Nchain.getNumCols() ? 1 : (Nt_inner - 1) / Nt_inner;

            //Use List<Integer> instead of JLineMatrix to store the index to avoid type transformation
            Map<Integer, List<Integer>> nnzclasses_eprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            Map<Integer, List<Integer>> nnzclasses_hprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            Map<Integer, List<Integer>> nnzclasses_ehprio = new HashMap<Integer, List<Integer>>(nnzclasses_inner.size());
            for(Integer r : nnzclasses_inner) {
                double prio = sn.classprio.get(0, r);
                List<Integer> eprio_list = new ArrayList<Integer>();
                List<Integer> hprio_list = new ArrayList<Integer>();
                for(int i = 0; i < sn.classprio.getNumCols(); i++) {
                    if (Double.compare(prio, sn.classprio.get(0,i)) == 0)
                        eprio_list.add(i);
                    else if (Double.compare(prio, sn.classprio.get(0,i)) < 0)
                        hprio_list.add(i);
                }
                List<Integer> eprio_common = new ArrayList<Integer>(nnzclasses_inner);
                List<Integer> hprio_common = new ArrayList<Integer>(nnzclasses_inner);
                eprio_common.retainAll(eprio_list);
                hprio_common.retainAll(hprio_list);
                nnzclasses_eprio.put(r, eprio_common);
                nnzclasses_hprio.put(r, hprio_common);

                Set<Integer> ehprio_common = new LinkedHashSet<Integer>(eprio_common);
                ehprio_common.addAll(hprio_common);
                nnzclasses_ehprio.put(r, new ArrayList<Integer>(ehprio_common));
            }
            while ((inner_iter < 2 || Qchain_inner.sub(1, Qchain).elementMaxAbs() > tol) && (inner_iter <= options.iter_max)) {
                inner_iter++;

                Qchain_inner = Qchain.clone();
                Xchain_inner = Xchain.clone();
                Uchain_inner = Uchain.clone();

                Pair<Matrix, Matrix> ret = solver_amvald_forward(gamma, tau, Qchain_inner, Xchain_inner, Uchain_inner, STchain, Vchain, Nchain, SCVchain,
                        Nt_inner, delta, deltaclass, ocl, ccl, nnzclasses_inner, nnzclasses_eprio, nnzclasses_hprio, nnzclasses_ehprio, sn, options);
                Matrix Wchain = ret.getLeft();
                STeff = ret.getRight();

                for(Integer r : nnzclasses) {
                    if (Wchain.sumCols(r) == 0) {
                        Xchain.remove(0, r);
                    } else {
                        if (Double.isInfinite(Nchain.get(0, r))) {
                            double val = 0;
                            for(int i = 0; i < Vchain.getNumRows(); i++)
                                val += Vchain.get(i, r) * Wchain.get(i, r);
                            Cchain_s.set(0, r, val);
                        } else if (Nchain.get(0, r) == 0) {
                            Xchain.remove(0, r);
                            Cchain_s.remove(0, r);
                        } else {
                            double val = 0;
                            for(int i = 0; i < Vchain.getNumRows(); i++)
                                val += Vchain.get(i, r) * Wchain.get(i, r);
                            Cchain_s.set(0, r, val);
                            Xchain.set(0, r, omicron * Nchain.get(0, r) / Cchain_s.get(0, r) + (1-omicron) * Xchain_inner.get(0, r));
                        }
                    }
                    for(int k = 0; k < M; k++) {
                        //Rchain(k,r) = Vchain(k,r) * Wchain(k,r); NOT USED
                        Qchain.set(k, r, omicron * Xchain.get(0, r) * Vchain.get(k, r) * Wchain.get(k, r) + (1-omicron) * Qchain_inner.get(k, r));
                        Tchain.set(k, r, Xchain.get(0, r) * Vchain.get(k, r));
                        Uchain.set(k, r, omicron * Vchain.get(k, r) * STeff.get(k, r) * Xchain.get(0, r) + (1-omicron) * Uchain_inner.get(k, r));
                    }
                }
            }
        }

        for(int k = 0; k < M; k++) {
            for(int r = 0; r < K; r++) {
                if (Vchain.get(k,r) * STeff.get(k,r) > 0) {
                    switch (sn.sched.get(sn.stations.get(k))) {
                        case FCFS:
                        case SIRO:
                        case PS:
                        case LCFSPR:
                        case DPS:
                        case HOL:
                            if (Uchain.sumRows(k) > 1) {
                                double sum_vchain_steff_xchain_k = 0;
                                for(int i = 0; i < Vchain.getNumCols(); i++)
                                    sum_vchain_steff_xchain_k += Vchain.get(k,i) * STeff.get(k, i) * Xchain.get(0, i);
                                Uchain.set(k, r, (Math.min(Uchain.sumRows(k), 1) * Vchain.get(k,r) * STeff.get(k,r) * Xchain.get(0, r)) / sum_vchain_steff_xchain_k);
                            }
                            break;
                        default:
                            break;
                    }
                }
            }
        }

        for(int k = 0; k < M; k++) {
            for(int r = 0; r < K; r++)
                Rchain.set(k, r, Qchain.get(k, r) / Tchain.get(k, r));
        }
        Xchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        Xchain.apply(Double.NaN, 0, "equal");
        Uchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        Uchain.apply(Double.NaN, 0, "equal");
        Rchain.apply(Double.POSITIVE_INFINITY, 0, "equal");
        Rchain.apply(Double.NaN, 0, "equal");

        for(int col = 0; col < K; col++) {
            if (Nchain.get(0, col) == 0) {
                Xchain.remove(0, col);
                for(int row = 0; row < M; row++) {
                    Uchain.remove(row, col);
                    Rchain.remove(row, col);
                    Tchain.remove(row, col);
                }
            }
        }

        SN.snDeaggregateChainResultsReturn ret = null;
        if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && (sn.cdscaling == null || sn.cdscaling.size() == 0))
            ret = SN.snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
        else
            ret = SN.snDeaggregateChainResults(sn, Lchain, null, STchain, Vchain, alpha, null, Uchain, Rchain, Tchain, null, Xchain);

        List<Integer> ccl = new ArrayList<Integer>();
        for(int i = 0; i < Nchain.getNumCols(); i++) {
            if (Double.isFinite(Nchain.get(0, i)))
                ccl.add(i);
        }
        Matrix Nclosed = new Matrix(1, ccl.size());
        Matrix Xclosed = new Matrix(1, ccl.size());
        for(int i = 0; i < ccl.size(); i++) {
            Nclosed.set(0, i, Nchain.get(0, ccl.get(i)));
            Xclosed.set(0, i, Xchain.get(0, ccl.get(i)));
        }
        double lG = 0;
        for(int i = 0; i < ccl.size(); i++) {
            if (Xclosed.get(0, i) > options.tol)
                lG += -Nclosed.get(0, i) * Math.log(Xclosed.get(0, i));
        }

        long endTime = System.currentTimeMillis();
        long runTime = endTime - startTime;

        SolverMVAResult result = new SolverMVAResult();
        result.method = options.method;
        result.QN = ret.Q;
        result.RN = ret.R;
        result.XN = ret.X;
        result.UN = ret.U;
        result.TN = ret.T;
        result.CN = ret.C;
        result.runtime = runTime/1000.0;
        result.logNormConstAggr = lG;
        return result;
    }

    public Pair<Matrix, Matrix> solver_amvald_forward(List<Matrix> gamma, Matrix tau, Matrix Qchain_in, Matrix Xchain_in, Matrix Uchain_in,
                                                      Matrix STchain_in, Matrix Vchain_in, Matrix Nchain_in, Matrix SCVchain_in,
                                                      double Nt, double delta, Matrix deltaclass, List<Integer> ocl, List<Integer> ccl, List<Integer> nnzclasses,
                                                      Map<Integer, List<Integer>> nnzclasses_eprio, Map<Integer, List<Integer>> nnzclasses_hprio, Map<Integer, List<Integer>> nnzclasses_ehprio, NetworkStruct sn, SolverOptions options) {

        int M = sn.nstations;
        int K = sn.nchains;
        Matrix nservers = sn.nservers;
        Matrix schedparam = sn.schedparam;
        Matrix lldscaling = null;
        Map<Station, SerializableFunction<Matrix, Double>> cdscaling = sn.cdscaling;

        if (gamma == null || gamma.size() == 0) {
            gamma = new ArrayList<Matrix>();
            gamma.add(new Matrix(K, M));
        }

        /* Evaluate lld and cd correction factors */
        Matrix totArvlQlenSeenByOpen = new Matrix(K, M);
        Matrix interpTotArvlQlen = new Matrix(M, 1);
        Matrix totArvlQlenSeenByClosed = new Matrix(M, K);
        Matrix stationaryQlen = new Matrix(M, K);
        Matrix selfArvlQlenSeenByClosed = new Matrix(M, K);
        for(int k = 0; k < M; k++) {
            double sum_qchain_k_nnz = 0;
            for(Integer r : nnzclasses)
                sum_qchain_k_nnz += Qchain_in.get(k, r);

            interpTotArvlQlen.set(k, 0, delta * sum_qchain_k_nnz);
            for(Integer r : nnzclasses) {
                selfArvlQlenSeenByClosed.set(k, r, deltaclass.get(r) * Qchain_in.get(k,r));
                if (Objects.requireNonNull(sn.sched.get(sn.stations.get(k))) == SchedStrategy.HOL) {
                    double sum_qchain_k_ehprio = 0;
                    for (Integer i : nnzclasses_ehprio.get(r))
                        sum_qchain_k_ehprio += Qchain_in.get(k, i);
                    totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_ehprio);
                    totArvlQlenSeenByClosed.set(k, r, deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_ehprio - Qchain_in.get(k, r));
                } else {
                    totArvlQlenSeenByOpen.set(r, k, sum_qchain_k_nnz);
                    totArvlQlenSeenByClosed.set(k, r, deltaclass.get(r) * Qchain_in.get(k, r) + sum_qchain_k_nnz - Qchain_in.get(k, r));
                }
                stationaryQlen.set(k, r, Qchain_in.get(k,r));
            }
        }

        if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
            //lldscaling = new Matrix(M, (int) Math.ceil(Nt)); // fails on LayeredExamples.ex3
            lldscaling = new Matrix(M, (int) 1);
            lldscaling.fill(1);
        } else {
            lldscaling = sn.lldscaling.clone();
        }
        Matrix lldterm = PFQN.pfqn_lldfun(interpTotArvlQlen.elementIncrease(1), lldscaling, null);

        Matrix cdterm = new Matrix(M, K);
        cdterm.fill(1);
        for(Integer r : nnzclasses) {
            if (!(cdscaling == null || cdscaling.size() == 0)) {
                if (Double.isFinite(Nchain_in.get(0,r)))
                    Matrix.extract(PFQN.pfqn_cdfun(selfArvlQlenSeenByClosed.elementIncrease(1.0), cdscaling, sn.stations), 0, M, 0, 1, cdterm, 0, r);
                else
                    Matrix.extract(PFQN.pfqn_cdfun(stationaryQlen.elementIncrease(1.0), cdscaling, sn.stations), 0, M, 0, 1, cdterm, 0, r);
            }
        }

        Matrix msterm = new Matrix(0,0);
        switch (options.config.multiserver){
            case "softmin":
                switch (options.method) {
                    case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                        Matrix g = new Matrix(ccl.size(), M);
                        for(Integer r : ccl) {
                            double param = ((Nt-1.0) / Nt) * Nchain_in.get(0,r);
                            Matrix gamma_r = gamma.get(r);
                            for(int i = 0; i < ccl.size(); i++) {
                                int idx = ccl.get(i);
                                for(int j = 0; j < M; j++)
                                    g.set(i, j, g.get(i,j) + param * gamma_r.get(idx, j));
                            }
                        }
                        Matrix input = interpTotArvlQlen.add(1, g.meanCol().transpose());
                        msterm = PFQN.pfqn_lldfun(input.elementIncrease(1), new Matrix(0,0), nservers);
                        break;
                    default:
                        g = new Matrix(ccl.size(), M);
                        for(int i = 0; i < ccl.size(); i++) {
                            int idx = ccl.get(i);
                            for(int j = 0; j < M; j++)
                                g.set(i, j, g.get(i,j) + (Nt - 1.0) * gamma.get(0).get(idx, j));
                        }
                        input = interpTotArvlQlen.add(1, g.meanCol().transpose());
                        //CommonOps_DSCC.add(1, interpTotArvlQlen, 1, CommonOps_DSCC.transpose(g.meanCol(), null, null), input, null, null);
                        msterm = PFQN.pfqn_lldfun(interpTotArvlQlen.elementIncrease(1), new Matrix(0,0), nservers);
                }
                break;
            case "seidmann":
                nservers.divide(1.0, msterm, false);
                msterm.apply(0, 1, "equal");
                break;
            case "default":
                switch (options.method) {
                    case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                        Matrix g = new Matrix(ccl.size(), M);
                        for(Integer r : ccl) {
                            double param = ((Nt-1.0) / Nt) * Nchain_in.get(0,r);
                            Matrix gamma_r = gamma.get(r);
                            for(int i = 0; i < ccl.size(); i++) {
                                int idx = ccl.get(i);
                                for(int j = 0; j < M; j++)
                                    g.set(i, j, g.get(i,j) + param * gamma_r.get(idx, j));
                            }
                        }
                        Matrix input = interpTotArvlQlen.add(1, g.meanCol().transpose());
                        //CommonOps_DSCC.add(1, interpTotArvlQlen, 1, CommonOps_DSCC.transpose(g.meanCol(), null, null), input, null, null);
                        msterm = PFQN.pfqn_lldfun(input.elementIncrease(1), new Matrix(0,0), nservers);
                        break;
                    default:
                        g = new Matrix(1, M);
                        for(Integer r : ccl) {
                            for(int j = 0; j < M; j++) {
                                g.set(0, j, g.get(0,j) + (Nt - 1.0) * gamma.get(0).get(r, j));
                            }
                        }
                        msterm = PFQN.pfqn_lldfun(interpTotArvlQlen.elementIncrease(1 + g.meanRow().get(0,0)), new Matrix(0,0), nservers);
                }
                for(int i = 0; i < M; i++) {
                    SchedStrategy schedStrategy = sn.sched.get(sn.stations.get(i));
                    if(schedStrategy.equals(SchedStrategy.FCFS) ||  schedStrategy.equals(SchedStrategy.SIRO) || schedStrategy.equals(SchedStrategy.LCFSPR))
                        msterm.set(i, 0, 1.0 / nservers.get(i, 0));
                }
                break;
            default:
                throw new RuntimeException("nrecognize multiserver approximation method");
        }

        Matrix Wchain = new Matrix(M, K);
        Matrix STeff = new Matrix(STchain_in.getNumRows(), STchain_in.getNumCols());
        lldterm = lldterm.repmat(1, K);
        for(Integer r : nnzclasses) {
            for(int k = 0; k < M; k++) {
                STeff.set(k, r, STchain_in.get(k, r) * lldterm.get(k, r) * msterm.get(k, 0) * cdterm.get(k,r));
            }
        }

        /* if amva.qli or amva.fli, update now totArvlQlenSeenByClosed with STeff */
        switch (options.method) {
            case "amva_qli": case "qli":
                List<Integer> infset = new ArrayList<Integer>();
                for(int i = 0; i < M; i++) {
                    if (sn.sched.get(sn.stations.get(i)).equals(SchedStrategy.INF))
                        infset.add(i);
                }

                for(int k = 0; k < M; k++) {
                    if (sn.sched.get(sn.stations.get(k)).equals(SchedStrategy.HOL)) {
                        for(Integer r : nnzclasses) {
                            //Calculate sum_Qchain_k_r_ephrio;
                            double sum_Qchain_k_r_ephrio = 0;
                            for(Integer i : nnzclasses_ehprio.get(r))
                                sum_Qchain_k_r_ephrio += Qchain_in.get(k, i);

                            if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20){
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k,r));
                            } else {
                                double qlinum = STeff.get(k,r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k,r));
                                double qliden = 0;
                                for(Integer i : infset)
                                    qliden += STeff.get(i, r);
                                for(int m = 0; m < M; m++) {
                                    double sum_Qchain_m_r_ephrio = 0;
                                    for(Integer i : nnzclasses_ehprio.get(r))
                                        sum_Qchain_m_r_ephrio += Qchain_in.get(m, i);
                                    qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m,r));
                                }
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - (1/(Nchain_in.get(0,r)-1))*(Qchain_in.get(k,r) - qlinum/qliden));
                            }
                        }
                    } else {
                        for(Integer r : nnzclasses) {
                            //Calculate sum_Qchain_k_r_ephrio;
                            double sum_Qchain_k_r_nnzclasses = 0;
                            for(Integer i : nnzclasses)
                                sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, i);

                            if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20){
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k,r));
                            } else {
                                double qlinum = STeff.get(k,r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k,r));
                                double qliden = 0;
                                for(Integer i : infset)
                                    qliden += STeff.get(i, r);
                                for(int m = 0; m < M; m++) {
                                    double sum_Qchain_m_r_nnzclasses = 0;
                                    for(Integer i : nnzclasses)
                                        sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, i);
                                    qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m,r));
                                }
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - (1/(Nchain_in.get(0,r)-1))*(Qchain_in.get(k,r) - qlinum/qliden));
                            }
                        }
                    }
                }
                break;
            case "amva_fli": case "fli":
                infset = new ArrayList<Integer>();
                for(int i = 0; i < M; i++) {
                    if (sn.sched.get(sn.stations.get(i)).equals(SchedStrategy.INF))
                        infset.add(i);
                }

                for(int k = 0; k < M; k++) {
                    if (sn.sched.get(sn.stations.get(k)).equals(SchedStrategy.HOL)) {
                        for(Integer r : nnzclasses) {
                            //Calculate sum_Qchain_k_r_ephrio;
                            double sum_Qchain_k_r_ephrio = 0;
                            for(Integer i : nnzclasses_ehprio.get(r))
                                sum_Qchain_k_r_ephrio += Qchain_in.get(k, i);

                            if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20){
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - Qchain_in.get(k,r));
                            } else {
                                double qlinum = STeff.get(k,r) * (1 + sum_Qchain_k_r_ephrio - Qchain_in.get(k,r));
                                double qliden = 0;
                                for(Integer i : infset)
                                    qliden += STeff.get(i, r);
                                for(int m = 0; m < M; m++) {
                                    double sum_Qchain_m_r_ephrio = 0;
                                    for(Integer i : nnzclasses_ehprio.get(r))
                                        sum_Qchain_m_r_ephrio += Qchain_in.get(m, i);
                                    qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_ephrio - Qchain_in.get(m,r));
                                }
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_ephrio - (2/Nchain_in.get(0,r))*(Qchain_in.get(k,r) + qlinum/qliden));
                            }
                        }
                    } else {
                        for(Integer r : nnzclasses) {
                            //Calculate sum_Qchain_k_r_ephrio;
                            double sum_Qchain_k_r_nnzclasses = 0;
                            for(Integer i : nnzclasses)
                                sum_Qchain_k_r_nnzclasses += Qchain_in.get(k, i);

                            if (Math.abs(Nchain_in.get(0, r) - 1) < 1e-20){
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - Qchain_in.get(k,r));
                            } else {
                                double qlinum = STeff.get(k,r) * (1 + sum_Qchain_k_r_nnzclasses - Qchain_in.get(k,r));
                                double qliden = 0;
                                for(Integer i : infset)
                                    qliden += STeff.get(i, r);
                                for(int m = 0; m < M; m++) {
                                    double sum_Qchain_m_r_nnzclasses = 0;
                                    for(Integer i : nnzclasses)
                                        sum_Qchain_m_r_nnzclasses += Qchain_in.get(m, i);
                                    qliden += STeff.get(m, r) * (1 + sum_Qchain_m_r_nnzclasses - Qchain_in.get(m,r));
                                }
                                totArvlQlenSeenByClosed.set(k, r, sum_Qchain_k_r_nnzclasses - (2/Nchain_in.get(0,r))*(Qchain_in.get(k,r) + qlinum/qliden));
                            }
                        }
                    }
                }
            default:
                break;
        }

        /* Compute response time */
        for(Integer r : nnzclasses) {
            List<Integer> sd = new ArrayList<Integer>(nnzclasses);
            List<Integer> sdprio = new ArrayList<Integer>(nnzclasses_ehprio.get(r));
            sd.remove(r);
            sdprio.remove(r);

            for(int k = 0; k < M; k++) {
                switch (sn.sched.get(sn.stations.get(k))){
                    case INF:
                        Wchain.set(k, r, STeff.get(k,r));
                        break;
                    case PS:
                        switch (options.method) {
                            case "def": case "amva": case "amva_qd": case "amva_qdamva": case "qd": case "qdamva": case "lin": case "qdlin":
                                if (options.config.multiserver.equals("seidmann")) {
                                    double val = STeff.get(k,r) * (nservers.get(k, 0) - 1);
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, val + STeff.get(k,r) * (1 + totArvlQlenSeenByOpen.get(r,k)));
                                    } else {
                                        switch (options.method) {
                                            case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                                                //Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1)
                                                double tmp = 0;
                                                for(Integer c : ccl)
                                                    tmp += Nchain_in.get(0, c) * gamma.get(c).get(r,k);
                                                Wchain.set(k, r, val + STeff.get(k,r) * (1 + interpTotArvlQlen.get(k,0) + tmp - gamma.get(r).get(r,k)));
                                                break;
                                            default:
                                                Wchain.set(k, r, val + STeff.get(k, r) * (1 + interpTotArvlQlen.get(k,0) + (Nt-1)*gamma.get(0).get(r,k)));
                                        }
                                    }
                                } else {
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, STeff.get(k,r) * (1 + totArvlQlenSeenByOpen.get(r,k)));
                                    } else {
                                        switch (options.method) {
                                            case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                                                double tmp = 0;
                                                for(Integer c : ccl)
                                                    tmp += Nchain_in.get(0, c) * gamma.get(c).get(r,k);
                                                Wchain.set(k, r, Wchain.get(k, r) + STeff.get(k,r) * (1 + interpTotArvlQlen.get(k,0) + tmp - gamma.get(r).get(r,k)));
                                                break;
                                            default:
                                                Wchain.set(k, r, STeff.get(k, r) * (1 + interpTotArvlQlen.get(k,0) + (Nt-1)*gamma.get(0).get(r,k)));
                                        }
                                    }
                                }
                                break;
                            default:
                                if (options.config.multiserver.equals("seidmann")) {
                                    double val = STeff.get(k,r) * (nservers.get(k, 0) - 1);
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, val + STeff.get(k,r) * (1 + totArvlQlenSeenByOpen.get(r,k)));
                                    } else {
                                        switch (options.method) {
                                            case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                                                double tmp = 0;
                                                for(Integer c : ccl)
                                                    tmp += Nchain_in.get(0, c) * gamma.get(c).get(r,k);
                                                Wchain.set(k, r, val + STeff.get(k,r) * (1 + totArvlQlenSeenByClosed.get(k,0) + tmp - gamma.get(r).get(r,k)));
                                                break;
                                            default:
                                                Wchain.set(k, r, val + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,0) + (Nt-1)*gamma.get(0).get(r,k)));
                                        }
                                    }
                                } else {
                                    double val = Wchain.get(k,r);
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, val + STeff.get(k,r) * (1 + totArvlQlenSeenByOpen.get(r, k)));
                                    } else {
                                        switch (options.method) {
                                            case "default": case "amva_lin": case "lin": case "amva_qdlin": case "qdlin":
                                                double tmp = 0;
                                                for(Integer c : ccl)
                                                    tmp += Nchain_in.get(0, c) * gamma.get(c).get(r,k);
                                                Wchain.set(k, r, val + STeff.get(k,r) * (1 + totArvlQlenSeenByClosed.get(k,0) + tmp - gamma.get(r).get(r,k)));
                                                break;
                                            default:
                                                Wchain.set(k, r, val + STeff.get(k, r) * (1 + totArvlQlenSeenByClosed.get(k,0) + (Nt-1)*gamma.get(0).get(r,k)));
                                        }
                                    }
                                }
                        }
                        break;
                    case DPS:
                        if (nservers.get(k,0) > 1)
                            throw new RuntimeException("Multi-server DPS not supported yet in AMVA solver");

                        Wchain.set(k, r, STeff.get(k, r) * (1 + selfArvlQlenSeenByClosed.get(k,r)));
                        for(Integer s : sd) {
                            if (schedparam.get(k,s) == schedparam.get(k,r))
                                Wchain.set(k, r, Wchain.get(k,r) + STeff.get(k,r) * stationaryQlen.get(k,s));
                            else if (schedparam.get(k,s) / schedparam.get(k,r) <= Double.POSITIVE_INFINITY)
                                Wchain.set(k, r, Wchain.get(k,r) + STeff.get(k,r) * stationaryQlen.get(k,s) * schedparam.get(k,s) / schedparam.get(k,r));
                        }
                        break;
                    case FCFS:
                    case SIRO:
                    case LCFSPR:
                        if (STeff.get(k,r) <= 0)
                            break;

                        //Uchain_r = Uchain_in ./ repmat(Xchain_in,M,1) .* (repmat(Xchain_in,M,1) + repmat(tau(r,:),M,1));
                        Matrix Uchain_r = new Matrix(Uchain_in.getNumRows(), Uchain_in.getNumCols());
                        for(int i = 0; i < Uchain_in.getNumRows(); i++) {
                            for(int j = 0; j < Uchain_in.getNumCols(); j++)
                                Uchain_r.set(i, j, (Uchain_in.get(i,j) / Xchain_in.get(0, j)) * (Xchain_in.get(0,j) + tau.get(r,j)));
                        }

                        Matrix Bk = new Matrix(1, K);
                        if (nservers.get(k,0) > 1) {
                            Matrix deltaclass_r = new Matrix(Xchain_in.getNumRows(), Xchain_in.getNumCols());
                            deltaclass_r.fill(1.0);
                            deltaclass_r.set(0, r, deltaclass.get(0, r));
                            //Compute load: deltaclass_r .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)
                            Matrix BK_tmp = new Matrix(1, K);
                            for(int i = 0; i < K; i++)
                                BK_tmp.set(0, i, deltaclass_r.get(0, i) * Xchain_in.get(0, i) * Vchain_in.get(k, i) * STeff.get(k, i));
                            if (BK_tmp.elementSum() < 0.75) {
                                Bk = BK_tmp.clone();
                            } else {
                                Bk = BK_tmp.power(nservers.get(k, 0) - 1);
                            }
                        } else {
                            Bk.fill(1);
                        }


                        if (nservers.get(k, 0) == 1 && (((lldscaling != null) && !lldscaling.isEmpty()) || ((cdscaling != null) && (cdscaling.size() != 0)))) {
                            if (options.config.highvar.equals("hvmva")) {
                                double sum_uchain_r = 0, val = 0;
                                for (Integer s : ccl) {
                                    sum_uchain_r += Uchain_r.get(k, s);
                                    val += STeff.get(k, s) * Uchain_r.get(k, s) * (1.0 + SCVchain_in.get(k, s)) / 2.0;
                                }
                                Wchain.set(k, r, val + STeff.get(k, r) * (1 - sum_uchain_r));
                            } else {
                                Wchain.set(k, r, STeff.get(k, r));
                            }

                            double steff_mult_stationaryQlen = 0;
                            for(Integer s : sd)
                                steff_mult_stationaryQlen += STeff.get(k, s) * stationaryQlen.get(k, s);

                            if (ocl.contains(r)) {
                                Wchain.set(k, r, Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k, r) + steff_mult_stationaryQlen));
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                Wchain.set(k, r, Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) + steff_mult_stationaryQlen));
                            }
                        } else {
                            double steff_mult_stationarQlen_mult_Bk = 0;
                            for(Integer s : sd)
                                steff_mult_stationarQlen_mult_Bk += STeff.get(k, s) * stationaryQlen.get(k, s) * Bk.get(0, s);

                            if (options.config.multiserver.equals("softmin")) {
                                if (ocl.contains(r)) {
                                    Wchain.set(k, r, STeff.get(k, r) + STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk);
                                } else {
                                    //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                    //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sd) * (stationaryQlen(k,sd) .* Bk(sd))' + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                    Wchain.set(k, r, STeff.get(k, r) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk);
                                }
                            } else {
                                if (ocl.contains(r)) {
                                    Wchain.set(k, r, STeff.get(k, r) * (nservers.get(k, 0) - 1) + STeff.get(k, r) + (STeff.get(k,r) * deltaclass.get(0, r) * stationaryQlen.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk));
                                } else {
                                    //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                    //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sd).*Bk(sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                    Wchain.set(k, r, STeff.get(k,r) * (nservers.get(k, 0) - 1) + STeff.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) + steff_mult_stationarQlen_mult_Bk));
                                }
                            }
                        }
                        break;
                    case HOL:
                        if (STeff.get(k,r) <= 0)
                            break;

                        Uchain_r = new Matrix(Uchain_in.getNumRows(), Uchain_in.getNumCols());
                        for(int i = 0; i < Uchain_in.getNumRows(); i++) {
                            for(int j = 0; j < Uchain_in.getNumCols(); j++)
                                Uchain_r.set(i, j, (Uchain_in.get(i,j) / Xchain_in.get(0, j)) * (Xchain_in.get(0,j) + tau.get(r,j)));
                        }

                        double prioScaling = 0;
                        switch (options.config.np_priority) {
                            case "default":
                            case "cl":
                                double UHigherPrio = 0;
                                for(Integer h : nnzclasses_hprio.get(r))
                                    UHigherPrio += Vchain_in.get(k,h)*STeff.get(k,h)*(Xchain_in.get(0,h) - Qchain_in.get(k,h)*tau.get(h));
                                prioScaling = Math.min(Math.max(options.tol, 1 - UHigherPrio), 1 - options.tol);
                                break;
                            case "shadow":
                                UHigherPrio = 0;
                                for(Integer h : nnzclasses_hprio.get(r))
                                    UHigherPrio += Vchain_in.get(k,h)*STeff.get(k,h)*Xchain_in.get(0,h);
                                prioScaling = Math.min(Math.max(options.tol, 1 - UHigherPrio), 1 - options.tol);
                        }

                        Bk = new Matrix(1,K);
                        if (nservers.get(k, 0) > 1) {
                            //Compute load: deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)
                            Matrix BK_tmp = new Matrix(1, K);
                            for(int i = 0; i < K; i++)
                                BK_tmp.set(0, i, deltaclass.get(0, i) * Xchain_in.get(0, i) * Vchain_in.get(k, i) * STeff.get(k, i));

                            if (BK_tmp.elementSum() < 0.75) {
                                switch (options.config.multiserver) {
                                    case "softmin":
                                        Bk = BK_tmp.clone();
                                        break;
                                    case "default":
                                    case "seidmann":
                                        BK_tmp.divide(nservers.get(k,0), Bk, true);
                                }
                            } else {
                                switch (options.config.multiserver) {
                                    case "softmin":
                                        Bk = BK_tmp.power(nservers.get(k,0));
                                        break;
                                    case "default":
                                    case "seidmann":
                                        BK_tmp.divide(nservers.get(k,0), BK_tmp, true);
                                        Bk = BK_tmp.power(nservers.get(k,0));
                                }
                            }
                        } else {
                            Bk.fill(1);
                        }

                        if (nservers.get(k, 0) == 1 && (((lldscaling != null) && !lldscaling.isEmpty()) || ((cdscaling != null) && (cdscaling.size() != 0)))) {
                            if (options.config.highvar.equals("hvmva")) {
                                double sum_uchain_r = 0;
                                for (Integer s : ccl)
                                    sum_uchain_r += Uchain_r.get(k, s);
                                Wchain.set(k, r, (STeff.get(k, r) / prioScaling) * (1 - sum_uchain_r));
                                for (Integer s : ccl) {
                                    double UHigherPrio_s = 0;
                                    for (Integer h : nnzclasses_hprio.get(s))
                                        UHigherPrio_s += Vchain_in.get(k, h) * STeff.get(k, h) * (Xchain_in.get(0, h) - Qchain_in.get(k, h) * tau.get(h));
                                    double prioScaling_s = Math.min(Math.max(options.tol, 1 - UHigherPrio_s), 1 - options.tol);
                                    Wchain.set(k, r, Wchain.get(k, r) + (STeff.get(k, s) / prioScaling_s) * Uchain_r.get(k, s) * (1 + SCVchain_in.get(k, s)) / 2);
                                }
                            } else {
                                Wchain.set(k, r, STeff.get(k, r) / prioScaling);
                            }

                            if (ocl.contains(r)) {
                                Wchain.set(k, r, Wchain.get(k, r) + (STeff.get(k, r) * stationaryQlen.get(k, r)) / prioScaling);
                            } else {
                                //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                //    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                //    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                                Wchain.set(k, r, Wchain.get(k, r) + (STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r)) / prioScaling);
                            }
                        } else {
                            switch (options.config.multiserver) {
                                case "softmin":
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * stationaryQlen.get(k, r) * Bk.get(0, r) / prioScaling);
                                    } else {
                                        //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                        //    %Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sdprio) * (stationaryQlen(k,sdprio) .* Bk(sdprio))' + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                        //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) / prioScaling + (STeff(k,[r]).*Nchain([r])*permute(gamma(r,k,[r]),3:-1:1) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                                        Wchain.set(k, r, (STeff.get(k, r) / prioScaling) + STeff.get(k, r) * selfArvlQlenSeenByClosed.get(k, r) * Bk.get(0, r) / prioScaling);
                                    }
                                    break;
                                case "seidmann":
                                case "default":
                                    if (ocl.contains(r)) {
                                        Wchain.set(k, r, (STeff.get(k,r) * (nservers.get(k,0) - 1) / prioScaling) + (STeff.get(k,r) / prioScaling) + (STeff.get(k,r) * stationaryQlen.get(k,r) * Bk.get(0, r)) / prioScaling);
                                    } else {
                                        //case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                        //    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sdprio).*Bk(sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                        //    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r)/prioScaling + (STeff(k,r).*Nchain(r)*permute(gamma(r,k,r),3:-1:1) - STeff(k,r)*gamma(r,k,r))/prioScaling;
                                        Wchain.set(k, r, (STeff.get(k,r) * (nservers.get(k,0) - 1) / prioScaling) + (STeff.get(k,r) / prioScaling) + (STeff.get(k,r) * selfArvlQlenSeenByClosed.get(k,r) * Bk.get(0, r)) / prioScaling);
                                    }
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        return new Pair<Matrix, Matrix>(Wchain, STeff);
    }
}
