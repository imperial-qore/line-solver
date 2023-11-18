package jline.solvers.mva.analyzers;

import jline.util.Maths;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.util.Matrix;

import static jline.api.PFQN.*;

/**
 * MVA Analyzer class for bounding methods
 */
public class SolverMVABoundAnalyzer implements MVAAnalyzer{

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        long startTime = System.currentTimeMillis();
        long endTime = startTime;
        int iter = 1;
        Matrix QN = new Matrix(0,0);
        Matrix UN = new Matrix(0,0);
        Matrix RN = new Matrix(0,0);
        Matrix TN = new Matrix(0,0);
        Matrix CN = new Matrix(0,0);
        Matrix WN = new Matrix(0,0);
        Matrix XN = new Matrix(0,0);
        double lG = Double.NaN;
        switch(options.method){
            case "aba.upper":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    double Dsum = D.elementSum();
                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + N * Dsum);
                    XN = new Matrix(1,1);
                    XN.set(0,0, Maths.min(1/Dmax, N / (Z + Dsum)));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            RN.set(i, 0, 1 / sn.rates.get(i));
                        } else {
                            RN.set(i,0, 1/sn.rates.get(i) * N);
                        }
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                           UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "aba.lower":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dsum = D.elementSum();
                    XN = new Matrix(1,1);
                    XN.set(0,0, N / (Z + N * Dsum));
                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + Dsum);
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        RN.set(i, 0, 1 / sn.rates.get(i));
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "bjb.upper":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    double Dsum = D.elementSum();

                    double Xaba_upper_1 = Maths.min(1/Dmax, (N - 1)/(Z + Dsum));
                    double Xaba_lower_1 = (N - 1)/(Z + (N - 1) * Dsum);

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1));
                    XN = new Matrix(1,1);
                    XN.set(0,0, Maths.min(1/Dmax, N / (Z + Dsum + D.meanCol().get(0,0) * (N - 1 - Z * Xaba_upper_1))));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA upper
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            RN.set(i, 0, 1 / sn.rates.get(i));
                        } else {
                            RN.set(i,0, 1/sn.rates.get(i) * N);
                        }
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "bjb.lower":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    double Dsum = D.elementSum();

                    double Xaba_upper_1 = Maths.min(1/Dmax, (N - 1)/(Z + Dsum));
                    double Xaba_lower_1 = (N - 1)/(Z + (N - 1) * Dsum);

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + Dsum + D.meanCol().get(0,0) * (N - 1 - Z * Xaba_upper_1));
                    XN = new Matrix(1,1);
                    XN.set(0,0,N / (Z + Dsum + Dmax * (N - 1 - Z * Xaba_lower_1)));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA lower
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        RN.set(i, 0, 1 / sn.rates.get(i));
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "pb.upper":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    double Dsum = D.elementSum();

                    double Xaba_upper_1 = Maths.min(1/Dmax, (N - 1)/(Z + Dsum));
                    double Xaba_lower_1 = (N - 1)/(Z + (N - 1) * Dsum);

                    double Dpb2 = D.power(2).elementSum() / Dsum;
                    Double DpbN = D.power(N).elementSum() / D.power(N - 1).elementSum();

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1));
                    XN = new Matrix(1,1);
                    XN.set(0,0, Maths.min(1/Dmax, N / (Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1))));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA upper
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            RN.set(i, 0, 1 / sn.rates.get(i));
                        } else {
                            RN.set(i,0, 1/sn.rates.get(i) * N);
                        }
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "pb.lower":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    double Dsum = D.elementSum();

                    double Xaba_upper_1 = Maths.min(1/Dmax, (N - 1)/(Z + Dsum));
                    double Xaba_lower_1 = (N - 1)/(Z + (N - 1) * Dsum);

                    double Dpb2 = D.power(2).elementSum() / Dsum;
                    Double DpbN = D.power(N).elementSum() / D.power(N - 1).elementSum();

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + Dsum + Dpb2 * (N - 1 - Z * Xaba_upper_1));
                    XN = new Matrix(1,1);
                    XN.set(0,0,N / (Z + Dsum + DpbN * (N - 1 - Z * Xaba_lower_1)));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA lower
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        RN.set(i, 0, 1 / sn.rates.get(i));
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "sb.upper":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF)
                            throw new RuntimeException("Unsupported method for a model with infinite-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();

                    double A3 = D.power(3).elementSum();
                    double A2 = D.power(2).elementSum();
                    double A1 = D.power(1).elementSum();

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + A1 + (N - 1) * (A1 * A2 + A3)/(A1 * A1 + A2));
                    XN = new Matrix(1,1);
                    XN.set(0,0, Maths.min(1/Dmax, N / CN.get(0,0)));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA lower
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        RN.set(i, 0, 1 / sn.rates.get(i));
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "sb.lower":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF)
                            throw new RuntimeException("Unsupported method for a model with infinite-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();

                    double AN = D.power(N).elementSum();
                    double A1 = D.power(1).elementSum();

                    CN = new Matrix(1,1);
                    CN.set(0,0,Z + A1 + (N - 1) * Math.pow(AN/A1, 1/(N-1)));
                    XN = new Matrix(1,1);
                    XN.set(0,0,N / CN.get(0,0));
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    // RN undefined in the literature so we use ABA lower
                    RN = new Matrix(sn.rates.getNumRows(), 1);
                    for(int i = 0; i < RN.getNumRows(); i++){
                        RN.set(i, 0, 1 / sn.rates.get(i));
                    }
                    QN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        QN.set(i, 0, TN.get(i, 0) * RN.get(i, 0));
                    }
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < TN.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        } else {
                            UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "gb.upper":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    double Dmax = D.elementMax();
                    XN = new Matrix(1,1);
                    XN.set(0,0, Maths.min(1/Dmax, pfqn_xzgsbup(D, N, Z)));
                    double ret = pfqn_xzgsblow(D,N,Z);
                    CN = new Matrix(1,1);
                    CN.set(0,0, N / ret);
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    double XNlow = ret;
                    int k = 0;
                    RN = new Matrix(sn.sched.size(), 1);
                    QN = new Matrix(sn.sched.size(), 1);
                    for(int i = 0; i < sn.sched.size(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            RN.set(i, 0, 1 / sn.rates.get(i));
                            QN.set(i, 0, XN.get(0,0) * RN.get(i, 0));
                        } else {
                            QN.set(i, 0, pfqn_qzgbup(D, N, Z, k));
                            RN.set(i, 0, QN.get(i, 0) / XNlow / V.get(i));
                            k++; // increment after setting QN and RN because of 0-indexing
                        }
                    }
                    /* RN(sn.schedid == SchedStrategy.ID_INF,1) = 1 ./ sn.rates(sn.schedid == SchedStrategy.ID_INF,1);
                    *   is redundant in LINE, RN is already set in the previous for-loop*/
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < UN.getNumRows(); i++){
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    for(int i = 0; i < sn.sched.size(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
            case "gb.lower":
                if(sn.nclasses == 1 && sn.nclosedjobs > 0){
                    // Closed single-class queueing network
                    for(int i = 0; i < sn.nservers.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF && sn.nservers.get(i) > 1)
                            throw new RuntimeException("Unsupported method for a model with multi-server stations.");
                    }
                    Matrix V = sn.visits.get(0);
                    double Z = 0;
                    int nondelays = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            Z += V.get(i) / sn.rates.get(i);
                        } else {
                            nondelays++;
                        }
                    }
                    Matrix D = new Matrix(nondelays, 1);
                    int idx = 0;
                    for(int i = 0; i < V.getNumRows(); i++){
                        if(sn.sched.get(sn.stations.get(i)) != SchedStrategy.INF){
                            D.set(idx, 0, V.get(i) / sn.rates.get(i));
                            idx++;
                        }
                    }
                    double N = sn.nclosedjobs;
                    XN = new Matrix(1,1);
                    XN.set(0,0, pfqn_xzgsblow(D, N, Z));
                    double ret = pfqn_xzgsbup(D,N,Z);
                    CN = new Matrix(1,1);
                    CN.set(0,0, N / ret);
                    TN = new Matrix(V.getNumRows(), 1);
                    for(int i = 0; i < V.getNumRows(); i++){
                        TN.set(i, 0, V.get(i, 0) * XN.get(0,0));
                    }
                    double XNup = ret;
                    int k = 0;
                    RN = new Matrix(sn.sched.size(), 1);
                    QN = new Matrix(sn.sched.size(), 1);
                    for(int i = 0; i < sn.sched.size(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            RN.set(i, 0, 1 / sn.rates.get(i));
                            QN.set(i, 0, XN.get(0,0) * RN.get(i, 0));
                        } else {
                            QN.set(i, 0, pfqn_qzgblow(D, N, Z, k));
                            RN.set(i, 0, QN.get(i, 0) / XNup / V.get(i));
                            k++; // increment after setting QN and RN because of 0-indexing
                        }
                    }
                    /* RN(sn.schedid == SchedStrategy.ID_INF,1) = 1 ./ sn.rates(sn.schedid == SchedStrategy.ID_INF,1);
                     *   is redundant in LINE, RN is already set in the previous for-loop*/
                    UN = new Matrix(TN.getNumRows(), 1);
                    for(int i = 0; i < UN.getNumRows(); i++){
                        UN.set(i, 0, TN.get(i, 0) / sn.rates.get(i));
                    }
                    for(int i = 0; i < sn.sched.size(); i++){
                        if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF){
                            UN.set(i, 0, QN.get(i, 0));
                        }
                    }
                    lG = -N * Math.log(XN.get(0,0));
                }
                endTime = System.currentTimeMillis();
                break;
        }
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = new Matrix(0,0);
        res.WN = WN;
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000.0;
        res.iter = iter;
    }
}
