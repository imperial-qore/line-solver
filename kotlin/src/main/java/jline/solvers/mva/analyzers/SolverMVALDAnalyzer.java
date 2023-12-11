package jline.solvers.mva.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.solvers.mva.handlers.AMVAHandler;
import jline.solvers.mva.handlers.MVALDHandler;

/**
 * MVALD Analyzer
 */
public class SolverMVALDAnalyzer implements MVAAnalyzer{
    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        long startTime = System.currentTimeMillis();
        SolverMVAResult ret = null;
        switch(options.method){
            case "exact": case "mva":
                if(sn.cdscaling != null && !sn.cdscaling.isEmpty()){
                    throw new RuntimeException("Exact class-dependent solver not available in MVA.");
                }
                ret = new MVALDHandler().solve(sn, options);
                break;
            case "default": case "amva": case "qd": case "aql": case "qdaql": case "lin": case "qdlin":
                ret = new AMVAHandler().solve(sn, options);
                break;
            default:
                throw new RuntimeException("The " + options.method + " method is not supported by the load-dependent MVA solver.");
        }
        long endTime = System.currentTimeMillis();
        res.QN = ret.QN;
        res.UN = ret.UN;
        res.RN = ret.RN;
        res.TN = ret.TN;
        res.CN = ret.CN;
        res.XN = ret.XN;
        res.AN = ret.AN;
        res.WN = ret.WN;
        res.logNormConstAggr = ret.logNormConstAggr;
        res.runtime = (endTime - startTime) / 1000.0;
        res.iter = ret.iter;
    }
}
