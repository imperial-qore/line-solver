package jline.solvers.mva.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;
import jline.solvers.mva.handlers.AMVAHandler;
import jline.solvers.mva.handlers.MVAHandler;
import jline.solvers.mva.handlers.QNAHandler;

/**
 * MVA Analyzer class
 */
public class SolverMVAAnalyzer implements MVAAnalyzer{

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res) {
        int iter = 0;
        long startTime = System.currentTimeMillis();
        SolverMVAResult ret = null;
        switch(options.method){
            case "exact": case "mva":
                ret = new MVAHandler().solve(sn, options);
                ret.iter = 0;
                break;
            case "qna":
                ret = new QNAHandler().solve(sn, options);
                break;
            case "default":
                ret = new AMVAHandler().solve(sn, options);
                break;
            case "amva": case "bs": case "qd": case "qli": case "fli": case "aql": case "qdaql": case "lin":
            case "qdlin": case "sqni": case "gflin": case "egflin":
                ret = new AMVAHandler().solve(sn, options);
                break;
            default:
                throw new RuntimeException("Unsupported SolverMVA method.");
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
        res.iter = ret.iter;
        res.runtime = (endTime - startTime) / 1000.0;
    }
}
