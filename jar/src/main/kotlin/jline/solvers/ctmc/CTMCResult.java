package jline.solvers.ctmc;

import jline.lang.nodes.StatefulNode;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Map;

public class CTMCResult extends SolverResult {
    public String solver;
    public CTMCResult.Prob prob;
    
    public CTMCResult() {
        super();
        this.prob = new Prob();
        this.Tran = new TRAN();
        this.Tran.Avg = new AVG();
        this.tranProb = new TranProbResult();
        this.tranProbAggr = new TranProbAggrResult();
        this.tranProbSys = new TranProbSysResult();
        this.tranProbSysAggr = new TranProbSysAggrResult();
    }
    public Matrix space;
    public Matrix infGen;
    public Matrix spaceAggr;
    public Map<StatefulNode, Matrix> nodeSpace;
    public MatrixCell eventFilt;
    public TRAN Tran;
    public Matrix solverSpecific;
    
    // Transient probability results
    public TranProbResult tranProb;
    public TranProbAggrResult tranProbAggr;
    public TranProbSysResult tranProbSys;
    public TranProbSysAggrResult tranProbSysAggr;

    public class Prob {
        public Double logNormConstAggr;
        public Matrix marginal;
        public Matrix joint;
    }

    public class TRAN {
        public AVG Avg;
    }

    public class AVG {
        public Map<Integer, Map<Integer, Matrix>> Q;
        public Map<Integer, Map<Integer, Matrix>> U;
        public Map<Integer, Map<Integer, Matrix>> T;
    }
    
    public class TranProbResult {
        public Matrix t;
        public Matrix pit;
        public Matrix stateSpace;
    }
    
    public class TranProbAggrResult {
        public Matrix t;
        public Matrix pit;
        public int node;
    }
    
    public class TranProbSysResult {
        public Matrix t;
        public Matrix pit;
        public Matrix stateSpace;
    }
    
    public class TranProbSysAggrResult {
        public Matrix t;
        public Matrix pit;
        public Matrix stateSpaceAggr;
    }
}
