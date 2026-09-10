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

    @Override
    public void reset() {
        super.reset();
        this.Tran = new TRAN();
        this.Tran.Avg = new AVG();
        this.tranProb = new TranProbResult();
        this.tranProbAggr = new TranProbAggrResult();
        this.tranProbSys = new TranProbSysResult();
        this.tranProbSysAggr = new TranProbSysAggrResult();
        this.space = null;
        this.infGen = null;
        this.spaceAggr = null;
        this.pi = null;
        this.spaceWork = null;
        this.spaceAggrWork = null;
        this.infGenWork = null;
        this.nodeSpace = null;
        this.eventFilt = null;
        this.startFilt = null;
        this.preemptFilt = null;
        this.startRate = null;
        this.preemptRate = null;
        this.solverSpecific = null;
        this.cftpSamples = null;
        this.cftpHorizon = null;
    }
    public Matrix space;
    public Matrix infGen;
    public Matrix spaceAggr;
    /**
     * Stationary distribution over the state space the analyzer solved on, and
     * the generator and state spaces it is indexed by. On a reducible chain the
     * analyzer restricts to the component supporting pi, so infGenWork/spaceWork/
     * spaceAggrWork are a subset of the rows of infGen/space/spaceAggr; on an
     * irreducible chain they are the same rows. They are one triple: a consumer
     * adopts all three or none, never a mix, or pi and the generator disagree on
     * what a row index means.
     */
    public Matrix pi;
    public Matrix spaceWork;
    public Matrix spaceAggrWork;
    public Matrix infGenWork;
    public Map<StatefulNode, Matrix> nodeSpace;
    public MatrixCell eventFilt;
    /**
     * Derived START/PREEMPT filtrations, indexed [station][class]. Kept apart
     * from eventFilt, which pairs one-to-one with sn.sync and is summed as D1:
     * a START rides on an arc eventFilt already carries.
     */
    public Matrix[][] startFilt;
    public Matrix[][] preemptFilt;
    /** (stations x classes) rates the two filtrations above reduce to. */
    public Matrix startRate;
    public Matrix preemptRate;
    public TRAN Tran;
    public Matrix solverSpecific;
    /**
     * Perfect-sampling output of the {@code cftp} method: the drawn states, one
     * per row (samples x stations), and the per-sample coalescence horizon
     * ({@code cftp}) or number of mixing steps ({@code cftp.approx}). Null under
     * every enumerating method. Mirrors MATLAB {@code result.cftp.samples} /
     * {@code result.cftp.horizon} and Python {@code _CFTPResult.cftpSamples} /
     * {@code cftpHorizon}; the empirical probabilities of the distinct states
     * are in {@code pi}, indexed by {@code spaceAggr}.
     */
    public Matrix cftpSamples;
    public Matrix cftpHorizon;
    
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
