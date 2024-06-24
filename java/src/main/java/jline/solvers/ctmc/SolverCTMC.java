package jline.solvers.ctmc;

import jline.lang.*;
import jline.lang.constant.SolverType;
import jline.lang.nodes.StatefulNode;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;

import jline.util.Maths;
import jline.util.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.Map;

public class SolverCTMC extends NetworkSolver {

  public SolverCTMC(Network model) {
    this(model, new SolverOptions(SolverType.CTMC));
  }

  public SolverCTMC(Network model, SolverOptions options) {
    super(model, "SolverCTMC", options);
    int ist = 0;
    //    getProbAggr(ist);
  }

  @Override
  protected void runAnalyzer()
      throws IllegalAccessException, ParserConfigurationException, IOException {
    long T0 = System.currentTimeMillis();
    if (options == null) {
      options = new SolverOptions(SolverType.CTMC);
    }

    if (!Double.isInfinite(options.timespan[0]) && options.timespan[0] == options.timespan[1]) {
      System.err.printf(
          "%s: timespan is a single point, spacing by options.tol (%e).\n",
          this.getClass().getSimpleName(), options.tol);
      options.timespan[1] = options.timespan[0] + options.tol;
    }

//run analyzercheck - to fix
    this.resetRandomGeneratorSeed(options.seed);

    if (this.enableChecks && !supports(this.model).bool){
      System.err.printf(
          "%s: This model contains features not supported by the solver.\n",
          this.getClass().getSimpleName());
    }

    sn = getStruct(this);
    int M = sn.nstations;
    int K = sn.nclasses;
    Matrix NK = sn.njobs;
    double sizeEstimator = 0;

//    worst-case estimate of the state space
    for (int k = 0; k < K; k++) {
      sizeEstimator += sizeEstimator + Maths.factln(1 + NK.get(k, 0) + M - 1)
      - Maths.factln(1 + M - 1) - Maths.factln(1 + NK.get(k, 0));
    }

    if (containsInfinity(sn.njobs)){
      if (Double.isInfinite(options.cutoff)){
        System.err.printf(
            "%s: The model has open chains, it is recommended to specify a finite cutoff value, e.g., SolverCTMC(model,''cutoff'',1).\n",
            this.getClass().getSimpleName());
      }
//      need to fix
      this.options.cutoff = Math.ceil(Math.pow(7000, 1.0 / (M * K)));
      options.cutoff = Math.ceil(Math.pow(7000, 1.0 / (M * K)));
      System.err.printf(
              "%s: Setting cutoff=%d.\n",
              this.getClass().getSimpleName(), (int) this.options.cutoff);

    }

    if (Double.isInfinite(options.timespan[0])){
      Map<StatefulNode, Matrix>  s0 = sn.state;
      Map<StatefulNode, Matrix> s0prior = sn.stateprior;
      for (int ind = 1; ind < sn.nnodes; ind++){
        if (sn.isstateful.get(ind, 0) == 1){
          int isf = (int) sn.nodeToStateful.get(ind, 0);
          Matrix initialStatePrior = s0prior.get(sn.stateful.get(0));
          int maxPos = maxpos(initialStatePrior);
          Matrix initialState = s0.get(sn.stateful.get(isf)).getRowsFrom(maxPos);
          sn.state.put(sn.stateful.get(isf), initialState);
        }
      }
//      call solver_ctmc_analyzer(sn, options)

    }
  }

  private int maxpos(Matrix matrix) {
    double max = Double.NEGATIVE_INFINITY;
    int maxPos = 0;
    for (int i = 0; i < matrix.getNumCols(); i++){
      double value = matrix.get(0, i);
      if (value > max){
        max = value;
        maxPos = i;
      }
    }
    return maxPos;
  }

  private boolean containsInfinity(Matrix matrix) {
    for (int i = 0; i < matrix.getNumRows(); i++){
      for (int j = 0; j < matrix.getNumCols(); j++){
        if (Double.isInfinite(matrix.get(i, j))){
          return true;
        }
      }
    }
    return false;
  }






  private NetworkStruct getStruct(SolverCTMC solverCTMC) {
//    return new NetworkStruct();
    return this.model.getStruct(true);
  }

  public static FeatureSet getFeatureSet(){
    FeatureSet featSupported = new FeatureSet();
    featSupported.setTrue(new String[] {
            "Source", "Sink",
            "ClassSwitch", "DelayStation", "Queue",
            "MAP", "APH", "MMPP2", "PH", "Coxian", "Erlang", "Exp", "HyperExp",
            "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
            "Cache", "CacheClassSwitcher",
            "Server", "JobSink", "RandomSource", "ServiceTunnel",
            "SchedStrategy_INF", "SchedStrategy_PS",
            "SchedStrategy_DPS", "SchedStrategy_GPS",
            "SchedStrategy_SIRO", "SchedStrategy_SEPT",
            "SchedStrategy_LEPT", "SchedStrategy_FCFS",
            "SchedStrategy_HOL", "SchedStrategy_LCFS",
            "SchedStrategy_LCFSPR",
            "RoutingStrategy_RROBIN",
            "RoutingStrategy_PROB", "RoutingStrategy_RAND",
            "ClosedClass", "OpenClass", "Replayer"
    });
    return featSupported;
  }

  public static SupportResult supports(Network model) {
    FeatureSet featUsed = model.getUsedLangFeatures();
    FeatureSet featSupported = SolverCTMC.getFeatureSet();
    boolean bool = SolverFeatureSet.supports(featSupported, featUsed);
    return new SupportResult(bool, featSupported, featUsed);
  }

  public static class SupportResult {
    public boolean bool;
    public FeatureSet featSupported;
    public FeatureSet featUsed;
    public SupportResult(boolean bool, FeatureSet featSupported, FeatureSet featUsed) {
      this.bool = bool;
      this.featSupported = featSupported;
      this.featUsed = featUsed;
    }
  }
}

//  private void getProbAggr(int ist) {
//    if (ist < 0) {
//      throw new IllegalArgumentException("getProb requires a non-negative parameter for the
// station of interest.");
//    }
//
//    if (ist > this.model.getNumberOfStations()) {
//      throw new RuntimeException("Station number exceeds the number of stations in the model.");
//    }
//
//    long startTimeMillis = System.currentTimeMillis();
//    NetworkStruct qn = getStruct();
//
//    qn.state.putAll(this.model.getState());
//
//
//
//
//    long endTimeMillis = System.currentTimeMillis();
//    this.result.runtime = (endTimeMillis - startTimeMillis) / 1000.0;
//  }
//
//
//  private NetworkStruct getStruct() {
//    if (this.sn == null)
//      this.sn = this.model.getStruct(false);
//    return this.sn;
//  }
// }

//
//
// public Map<Node, Map<JobClass, Double>> cutoffMatrix;
// public Map<Node, Double> nodeCutoffMatrix;
// public EventStack eventStack;
// private final SolverCTMCResult ctmcResult;
//
// private Map<SSAStateMatrix,Integer> indexMap = null;
// private int stateN = 0;
// private SSAStateMatrix initialNetworkState;
//
// public SolverCTMC(Network model) {
//    this(model, new SolverOptions(SolverType.CTMC));
// }
//
// public SolverCTMC(Network model, SolverOptions options) {
//    super(model, "SolverCTMC", options);
//    this.ctmcResult = new SolverCTMCResult();
//    computeInitialStateMatrix();
//
//    this.eventStack = new EventStack();
//
//    // loop through each node and add active events to the eventStack
//    ListIterator<Node> nodeIter = model.getNodes().listIterator();
//    int nodeIdx = -1;
//    while (nodeIter.hasNext()) {
//        Node node = nodeIter.next();
//        if (!(node instanceof StatefulNode)) {
//            continue;
//        }
//
//        nodeIdx++;
//        Iterator<JobClass> jobClassIter = model.getClasses().listIterator();
//
//        while (jobClassIter.hasNext()) {
//            JobClass jobClass = jobClassIter.next();
//            int jobClassIdx = jobClass.getJobClassIdx();
////                if (network.getClassLinks(node, jobClass) == 0) {
////                    this.simStruct.classcap[nodeIdx][jobClassIdx] = 0;
////                } else {
////                    double jobCap = jobClass.getNumberOfJobs();
////                    jobCap = Math.min(jobCap, node.getClassCap(jobClass));
////                    if ((jobCap == Double.POSITIVE_INFINITY) || (node.getDropStrategy() ==
// DropStrategy.WaitingQueue)) {
////                        this.simStruct.classcap[nodeIdx][jobClassIdx] = Integer.MAX_VALUE;
////                    } else {
////                        this.simStruct.classcap[nodeIdx][jobClassIdx] = (int) jobCap;
////                    }
////                }
//            Event dEvent = DepartureEvent.fromNodeAndClass(node, jobClass);
//            this.eventStack.addEvent(dEvent);
//            if (dEvent instanceof DepartureEvent) {
//                if (((DepartureEvent) dEvent).getPhaseEvent() != null) {
//                    this.eventStack.addEvent(((DepartureEvent) dEvent).getPhaseEvent());
//                }
//            }
//        }
//
////            double nodeCap = node.getCap();
////            if (nodeCap == Double.POSITIVE_INFINITY) {
////                this.simStruct.cap[nodeIdx] = Integer.MAX_VALUE;
////            } else {
////                this.simStruct.cap[nodeIdx] = (int) nodeCap;
////            }
//    }
// }
//
// public void computeInitialStateMatrix() {
//    SSAStateMatrix networkState = new SSAStateMatrix(this.sn,this.random);
//    for (JobClass jobClass : this.model.getClasses()) {
//        if (jobClass instanceof ClosedClass) {
//            int classIdx = this.model.getJobClassIndex(jobClass);
//            ClosedClass cClass = (ClosedClass) jobClass;
//            int stationIdx = this.model.getStatefulNodeIndex(cClass.getRefstat());
//            networkState.setState(stationIdx, classIdx, (int)cClass.getPopulation());
//            for (int i = 0; i < cClass.getPopulation(); i++) {
//                networkState.addToBuffer(stationIdx, classIdx);
//            }
//        }
//    }
//    this.initialNetworkState = networkState;
// }
//
// public SolverOptions setOptions() {
//    return this.options;
// }
//
// public ArrayList<SSAStateMatrix> getStateSpace() {
//    if(ctmcResult.stateSpace == null) {
//        solver_ctmc();
//    }
////        System.out.println("\nStateSpace =");
////        for(SSAStateMatrix state : ctmcResult.stateSpace) {
////            state.printStateVector();
////        }
//    return ctmcResult.stateSpace;
// }
//
// public Matrix getGenerator() {
//    if(ctmcResult.infGen == null) {
//        solver_ctmc();
//    }
////        System.out.println("\nInfGen =");
////        for (int i = 0; i < ctmcResult.infGen.getNumRows(); i++) {
////            for (int j = 0; j < ctmcResult.infGen.getNumCols(); j++) {
////                System.out.print(ctmcResult.infGen.get(i, j) + "  ");
////            }
////            System.out.println();
////        }
//    return ctmcResult.infGen;
// }
//
// public void getProbabilityVector() {
//    if(ctmcResult.piVector == null) {
//        solver_ctmc_analyzer();
//    }
//    System.out.println("\nProbability Vector =");
//    ctmcResult.piVector.print();
// }
//
//
// @Override
// public void runAnalyzer(){
//    if (this.model == null)
//        throw new RuntimeException("Model is not provided");
//    if (this.sn == null)
//        this.sn = this.model.getStruct(false);
//
////        simCache.applyCutoff(simOptions, model);
//    // Add ClosedClass instances to the reference station
////        if(this.options.timespan[1] == 0) {
//    solver_ctmc_analyzer();
////        }
////        else {
//    //TODO transient analysis
////        }
// }
//
// public void applyCutoff(double cutoff){
//    this.cutoffMatrix = new HashMap<Node, Map<JobClass, Double>>();
//
//    this.options.cutoff = cutoff;
//    for (int i = 0; i < this.sn.nstateful; i++) {
//        Node nodeIter = model.getStatefulNodeFromIndex(i);
//        for (int j = 0; j < this.sn.nclasses; j++) {
//            JobClass jobClassIter = model.getJobClassFromIndex(j);
//            double cutoffVal = cutoff;
//            if(this.cutoffMatrix.get(nodeIter) != null &&
// this.cutoffMatrix.get(nodeIter).get(jobClassIter) != null) {
//                cutoffVal = Math.min(cutoffVal,
// this.cutoffMatrix.get(nodeIter).get(jobClassIter));
//            }
//            if (cutoffVal != Double.POSITIVE_INFINITY) {
//                this.sn.classcap.set(Math.min((int)this.sn.classcap.get(i,j), (int) cutoff),i,j);
//            }
//        }
//        double nodeCutoff = cutoff;
//        if(this.cutoffMatrix.get(nodeIter) != null) {
//            nodeCutoff = Math.min(nodeCutoff, this.nodeCutoffMatrix.get(nodeIter));
//        }
//        if (nodeCutoff != Double.POSITIVE_INFINITY) {
//            this.sn.cap.set(Math.min((int)this.sn.cap.get(i), (int) nodeCutoff),i);
//        }
//    }
// }
//
// public void solver_ctmc_analyzer() {
//    long startTime = System.currentTimeMillis();
//    solver_ctmc();
//
//    // Compute probability vector
//
//    ctmcResult.piVector = CTMC.ctmc_solve(ctmcResult.infGen);
//    this.stateN = this.initialNetworkState.state.length;
//
//    int nClasses = model.getNumberOfClasses();
//
//    double[][] ssprobabilities = ctmcResult.piVector.toArray2D();
//    double [][][] arrRates = new double[stateN][nClasses][ctmcResult.stateSpace.size()];
//    Matrix UN = new Matrix(stateN, nClasses);
//    Matrix QN = new Matrix(stateN, nClasses);
//    Matrix TN = new Matrix(stateN, nClasses);
//    Matrix RN = new Matrix(stateN, nClasses);
//    Matrix XN = new Matrix(1, nClasses);
//    Matrix CN = new Matrix(1, nClasses);
//
//
//
//    // q length
//    for (int i=0;i<ctmcResult.stateSpace.size();i++){
//        if(ssprobabilities[0][i]>0){
//            SSAStateMatrix currState = ctmcResult.stateSpace.get(i);
//            int[][] state = currState.state;
//            for(int j = 0;j<state.length;j++){
//                for(int k=0;k<state[j].length;k++){
//                    QN.set(j, k, QN.get(j, k) + ssprobabilities[0][i]*state[j][k]);
//                }
//            }
//        }
//    }
//
//    // arrival rates, throughput
//    for(EventData eventData : ctmcResult.eventSpace){
//        int eventIdx = eventData.getValue0().getNode().getNodeIdx();
//        int event2Idx;
//        int classIdx = eventData.getValue1().getLeft().getClassIdx();
//
//
//        //departure event
//        if(!eventData.getValue1().getLeft().isDummy()) {
//            event2Idx = eventData.getValue1().getLeft().getNode().getNodeIdx();
//            SSAStateMatrix depState = eventData.getValue2();
//            double rate = eventData.getValue1().getRight();
//            double prob = ctmcResult.piVector.get(indexMap.get(depState));
//            if(!(eventData.getValue0() instanceof ErlangPhaseEvent)) {
//                rate *= eventData.getValue0().getRate(eventData.getValue2());
//            }
//            else if(!eventData.getValue1().getLeft().isDummy()) {
//                if(eventData.getValue0() instanceof ErlangPhaseEvent) {
//                    rate *= ((ErlangPhaseEvent)
// eventData.getValue0()).getDepartureRate(eventData.getValue2());
//                }
//            }
//            if(eventData.getValue1().getLeft().getNode().isStateful()) {
//                arrRates[event2Idx][classIdx][indexMap.get(eventData.getValue2())] += rate;
//            }
//            TN.set(eventIdx, classIdx, TN.get(eventIdx, classIdx) + rate * prob);
//            if (sn.refstat.get(classIdx, 0) == event2Idx) {
//                XN.set(0, classIdx, XN.get(0, classIdx) + prob * rate);
//            }
//        }
//    }
//
//    for(int k = 0; k < nClasses; k++) {
//        CN.set(0, k, sn.njobs.get(0, k) / (k + 1));
//    }
//
//    // utilisation
//    for(int i = 0; i < stateN; i++) {
//        Node node = model.getStatefulNodeFromIndex(i);
//        SchedStrategy strategy = this.sn.sched.get(this.sn.stations.get(i));
//        if(!(node instanceof Source)) {
//            StatefulNode statefulNode;
//            if(node instanceof StatefulNode) {
//                statefulNode = (StatefulNode) node;
//            }
//            else {
//                continue;
//            }
//            switch (strategy) {
//                case INF:
//                    for(int k = 0; k < nClasses; k++) {
//                        UN.set(i, k, QN.get(i, k));
//                    }
//                    break;
//                case GPS:
//                case DPS:
//                case PS:
//                    if(this.model.getLimitedLoadDependence().isEmpty() &&
// (this.model.getLimitedClassDependence() == null ||
// this.model.getLimitedClassDependence().isEmpty())) {
//                        for(int k = 0; k < nClasses; k++) {
//                            Distribution dist =
// node.getServer().getServiceDistribution(this.model.getJobClassFromIndex(k));
//                            double arrEstimate = 0.0;
//                            for(int j = 0; j < ctmcResult.stateSpace.size(); j++) {
//                                arrEstimate += ctmcResult.piVector.get(j) * arrRates[i][k][j];
//                            }
//                            arrEstimate = arrEstimate * dist.getMean() /
// statefulNode.getNumberOfServers();
//                            double depEstimate = TN.get(i, k) * dist.getMean() /
// statefulNode.getNumberOfServers();
//                            UN.set(i, k, UN.get(i, k) + Math.max(arrEstimate, depEstimate));
//                        }
//                    }
//                    else {
//                        //TODO lld/cd cases not finished
//                        double arrEstimate, depEstimate;
//                        for(int j = 0; j < ctmcResult.stateSpace.size(); j++) {
//                            State.StateMarginalStatistics stats = State.toMarginal(sn, i,
// sn.state.get(sn.stations.get(i)), null, null, null, null, null);
//                            Matrix ni = stats.ni;
//                            Matrix nir = stats.nir;
//                            if(ni.get(0,0) > 0) {
//                                for(int k = 0; k < nClasses; k++) {
//                                }
//                            }
//                        }
//                    }
//                    break;
//                default:
//                    if(this.model.getLimitedLoadDependence().isEmpty() &&
// (this.model.getLimitedClassDependence() == null ||
// this.model.getLimitedClassDependence().isEmpty())) {
//                        for(int k = 0; k < nClasses; k++) {
//                            Distribution dist =
// node.getServer().getServiceDistribution(this.model.getJobClassFromIndex(k));
//                            double arrEstimate = 0.0;
//                            for(int j = 0; j < ctmcResult.stateSpace.size(); j++) {
//                                arrEstimate += ctmcResult.piVector.get(j) * arrRates[i][k][j];
//                            }
//                            arrEstimate = arrEstimate * dist.getMean() /
// statefulNode.getNumberOfServers();
//                            double depEstimate = TN.get(i, k) * dist.getMean() /
// statefulNode.getNumberOfServers();
//                            UN.set(i, k, UN.get(i, k) + Math.max(arrEstimate, depEstimate));
//                        }
//                    }
//                    else {
//                        //TODO lld/cd cases
//                    }
//                    break;
//            }
//        }
//    }
//    for(int k = 0; k < nClasses; k++) {
//        for(int i = 0; i < stateN; i++) {
//            if(TN.get(i, k) > 0) {
//                RN.set(i, k, QN.get(i, k) / TN.get(i ,k));
//            }
//            else {
//                RN.set(i, k, 0);
//            }
//        }
//    }
//    long endTime = System.currentTimeMillis();
//    long runTime = endTime - startTime;
//    result.QN = QN;
//    result.UN = UN;
//    result.TN = TN;
//    result.RN = RN;
//    result.AN = new Matrix(stateN,nClasses);
//    result.XN = XN;
//    result.CN = CN;
//    result.runtime = runTime/1000.0;
//    System.out.println("CTMC Analysis completed. Runtime: " + result.runtime + " seconds");
// }
//
// public void solver_ctmc() {
//    ArrayList<SSAStateMatrix> stateSpace = new ArrayList<>();
//    ArrayList<EventData>  eventSpace = new ArrayList<>();
//    Queue<SSAStateMatrix> queue = new LinkedList<>();
//    Set<EventData> eventSet = new HashSet<>();
//    // Compute state space
//    if(this.options.cutoff != -1) {
//        initialNetworkState.setCutoff(this.options.cutoff);
//    }
//    stateSpace.add(this.initialNetworkState);
//    queue.add(this.initialNetworkState);
//    Set<SSAStateMatrix> stateSet = new HashSet<>();
//    stateSet.add(initialNetworkState);
//    this.eventStack.updateStateSpace(stateSpace,queue, stateSet);
//    ctmcResult.stateSpace = stateSpace;
//    // Compute event space
//    eventSet.add(new EventData(null,null,this.initialNetworkState,null));
//    eventSpace.add(new EventData(null,null,this.initialNetworkState,null));
//    queue.add(this.initialNetworkState);
//
//    this.eventStack.updateEventSpace(eventSpace,queue, eventSet);
//
//    eventSpace.remove(0);
//    ctmcResult.eventSpace = eventSpace;
//    int size = stateSpace.size();
//    Map<SSAStateMatrix,Integer> indexMap = new HashMap<>();
//    for(int i = 0; i< stateSpace.size(); i++){
//        indexMap.put(stateSpace.get(i),i);
//    }
//    this.indexMap = indexMap;
//
//    // Compute inf generator
//    Matrix infGen = new Matrix(size, size);
//    for(int i =0;i<size;i++){
//        for(int j=0;j<size;j++){
//            infGen.set(i, j, 0.0);
//        }
//    }
//    for(EventData eventData : eventSpace){
//        double rate = eventData.getValue1().getRight();
//        if(!(eventData.getValue0() instanceof ErlangPhaseEvent)) {
//            rate *= eventData.getValue0().getRate(eventData.getValue2());
//        }
//        else{
//            if(!eventData.getValue1().getLeft().isDummy()) {
//                if(eventData.getValue0() instanceof ErlangPhaseEvent) {
//                    rate *= ((ErlangPhaseEvent)
// eventData.getValue0()).getDepartureRate(eventData.getValue2());
//                }
//            }
//        }
//        infGen.set(indexMap.get(eventData.getValue2()), indexMap.get(eventData.getValue3()),
// infGen.get(indexMap.get(eventData.getValue2()), indexMap.get(eventData.getValue3())) + rate);
//    }
//
//    for(int i = 0; i<size;i++){
//        double sum = 0;
//        for(int j=0;j<size;j++){
//            if(j!=i){
//                sum+=infGen.get(i, j);
//            }
//        }
//        infGen.set(i, i, - sum);
//    }
//    ctmcResult.infGen = infGen;
// }
// public static SolverOptions defaultOptions() {
//    return new SolverOptions(SolverType.CTMC);
// }
