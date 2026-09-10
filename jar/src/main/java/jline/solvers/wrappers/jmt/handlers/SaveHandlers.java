package jline.solvers.wrappers.jmt.handlers;

import static jline.GlobalConstants.Inf;

import jline.api.mam.*;
import jline.io.DocumentSectionPair;
import jline.io.ElementDocumentPair;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodeparam.*;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import jline.lang.sections.PollingServer;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;
import jline.lang.nodes.Transition;
import jline.solvers.AvgHandle;
import jline.solvers.SolverAvgHandles;
import jline.util.Maths;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Handles the generation and serialization of JMT (Java Modelling Tools) simulation models.
 * 
 * <p>This class is responsible for converting LINE network models into JMT XML format
 * for discrete-event simulation. It provides comprehensive support for translating 
 * various network components including nodes, routing strategies, service disciplines,
 * and performance metrics.
 * 
 * <p>The class supports both JSIMg (JMT Simulation Graph) and JSIM (JMT Simulation)
 * XML formats, with JSIMg being the default format that includes additional simulation
 * parameters such as confidence intervals, maximum relative error, and stopping criteria.
 * 
 * <h3>Key Features:</h3>
 * <ul>
 * <li>Conversion of LINE Network models to JMT XML format</li>
 * <li>Support for complex routing strategies and class switching</li>
 * <li>Handling of various node types (Sources, Queues, Delays, Sinks, etc.)</li>
 * <li>Translation of service disciplines and scheduling strategies</li>
 * <li>Cache modeling with different replacement strategies</li>
 * <li>Fork-join synchronization patterns</li>
 * <li>Performance metric collection and measurement configuration</li>
 * </ul>
 * 
 * <h3>Supported Node Types:</h3>
 * <ul>
 * <li>Sources with various arrival processes</li>
 * <li>Queues with different scheduling disciplines</li>
 * <li>Delay stations for pure delays</li>
 * <li>Sinks for job termination</li>
 * <li>Routers for probabilistic routing</li>
 * <li>Forks and Joins for synchronization</li>
 * <li>Caches with replacement strategies</li>
 * <li>Transitions for Petri net models</li>
 * </ul>
 * 
 * @see jline.lang.Network
 * @see jline.solvers.wrappers.jmt.SolverJMT
 * @see jline.lang.NetworkStruct
 * @see jline.solvers.SolverAvgHandles
 */
public class SaveHandlers {
    private NetworkStruct sn;
    private final double simMaxRelErr;
    private final double simConfInt;
    private final Network simModel;
    private final SolverAvgHandles avgHandles;
    private final long seed;
    private final String simFileName;
    private final long maxEvents;
    private final long maxSamples;
    private final double maxSimulatedTime;

    /**
     * Constructs a SaveHandlers instance with full simulation parameter control.
     * 
     * <p>This constructor allows complete customization of all simulation parameters
     * for JMT model generation. It is typically used when specific simulation
     * configuration is required beyond the default settings.
     * 
     * @param simModel The LINE network model to be converted to JMT format
     * @param simMaxRelErr Maximum relative error for simulation stopping criterion (e.g., 0.03 for 3%)
     * @param simConfInt Confidence interval for simulation results (e.g., 0.99 for 99%)
     * @param avgHandles Collection of performance metric handlers to be included in the simulation
     * @param seed Random number generator seed for reproducible simulation results
     * @param simFileName Output filename for the simulation model (empty string for default naming)
     * @param maxEvents Maximum number of events to simulate (0 for unlimited)
     * @param maxSamples Maximum number of samples to collect for each metric
     * @param maxSimulatedTime Maximum simulation time (use GlobalConstants.Inf for unlimited)
     */
    public SaveHandlers(Network simModel, double simMaxRelErr, double simConfInt,
                        SolverAvgHandles avgHandles, long seed,
                        String simFileName, long maxEvents, long maxSamples, double maxSimulatedTime) {
        this.simModel = simModel;
        this.sn = simModel.getStruct();
        this.simMaxRelErr = simMaxRelErr;
        this.simConfInt = simConfInt;
        this.avgHandles = avgHandles;
        this.seed = seed;
        this.simFileName = simFileName;
        this.maxEvents = maxEvents;
        this.maxSamples = maxSamples;
        this.maxSimulatedTime = maxSimulatedTime;
    }

    /**
     * Constructs a SaveHandlers instance with default simulation parameters.
     * 
     * <p>This convenience constructor uses SolverJMT default values for all simulation
     * parameters, making it suitable for most standard simulation scenarios. The
     * default configuration provides reasonable accuracy and performance for typical
     * queueing network analysis.
     * 
     * <p><strong>Default Parameters:</strong>
     * <ul>
     * <li>Maximum relative error: 3%</li>
     * <li>Confidence interval: 99%</li>
     * <li>Random seed: 23000</li>
     * <li>Maximum samples: 10,000</li>
     * <li>Maximum events: unlimited</li>
     * <li>Maximum simulation time: unlimited</li>
     * </ul>
     * 
     * @param simModel The LINE network model to be converted to JMT format
     */
    public SaveHandlers(Network simModel) {
        this(simModel, 
             0.03,                           // simMaxRelErr (SolverJMT default)
             0.99,                           // simConfInt (SolverJMT default)
             simModel.getAvgHandles(),       // avgHandles (from model)
             23000L,                         // seed (specified default)
             "",                             // simFileName (empty default)
             0L,                             // maxEvents (unlimited)
             10000L,                         // maxSamples (SolverOptions default)
             Inf        // maxSimulatedTime (unlimited)
        );
    }
    
    /**
     * Updates the internal network structure used for model conversion.
     * 
     * <p>This method allows updating the network structure after the SaveHandlers
     * instance has been created. This is useful when the network model has been
     * modified or when working with multiple network configurations.
     * 
     * @param sn The new network structure to use for subsequent conversions
     */
    public void updateNetworkStruct(NetworkStruct sn) {
        this.sn = sn;
    }

    /**
     * Returns a boolean array indicating which classes should be exported to JMT.
     * Classes with 0 customers are normally skipped unless they are used as
     * cache hit/miss classes (since jobs will switch into them) or can receive
     * jobs via class-switching from another class in the same chain.
     *
     * @return boolean array of length nclasses, true for classes to export
     */
    public boolean[] getExportableClasses() {
        int numOfClasses = sn.nclasses;

        // Build set of classes used as cache hit/miss classes
        Set<Integer> cacheClasses = new HashSet<Integer>();
        for (Map.Entry<Node, NodeParam> entry : sn.nodeparam.entrySet()) {
            NodeParam param = entry.getValue();
            if (param instanceof CacheNodeParam) {
                CacheNodeParam cacheParam = (CacheNodeParam) param;
                // hitclass/missclass are 0-indexed in JAR (-1 = invalid)
                if (cacheParam.hitclass != null) {
                    for (int k = 0; k < cacheParam.hitclass.length(); k++) {
                        int h = (int) cacheParam.hitclass.get(k);
                        if (h >= 0) {
                            cacheClasses.add(h);
                        }
                    }
                }
                if (cacheParam.missclass != null) {
                    for (int k = 0; k < cacheParam.missclass.length(); k++) {
                        int m = (int) cacheParam.missclass.get(k);
                        if (m >= 0) {
                            cacheClasses.add(m);
                        }
                    }
                }
            }
        }

        // Build set of classes that can receive jobs via class-switching
        Set<Integer> classSwitchClasses = new HashSet<Integer>();
        for (int c = 0; c < sn.nchains; c++) {
            List<Integer> classesInChain = new ArrayList<Integer>();
            for (int r = 0; r < numOfClasses; r++) {
                if (sn.chains.get(c, r) > 0) {
                    classesInChain.add(r);
                }
            }
            if (classesInChain.size() > 1) {
                boolean chainHasJobs = false;
                for (int r : classesInChain) {
                    if (sn.njobs.get(r) > 0) {
                        chainHasJobs = true;
                        break;
                    }
                }
                if (chainHasJobs) {
                    classSwitchClasses.addAll(classesInChain);
                }
            }
        }

        // If the model has ClassSwitch nodes, all classes must be exported because
        // the ClassSwitch matrix must include all classes (even 0-job classes that
        // pass through with identity routing)
        boolean hasClassSwitch = false;
        for (NodeType t : sn.nodetype) {
            if (t == NodeType.ClassSwitch) {
                hasClassSwitch = true;
                break;
            }
        }

        // Determine which classes to export
        boolean[] exportClasses = new boolean[numOfClasses];
        for (int r = 0; r < numOfClasses; r++) {
            // Skip closed classes with 0 customers unless they are cache classes,
            // can receive jobs via class-switching, or the model has ClassSwitch nodes
            if (!Double.isInfinite(sn.njobs.get(r)) && sn.njobs.get(r) == 0
                    && !cacheClasses.contains(r) && !classSwitchClasses.contains(r)
                    && !hasClassSwitch) {
                exportClasses[r] = false;
            } else {
                exportClasses[r] = true;
            }
        }

        return exportClasses;
    }

    private String cacheStrategyMap(ReplacementStrategy method) {
        Map<ReplacementStrategy, String> cacheStrategyMap = new HashMap<>();
        cacheStrategyMap.put(ReplacementStrategy.FIFO, "jmt.engine.NetStrategies.CacheStrategies.FIFOCache");
        cacheStrategyMap.put(ReplacementStrategy.RR, "jmt.engine.NetStrategies.CacheStrategies.RandomCache");
        cacheStrategyMap.put(ReplacementStrategy.LRU, "jmt.engine.NetStrategies.CacheStrategies.LRUCache");
        //cacheStrategyMap.put("LFU", "jmt.engine.NetStrategies.CacheStrategies.LFUCache");
        //cacheStrategyMap.put(ReplacementStrategy.LRUM, "jmt.engine.NetStrategies.CacheStrategies.LRUCache");
        //cacheStrategyMap.put(ReplacementStrategy.HLRU, "jmt.engine.NetStrategies.CacheStrategies.HLRUCache");
        //cacheStrategyMap.put("TTL", "jmt.engine.NetStrategies.CacheStrategies.TTLCache");

        return cacheStrategyMap.get(method);
    }

    public DocumentSectionPair saveArrivalStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        strategyNode.setAttribute("name", "ServiceStrategy");
        int numOfClasses = sn.nclasses;
        int ist = (int) sn.nodeToStation.get(ind);
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode2);
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

            if (!Utils.isInf(sn.njobs.get(r))) {   // Closed
                Element subParValue = simDoc.createElement("value");
                subParValue.appendChild(simDoc.createTextNode("null"));
                serviceTimeStrategyNode.appendChild(subParValue);
            } else { // Open
                Station istStation = sn.stations.get(ist);
                JobClass rstJobClass = sn.jobclasses.get(r);
                if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.DISABLED) {
                    Element subParValue = simDoc.createElement("value");
                    subParValue.appendChild(simDoc.createTextNode("null"));
                    serviceTimeStrategyNode.appendChild(subParValue);
                } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.IMMEDIATE) {
                    serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                    serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
                } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.PH
                        || sn.procid.get(istStation).get(rstJobClass) == ProcessType.APH
                        || sn.procid.get(istStation).get(rstJobClass) == ProcessType.COX2
                        || sn.procid.get(istStation).get(rstJobClass) == ProcessType.COXIAN
                        || sn.procid.get(istStation).get(rstJobClass) == ProcessType.HYPEREXP) {
                    Element distributionNode = simDoc.createElement("subParameter");
                    distributionNode.setAttribute("classPath", "jmt.engine.random.PhaseTypeDistr");
                    distributionNode.setAttribute("name", "Phase-Type");
                    Element distrParNode = simDoc.createElement("subParameter");
                    distrParNode.setAttribute("classPath", "jmt.engine.random.PhaseTypePar");
                    distrParNode.setAttribute("name", "distrPar");

                    Element subParNodeAlpha = simDoc.createElement("subParameter");
                    subParNodeAlpha.setAttribute("array", "true");
                    subParNodeAlpha.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlpha.setAttribute("name", "alpha");
                    Element subParNodeAlphaVec = simDoc.createElement("subParameter");
                    subParNodeAlphaVec.setAttribute("array", "true");
                    subParNodeAlphaVec.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlphaVec.setAttribute("name", "vector");
                    MatrixCell PH = sn.proc.get(istStation).get(rstJobClass);
                    Matrix alpha = sn.pie.get(istStation).get(rstJobClass);
                    alpha.absEq();
                    for (int k = 0; k < sn.phases.get(ist, r); k++) {
                        Element subParNodeAlphaElem = simDoc.createElement("subParameter");
                        subParNodeAlphaElem.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlphaElem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha.get(k))));
                        subParNodeAlphaElem.appendChild(subParValue);
                        subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
                    }

                    Element subParNodeT = simDoc.createElement("subParameter");
                    subParNodeT.setAttribute("array", "true");
                    subParNodeT.setAttribute("classPath", "java.lang.Object");
                    subParNodeT.setAttribute("name", "T");
                    Matrix T = PH.get(0);
                    for (int k = 0; k < sn.phases.get(ist, r); k++) {
                        Element subParNodeTvec = simDoc.createElement("subParameter");
                        subParNodeTvec.setAttribute("array", "true");
                        subParNodeTvec.setAttribute("classPath", "java.lang.Object");
                        subParNodeTvec.setAttribute("name", "vector");
                        for (int j = 0; j < sn.phases.get(ist, r); j++) {
                            Element subParNodeTElem = simDoc.createElement("subParameter");
                            subParNodeTElem.setAttribute("classPath", "java.lang.Double");
                            subParNodeTElem.setAttribute("name", "entry");
                            Element subParValue = simDoc.createElement("value");
                            if (k == j) {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * FastMath.abs(T.get(k, j)))));
                            } else {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", FastMath.abs(T.get(k, j)))));
                            }
                            subParNodeTElem.appendChild(subParValue);
                            subParNodeTvec.appendChild(subParNodeTElem);
                        }
                        subParNodeT.appendChild(subParNodeTvec);
                    }
                    subParNodeAlpha.appendChild(subParNodeAlphaVec);
                    distrParNode.appendChild(subParNodeAlpha);
                    distrParNode.appendChild(subParNodeT);
                    serviceTimeStrategyNode.appendChild(distributionNode);
                    serviceTimeStrategyNode.appendChild(distrParNode);
                } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.MAP
                        || sn.procid.get(istStation).get(rstJobClass) == ProcessType.MMPP2) {
                    Element distributionNode = simDoc.createElement("subParameter");
                    distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                    distributionNode.setAttribute("name", "Burst (MAP)");
                    Element distrParNode = simDoc.createElement("subParameter");
                    distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                    distrParNode.setAttribute("name", "distrPar");

                    MatrixCell MAP = sn.proc.get(istStation).get(rstJobClass);

                    Element subParNodeD0 = simDoc.createElement("subParameter");
                    subParNodeD0.setAttribute("array", "true");
                    subParNodeD0.setAttribute("classPath", "java.lang.Object");
                    subParNodeD0.setAttribute("name", "D0");
                    Matrix D0 = MAP.get(0);

                    for (int k = 0; k < sn.phases.get(ist, r); k++) {
                        Element subParNodeD0vec = simDoc.createElement("subParameter");
                        subParNodeD0vec.setAttribute("array", "true");
                        subParNodeD0vec.setAttribute("classPath", "java.lang.Object");
                        subParNodeD0vec.setAttribute("name", "vector");
                        for (int j = 0; j < sn.phases.get(ist, r); j++) {
                            Element subParNodeD0Elem = simDoc.createElement("subParameter");
                            subParNodeD0Elem.setAttribute("classPath", "java.lang.Double");
                            subParNodeD0Elem.setAttribute("name", "entry");
                            Element subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D0.get(k, j))));
                            subParNodeD0Elem.appendChild(subParValue);
                            subParNodeD0vec.appendChild(subParNodeD0Elem);
                        }
                        subParNodeD0.appendChild(subParNodeD0vec);
                    }
                    distrParNode.appendChild(subParNodeD0);

                    Element subParNodeD1 = simDoc.createElement("subParameter");
                    subParNodeD1.setAttribute("array", "true");
                    subParNodeD1.setAttribute("classPath", "java.lang.Object");
                    subParNodeD1.setAttribute("name", "D1");
                    Matrix D1 = MAP.get(1);
                    for (int k = 0; k < sn.phases.get(ist, r); k++) {
                        Element subParNodeD1vec = simDoc.createElement("subParameter");
                        subParNodeD1vec.setAttribute("array", "true");
                        subParNodeD1vec.setAttribute("classPath", "java.lang.Object");
                        subParNodeD1vec.setAttribute("name", "vector");
                        for (int j = 0; j < sn.phases.get(ist, r); j++) {
                            Element subParNodeD1Elem = simDoc.createElement("subParameter");
                            subParNodeD1Elem.setAttribute("classPath", "java.lang.Double");
                            subParNodeD1Elem.setAttribute("name", "entry");
                            Element subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D1.get(k, j))));
                            subParNodeD1Elem.appendChild(subParValue);
                            subParNodeD1vec.appendChild(subParNodeD1Elem);
                        }
                        subParNodeD1.appendChild(subParNodeD1vec);
                    }
                    distrParNode.appendChild(subParNodeD1);
                    serviceTimeStrategyNode.appendChild(distributionNode);
                    serviceTimeStrategyNode.appendChild(distrParNode);
                } else {
                    Element distributionNode = simDoc.createElement("subParameter");
                    String javaClass = "";
                    String javaParClass = "";
                    switch (sn.procid.get(istStation).get(rstJobClass)) {
                        case DET:
                            javaClass = "jmt.engine.random.DeterministicDistr";
                            javaParClass = "jmt.engine.random.DeterministicDistrPar";
                            break;
                        case COX2:
                        case COXIAN:
                            javaClass = "jmt.engine.random.CoxianDistr";
                            javaParClass = "jmt.engine.random.CoxianPar";
                            break;
                        case ERLANG:
                            javaClass = "jmt.engine.random.Erlang";
                            javaParClass = "jmt.engine.random.ErlangPar";
                            break;
                        case EXP:
                            javaClass = "jmt.engine.random.Exponential";
                            javaParClass = "jmt.engine.random.ExponentialPar";
                            break;
                        case GAMMA:
                            javaClass = "jmt.engine.random.GammaDistr";
                            javaParClass = "jmt.engine.random.GammaDistrPar";
                            break;
                        case HYPEREXP:
                            javaClass = "jmt.engine.random.HyperExp";
                            javaParClass = "jmt.engine.random.HyperExpPar";
                            break;
                        case PARETO:
                            javaClass = "jmt.engine.random.Pareto";
                            javaParClass = "jmt.engine.random.ParetoPar";
                            break;
                        case WEIBULL:
                            javaClass = "jmt.engine.random.Weibull";
                            javaParClass = "jmt.engine.random.WeibullPar";
                            break;
                        case LOGNORMAL:
                            javaClass = "jmt.engine.random.Lognormal";
                            javaParClass = "jmt.engine.random.LognormalPar";
                            break;
                        case UNIFORM:
                            javaClass = "jmt.engine.random.Uniform";
                            javaParClass = "jmt.engine.random.UniformPar";
                            break;
                        case MMPP2:
                            javaClass = "jmt.engine.random.MMPP2Distr";
                            javaParClass = "jmt.engine.random.MMPP2Par";
                            break;
                        case REPLAYER:
                        case TRACE:
                            javaClass = "jmt.engine.random.Replayer";
                            javaParClass = "jmt.engine.random.ReplayerPar";
                            break;
                    }
                    distributionNode.setAttribute("classPath", javaClass);
                    switch (sn.procid.get(istStation).get(rstJobClass)) {
                        case REPLAYER:
                        case TRACE:
                            distributionNode.setAttribute("name", "Replayer");
                            break;
                        case EXP:
                            distributionNode.setAttribute("name", "Exponential");
                            break;
                        case HYPEREXP:
                            distributionNode.setAttribute("name", "Hyperexponential");
                            break;
                        default:
                            distributionNode.setAttribute("name", ProcessType.toText(sn.procid.get(istStation).get(rstJobClass)));
                            break;
                    }
                    serviceTimeStrategyNode.appendChild(distributionNode);

                    Element distrParNode = simDoc.createElement("subParameter");
                    distrParNode.setAttribute("classPath", javaParClass);
                    distrParNode.setAttribute("name", "distrPar");

                    Element subParNodeAlpha = simDoc.createElement("subParameter");
                    Element subParValue = simDoc.createElement("value");
                    double c = FastMath.sqrt(sn.scv.get(ist, r));

                    switch (sn.procid.get(istStation).get(rstJobClass)) {
                        case DET:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "t");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / sn.rates.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case EXP:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "lambda");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case HYPEREXP:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "p");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.pie.get(istStation).get(rstJobClass).get(0))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "lambda1");
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).value())));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "lambda2");
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).get(1, 1))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case ERLANG:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "alpha");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(ist, r) * sn.phases.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Long");
                            subParNodeAlpha.setAttribute("name", "r");
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%d", (int) sn.phases.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case GAMMA:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "alpha");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / sn.scv.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "beta");
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.scv.get(ist, r) / sn.rates.get(ist, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case PARETO:
                            double shape = FastMath.sqrt(1 + 1.0 / sn.scv.get(ist, r)) + 1;
                            double scale = 1.0 / sn.rates.get(ist, r) * (shape - 1) / shape;
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "alpha"); //shape
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", shape)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "k"); //scale
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", scale)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case WEIBULL:
                            double rval = FastMath.pow(c, -1.086); //Justus approximation (1976)
                            double alpha = 1.0 / sn.rates.get(ist, r) / Maths.gammaFunction(1 + 1.0 / rval);
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "alpha"); //shape
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "r"); //scale
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", rval)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case LOGNORMAL:
                            double mu = FastMath.log(1.0 / sn.rates.get(ist, r) / FastMath.sqrt(c * c + 1));
                            double sigma = FastMath.sqrt(Math.log(c * c + 1));
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "mu"); //shape
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", mu)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "sigma"); //scale
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sigma)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case UNIFORM:
                            double maxVal = (Math.sqrt(12 * sn.scv.get(ist, r) / FastMath.pow(sn.rates.get(ist, r), 2)) + 2 / sn.rates.get(ist, r)) / 2;
                            double minVal = 2 / sn.rates.get(ist, r) - maxVal;
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "min"); //shape
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", minVal)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            subParNodeAlpha = simDoc.createElement("subParameter");
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "max"); //scale
                            subParValue = simDoc.createElement("value");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", maxVal)));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case REPLAYER:
                        case TRACE:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.String");
                            subParNodeAlpha.setAttribute("name", "fileName");
                            String istFileName = ((ServiceNodeParam) sn.nodeparam.get(simModel.getNodes().get(ind))).fileName.get(r);
                            subParValue.appendChild(simDoc.createTextNode(istFileName));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                    }
                    serviceTimeStrategyNode.appendChild(distrParNode);
                }
            }
            strategyNode.appendChild(serviceTimeStrategyNode);
            section.appendChild(strategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves buffer capacity for JMT XML.
     * LINE uses Kendall notation where cap = K = total system capacity (queue + in-service).
     * JMT's "size" parameter also represents total capacity K.
     */
    public DocumentSectionPair saveBufferCapacity(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "size");
        Element valueNode = simDoc.createElement("value");
        int istStation = (int) sn.nodeToStation.get(ind);
        int cap;
        int numServers;
        if (istStation >= 0) {
            cap = (int) sn.cap.get(istStation);
            numServers = (int) sn.nservers.get(istStation);
        } else {
            cap = Integer.MAX_VALUE;
            numServers = 1;
        }
        // jmtCapIsUnbounded is the same question SolverJMT.jmtBufferCapacityRefusal
        // asks, so the gate and the writer cannot disagree about which stations
        // carry a buffer at all. It used to be Utils.isInf on the int-NARROWED
        // value, which saturated a derived sum of sentinels back onto
        // Integer.MAX_VALUE and so answered correctly by accident.
        if (sn.isstation.get(ind) == 0 || jmtCapIsUnbounded(sn, istStation)) {
            valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
        } else {
            if (cap >= jmtReachablePopulation(istStation)) {
                // A capacity the population cannot reach is unbounded, and the
                // test is >= and not ==: refreshCapacity DERIVES sn.cap for a
                // station the user never capped, as sum over the classes served
                // there of the chain population, so a multi-class station gets
                // (#classes) x N -- 8 on a two-class model of 4 jobs. Under ==
                // only the single-class case matched, and every multi-class one
                // fell through to jmtStationCapAssert and was refused as a
                // "binding" buffer nobody declared.
                //
                // THE POPULATION COMPARED AGAINST IS THE ONE THAT CAN REACH
                // istStation, not the model total. A class that never visits
                // this station cannot fill it, so counting its jobs makes a
                // capacity that is exactly the reachable population look like a
                // buffer -- which is what a SELF-LOOPING CLASS does. See the
                // MATLAB twin in JMTIO/saveBufferCapacity.m.
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else if (Utils.isInf(numServers)) {
                // Infinite servers (delay node) - no buffer constraint
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else {
                jmtStationCapAssert(istStation);
                // Send LINE's total capacity K directly to JMT
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(cap)));
            }
        }
        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * The most jobs that can be present at station {@code ist}, read the way
     * refreshCapacity derives the capacity itself: per CHAIN, because a chain's
     * whole population can reach a station that serves any one of its classes
     * (class switching moves jobs between them), and a chain none of whose
     * classes is served there cannot put a single job on it.
     * <p>
     * Deliberately NOT read off {@code sn.classcap}, which refreshCapacity has
     * already clamped by the station's own cap: comparing a capacity against a
     * quantity derived from it would make every user-declared buffer look
     * non-binding. Infinite when an open chain is served here, which is what the
     * model total gave before and which sends the station to
     * {@link #jmtStationCapAssert}, where the open classes are skipped by name.
     * </p>
     *
     * @param ist station index
     * @return the reachable population, possibly infinite
     */
    private double jmtReachablePopulation(int ist) {
        return jmtReachablePopulation(sn, ist);
    }

    /**
     * {@link #jmtReachablePopulation(int)} against an explicit struct, so that
     * SolverJMT.jmtMethodRefusal can apply the writer's own BINDING test without
     * building a writer. sn.cap is DERIVED for a station the user never capped, so
     * "the cap is finite" is not the question -- "the cap is below what can reach
     * the station" is.
     *
     * @param sn  the model struct
     * @param ist station index
     * @return the reachable population, possibly infinite
     */
    public static double jmtReachablePopulation(NetworkStruct sn, int ist) {
        double n = 0;
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            if (inchain == null) {
                continue;
            }
            boolean served = false;
            double chainJobs = 0;
            for (int k = 0; k < inchain.length(); k++) {
                int r = (int) inchain.get(k);
                if (!Double.isNaN(sn.rates.get(ist, r))) {
                    served = true;
                }
                chainJobs += sn.njobs.get(r);
            }
            if (served) {
                n += chainJobs;
            }
        }
        return n;
    }

    /**
     * Refuses a binding station capacity JMT cannot express, on two counts.
     * <p>
     * (1) THE RULE IS ONE JMT CANNOT READ -- BBS, RSRD or retrial-with-limit; see
     * {@link #jmtReadsDropStrategy}. Such a value is not approximated, it is IGNORED,
     * so the capacity stops being enforced and JMT returns the unconstrained answer.
     * </p><p>
     * (2) THE RULE IS WAITQ AND A CLOSED CLASS CAN REACH THE LIMIT, for the same
     * reason {@code jmtClassCapAssert} in {@code saveRegions} refuses the per-class
     * one: JMT cannot hold a blocked closed job at its upstream station. Note this is
     * the case where NO blocking rule is declared. That the limit CAN be reached is
     * the caller's to establish and not retested here: {@code saveBufferCapacity}
     * reaches this method only for a capacity strictly below the total population,
     * which is the one thing that makes a buffer a buffer. A model that does declare BAS is
     * exported as JMT "BAS blocking", which is the same queueing model, under either
     * declaration form -- see {@link #jmtIsBasDestination}.
     * </p><p>
     * That entry advised expressing the limit as the STATION capacity instead.
     * Measured on 2026-08-19, the advice was wrong, and neither of the two strategies
     * a WAITQ station maps onto reproduces the UNDECLARED case:
     * </p><ul>
     * <li>{@code waiting queue} does not enforce {@code size} at all. On a closed
     * 3-queue tandem, N=6, Exp(1) FCFS, cap 2 at Q2, JMT returned the
     * UNCONSTRAINED [2.03 1.99 1.98], X = 0.750, against the exact
     * [3.6090 0.9711 1.4199], X = 0.6522.</li>
     * <li>{@code BAS blocking} enforces it, but completes the service BEFORE
     * blocking, so the blocked job moves the instant room frees. That is a
     * different queueing model, not a rounding: same fixture, [2.871 1.373 1.756],
     * X = 0.7126.</li>
     * </ul><p>
     * With no rule declared LINE instead disables the upstream departure while the
     * destination is full, which for exponential service is repetitive service (RS)
     * and is what SolverCTMC, SolverSSA and SolverLDES all agree on. So THAT model is
     * refused rather than exported as either of the two things JMT can say. See
     * BUG-81. A declared-BAS model is a different model and is exported, not refused:
     * blocking after service is precisely what JMT's "BAS blocking" does.
     * </p>
     */
    private void jmtStationCapAssert(int ist) {
        String reason = jmtStationCapRefusal(sn, ist);
        if (!reason.isEmpty()) {
            throw new RuntimeException(reason);
        }
    }

    /**
     * Is the buffer at station {@code ist} UNBOUNDED, as the JSIM writer reads it?
     *
     * <p>THIS PORT DOES NOT CARRY Inf ON sn.cap FOR AN UNBOUNDED STATION, and
     * that is what makes the question worth a named predicate. {@code
     * Station.cap} is an {@code int}, so it cannot hold Inf: its "no bound"
     * value is {@code Integer.MAX_VALUE} ({@link
     * jline.lang.nodes.Station#hasFiniteCap} is the predicate for the declared
     * one). refreshCapacity then derives a station's capacity as
     * {@code min(sum(chaincap row), sum(classcap row))}, and it SUMS that
     * sentinel across the classes served there -- so an unbounded MIXED station
     * comes out as 2147483647 + 2 = 2147483649 for one open and one closed
     * class, and a Source serving two open classes as 2 x 2147483647 =
     * 4294967294. {@link Utils#isInf} recognises the sentinel EXACTLY and so
     * does not see a sum of them.
     *
     * <p>{@code >= Integer.MAX_VALUE} is refreshCapacity's own idiom for the
     * same question ({@code chainCap >= Integer.MAX_VALUE},
     * {@code station.getCap() >= Integer.MAX_VALUE}), and it is safe because
     * setCapacity takes an {@code int}: no capacity a caller can declare is
     * larger. MATLAB, native python and C++ all carry a genuine Inf here, so
     * this narrowing is the JAR's alone.
     *
     * <p>saveBufferCapacity used to ask the question by NARROWING to int first,
     * which saturates the sum back onto Integer.MAX_VALUE and made Utils.isInf
     * answer correctly by accident. Naming it here is what lets the gate ask the
     * same question as the writer instead of reading the raw double and
     * concluding that every mixed model has a binding buffer.
     *
     * @param sn  the model struct
     * @param ist station index
     * @return true when the station has no bound a job can reach
     */
    public static boolean jmtCapIsUnbounded(NetworkStruct sn, int ist) {
        if (ist < 0 || sn.cap == null || ist >= sn.cap.length()) {
            return true;
        }
        double cap = sn.cap.get(ist);
        return Utils.isInf(cap) || cap >= Integer.MAX_VALUE;
    }

    /**
     * {@link #jmtStationCapAssert(int)} as a SENTENCE rather than an exception,
     * against an explicit struct; empty when the buffer at {@code ist} is
     * exportable.
     *
     * <p>ONE PREDICATE, TWO CALLERS. saveBufferCapacity raises it while writing
     * the JSIM document, and SolverJMT.jmtMethodRefusal returns it so that
     * findSolver and SolverAUTO never offer jmt.jsim on a model the writer will
     * refuse. It used to be reachable only from inside the writer, which is why
     * the gate could not see it.
     *
     * @param sn  the model struct
     * @param ist station index
     * @return empty string when exportable, otherwise the refusal
     */
    public static String jmtStationCapRefusal(NetworkStruct sn, int ist) {
        if (jmtCapIsUnbounded(sn, ist)) {
            return "";
        }
        for (int r = 0; r < sn.nclasses; r++) {
            if (Double.isNaN(sn.rates.get(ist, r))) {
                continue;   // class r is not served here
            }
            Map<JobClass, DropStrategy> byClass = sn.droprule.get(sn.stations.get(ist));
            DropStrategy dr = byClass == null ? null : byClass.get(sn.jobclasses.get(r));
            if (dr != null && !jmtReadsDropStrategy(dr)) {
                // Unmappable for EITHER class type, so this test precedes the open-class skip
                return "SolverJMT: station '"
                        + sn.nodenames.get((int) sn.stationToNode.get(ist))
                        + "' applies drop strategy \"" + DropStrategy.toText(dr)
                        + "\" to class '" + sn.classnames.get(r) + "' and carries a finite capacity "
                        + (long) sn.cap.get(ist) + " it can reach. JMT's queue section reads only"
                        + " \"drop\", \"waiting queue\", \"BAS blocking\" and \"retrial\"; it does not"
                        + " approximate anything else, it ignores it, so the capacity would stop being"
                        + " enforced and the run would return the unconstrained answer."
                        + " Use SolverCTMC, SolverSSA or SolverLDES.";
            }
            if (Utils.isInf(sn.njobs.get(r))) {
                continue;   // open class: JMT loses its arrivals, as LINE does
            }
            if (dr != DropStrategy.WaitingQueue) {
                continue;   // a mappable declared blocking rule is exported as itself
            }
            if (jmtIsBasDestination(sn, ist, r)) {
                // BAS declared on the UPSTREAM station: jmtDropStrategyText moves it
                // onto this one, which is where JMT reads it
                continue;
            }
            return "SolverJMT: station '" + sn.nodenames.get((int) sn.stationToNode.get(ist))
                    + "' carries a finite capacity " + (long) sn.cap.get(ist)
                    + " that binds for the closed class '" + sn.classnames.get(r)
                    + "'. LINE blocks a closed job that finds no room -- the upstream departure is disabled"
                    + " and the job stays where it is -- and no JMT drop strategy reproduces that:"
                    + " \"waiting queue\" does not enforce the size at all, and \"BAS blocking\" completes"
                    + " the service before blocking, which is a different queueing model."
                    + " Use SolverCTMC, SolverSSA or SolverLDES, or declare DropStrategy.BAS"
                    + " if blocking after service is the model you want, which SolverJMT"
                    + " does export.";
        }
        return "";
    }

    public ElementDocumentPair saveCache(ElementDocumentPair elementDocumentPair) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;


        //add measure for cache hit/ miss rate
        int numOfNodes = sn.nnodes;
        Cache cacheNode = null;
        for (int r = 0; r < numOfNodes; r++) {
            Node node = simModel.getNodes().get(r);
            if (node instanceof Cache) {
                cacheNode = (Cache) simModel.getNodes().get(r);
                int numOfCacheClasses = cacheNode.getMissClass().getNumCols();
                for (int i = 0; i < numOfCacheClasses; i++) {
                    if (cacheNode.getHitClass().get(i) >= 0) {
                        // Add Throughput measure for InitClass
                        Element cacheOriginNode = simDoc.createElement("measure");
                        cacheOriginNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                        cacheOriginNode.setAttribute("name", "cacheHit");
                        cacheOriginNode.setAttribute("nodeType", "section");
                        cacheOriginNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                        cacheOriginNode.setAttribute("referenceNode", cacheNode.getName());
                        cacheOriginNode.setAttribute("referenceUserClass", sn.classnames.get(i));
                        cacheOriginNode.setAttribute("type", "Throughput");
                        cacheOriginNode.setAttribute("verbose", "false");
                        simElem.appendChild(cacheOriginNode);

                        // Add Arrival Rate measure for InitClass
                        Element cacheArrivalNode = simDoc.createElement("measure");
                        cacheArrivalNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                        cacheArrivalNode.setAttribute("name", "cacheArrival");
                        cacheArrivalNode.setAttribute("nodeType", "station");
                        cacheArrivalNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                        cacheArrivalNode.setAttribute("referenceNode", cacheNode.getName());
                        cacheArrivalNode.setAttribute("referenceUserClass", sn.classnames.get(i));
                        cacheArrivalNode.setAttribute("type", "Arrival Rate");
                        cacheArrivalNode.setAttribute("verbose", "false");
                        simElem.appendChild(cacheArrivalNode);

                        // Add Throughput measure for HitClass
                        Element cacheHitNode = simDoc.createElement("measure");
                        cacheHitNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                        cacheHitNode.setAttribute("name", "cacheHit");
                        cacheHitNode.setAttribute("nodeType", "section");
                        cacheHitNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                        cacheHitNode.setAttribute("referenceNode", cacheNode.getName());
                        cacheHitNode.setAttribute("referenceUserClass", sn.classnames.get((int) (cacheNode.getHitClass().get(i))));
                        cacheHitNode.setAttribute("type", "Throughput");
                        cacheHitNode.setAttribute("verbose", "false");
                        simElem.appendChild(cacheHitNode);
                    }

                    if (cacheNode.getMissClass().get(i) >= 0) {
                        Element cacheMissNode = simDoc.createElement("measure");
                        cacheMissNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                        cacheMissNode.setAttribute("name", "cacheHit");
                        cacheMissNode.setAttribute("nodeType", "section");
                        cacheMissNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                        cacheMissNode.setAttribute("referenceNode", cacheNode.getName());
                        //cacheMissNode.setAttribute("referenceUserClass", sn.classnames.get(i));
                        //cacheMissNode.setAttribute("referenceUserClass", sn.classnames.get((int) (cacheNode.getHitClass().get(i))));
                        cacheMissNode.setAttribute("referenceUserClass", sn.classnames.get((int) (cacheNode.getMissClass().get(i))));
                        cacheMissNode.setAttribute("type", "Throughput");
                        cacheMissNode.setAttribute("verbose", "false");
                        simElem.appendChild(cacheMissNode);
                    }


                }


            }
        }

        return new ElementDocumentPair(simElem, simDoc);
    }

    public DocumentSectionPair saveCacheStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;


        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        Cache cacheNode = (Cache) simModel.getNodes().get(ind);
        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractRows(sn.connmatrix, ind, ind + 1, conn_i);
        Matrix conn_i_find = conn_i.find();
        int j = (int) conn_i_find.get(0);


        Element paramNode = simDoc.createElement("parameter");
        paramNode.setAttribute("classPath", "java.lang.Integer");
        paramNode.setAttribute("name", "maxItems");
        Element valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(String.format("%d", cacheNode.getItems().getNumberOfItems())));
        paramNode.appendChild(valueNode);
        section.appendChild(paramNode);


        Element cacheCapacityParamNode = simDoc.createElement("parameter");
        cacheCapacityParamNode.setAttribute("array", "true");
        cacheCapacityParamNode.setAttribute("classPath", "java.lang.Integer");
        cacheCapacityParamNode.setAttribute("name", "cacheCapacity");

        for (int r = 0; r < cacheNode.getItemLevelCap().getNumRows(); r++) {
            for (int c = 0; c < cacheNode.getItemLevelCap().getNumCols(); c++) {
                Element subcacheCapacityParamNode = simDoc.createElement("subParameter");
                subcacheCapacityParamNode.setAttribute("classPath", "java.lang.Integer");
                subcacheCapacityParamNode.setAttribute("name", "cache");
                Element cacheCapacityvalueNode = simDoc.createElement("value");
                cacheCapacityvalueNode.appendChild(simDoc.createTextNode(String.format("%d", (int) cacheNode.getItemLevelCap().get(r, c))));
                subcacheCapacityParamNode.appendChild(cacheCapacityvalueNode);
                cacheCapacityParamNode.appendChild(subcacheCapacityParamNode);
            }
        }

        section.appendChild(cacheCapacityParamNode);

        Element accessProbNodeMatrix = simDoc.createElement("parameter");
        accessProbNodeMatrix.setAttribute("array", "true");
        accessProbNodeMatrix.setAttribute("classPath", "java.lang.Object");
        accessProbNodeMatrix.setAttribute("name", "matrix");
        for (int r = 0; r < cacheNode.getItemLevelCap().getNumCols(); r++) {
            Element accessProbNodeRows = simDoc.createElement("subParameter");
            accessProbNodeRows.setAttribute("array", "true");
            accessProbNodeRows.setAttribute("classPath", "java.lang.Float");
            accessProbNodeRows.setAttribute("name", "row");
            for (int c = 0; c < cacheNode.getItemLevelCap().getNumCols(); c++) {
                Element accessProbNodeCell = simDoc.createElement("subParameter");
                accessProbNodeCell.setAttribute("classPath", "java.lang.Float");
                accessProbNodeCell.setAttribute("name", "cell");
                Element valNode = simDoc.createElement("value");
                double prob = 0.00;
                if ((r == c - 1) || ((r == c) && (c == cacheNode.getItemLevelCap().getNumCols() - 1))) {
                    prob = 1.00;
                }
                valNode.appendChild(simDoc.createTextNode(String.format("%12.12f", prob)));
                accessProbNodeCell.appendChild(valNode);
                accessProbNodeRows.appendChild(accessProbNodeCell);
            }
            accessProbNodeMatrix.appendChild(accessProbNodeRows);
        }
        section.appendChild(accessProbNodeMatrix);


        // save hit class and miss class
        ArrayList<JobClass> initialClass = new ArrayList<JobClass>();
        ArrayList<JobClass> missClass = new ArrayList<JobClass>();
        ArrayList<JobClass> hitClass = new ArrayList<JobClass>();
        for (int r = 0; r < cacheNode.getMissClass().getNumElements(); r++) {
            JobClass tempjobClass = sn.jobclasses.get(r);
            if (cacheNode.getMissClass().get(r) >= 0 && cacheNode.getHitClass().get(r) >= 0) {
                initialClass.add(tempjobClass);
                missClass.add(sn.jobclasses.get((int) cacheNode.getMissClass().get(r)));
                hitClass.add(sn.jobclasses.get((int) cacheNode.getHitClass().get(r)));
            }
        }


        Element cacheInitNode = simDoc.createElement("parameter");
        cacheInitNode.setAttribute("array", "true");
        cacheInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
        cacheInitNode.setAttribute("name", "jobClasses");


        for (int r = 0; r < initialClass.size(); r++) {
            Element subCacheInitNode = simDoc.createElement("subParameter");
            subCacheInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
            subCacheInitNode.setAttribute("name", "initClass");
            Element subsubCacheInitNode = simDoc.createElement("subParameter");
            subsubCacheInitNode.setAttribute("classPath", "java.lang.String");
            subsubCacheInitNode.setAttribute("name", "name");
            Element valNode = simDoc.createElement("value");
            valNode.appendChild(simDoc.createTextNode(initialClass.get(r).getName()));
            subsubCacheInitNode.appendChild(valNode);
            subCacheInitNode.appendChild(subsubCacheInitNode);
            cacheInitNode.appendChild(subCacheInitNode);
        }
        section.appendChild(cacheInitNode);


        Element hitClassInitNode = simDoc.createElement("parameter");
        hitClassInitNode.setAttribute("array", "true");
        hitClassInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
        hitClassInitNode.setAttribute("name", "hitClasses");


        for (int r = 0; r < hitClass.size(); r++) {
            Element subHitInitNode = simDoc.createElement("subParameter");
            subHitInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
            subHitInitNode.setAttribute("name", "hitClass");
            Element subsubHitInitNode = simDoc.createElement("subParameter");
            subsubHitInitNode.setAttribute("classPath", "java.lang.String");
            subsubHitInitNode.setAttribute("name", "name");
            Element valNode = simDoc.createElement("value");
            valNode.appendChild(simDoc.createTextNode(hitClass.get(r).getName()));
            subsubHitInitNode.appendChild(valNode);
            subHitInitNode.appendChild(subsubHitInitNode);
            hitClassInitNode.appendChild(subHitInitNode);
        }
        section.appendChild(hitClassInitNode);

        Element missClassInitNode = simDoc.createElement("parameter");
        missClassInitNode.setAttribute("array", "true");
        missClassInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
        missClassInitNode.setAttribute("name", "missClasses");


        for (int r = 0; r < missClass.size(); r++) {
            Element subMissInitNode = simDoc.createElement("subParameter");
            subMissInitNode.setAttribute("classPath", "jmt.engine.QueueNet.JobClass");
            subMissInitNode.setAttribute("name", "missClass");
            Element subsubMissInitNode = simDoc.createElement("subParameter");
            subsubMissInitNode.setAttribute("classPath", "java.lang.String");
            subsubMissInitNode.setAttribute("name", "name");
            Element valNode = simDoc.createElement("value");
            valNode.appendChild(simDoc.createTextNode(missClass.get(r).getName()));
            subsubMissInitNode.appendChild(valNode);
            subMissInitNode.appendChild(subsubMissInitNode);
            missClassInitNode.appendChild(subMissInitNode);
        }
        section.appendChild(missClassInitNode);


        Element cacheStrategyParamNode = simDoc.createElement("parameter");
        String strategyPath = cacheStrategyMap(cacheNode.getReplacementStrategy());
        cacheStrategyParamNode.setAttribute("classPath", strategyPath);
        cacheStrategyParamNode.setAttribute("name", "replacePolicy");
        section.appendChild(cacheStrategyParamNode);

        Element cachePopularityParamNode = simDoc.createElement("parameter");
        cachePopularityParamNode.setAttribute("array", "true");
        cachePopularityParamNode.setAttribute("classPath", "jmt.engine.random.discrete.DiscreteDistribution");
        cachePopularityParamNode.setAttribute("name", "popularity");


        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            cachePopularityParamNode.appendChild(refClassNode);

            Element subParNodeRow = null;


            if (cacheNode.popularityGet(0, r) == null) {
                subParNodeRow = simDoc.createElement("subParameter");
                subParNodeRow.setAttribute("array", "true");
                subParNodeRow.setAttribute("classPath", "jmt.engine.random.discrete.DiscreteDistribution");
                subParNodeRow.setAttribute("name", "popularity");

                Element valueNodeCell = simDoc.createElement("value");
                valueNodeCell.appendChild(simDoc.createTextNode("null"));
                subParNodeRow.appendChild(valueNodeCell);
            } else if (cacheNode.popularityGet(0, r).getName().equals("Zipf")) {

                subParNodeRow = simDoc.createElement("subParameter");
                // subParNodeRow.setAttribute("array", "true");
                subParNodeRow.setAttribute("classPath", "jmt.engine.random.discrete.Zipf");
                subParNodeRow.setAttribute("name", "popularity");

                Element subParNodeCell1 = simDoc.createElement("subParameter");
                subParNodeCell1.setAttribute("classPath", "java.lang.Double");
                subParNodeCell1.setAttribute("name", "alpha");
                Element valNode1 = simDoc.createElement("value");
                valNode1.appendChild(simDoc.createTextNode(String.format("%12.12f", cacheNode.popularityGet(0, r).getParam(3).getValue())));


                subParNodeCell1.appendChild(valNode1);
                subParNodeRow.appendChild(subParNodeCell1);

                Element subParNodeCell2 = simDoc.createElement("subParameter");
                subParNodeCell2.setAttribute("classPath", "java.lang.Integer");
                subParNodeCell2.setAttribute("name", "numberOfElements");
                Element valNode2 = simDoc.createElement("value");
                valNode2.appendChild(simDoc.createTextNode(
                        String.format("%d", cacheNode.popularityGet(0, r).getParam(2).getValue())));
                subParNodeCell2.appendChild(valNode2);
                subParNodeRow.appendChild(subParNodeCell2);
            } else if (cacheNode.popularityGet(0, r).getName().equals("DiscreteSampler")) {
                subParNodeRow = simDoc.createElement("subParameter");
                // subParNodeRow.setAttribute("array", "true");
                subParNodeRow.setAttribute("classPath", "jmt.engine.random.discrete.Uniform");
                subParNodeRow.setAttribute("name", "popularity");

                Element subParNodeCell1 = simDoc.createElement("subParameter");
                subParNodeCell1.setAttribute("classPath", "java.lang.Integer");
                subParNodeCell1.setAttribute("name", "min");
                Element valNode1 = simDoc.createElement("value");
                valNode1.appendChild(simDoc.createTextNode(String.format("%d", cacheNode.popularityGet(0, r).getParam(4).getValue())));

                subParNodeCell1.appendChild(valNode1);
                subParNodeRow.appendChild(subParNodeCell1);

                Element subParNodeCell2 = simDoc.createElement("subParameter");
                subParNodeCell2.setAttribute("classPath", "java.lang.Integer");
                subParNodeCell2.setAttribute("name", "max");
                Element valNode2 = simDoc.createElement("value");
                valNode2.appendChild(simDoc.createTextNode(String.format("%d", cacheNode.popularityGet(0, r).getParam(5).getValue())));
                subParNodeCell2.appendChild(valNode2);
                subParNodeRow.appendChild(subParNodeCell2);
            } else {
                subParNodeRow = simDoc.createElement("subParameter");
                subParNodeRow.setAttribute("array", "true");
                subParNodeRow.setAttribute("classPath", "jmt.engine.random.discrete.DiscreteDistribution");
                subParNodeRow.setAttribute("name", "popularity");

                Element valueNodeCell = simDoc.createElement("value");
                valueNodeCell.appendChild(simDoc.createTextNode("null"));
                subParNodeRow.appendChild(valueNodeCell);

            }


            if (subParNodeRow != null) {
                cachePopularityParamNode.appendChild(subParNodeRow);
            }

        }

        section.appendChild(cachePopularityParamNode);

        Element paramNodeMatrix = simDoc.createElement("parameter");
        paramNodeMatrix.setAttribute("array", "true");
        paramNodeMatrix.setAttribute("classPath", "java.lang.Object");
        paramNodeMatrix.setAttribute("name", "matrix");

        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            paramNodeMatrix.appendChild(refClassNode);

            Element subParNodeRow = simDoc.createElement("subParameter");
            subParNodeRow.setAttribute("array", "true");
            subParNodeRow.setAttribute("classPath", "java.lang.Float");
            subParNodeRow.setAttribute("name", "row");
            for (int s = 0; s < numOfClasses; s++) {
                refClassNode = simDoc.createElement("refClass");
                refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(s)));
                subParNodeRow.appendChild(refClassNode);

                Element subParNodeCell = simDoc.createElement("subParameter");
                subParNodeCell.setAttribute("classPath", "java.lang.Float");
                subParNodeCell.setAttribute("name", "cell");
                Element valNode = simDoc.createElement("value");
                valNode.appendChild(simDoc.createTextNode(String.format("%12.12f", sn.rtnodes.get(ind * numOfClasses + r, j * numOfClasses + s))));
                subParNodeCell.appendChild(valNode);
                subParNodeRow.appendChild(subParNodeCell);
            }
            paramNodeMatrix.appendChild(subParNodeRow);
        }

        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveClassSwitchStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVECLASSSWITCHSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element paramNode = simDoc.createElement("parameter");
        paramNode.setAttribute("array", "true");
        paramNode.setAttribute("classPath", "java.lang.Object");
        paramNode.setAttribute("name", "matrix");

        int K = sn.nclasses;
        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractRows(sn.connmatrix, ind, ind + 1, conn_i);
        Matrix jset = conn_i.find();
        boolean[] exportClasses = getExportableClasses();
        for (int r = 0; r < K; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            paramNode.appendChild(refClassNode);

            Element subParNodeRow = simDoc.createElement("subParameter");
            subParNodeRow.setAttribute("array", "true");
            subParNodeRow.setAttribute("classPath", "java.lang.Float");
            subParNodeRow.setAttribute("name", "row");
            for (int s = 0; s < K; s++) {
                if (!exportClasses[s]) continue;
                refClassNode = simDoc.createElement("refClass");
                refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(s)));
                subParNodeRow.appendChild(refClassNode);

                Element subParNodeCell = simDoc.createElement("subParameter");
                subParNodeCell.setAttribute("classPath", "java.lang.Float");
                subParNodeCell.setAttribute("name", "cell");
                Element valNode = simDoc.createElement("value");
                double val = 0.0;
                for (int ij = 0; ij < jset.length(); ij++) {
                    val += sn.rtnodes.get(ind * K + r, (int) (jset.get(ij) * K + s));
                }
                valNode.appendChild(simDoc.createTextNode(String.format("%12.12f", val)));
                subParNodeCell.appendChild(valNode);
                subParNodeRow.appendChild(subParNodeCell);
            }
            paramNode.appendChild(subParNodeRow);
        }
        section.appendChild(paramNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public ElementDocumentPair saveClasses(ElementDocumentPair elementDocumentPair) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;

        // [SIMELEM, SIMDOC] = SAVECLASSES(SIMELEM, SIMDOC)
        int numOfClasses = sn.nclasses;
        int numOfNodes = sn.nnodes;

        Cache cacheNode = null;

        for (int i = 0; i < numOfNodes; i++) {
            Node node = simModel.getNodes().get(i);

            if (node instanceof Cache) {
                cacheNode = (Cache) node;
            }
        }

        // JMT uses higher priority value = higher priority, while LINE uses lower value = higher priority
        // We need to invert priorities when writing to JMT
        boolean[] exportClasses = getExportableClasses();
        int maxPrio = 0;
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            int prio = (int) sn.classprio.get(r);
            if (prio > maxPrio) {
                maxPrio = prio;
            }
        }

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element userClass = simDoc.createElement("userClass");
            userClass.setAttribute("name", sn.classnames.get(r));
            if (sn.jobclasses.get(r) instanceof OpenClass) {
                userClass.setAttribute("type", "open");
            } else {
                userClass.setAttribute("type", "closed");
            }
            // Set soft deadline from class deadline (0.0 if no deadline)
            double classDeadline = sn.jobclasses.get(r).getDeadline();
            if (Double.isFinite(classDeadline)) {
                userClass.setAttribute("softDeadline", String.valueOf(classDeadline));
            } else {
                userClass.setAttribute("softDeadline", "0.0");
            }
            // Invert priority: LINE uses lower=higher, JMT uses higher=higher
            int jmtPrio = maxPrio - (int) sn.classprio.get(r);
            userClass.setAttribute("priority", String.valueOf(jmtPrio));

            int refStatIndex = (int) sn.refstat.get(r);
            MatrixCell integerMatrixMap = sn.proc.get(sn.stations.get(refStatIndex)).get(sn.jobclasses.get(r));

            if (!integerMatrixMap.isEmpty()) {
                if (cacheNode != null) {
                    userClass.setAttribute("customers", String.valueOf((int) sn.njobs.get(r)));
                    userClass.setAttribute("referenceSource", (simModel.getClasses().get(r)).getReferenceStation().getName());
                } else if (!Utils.isInf(sn.njobs.get(r)) && !Double.isNaN(sn.njobs.get(r))) {
                    userClass.setAttribute("customers", String.valueOf((int) sn.njobs.get(r)));
                    userClass.setAttribute("referenceSource", sn.nodenames.get((int) sn.stationToNode.get(refStatIndex)));
                } else if (integerMatrixMap.get(0).hasNaN()) { // open disabled in source
                    userClass.setAttribute("referenceSource", "ClassSwitch");
                } else { // if other open
                    userClass.setAttribute("referenceSource", sn.nodenames.get((int) sn.stationToNode.get((int) sn.refstat.get(r))));
                }
            } else {
                userClass.setAttribute("referenceSource", sn.nodenames.get((int) sn.stationToNode.get((int) sn.refstat.get(r))));
            }

            simElem.appendChild(userClass);
        }
        return new ElementDocumentPair(simElem, simDoc);
    }

    /**
     * True when station {@code ist} is the RECEIVING side of a true-BAS relation for
     * class {@code r}, i.e. an arrival of {@code r} that finds {@code ist} full must
     * block an upstream station rather than be lost.
     * <p>
     * LINE accepts the BAS declaration in two places -- on the blocking (upstream)
     * station, as {@code cqn_bas_blocking} does, or on the full destination, as a model
     * read back from JMT does -- and {@code Network.refreshLocalVars} resolves both into
     * {@code sn.isbasdestination} (BUG-83). Reading {@code sn.droprule} at the capped
     * station sees only the second form, which is what made SolverJMT refuse the first.
     * </p>
     */
    private boolean jmtIsBasDestination(int ist, int r) {
        return jmtIsBasDestination(sn, ist, r);
    }

    /** {@link #jmtIsBasDestination(int, int)} against an explicit struct. */
    private static boolean jmtIsBasDestination(NetworkStruct sn, int ist, int r) {
        if (ist < 0 || sn.isbasdestination == null) {
            return false;
        }
        if (ist >= sn.isbasdestination.getNumRows() || r >= sn.isbasdestination.getNumCols()) {
            return false;
        }
        return sn.isbasdestination.get(ist, r) == 1.0;
    }

    /**
     * Whether JMT's queue section can read this drop strategy at all.
     * <p>
     * It recognizes exactly four {@code dropStrategies} strings -- 'drop', 'BAS
     * blocking', 'waiting queue', 'retrial' (a lookupswitch on String.hashCode in
     * {@code jmt/engine/NodeSections/Queue.class}; the Storage section of a Place is
     * even narrower and drops 'retrial'). An unrecognized value falls through the
     * default arm with NO flag set, so BBS, RSRD and retrial-with-limit are not
     * approximated, they are IGNORED.
     * </p>
     */
    private static boolean jmtReadsDropStrategy(DropStrategy dr) {
        return dr == DropStrategy.Drop || dr == DropStrategy.BlockingAfterService
                || dr == DropStrategy.WaitingQueue || dr == DropStrategy.Retrial;
    }

    /**
     * The JMT dropStrategy/dropRule string for station {@code ist}, class {@code r}.
     * <p>
     * Beyond {@code DropStrategy.toText} this does two things.
     * </p><p>
     * It resolves the two ways LINE can declare BAS blocking onto the one way JMT can
     * read it. JMT's queue section says what happens to an arrival that finds THIS buffer
     * full, so it only understands the rule on the destination; a WAITQ slot that
     * {@link #jmtIsBasDestination} marks is therefore written out as 'BAS blocking'.
     * </p><p>
     * And it keeps the written file VALID. JMT recognizes exactly four strings --
     * 'drop', 'BAS blocking', 'waiting queue', 'retrial' -- so BBS, RSRD and
     * retrial-with-limit are spelled 'waiting queue', JMT's own no-limit default. That
     * substitution is only ever reached where the rule cannot be consulted (infinite
     * size, or a closed capacity equal to the population): a buffer that can actually
     * fill under one of those three is refused outright by {@link #jmtStationCapAssert}.
     * </p>
     */
    private String jmtDropStrategyText(int ist, int r) {
        if (ist < 0) {
            return "drop";     // JMT sets the field to 'drop' for nodes without a buffer
        }
        Map<JobClass, DropStrategy> byClass = sn.droprule.get(sn.stations.get(ist));
        DropStrategy dr = byClass == null ? null : byClass.get(sn.jobclasses.get(r));
        if (dr == null || dr.getID() == 0) {
            return "drop";     // the "no rule recorded" slot, which JMT also spells 'drop'
        }
        if (dr == DropStrategy.WaitingQueue && jmtIsBasDestination(ist, r)) {
            return DropStrategy.toText(DropStrategy.BlockingAfterService);
        }
        if (!jmtReadsDropStrategy(dr)) {
            return DropStrategy.toText(DropStrategy.WaitingQueue);
        }
        return DropStrategy.toText(dr);
    }

    public DocumentSectionPair saveDropRule(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element schedStrategyNode = simDoc.createElement("parameter");
        schedStrategyNode.setAttribute("array", "true");
        schedStrategyNode.setAttribute("classPath", "java.lang.String");
        schedStrategyNode.setAttribute("name", "dropRules");

        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            schedStrategyNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.String");
            subParameterNode.setAttribute("name", "dropRule");

            Element valueNode2 = simDoc.createElement("value");

            valueNode2.appendChild(simDoc.createTextNode(jmtDropStrategyText(i, r)));
            subParameterNode.appendChild(valueNode2);
            schedStrategyNode.appendChild(subParameterNode);
            section.appendChild(schedStrategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveDropStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();

        Element schedStrategyNode = simDoc.createElement("parameter");
        schedStrategyNode.setAttribute("array", "true");
        schedStrategyNode.setAttribute("classPath", "java.lang.String");
        schedStrategyNode.setAttribute("name", "dropStrategies");
        double i = sn.nodeToStation.get(ind);
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            schedStrategyNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.String");
            subParameterNode.setAttribute("name", "dropStrategy");

            Element valueNode2 = simDoc.createElement("value");
            valueNode2.appendChild(simDoc.createTextNode(
                    Double.isNaN(i) ? "drop" : jmtDropStrategyText((int) i, r)));
            subParameterNode.appendChild(valueNode2);
            schedStrategyNode.appendChild(subParameterNode);
            section.appendChild(schedStrategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveEnablingConditions(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element enablingNode = simDoc.createElement("parameter");
        enablingNode.setAttribute("array", "true");
        enablingNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
        enablingNode.setAttribute("name", "enablingConditions");

        int numOfNodes = sn.nnodes;
        int numOfClasses = sn.nclasses;
        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subEnablingConditionNode = simDoc.createElement("subParameter");
            subEnablingConditionNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            subEnablingConditionNode.setAttribute("name", "enablingCondition");

            Element subEnablingVectorsNode = simDoc.createElement("subParameter");
            subEnablingVectorsNode.setAttribute("array", "true");
            subEnablingVectorsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
            subEnablingVectorsNode.setAttribute("name", "enablingVectors");

            for (int k = 0; k < numOfNodes; k++) {
                boolean hasRelevantEntry = false;

                for (int r = 0; r < numOfClasses; r++) {
                    double en_val = ((TransitionNodeParam) sn.nodeparam.get(istNode)).enabling.get(m).get(k, r);
                    double in_val = ((TransitionNodeParam) sn.nodeparam.get(istNode)).inhibiting.get(m).get(k, r);
                    if ((!Double.isInfinite(en_val) && en_val > 0) || (!Double.isInfinite(in_val) && in_val > 0)) {
                        hasRelevantEntry = true;
                        break;
                    }
                }

                if (!hasRelevantEntry) continue; // Skip this node if no relevant enabling entries

                Element subEnablingVectorNode = simDoc.createElement("subParameter");
                subEnablingVectorNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                subEnablingVectorNode.setAttribute("name", "enablingVector");

                Element subStationNameNode = simDoc.createElement("subParameter");
                subStationNameNode.setAttribute("classPath", "java.lang.String");
                subStationNameNode.setAttribute("name", "stationName");

                Element placeNameValueNode = simDoc.createElement("value");
                placeNameValueNode.appendChild(simDoc.createTextNode(sn.nodenames.get(k)));
                subStationNameNode.appendChild(placeNameValueNode);

                subEnablingVectorNode.appendChild(subStationNameNode);

                Element subEnablingEntriesNode = simDoc.createElement("subParameter");
                subEnablingEntriesNode.setAttribute("array", "true");
                subEnablingEntriesNode.setAttribute("classPath", "java.lang.Integer");
                subEnablingEntriesNode.setAttribute("name", "enablingEntries");

                for (int r = 0; r < numOfClasses; r++) {
                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    subEnablingEntriesNode.appendChild(refClassNode);

                    Element subParameterNode = simDoc.createElement("subParameter");
                    subParameterNode.setAttribute("classPath", "java.lang.Integer");
                    subParameterNode.setAttribute("name", "enablingEntry");

                    Element valueNode = simDoc.createElement("value");
                    double val = ((TransitionNodeParam) sn.nodeparam.get(istNode)).enabling.get(m).get(k, r);
                    if (Double.isInfinite(val)) {
                        valueNode.appendChild(simDoc.createTextNode("-1"));
                    } else {
                        valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) val)));
                    }

                    subParameterNode.appendChild(valueNode);
                    subEnablingEntriesNode.appendChild(subParameterNode);
                }

                subEnablingVectorNode.appendChild(subEnablingEntriesNode);
                subEnablingVectorsNode.appendChild(subEnablingVectorNode);
            }

            subEnablingConditionNode.appendChild(subEnablingVectorsNode);
            enablingNode.appendChild(subEnablingConditionNode);
        }

        section.appendChild(enablingNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveFiringOutcomes(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element firingOutcomesNode = simDoc.createElement("parameter");
        firingOutcomesNode.setAttribute("array", "true");
        firingOutcomesNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
        firingOutcomesNode.setAttribute("name", "firingOutcomes");


        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractRows(sn.connmatrix, ind, ind + 1, conn_i);
        Matrix outputs = conn_i.find();
        List<String> connections = new ArrayList<String>();
        for (int idx = 0; idx < outputs.length(); idx++) {
            connections.add(sn.nodenames.get((int) outputs.get(idx)));
        }
        int numOfOutput = connections.size();
        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;
        int numOfClasses = simModel.getClasses().size();

        for (int m = 0; m < numOfModes; m++) {
            Element subFiringOutcomeNode = simDoc.createElement("subParameter");
            subFiringOutcomeNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            subFiringOutcomeNode.setAttribute("name", "firingOutcome");

            Element subFiringVectorsNode = simDoc.createElement("subParameter");
            subFiringVectorsNode.setAttribute("array", "true");
            subFiringVectorsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
            subFiringVectorsNode.setAttribute("name", "firingVectors");

            for (int k = 0; k < numOfOutput; k++) {
                Element subFiringVectorNode = simDoc.createElement("subParameter");
                subFiringVectorNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                subFiringVectorNode.setAttribute("name", "firingVector");

                Element subStationNameNode = simDoc.createElement("subParameter");
                subStationNameNode.setAttribute("classPath", "java.lang.String");
                subStationNameNode.setAttribute("name", "stationName");

                Element placeNameValueNode = simDoc.createElement("value");
                placeNameValueNode.appendChild(simDoc.createTextNode(connections.get(k)));
                subStationNameNode.appendChild(placeNameValueNode);

                subFiringVectorNode.appendChild(subStationNameNode);

                Element subFiringEntriesNode = simDoc.createElement("subParameter");
                subFiringEntriesNode.setAttribute("array", "true");
                subFiringEntriesNode.setAttribute("classPath", "java.lang.Integer");
                subFiringEntriesNode.setAttribute("name", "firingEntries");

                for (int j = 0; j < numOfClasses; j++) {
                    JobClass currentClass = simModel.getClasses().get(j);

                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(currentClass.getName()));
                    subFiringEntriesNode.appendChild(refClassNode);

                    Element subParameterNode = simDoc.createElement("subParameter");
                    subParameterNode.setAttribute("classPath", "java.lang.Integer");
                    subParameterNode.setAttribute("name", "firingEntry");

                    Element valueNode2 = simDoc.createElement("value");
                    valueNode2.appendChild(simDoc.createTextNode(String.valueOf((int) ((TransitionNodeParam) sn.nodeparam.get(istNode)).firing.get(m).get((int) outputs.get(k), j))));

                    subParameterNode.appendChild(valueNode2);
                    subFiringEntriesNode.appendChild(subParameterNode);
                    subFiringVectorNode.appendChild(subFiringEntriesNode);
                }
                subFiringVectorsNode.appendChild(subFiringVectorNode);
            }
            subFiringOutcomeNode.appendChild(subFiringVectorsNode);
            firingOutcomesNode.appendChild(subFiringOutcomeNode);
        }
        section.appendChild(firingOutcomesNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveFiringPriorities(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element firingPrioritiesNode = simDoc.createElement("parameter");
        firingPrioritiesNode.setAttribute("classPath", "java.lang.Integer");
        firingPrioritiesNode.setAttribute("name", "firingPriorities");
        firingPrioritiesNode.setAttribute("array", "true");

        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subFiringPriorityNode = simDoc.createElement("subParameter");
            subFiringPriorityNode.setAttribute("classPath", "java.lang.Integer");
            subFiringPriorityNode.setAttribute("name", "firingPriority");

            Element valueNode = simDoc.createElement("value");
            double firingPrio = ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprio.get(m);
            if (Double.isInfinite(firingPrio) || firingPrio == Integer.MAX_VALUE) {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) firingPrio)));
            }
            subFiringPriorityNode.appendChild(valueNode);
            firingPrioritiesNode.appendChild(subFiringPriorityNode);
        }
        section.appendChild(firingPrioritiesNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveFiringWeights(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element firingWeightsNode = simDoc.createElement("parameter");
        firingWeightsNode.setAttribute("classPath", "java.lang.Double");
        firingWeightsNode.setAttribute("name", "firingWeights");
        firingWeightsNode.setAttribute("array", "true");

        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subFiringWeightNode = simDoc.createElement("subParameter");
            subFiringWeightNode.setAttribute("classPath", "java.lang.Double");
            subFiringWeightNode.setAttribute("name", "firingWeight");

            Element valueNode = simDoc.createElement("value");
            double firingWeights = ((TransitionNodeParam) sn.nodeparam.get(istNode)).fireweight.get(m);

            if (Double.isInfinite(firingWeights) || firingWeights == Integer.MAX_VALUE) {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(firingWeights)));
            }
            subFiringWeightNode.appendChild(valueNode);
            firingWeightsNode.appendChild(subFiringWeightNode);
        }
        section.appendChild(firingWeightsNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveForkStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEFORKSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Node istNode = simModel.getNodes().get(ind);

        Element jplNode = simDoc.createElement("parameter");
        jplNode.setAttribute("classPath", "java.lang.Integer");
        jplNode.setAttribute("name", "jobsPerLink");
        Element valueNode = simDoc.createElement("value");
        
        // Check if the node parameter exists and is not null
        NodeParam nodeParam = sn.nodeparam.get(istNode);
        int fanOut = 1; // default value
        Matrix fanOutLink = null, fanOutProb = null;
        jline.lang.processes.DiscreteSampler[][] fanOutDist = null;
        if (nodeParam instanceof ForkNodeParam) {
            ForkNodeParam forkParam = (ForkNodeParam) nodeParam;
            if (!Double.isNaN(forkParam.fanOut)) {
                fanOut = (int) forkParam.fanOut;
            }
            fanOutLink = forkParam.fanOutLink;
            fanOutProb = forkParam.fanOutProb;
            fanOutDist = forkParam.fanOutDist;
        }

        valueNode.appendChild(simDoc.createTextNode(String.valueOf(fanOut)));
        jplNode.appendChild(valueNode);
        section.appendChild(jplNode);

        Element blockNode = simDoc.createElement("parameter");
        blockNode.setAttribute("classPath", "java.lang.Integer");
        blockNode.setAttribute("name", "block");
        valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode("-1"));
        blockNode.appendChild(valueNode);
        section.appendChild(blockNode);

        // isSimplifiedFork lets JMT ignore the branch list and send one job down
        // every link. That is only the same model when every branch is certain
        // and carries the same number of tasks, so a variable forking level
        // switches it off and makes JMT read the per-branch entries below.
        boolean isSimplified = true;
        if (fanOutLink != null) {
            for (int k = 0; k < fanOutLink.getNumRows() && isSimplified; k++)
                for (int r = 0; r < fanOutLink.getNumCols() && isSimplified; r++) {
                    if (fanOutProb.get(k, r) == 0.0) continue;  // link not taken
                    if (fanOutProb.get(k, r) != 1.0) isSimplified = false;
                    if (fanOutLink.get(k, r) != fanOut) isSimplified = false;
                    if (fanOutDist != null && fanOutDist[k][r] != null) isSimplified = false;
                }
        }

        Element issimplNode = simDoc.createElement("parameter");
        issimplNode.setAttribute("classPath", "java.lang.Boolean");
        issimplNode.setAttribute("name", "isSimplifiedFork");
        valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(isSimplified ? "true" : "false"));
        issimplNode.appendChild(valueNode);
        section.appendChild(issimplNode);

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategy");
        strategyNode.setAttribute("name", "ForkStrategy");

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            JobClass rstJobClass = sn.jobclasses.get(r);
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode);

            Element classStratNode = simDoc.createElement("subParameter");
            classStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork");
            classStratNode.setAttribute("name", "Branch Probabilities");
            Element classStratNode2 = simDoc.createElement("subParameter");
            classStratNode2.setAttribute("array", "true");
            classStratNode2.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategies.OutPath");
            classStratNode2.setAttribute("name", "EmpiricalEntryArray");
            switch (sn.routing.get(istNode).get(rstJobClass)) {
                case PROB:
                    Matrix conn_i = new Matrix(0, 0);
                    Matrix.extractRows(sn.connmatrix, ind, ind + 1, conn_i);
                    Matrix conn_i_find = conn_i.find();

                    // One OutPathEntry per outgoing link. An earlier version
                    // built the entry inside this loop but appended it outside,
                    // so only the last connected node survived; that was
                    // invisible because isSimplifiedFork makes JMT send one job
                    // down every link and ignore the branch list.
                    for (int idx = 0; idx < conn_i_find.length(); idx++) {
                        int k = (int) conn_i_find.get(idx);
                        Element classStratNode3 = simDoc.createElement("subParameter");
                        classStratNode3.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategies.OutPath");
                        classStratNode3.setAttribute("name", "OutPathEntry");

                        Element classStratNode4 = simDoc.createElement("subParameter");
                        classStratNode4.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                        classStratNode4.setAttribute("name", "outUnitProbability");
                        Element classStratNode4Station = simDoc.createElement("subParameter");
                        classStratNode4Station.setAttribute("classPath", "java.lang.String");
                        classStratNode4Station.setAttribute("name", "stationName");
                        Element classStratNode4StationValueNode = simDoc.createElement("value");
                        classStratNode4StationValueNode.appendChild(simDoc.createTextNode(String.format("%s", sn.nodenames.get(k))));
                        classStratNode4Station.appendChild(classStratNode4StationValueNode);
                        classStratNode4.appendChild(classStratNode4Station);
                        // branch activation probability: JMT's outUnitProbability
                        double branchP = (fanOutProb == null) ? 1.0 : fanOutProb.get(k, r);
                        Element classStratNode4Probability = simDoc.createElement("subParameter");
                        classStratNode4Probability.setAttribute("classPath", "java.lang.Double");
                        classStratNode4Probability.setAttribute("name", "probability");
                        Element classStratNode4ProbabilityValueNode = simDoc.createElement("value");
                        classStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode(String.valueOf(branchP)));
                        classStratNode4Probability.appendChild(classStratNode4ProbabilityValueNode);
                        classStratNode4.appendChild(classStratNode4Probability);
                        classStratNode3.appendChild(classStratNode4);

                        // JobsPerLinkDis is an EmpiricalEntry ARRAY: one entry
                        // per point of the jobs-per-link distribution. A
                        // deterministic fork emits the single degenerate entry
                        // it always did.
                        double[] jplPoints;
                        double[] jplProbs;
                        jline.lang.processes.DiscreteSampler jplDist =
                                (fanOutDist == null) ? null : fanOutDist[k][r];
                        if (jplDist != null) {
                            Matrix pm = (Matrix) jplDist.getParam(1).getValue();
                            Matrix xm = (Matrix) jplDist.getParam(2).getValue();
                            jplPoints = new double[xm.length()];
                            jplProbs = new double[pm.length()];
                            double tot = 0.0;
                            for (int e = 0; e < pm.length(); e++) tot += pm.get(e);
                            for (int e = 0; e < pm.length(); e++) {
                                jplPoints[e] = xm.get(e);
                                jplProbs[e] = pm.get(e) / tot;
                            }
                        } else if (fanOutLink != null) {
                            jplPoints = new double[]{fanOutLink.get(k, r)};
                            jplProbs = new double[]{1.0};
                        } else {
                            jplPoints = new double[]{fanOut};
                            jplProbs = new double[]{1.0};
                        }

                        Element classStratNode4b = simDoc.createElement("subParameter");
                        classStratNode4b.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                        classStratNode4b.setAttribute("array", "true");
                        classStratNode4b.setAttribute("name", "JobsPerLinkDis");
                        for (int e = 0; e < jplPoints.length; e++) {
                            Element classStratNode5b = simDoc.createElement("subParameter");
                            classStratNode5b.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                            classStratNode5b.setAttribute("name", "EmpiricalEntry");
                            Element classStratNode5bStation = simDoc.createElement("subParameter");
                            classStratNode5bStation.setAttribute("classPath", "java.lang.String");
                            classStratNode5bStation.setAttribute("name", "numbers");
                            Element classStratNode5bStationValueNode = simDoc.createElement("value");
                            classStratNode5bStationValueNode.appendChild(
                                    simDoc.createTextNode(String.valueOf((int) Math.round(jplPoints[e]))));
                            classStratNode5bStation.appendChild(classStratNode5bStationValueNode);
                            classStratNode5b.appendChild(classStratNode5bStation);
                            Element classStratNode5bProbability = simDoc.createElement("subParameter");
                            classStratNode5bProbability.setAttribute("classPath", "java.lang.Double");
                            classStratNode5bProbability.setAttribute("name", "probability");
                            Element classStratNode5bProbabilityValueNode = simDoc.createElement("value");
                            classStratNode5bProbabilityValueNode.appendChild(
                                    simDoc.createTextNode(String.valueOf(jplProbs[e])));
                            classStratNode5bProbability.appendChild(classStratNode5bProbabilityValueNode);
                            classStratNode5b.appendChild(classStratNode5bProbability);
                            classStratNode4b.appendChild(classStratNode5b);
                        }
                        classStratNode3.appendChild(classStratNode4b);

                        classStratNode2.appendChild(classStratNode3);
                    }
            }
            classStratNode.appendChild(classStratNode2);
            strategyNode.appendChild(classStratNode);
        }
        section.appendChild(strategyNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveGetStrategy(DocumentSectionPair documentSectionPair) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        // the get strategy is always fcfs
        Element queueGetStrategyNode = simDoc.createElement("parameter");
        queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy");
        queueGetStrategyNode.setAttribute("name", "FCFSstrategy");
        section.appendChild(queueGetStrategyNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveGetStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        int ist = (int) sn.nodeToStation.get(ind);
        // For non-station nodes (e.g., Logger nodes have nodeToStation = -1),
        // use default FCFS strategy since they're not real queues
        if (ist < 0) {
            Element queueGetStrategyNode = simDoc.createElement("parameter");
            queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy");
            queueGetStrategyNode.setAttribute("name", "FCFSstrategy");
            section.appendChild(queueGetStrategyNode);
            return documentSectionPair;
        }
        Station istStation = sn.stations.get(ist);

        if (sn.nodetype.get(ind) == NodeType.Queue && sn.sched.get(istStation) == SchedStrategy.POLLING) {
            Element queueGetStrategyNode = simDoc.createElement("parameter");
            // Get the polling server from the queue to determine polling type
            Node node = simModel.getNodes().get(ind);
            if (node instanceof jline.lang.nodes.Queue) {
                jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                if (queue.getServer() instanceof jline.lang.sections.PollingServer) {
                    jline.lang.sections.PollingServer pollingServer = (jline.lang.sections.PollingServer) queue.getServer();
                    PollingType pollingType = pollingServer.getPollingType();
                    
                    switch (pollingType) {
                        case GATED:
                            queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.GatedPollingGetStrategy");
                            break;
                        case EXHAUSTIVE:
                            queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.ExhaustivePollingGetStrategy");
                            break;
                        case KLIMITED:
                            queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.LimitedPollingGetStrategy");
                            Element pollingKNode = simDoc.createElement("subParameter");
                            pollingKNode.setAttribute("classPath", "java.lang.Integer");
                            pollingKNode.setAttribute("name", "pollingKValue");
                            Element valueNode = simDoc.createElement("value");
                            valueNode.appendChild(simDoc.createTextNode("1")); // Default K value
                            pollingKNode.appendChild(valueNode);
                            queueGetStrategyNode.appendChild(pollingKNode);
                            break;
                        case DECREMENTING:
                            throw new RuntimeException("JMT does not support the decrementing "
                                    + "(semiexhaustive) polling discipline; use the LDES solver.");
                    }
                }
            }
            queueGetStrategyNode.setAttribute("name", "FCFSstrategy");
            section.appendChild(queueGetStrategyNode);
        } else {
            // the get strategy is always fcfs for queues unless a polling queue
            Element queueGetStrategyNode = simDoc.createElement("parameter");
            queueGetStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy");
            queueGetStrategyNode.setAttribute("name", "FCFSstrategy");
            section.appendChild(queueGetStrategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair setPollingServerClassName(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;
        
        Node node = simModel.getNodes().get(ind);
        
        if (node instanceof jline.lang.nodes.Queue) {
            jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
            if (queue.getServer() instanceof PollingServer) {
                PollingServer pollingServer = (PollingServer) queue.getServer();
                PollingType pollingType = pollingServer.getPollingType();
                
                switch (pollingType) {
                    case GATED:
                        section.setAttribute("className", "GatedPollingServer");
                        break;
                    case EXHAUSTIVE:
                        section.setAttribute("className", "ExhaustivePollingServer");
                        break;
                    case KLIMITED:
                        section.setAttribute("className", "LimitedPollingServer");
                        break;
                    case DECREMENTING:
                        throw new RuntimeException("JMT does not support the decrementing "
                                + "(semiexhaustive) polling discipline; use the LDES solver.");
                }
            }
        }
        
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveSwitchoverStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        int numOfClasses = sn.nclasses;
        int isf = (int) sn.nodeToStation.get(ind);
        Station isfStation = sn.stations.get(isf);

        if (sn.sched.get(isfStation) == SchedStrategy.POLLING) {
            Element paramNode = simDoc.createElement("parameter");
            paramNode.setAttribute("array", "true");
            paramNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
            paramNode.setAttribute("name", "SwitchoverStrategy");
            
            Node node = simModel.getNodes().get(ind);
            if (node instanceof jline.lang.nodes.Queue) {
                jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                if (queue.getServer() instanceof PollingServer) {
                    PollingServer pollingServer = (PollingServer) queue.getServer();
                    
                    for (int r = 0; r < numOfClasses; r++) {
                        Element refClassNode = simDoc.createElement("refClass");
                        refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                        paramNode.appendChild(refClassNode);

                        JobClass jobClass = sn.jobclasses.get(r);
                        jline.lang.processes.Distribution switchover = pollingServer.getSwitchover(jobClass);

                        Element serviceTimeStrategyNode = simDoc.createElement("subParameter");

                        // Check if switchover is zero/immediate (null, Immediate, or Exp with zero mean)
                        boolean isZeroSwitchover = (switchover == null) ||
                            (switchover instanceof jline.lang.processes.Immediate) ||
                            (switchover instanceof jline.lang.processes.Exp &&
                             ((jline.lang.processes.Exp) switchover).getMean() == 0);

                        if (isZeroSwitchover) {
                            // Use ZeroServiceTimeStrategy for zero switchover time
                            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                            serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
                        } else if (switchover instanceof jline.lang.processes.Exp) {
                            // Use exponential distribution for non-zero Exp switchover
                            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                            serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

                            Element distributionNode = simDoc.createElement("subParameter");
                            distributionNode.setAttribute("classPath", "jmt.engine.random.Exponential");
                            distributionNode.setAttribute("name", "Exponential");
                            serviceTimeStrategyNode.appendChild(distributionNode);

                            Element distrParNode = simDoc.createElement("subParameter");
                            distrParNode.setAttribute("classPath", "jmt.engine.random.ExponentialPar");
                            distrParNode.setAttribute("name", "distrPar");

                            Element subParNodeLambda = simDoc.createElement("subParameter");
                            subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                            subParNodeLambda.setAttribute("name", "lambda");
                            Element subParValue = simDoc.createElement("value");
                            double lambda = ((jline.lang.processes.Exp) switchover).getRate();
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", lambda)));
                            subParNodeLambda.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeLambda);
                            serviceTimeStrategyNode.appendChild(distrParNode);
                        } else {
                            // Default to ZeroServiceTimeStrategy for unsupported distributions
                            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                            serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
                        }
                        paramNode.appendChild(serviceTimeStrategyNode);
                    }
                }
            }
            section.appendChild(paramNode);
        } else {
            // For GatedPollingServer and other polling servers even when not using POLLING strategy
            // We need to check if this is actually a polling server being saved as GatedPollingServer
            Node node = simModel.getNodes().get(ind);
            boolean isPollingServer = false;
            
            if (node instanceof jline.lang.nodes.Queue) {
                jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                isPollingServer = queue.getServer() instanceof PollingServer;
            }
            
            if (isPollingServer) {
                // Generate proper switchover strategy for polling servers
                Element paramNode = simDoc.createElement("parameter");
                paramNode.setAttribute("array", "true");
                paramNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
                paramNode.setAttribute("name", "SwitchoverStrategy");
                
                jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                PollingServer pollingServer = (PollingServer) queue.getServer();
                
                for (int r = 0; r < numOfClasses; r++) {
                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    paramNode.appendChild(refClassNode);

                    JobClass jobClass = sn.jobclasses.get(r);
                    jline.lang.processes.Distribution switchover = pollingServer.getSwitchover(jobClass);

                    Element serviceTimeStrategyNode = simDoc.createElement("subParameter");

                    // Check if switchover is zero/immediate (null, Immediate, or Exp with zero mean)
                    boolean isZeroSwitchover = (switchover == null) ||
                        (switchover instanceof jline.lang.processes.Immediate) ||
                        (switchover instanceof jline.lang.processes.Exp &&
                         ((jline.lang.processes.Exp) switchover).getMean() == 0);

                    if (isZeroSwitchover) {
                        // Use ZeroServiceTimeStrategy for zero switchover time
                        serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                        serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
                        paramNode.appendChild(serviceTimeStrategyNode);
                        continue;
                    }

                    serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                    serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

                    if (switchover instanceof jline.lang.processes.Exp) {
                        Element distributionNode = simDoc.createElement("subParameter");
                        distributionNode.setAttribute("classPath", "jmt.engine.random.Exponential");
                        distributionNode.setAttribute("name", "Exponential");
                        serviceTimeStrategyNode.appendChild(distributionNode);
                        
                        Element distrParNode = simDoc.createElement("subParameter");
                        distrParNode.setAttribute("classPath", "jmt.engine.random.ExponentialPar");
                        distrParNode.setAttribute("name", "distrPar");
                        
                        Element subParNodeLambda = simDoc.createElement("subParameter");
                        subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda.setAttribute("name", "lambda");
                        Element subParValue = simDoc.createElement("value");
                        double lambda = ((jline.lang.processes.Exp) switchover).getRate();
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", lambda)));
                        subParNodeLambda.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeLambda);
                        serviceTimeStrategyNode.appendChild(distrParNode);
                    } else if (switchover instanceof jline.lang.processes.Erlang) {
                        Element distributionNode = simDoc.createElement("subParameter");
                        distributionNode.setAttribute("classPath", "jmt.engine.random.Erlang");
                        distributionNode.setAttribute("name", "Erlang");
                        serviceTimeStrategyNode.appendChild(distributionNode);
                        
                        Element distrParNode = simDoc.createElement("subParameter");
                        distrParNode.setAttribute("classPath", "jmt.engine.random.ErlangPar");
                        distrParNode.setAttribute("name", "distrPar");
                        
                        jline.lang.processes.Erlang erlang = (jline.lang.processes.Erlang) switchover;
                        
                        Element alphaNode = simDoc.createElement("subParameter");
                        alphaNode.setAttribute("classPath", "java.lang.Double");
                        alphaNode.setAttribute("name", "alpha");
                        Element alphaValue = simDoc.createElement("value");
                        double alpha = (double) erlang.getParam(1).getValue();
                        alphaValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha)));
                        alphaNode.appendChild(alphaValue);
                        distrParNode.appendChild(alphaNode);
                        
                        Element rNode = simDoc.createElement("subParameter");
                        rNode.setAttribute("classPath", "java.lang.Long");
                        rNode.setAttribute("name", "r");
                        Element rValue = simDoc.createElement("value");
                        int phases = (int) erlang.getParam(2).getValue();
                        rValue.appendChild(simDoc.createTextNode(String.valueOf((long)phases)));
                        rNode.appendChild(rValue);
                        distrParNode.appendChild(rNode);
                        
                        serviceTimeStrategyNode.appendChild(distrParNode);
                    } else {
                        // Default to Immediate (zero time) for other distributions
                        Element distributionNode = simDoc.createElement("subParameter");
                        distributionNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistr");
                        distributionNode.setAttribute("name", "Deterministic");
                        serviceTimeStrategyNode.appendChild(distributionNode);
                        
                        Element distrParNode = simDoc.createElement("subParameter");
                        distrParNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistrPar");
                        distrParNode.setAttribute("name", "distrPar");
                        
                        Element tNode = simDoc.createElement("subParameter");
                        tNode.setAttribute("classPath", "java.lang.Double");
                        tNode.setAttribute("name", "t");
                        Element tValue = simDoc.createElement("value");
                        tValue.appendChild(simDoc.createTextNode("0.0"));
                        tNode.appendChild(tValue);
                        distrParNode.appendChild(tNode);
                        serviceTimeStrategyNode.appendChild(distrParNode);
                    }
                    
                    paramNode.appendChild(serviceTimeStrategyNode);
                }
                section.appendChild(paramNode);
            } else {
                // Non-polling case - simplified implementation
                boolean[] exportClasses = getExportableClasses();
                Element paramNode = simDoc.createElement("parameter");
                paramNode.setAttribute("array", "true");
                paramNode.setAttribute("classPath", "java.lang.Object");
                paramNode.setAttribute("name", "SwitchoverStrategy");

                for (int r = 0; r < numOfClasses; r++) {
                    if (!exportClasses[r]) continue;
                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    paramNode.appendChild(refClassNode);
                }
                section.appendChild(paramNode);
            }
        }
        
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves delay-off and setup time strategies for Queue nodes.
     * This exports the switchoverStrategies, delayOffTime and setUpTime parameters to JMT XML format
     * when the queue has delay-off times enabled.
     *
     * <p>JMT's Server constructor expects parameters in order:
     * switchoverStrategies, delayOffStrategies, setUpStrategies
     *
     * @param documentSectionPair the document/section pair to append to
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveDelayOffStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Node node = simModel.getNodes().get(ind);

        // Only process Queue nodes with delay-off enabled
        if (!(node instanceof jline.lang.nodes.Queue)) {
            return new DocumentSectionPair(simDoc, section);
        }

        jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
        if (!queue.isDelayOffEnabled()) {
            return new DocumentSectionPair(simDoc, section);
        }

        int numOfClasses = sn.nclasses;

        // NOTE: Export a null SwitchoverStrategy for delayoff queues
        // JMT Server.java checks "if (switchoverStrategies != null)" at line 1123 and uses
        // the switchover code path instead of the delayoff code path. To properly use
        // delayoff/setup logic, switchoverStrategies must be null.
        // The Server constructor with delayoff requires 6 parameters:
        // (numberOfServers, numberOfVisits, serviceStrategies, switchoverStrategies, delayOffStrategies, setUpStrategies)
        // So we must explicitly pass null for switchoverStrategies.
        Element switchoverParamNode = simDoc.createElement("parameter");
        switchoverParamNode.setAttribute("array", "true");
        switchoverParamNode.setAttribute("classPath", "java.lang.Object");
        switchoverParamNode.setAttribute("name", "switchoverStrategies");
        Element nullValue = simDoc.createElement("value");
        nullValue.appendChild(simDoc.createTextNode("null"));
        switchoverParamNode.appendChild(nullValue);
        section.appendChild(switchoverParamNode);

        // Export delayOffTime parameter
        // JMT expects a 2D array structure:
        // <parameter name="delayOffTime" array="true" classPath="java.lang.Object">
        //   <refClass>ClassName</refClass>
        //   <subParameter name="delayOffTime" classPath="jmt.engine.NetStrategies.ServiceStrategy">
        //     <refClass>ClassName</refClass>  <!-- Inner refClass with same class name -->
        //     <subParameter name="ServiceTimeStrategy" classPath="...ServiceTimeStrategy">
        //       <!-- distribution content -->
        //     </subParameter>
        //   </subParameter>
        // </parameter>
        Element delayOffParamNode = simDoc.createElement("parameter");
        delayOffParamNode.setAttribute("array", "true");
        delayOffParamNode.setAttribute("classPath", "java.lang.Object");
        delayOffParamNode.setAttribute("name", "delayOffTime");

        for (int r = 0; r < numOfClasses; r++) {
            JobClass jobClass = sn.jobclasses.get(r);
            String className = sn.classnames.get(r);

            // Outer refClass element
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(className));
            delayOffParamNode.appendChild(refClassNode);

            // Intermediate subParameter (delayOffTime -> ServiceStrategy)
            // array="true" is required so SimLoader creates a ServiceStrategy[] array
            Element intermediateSubParam = simDoc.createElement("subParameter");
            intermediateSubParam.setAttribute("array", "true");
            intermediateSubParam.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
            intermediateSubParam.setAttribute("name", "delayOffTime");

            // Inner refClass element (same class name)
            Element innerRefClassNode = simDoc.createElement("refClass");
            innerRefClassNode.appendChild(simDoc.createTextNode(className));
            intermediateSubParam.appendChild(innerRefClassNode);

            // Inner ServiceTimeStrategy subParameter with actual distribution
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

            jline.lang.processes.Distribution delayOffDist = queue.getDelayOffTime(jobClass);
            appendDistributionXml(simDoc, serviceTimeStrategyNode, delayOffDist);

            intermediateSubParam.appendChild(serviceTimeStrategyNode);
            delayOffParamNode.appendChild(intermediateSubParam);
        }
        section.appendChild(delayOffParamNode);

        // Export setUpTime parameter (same 2D array structure as delayOffTime)
        Element setupParamNode = simDoc.createElement("parameter");
        setupParamNode.setAttribute("array", "true");
        setupParamNode.setAttribute("classPath", "java.lang.Object");
        setupParamNode.setAttribute("name", "setUpTime");

        for (int r = 0; r < numOfClasses; r++) {
            JobClass jobClass = sn.jobclasses.get(r);
            String className = sn.classnames.get(r);

            // Outer refClass element
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(className));
            setupParamNode.appendChild(refClassNode);

            // Intermediate subParameter (setUpTime -> ServiceStrategy)
            // array="true" is required so SimLoader creates a ServiceStrategy[] array
            Element intermediateSubParam = simDoc.createElement("subParameter");
            intermediateSubParam.setAttribute("array", "true");
            intermediateSubParam.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
            intermediateSubParam.setAttribute("name", "setUpTime");

            // Inner refClass element (same class name)
            Element innerRefClassNode = simDoc.createElement("refClass");
            innerRefClassNode.appendChild(simDoc.createTextNode(className));
            intermediateSubParam.appendChild(innerRefClassNode);

            // Inner ServiceTimeStrategy subParameter with actual distribution
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

            jline.lang.processes.Distribution setupDist = queue.getSetupTime(jobClass);
            appendDistributionXml(simDoc, serviceTimeStrategyNode, setupDist);

            intermediateSubParam.appendChild(serviceTimeStrategyNode);
            setupParamNode.appendChild(intermediateSubParam);
        }
        section.appendChild(setupParamNode);

        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Appends distribution XML elements to a parent node for JMT export.
     * Supports Exp, Erlang, Det, and Immediate distributions.
     *
     * @param simDoc the XML document
     * @param parentNode the parent element to append to
     * @param dist the distribution to serialize
     */
    private void appendDistributionXml(Document simDoc, Element parentNode, jline.lang.processes.Distribution dist) {
        if (dist == null || dist instanceof jline.lang.processes.Immediate) {
            // Zero/Immediate service time
            Element distributionNode = simDoc.createElement("subParameter");
            distributionNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistr");
            distributionNode.setAttribute("name", "Deterministic");
            parentNode.appendChild(distributionNode);

            Element distrParNode = simDoc.createElement("subParameter");
            distrParNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistrPar");
            distrParNode.setAttribute("name", "distrPar");

            Element tNode = simDoc.createElement("subParameter");
            tNode.setAttribute("classPath", "java.lang.Double");
            tNode.setAttribute("name", "t");
            Element tValue = simDoc.createElement("value");
            tValue.appendChild(simDoc.createTextNode("0.0"));
            tNode.appendChild(tValue);
            distrParNode.appendChild(tNode);
            parentNode.appendChild(distrParNode);

        } else if (dist instanceof jline.lang.processes.Exp) {
            jline.lang.processes.Exp expDist = (jline.lang.processes.Exp) dist;

            Element distributionNode = simDoc.createElement("subParameter");
            distributionNode.setAttribute("classPath", "jmt.engine.random.Exponential");
            distributionNode.setAttribute("name", "Exponential");
            parentNode.appendChild(distributionNode);

            Element distrParNode = simDoc.createElement("subParameter");
            distrParNode.setAttribute("classPath", "jmt.engine.random.ExponentialPar");
            distrParNode.setAttribute("name", "distrPar");

            Element lambdaNode = simDoc.createElement("subParameter");
            lambdaNode.setAttribute("classPath", "java.lang.Double");
            lambdaNode.setAttribute("name", "lambda");
            Element lambdaValue = simDoc.createElement("value");
            lambdaValue.appendChild(simDoc.createTextNode(String.format("%.12f", expDist.getRate())));
            lambdaNode.appendChild(lambdaValue);
            distrParNode.appendChild(lambdaNode);
            parentNode.appendChild(distrParNode);

        } else if (dist instanceof jline.lang.processes.Erlang) {
            jline.lang.processes.Erlang erlangDist = (jline.lang.processes.Erlang) dist;

            Element distributionNode = simDoc.createElement("subParameter");
            distributionNode.setAttribute("classPath", "jmt.engine.random.Erlang");
            distributionNode.setAttribute("name", "Erlang");
            parentNode.appendChild(distributionNode);

            Element distrParNode = simDoc.createElement("subParameter");
            distrParNode.setAttribute("classPath", "jmt.engine.random.ErlangPar");
            distrParNode.setAttribute("name", "distrPar");

            Element alphaNode = simDoc.createElement("subParameter");
            alphaNode.setAttribute("classPath", "java.lang.Double");
            alphaNode.setAttribute("name", "alpha");
            Element alphaValue = simDoc.createElement("value");
            double mean = erlangDist.getMean();
            int phases = (int) erlangDist.getParam(2).getValue();
            double alpha = phases / mean;
            alphaValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha)));
            alphaNode.appendChild(alphaValue);
            distrParNode.appendChild(alphaNode);

            Element rNode = simDoc.createElement("subParameter");
            rNode.setAttribute("classPath", "java.lang.Long");
            rNode.setAttribute("name", "r");
            Element rValue = simDoc.createElement("value");
            rValue.appendChild(simDoc.createTextNode(String.valueOf((long) phases)));
            rNode.appendChild(rValue);
            distrParNode.appendChild(rNode);
            parentNode.appendChild(distrParNode);

        } else if (dist instanceof jline.lang.processes.Det) {
            jline.lang.processes.Det detDist = (jline.lang.processes.Det) dist;

            Element distributionNode = simDoc.createElement("subParameter");
            distributionNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistr");
            distributionNode.setAttribute("name", "Deterministic");
            parentNode.appendChild(distributionNode);

            Element distrParNode = simDoc.createElement("subParameter");
            distrParNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistrPar");
            distrParNode.setAttribute("name", "distrPar");

            Element tNode = simDoc.createElement("subParameter");
            tNode.setAttribute("classPath", "java.lang.Double");
            tNode.setAttribute("name", "t");
            Element tValue = simDoc.createElement("value");
            tValue.appendChild(simDoc.createTextNode(String.format("%.12f", detDist.getMean())));
            tNode.appendChild(tValue);
            distrParNode.appendChild(tNode);
            parentNode.appendChild(distrParNode);

        } else {
            // Default fallback: use mean as deterministic value
            Element distributionNode = simDoc.createElement("subParameter");
            distributionNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistr");
            distributionNode.setAttribute("name", "Deterministic");
            parentNode.appendChild(distributionNode);

            Element distrParNode = simDoc.createElement("subParameter");
            distrParNode.setAttribute("classPath", "jmt.engine.random.DeterministicDistrPar");
            distrParNode.setAttribute("name", "distrPar");

            Element tNode = simDoc.createElement("subParameter");
            tNode.setAttribute("classPath", "java.lang.Double");
            tNode.setAttribute("name", "t");
            Element tValue = simDoc.createElement("value");
            tValue.appendChild(simDoc.createTextNode(String.format("%.12f", dist.getMean())));
            tNode.appendChild(tValue);
            distrParNode.appendChild(tNode);
            parentNode.appendChild(distrParNode);
        }
    }

    public DocumentSectionPair saveInhibitingConditions(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element inhibitingConditionsNode = simDoc.createElement("parameter");
        inhibitingConditionsNode.setAttribute("array", "true");
        inhibitingConditionsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
        inhibitingConditionsNode.setAttribute("name", "inhibitingConditions");


        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractColumn(sn.connmatrix, ind, conn_i);
        Matrix inputs = conn_i.find();
        List<String> connections = new ArrayList<String>();
        for (int idx = 0; idx < inputs.length(); idx++) {
            connections.add(sn.nodenames.get((int) inputs.get(idx)));
        }
        int numOfInputs = connections.size();

        int numOfClasses = sn.nclasses;
        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subInhibitingConditionNode = simDoc.createElement("subParameter");
            subInhibitingConditionNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            subInhibitingConditionNode.setAttribute("name", "inhibitingCondition");

            Element subInhibitingVectorsNode = simDoc.createElement("subParameter");
            subInhibitingVectorsNode.setAttribute("array", "true");
            subInhibitingVectorsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
            subInhibitingVectorsNode.setAttribute("name", "inhibitingVectors");

            for (int k = 0; k < numOfInputs; k++) {
                Element subInhibitingVectorNode = simDoc.createElement("subParameter");
                subInhibitingVectorNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
                subInhibitingVectorNode.setAttribute("name", "inhibitingVector");

                Element subStationNameNode = simDoc.createElement("subParameter");
                subStationNameNode.setAttribute("classPath", "java.lang.String");
                subStationNameNode.setAttribute("name", "stationName");

                Element placeNameValueNode = simDoc.createElement("value");
                placeNameValueNode.appendChild(simDoc.createTextNode(connections.get(k)));
                subStationNameNode.appendChild(placeNameValueNode);

                subInhibitingVectorNode.appendChild(subStationNameNode);

                Element subInhibitingEntriesNode = simDoc.createElement("subParameter");
                subInhibitingEntriesNode.setAttribute("array", "true");
                subInhibitingEntriesNode.setAttribute("classPath", "java.lang.Integer");
                subInhibitingEntriesNode.setAttribute("name", "inhibitingEntries");

                for (int r = 0; r < numOfClasses; r++) {
                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    subInhibitingEntriesNode.appendChild(refClassNode);

                    Element subParameterNode = simDoc.createElement("subParameter");
                    subParameterNode.setAttribute("classPath", "java.lang.Integer");
                    subParameterNode.setAttribute("name", "inhibitingEntry");

                    Element valueNode2 = simDoc.createElement("value");

                    if (Utils.isInf(((TransitionNodeParam) sn.nodeparam.get(istNode)).inhibiting.get(m).get((int) inputs.get(k), r))) {
                        valueNode2.appendChild(simDoc.createTextNode("0"));
                    } else {
                        valueNode2.appendChild(simDoc.createTextNode(String.valueOf((int) ((TransitionNodeParam) sn.nodeparam.get(istNode)).inhibiting.get(m).get((int) inputs.get(k), r))));
                    }

                    subParameterNode.appendChild(valueNode2);
                    subInhibitingEntriesNode.appendChild(subParameterNode);
                    subInhibitingVectorNode.appendChild(subInhibitingEntriesNode);
                }
                subInhibitingVectorsNode.appendChild(subInhibitingVectorNode);
            }
            subInhibitingConditionNode.appendChild(subInhibitingVectorsNode);
            inhibitingConditionsNode.appendChild(subInhibitingConditionNode);
        }
        section.appendChild(inhibitingConditionsNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveJoinStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEJOINSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.JoinStrategy");
        strategyNode.setAttribute("name", "JoinStrategy");

        int numOfClasses = sn.nclasses;

        Element refClassNode2 = simDoc.createElement("refClass");
        // Declare variables outside switch to avoid scope issues
        Element joinStrategyNode;
        Element reqNode;
        Element valueNode;
        
        for (int r = 0; r < numOfClasses; r++) {
            Node istNode = simModel.getNodes().get(ind);
            JobClass rstJobClass = sn.jobclasses.get(r);
            switch (((JoinNodeParam) sn.nodeparam.get(istNode)).joinStrategy.get(rstJobClass)) {
                case STD:
                    refClassNode2 = simDoc.createElement("refClass");
                    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    strategyNode.appendChild(refClassNode2);

                    joinStrategyNode = simDoc.createElement("subParameter");
                    joinStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.JoinStrategies.NormalJoin");
                    joinStrategyNode.setAttribute("name", "Standard Join");
                    reqNode = simDoc.createElement("subParameter");
                    reqNode.setAttribute("classPath", "java.lang.Integer");
                    reqNode.setAttribute("name", "numRequired");
                    valueNode = simDoc.createElement("value");
                    valueNode.appendChild(simDoc.createTextNode(String.valueOf(((JoinNodeParam) sn.nodeparam.get(istNode)).fanIn.get(rstJobClass).intValue())));
                    reqNode.appendChild(valueNode);
                    joinStrategyNode.appendChild(reqNode);
                    strategyNode.appendChild(joinStrategyNode);
                    section.appendChild(strategyNode);
                    break;
                case PARTIAL:
                    refClassNode2 = simDoc.createElement("refClass");
                    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    strategyNode.appendChild(refClassNode2);

                    joinStrategyNode = simDoc.createElement("subParameter");
                    joinStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.JoinStrategies.PartialJoin");
                    joinStrategyNode.setAttribute("name", "Quorum");
                    reqNode = simDoc.createElement("subParameter");
                    reqNode.setAttribute("classPath", "java.lang.Integer");
                    reqNode.setAttribute("name", "numRequired");
                    valueNode = simDoc.createElement("value");
                    valueNode.appendChild(simDoc.createTextNode(String.valueOf(((JoinNodeParam) sn.nodeparam.get(istNode)).joinRequired.get(rstJobClass))));
                    reqNode.appendChild(valueNode);
                    joinStrategyNode.appendChild(reqNode);
                    strategyNode.appendChild(joinStrategyNode);
                    section.appendChild(strategyNode);
                    break;

            }
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public ElementDocumentPair saveLinks(ElementDocumentPair elementDocumentPair) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;


        for (int j = 0; j < sn.connmatrix.getNumCols(); j++) {
            for (int i = 0; i < sn.connmatrix.getNumRows(); i++) {
                if (sn.connmatrix.get(i, j) != 0) {
                    Element connectionNode = simDoc.createElement("connection");
                    connectionNode.setAttribute("source", sn.nodenames.get(i));
                    connectionNode.setAttribute("target", sn.nodenames.get(j));
                    simElem.appendChild(connectionNode);
                }
            }
        }
        return new ElementDocumentPair(simElem, simDoc);
    }

    public DocumentSectionPair saveLogTunnel(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVELOGTUNNEL(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Node istNode = simModel.getNodes().get(ind);
        List<String> loggerNodesCP = new ArrayList<String>();
        loggerNodesCP.add("java.lang.String"); // index 0
        loggerNodesCP.add("java.lang.String"); // index 1
        for (int i = 2; i < 9; i++) { // indices 2 to 8
            loggerNodesCP.add("java.lang.Boolean");
        }
        loggerNodesCP.add("java.lang.Integer"); // index 9
        List<String> loggerNodesNames = Arrays.asList("logfileName", "logfilePath", "logExecTimestamp",
                "logLoggerName", "logTimeStamp", "logJobID",
                "logJobClass", "logTimeSameClass", "logTimeAnyClass",
                "numClasses");

        int numOfClasses = sn.nclasses;

        // logger specific path does not work in JMT at the moment
        if (!((LoggerNodeParam) sn.nodeparam.get(istNode)).filePath.endsWith(File.separator)) {
            ((LoggerNodeParam) sn.nodeparam.get(istNode)).filePath = ((LoggerNodeParam) sn.nodeparam.get(istNode)).filePath + File.separator;
        }

        List<String> loggerNodesValues = new ArrayList<String>();
        loggerNodesValues.add(((LoggerNodeParam) sn.nodeparam.get(istNode)).fileName.get(0));
        loggerNodesValues.add(((LoggerNodeParam) sn.nodeparam.get(istNode)).filePath);
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).startTime));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).loggerName));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).timestamp));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).jobID));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).jobClass));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).timeSameClass));
        loggerNodesValues.add(String.valueOf(((LoggerNodeParam) sn.nodeparam.get(istNode)).timeAnyClass));
        loggerNodesValues.add(String.valueOf(numOfClasses));

        for (int j = 0; j < loggerNodesValues.size(); j++) {
            Element loggerNode = simDoc.createElement("parameter");
            loggerNode.setAttribute("classPath", loggerNodesCP.get(j));
            loggerNode.setAttribute("name", loggerNodesNames.get(j));
            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(loggerNodesValues.get(j)));
            loggerNode.appendChild(valueNode);
            section.appendChild(loggerNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public ElementDocumentPair saveMetric(ElementDocumentPair elementDocumentPair, AvgHandle handles) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;

        boolean[] exportClasses = getExportableClasses();
        for (int i = 0; i < handles.keySet().size(); i++) {
            Station istStation = sn.stations.get(i);
            for (int r = 0; r < handles.get(istStation).size(); r++) {
                if (!exportClasses[r]) continue;
                JobClass rstJobClass = sn.jobclasses.get(r);
                Metric currentPerformanceIndex = handles.get(istStation).get(rstJobClass);
                if (!currentPerformanceIndex.isDisabled) {
                    Element performanceNode = simDoc.createElement("measure");
                    performanceNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                    performanceNode.setAttribute("name", "Performance_" + (i + 1));
                    // System-level metrics (station is null) use nodeType=""
                    if (currentPerformanceIndex.station == null) {
                        performanceNode.setAttribute("nodeType", "");
                        performanceNode.setAttribute("referenceNode", "");
                    } else {
                        performanceNode.setAttribute("nodeType", "station");
                        performanceNode.setAttribute("referenceNode", currentPerformanceIndex.station.getName());
                    }
                    performanceNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                    performanceNode.setAttribute("referenceUserClass", currentPerformanceIndex.jobClass.getName());
                    performanceNode.setAttribute("type", currentPerformanceIndex.type);
                    performanceNode.setAttribute("verbose", "false");
                    simElem.appendChild(performanceNode);
                }
            }
        }

        return new ElementDocumentPair(simElem, simDoc);
    }

    public ElementDocumentPair saveMetrics(ElementDocumentPair elementDocumentPair) {
        SolverAvgHandles handles = avgHandles;

        ElementDocumentPair res = saveMetric(elementDocumentPair, handles.Q);
        res = saveMetric(res, handles.U);
        res = saveMetric(res, handles.R);
        res = saveMetric(res, handles.T);
        res = saveMetric(res, handles.A);
        if (handles.Tard != null) {
            res = saveMetric(res, handles.Tard);
        }
        if (handles.SysTard != null) {
            res = saveMetric(res, handles.SysTard);
        }
        res = saveCache(res);
        res = saveFCRMetrics(res);

        // JMT ResidT is inconsistently defined with LINE"s on some
        // difficult class switching cases, hence we recompute it at the
        // level of the NetworkSolver class to preserve consistency.
        return res;
    }

    /**
     * Saves performance metrics for Finite Capacity Regions (blocking regions).
     * Requests QLen, RespT, ResidT, and Tput metrics for each FCR.
     */
    public ElementDocumentPair saveFCRMetrics(ElementDocumentPair elementDocumentPair) {
        Document simDoc = elementDocumentPair.simDoc;
        Element simElem = elementDocumentPair.simElem;

        if (simModel.getRegions() == null || simModel.getRegions().isEmpty()) {
            return elementDocumentPair;
        }

        // Metric types to request for FCRs. "FCR Capacity" is JMT's weighted
        // occupation (Total Weight); "FCR Memory" is the memory occupation.
        String[] metricTypes = {"Number of Customers", "Response Time", "Residence Time", "Throughput",
                "FCR Capacity", "FCR Memory"};
        int metricCounter = 0;

        for (int r = 0; r < simModel.getRegions().size(); r++) {
            String fcrName = "FCRegion" + (r + 1);

            for (String metricType : metricTypes) {
                Element performanceNode = simDoc.createElement("measure");
                performanceNode.setAttribute("alpha", String.format("%.2f", 1 - simConfInt));
                performanceNode.setAttribute("name", "FCR_" + fcrName + "_" + metricType.replace(" ", "") + "_" + metricCounter);
                performanceNode.setAttribute("nodeType", "region");
                performanceNode.setAttribute("precision", String.format("%.2f", simMaxRelErr));
                performanceNode.setAttribute("referenceNode", fcrName);
                performanceNode.setAttribute("referenceUserClass", "");  // FCR metrics are not class-specific
                performanceNode.setAttribute("type", metricType);
                performanceNode.setAttribute("verbose", "false");
                simElem.appendChild(performanceNode);
                metricCounter++;
            }
        }

        return new ElementDocumentPair(simElem, simDoc);
    }

    public DocumentSectionPair saveModeNames(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element modeNamesNode = simDoc.createElement("parameter");
        modeNamesNode.setAttribute("classPath", "java.lang.String");
        modeNamesNode.setAttribute("name", "modeNames");
        modeNamesNode.setAttribute("array", "true");

        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;
        for (int m = 0; m < numOfModes; m++) {
            Element subModeNameNode = simDoc.createElement("subParameter");
            subModeNameNode.setAttribute("classPath", "java.lang.String");
            subModeNameNode.setAttribute("name", "modeName");

            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(((TransitionNodeParam) sn.nodeparam.get(istNode)).modenames.get(m)));

            subModeNameNode.appendChild(valueNode);
            modeNamesNode.appendChild(subModeNameNode);
        }
        section.appendChild(modeNamesNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveNumberOfServers(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVENUMBEROFSERVERS(SIMDOC, SECTION, CURRENTNODE)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "maxJobs");

        int istStation = (int) sn.nodeToStation.get(ind);
        int maxJobs;

        // For LPS, use limit from schedparam, otherwise use nservers
        if (sn.sched.get(istStation) == SchedStrategy.LPS) {
            maxJobs = (int) sn.schedparam.get(istStation, 0);  // LPS limit stored in first column
        } else {
            // Use the maximum of nservers and lldscaling (for load-dependent models).
            // Load-dependent scaling with min(1:N,c) represents a c-server queue,
            // where max(lldscaling) = c.
            maxJobs = (int) sn.nservers.get(istStation);
            if (sn.lldscaling != null && sn.lldscaling.getNumRows() > istStation) {
                int effectiveServers = (int) sn.lldscaling.getRow(istStation).elementMax();
                maxJobs = Math.max(maxJobs, effectiveServers);
            }
        }

        Element valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(String.valueOf(maxJobs)));

        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);

        return new DocumentSectionPair(simDoc, section);
    }

    // ==================== Heterogeneous Server Save Methods ====================

    /**
     * Checks if a station has heterogeneous server configuration.
     *
     * @param stationIdx the station index
     * @return true if the station has heterogeneous servers
     */
    public boolean isHeterogeneousStation(int stationIdx) {
        jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(sn.stations.get(stationIdx));
        return snp != null && snp.nservertypes > 0;
    }

    /**
     * Server pools and job parallelism of a node, as the JMT Server section needs them.
     *
     * <p>JMT's Server section takes classParallelism, serverNames,
     * serversPerServerType, serverCompatibilities and schedulingPolicy as one
     * positional block of its constructor ({@code jmt.engine.NodeSections.Server}),
     * so the five are emitted together or not at all, and always after the service
     * strategies. A station declaring parallelism alone is therefore given one
     * synthetic pool holding all of its servers, since the pool counts, not
     * numberOfServers, size the server pool once any pool is declared.</p>
     */
    public static class ServerPools {
        /** Pool names, one per server type. */
        public final java.util.List<String> names;
        /** Servers held by each pool. */
        public final Matrix counts;
        /** Pool-class compatibility, (nTypes x K), 1.0 where the pool serves the class. */
        public final Matrix compat;
        /** Assignment policy across pools. */
        public final HeteroSchedPolicy policy;
        /** Servers seized at once by a job, (1 x K). */
        public final Matrix parallelism;

        ServerPools(java.util.List<String> names, Matrix counts, Matrix compat,
                    HeteroSchedPolicy policy, Matrix parallelism) {
            this.names = names;
            this.counts = counts;
            this.compat = compat;
            this.policy = policy;
            this.parallelism = parallelism;
        }
    }

    /**
     * Builds the server pools of node {@code ind}, or null when the station declares
     * neither server types nor job parallelism.
     *
     * @param ind the node index
     * @return the pools, or null when the Server section needs no pool block
     */
    public ServerPools serverPools(int ind) {
        int istStation = (int) sn.nodeToStation.get(ind);
        if (istStation < 0 || istStation >= sn.stations.size()) {
            return null;
        }
        Station station = sn.stations.get(istStation);
        jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(station);
        if (snp == null) {
            return null;
        }
        int K = sn.nclasses;
        boolean hasTypes = snp.nservertypes > 0 && snp.servertypenames != null
                && !snp.servertypenames.isEmpty() && snp.serverspertype != null
                && snp.servercompat != null;
        boolean hasParallelism = false;
        if (snp.serverparallelism != null) {
            for (int r = 0; r < Math.min(K, snp.serverparallelism.getNumCols()); r++) {
                if (snp.serverparallelism.get(0, r) > 1) {
                    hasParallelism = true;
                    break;
                }
            }
        }
        if (!hasTypes && !hasParallelism) {
            return null;
        }

        Matrix par = new Matrix(1, K, K);
        for (int r = 0; r < K; r++) {
            double n = 1;
            if (snp.serverparallelism != null && r < snp.serverparallelism.getNumCols()) {
                n = Math.max(1, snp.serverparallelism.get(0, r));
            }
            par.set(0, r, n);
        }

        if (hasTypes) {
            HeteroSchedPolicy policy = snp.heteroschedpolicy == null
                    ? HeteroSchedPolicy.ORDER : snp.heteroschedpolicy;
            return new ServerPools(snp.servertypenames, snp.serverspertype, snp.servercompat, policy, par);
        }

        java.util.List<String> names = new java.util.ArrayList<String>();
        names.add(sn.nodenames.get(ind) + " - Server Type 1");
        Matrix counts = new Matrix(1, 1, 1);
        counts.set(0, 0, sn.nservers.get(istStation));
        Matrix compat = new Matrix(1, K, K);
        for (int r = 0; r < K; r++) {
            compat.set(0, r, 1.0);
        }
        return new ServerPools(names, counts, compat, HeteroSchedPolicy.ORDER, par);
    }

    /**
     * Saves the per-class job parallelism to JMT XML.
     * Generates the classParallelism parameter array (Server.serverNumRequired).
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveClassParallelism(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        ServerPools pools = serverPools(ind);
        if (pools == null) {
            return documentSectionPair;
        }

        Element parNode = simDoc.createElement("parameter");
        parNode.setAttribute("classPath", "java.lang.Integer");
        parNode.setAttribute("name", "classParallelism");
        parNode.setAttribute("array", "true");

        for (int r = 0; r < sn.nclasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            parNode.appendChild(refClassNode);

            Element subNode = simDoc.createElement("subParameter");
            subNode.setAttribute("classPath", "java.lang.Integer");
            subNode.setAttribute("name", "serverParallelism");

            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) pools.parallelism.get(0, r))));
            subNode.appendChild(valueNode);
            parNode.appendChild(subNode);
        }

        section.appendChild(parNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves heterogeneous server type names to JMT XML.
     * Generates the serverNames parameter array.
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveServerTypeNames(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        ServerPools pools = serverPools(ind);
        if (pools == null) {
            return documentSectionPair;
        }

        Element serverNamesNode = simDoc.createElement("parameter");
        serverNamesNode.setAttribute("classPath", "java.lang.String");
        serverNamesNode.setAttribute("name", "serverNames");
        serverNamesNode.setAttribute("array", "true");

        for (String name : pools.names) {
            Element subNode = simDoc.createElement("subParameter");
            subNode.setAttribute("classPath", "java.lang.String");
            subNode.setAttribute("name", "serverTypesNames");

            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(name));
            subNode.appendChild(valueNode);
            serverNamesNode.appendChild(subNode);
        }

        section.appendChild(serverNamesNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves the number of servers per server type to JMT XML.
     * Generates the serversPerServerType parameter array.
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveServersPerType(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        ServerPools pools = serverPools(ind);
        if (pools == null) {
            return documentSectionPair;
        }

        Matrix serversPerType = pools.counts;

        Element serversPerTypeNode = simDoc.createElement("parameter");
        serversPerTypeNode.setAttribute("classPath", "java.lang.Integer");
        serversPerTypeNode.setAttribute("name", "serversPerServerType");
        serversPerTypeNode.setAttribute("array", "true");

        for (int t = 0; t < serversPerType.length(); t++) {
            Element subNode = simDoc.createElement("subParameter");
            subNode.setAttribute("classPath", "java.lang.Integer");
            subNode.setAttribute("name", "serverTypesNumOfServers");

            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) serversPerType.get(t))));
            subNode.appendChild(valueNode);
            serversPerTypeNode.appendChild(subNode);
        }

        section.appendChild(serversPerTypeNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves server-class compatibility matrix to JMT XML.
     * Generates the serverCompatibilities parameter array.
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveServerCompatibilities(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        ServerPools pools = serverPools(ind);
        if (pools == null) {
            return documentSectionPair;
        }

        Matrix compat = pools.compat;
        int nTypes = compat.getNumRows();
        int nClasses = compat.getNumCols();

        Element compatNode = simDoc.createElement("parameter");
        compatNode.setAttribute("classPath", "java.lang.Object");
        compatNode.setAttribute("name", "serverCompatibilities");
        compatNode.setAttribute("array", "true");

        // For each server type
        for (int t = 0; t < nTypes; t++) {
            Element typeNode = simDoc.createElement("subParameter");
            typeNode.setAttribute("classPath", "java.lang.Boolean");
            typeNode.setAttribute("name", "serverTypesCompatibilities");
            typeNode.setAttribute("array", "true");

            // For each class
            for (int r = 0; r < nClasses; r++) {
                Element classNode = simDoc.createElement("subParameter");
                classNode.setAttribute("classPath", "java.lang.Boolean");
                classNode.setAttribute("name", "compatibilities");

                Element valueNode = simDoc.createElement("value");
                boolean isCompatible = compat.get(t, r) > 0;
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(isCompatible)));
                classNode.appendChild(valueNode);
                typeNode.appendChild(classNode);
            }
            compatNode.appendChild(typeNode);
        }

        section.appendChild(compatNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Saves heterogeneous scheduling policy to JMT XML.
     * Generates the schedulingPolicy parameter.
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveHeteroSchedPolicy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        ServerPools pools = serverPools(ind);
        if (pools == null) {
            return documentSectionPair;
        }

        Element policyNode = simDoc.createElement("parameter");
        policyNode.setAttribute("classPath", "java.lang.String");
        policyNode.setAttribute("name", "schedulingPolicy");

        Element valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(HeteroSchedPolicy.toJMTText(pools.policy)));
        policyNode.appendChild(valueNode);

        section.appendChild(policyNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Warns that per-server-type service rates cannot reach the JMT engine.
     *
     * <p>JMT keys the ServiceStrategy array of a station by refClass, so its loader
     * ({@code jmt.engine.simEngine.SimLoader}) keeps one strategy per class however
     * many (type, class) entries are written, and every pool of a station ends up
     * serving at the class rate. Pool sizes, class compatibilities and the assignment
     * policy do cross; setService(class, type, distribution) rates do not.</p>
     *
     * @param ind the node index
     */
    public void warnHeteroRates(int ind) {
        int istRates = (int) sn.nodeToStation.get(ind);
        if (istRates < 0 || istRates >= sn.stations.size()) {
            return;
        }
        Station station = sn.stations.get(istRates);
        jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(station);
        if (snp == null || snp.heterorates == null || snp.heterorates.size() < 2) {
            return;
        }
        Double first = null;
        boolean distinct = false;
        for (Map<Integer, Double> typeRates : snp.heterorates.values()) {
            for (Double rate : typeRates.values()) {
                if (rate == null || rate <= 0) {
                    continue;
                }
                if (first == null) {
                    first = rate;
                } else if (Math.abs(rate - first) > 1e-12) {
                    distinct = true;
                }
            }
        }
        if (!distinct) {
            return;
        }
        System.err.format("Warning: JMT keys service strategies by job class, so the per-server-type "
                + "service rates of station %s cannot be exported; every pool will serve at the class "
                + "service rate. Use the LDES or CTMC solver for per-type rates.%n", station.getName());
    }

    /**
     * Saves all heterogeneous server configuration to JMT XML.
     * This is a convenience method that calls all hetero save methods.
     *
     * @param documentSectionPair the document/section pair
     * @param ind the node index
     * @return updated document/section pair
     */
    public DocumentSectionPair saveHeterogeneousServerConfig(DocumentSectionPair documentSectionPair, int ind) {
        documentSectionPair = saveClassParallelism(documentSectionPair, ind);
        documentSectionPair = saveServerTypeNames(documentSectionPair, ind);
        documentSectionPair = saveServersPerType(documentSectionPair, ind);
        documentSectionPair = saveServerCompatibilities(documentSectionPair, ind);
        documentSectionPair = saveHeteroSchedPolicy(documentSectionPair, ind);
        warnHeteroRates(ind);
        return documentSectionPair;
    }

    // ==================== End Heterogeneous Server Methods ====================

    public DocumentSectionPair saveNumbersOfServers(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element numbersOfServersNode = simDoc.createElement("parameter");
        numbersOfServersNode.setAttribute("classPath", "java.lang.Integer");
        numbersOfServersNode.setAttribute("name", "numbersOfServers");
        numbersOfServersNode.setAttribute("array", "true");

        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;
        for (int m = 0; m < numOfModes; m++) {
            Element subNumberOfServersNode = simDoc.createElement("subParameter");
            subNumberOfServersNode.setAttribute("classPath", "java.lang.Integer");
            subNumberOfServersNode.setAttribute("name", "numberOfServers");

            Element valueNode = simDoc.createElement("value");
            double nmodeservers = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodeservers.get(m);

            if (Utils.isInf(nmodeservers)) {
                valueNode.appendChild(simDoc.createTextNode("-1"));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) nmodeservers)));
            }
            subNumberOfServersNode.appendChild(valueNode);
            numbersOfServersNode.appendChild(subNumberOfServersNode);
        }
        section.appendChild(numbersOfServersNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePlaceCapacities(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element placeCapacityNode = simDoc.createElement("parameter");
        placeCapacityNode.setAttribute("array", "true");
        placeCapacityNode.setAttribute("classPath", "java.lang.Integer");
        placeCapacityNode.setAttribute("name", "capacities");
        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();
        double i = sn.nodeToStation.get(ind);
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            placeCapacityNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.Integer");
            subParameterNode.setAttribute("name", "capacity");

            Element valueNode2 = simDoc.createElement("value");
            if (Utils.isInf(sn.cap.get((int) sn.nodeToStation.get(ind)))) {
                valueNode2.appendChild(simDoc.createTextNode("-1"));
            } else if (Utils.isInf(sn.classcap.get((int) i, r))) {
                valueNode2.appendChild(simDoc.createTextNode("-1"));
            } else {
                valueNode2.appendChild(simDoc.createTextNode(String.valueOf((int) sn.classcap.get((int) i, r))));
            }
            subParameterNode.appendChild(valueNode2);
            placeCapacityNode.appendChild(subParameterNode);
        }
        section.appendChild(placeCapacityNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePreemptiveStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element visitsNode = simDoc.createElement("parameter");
        visitsNode.setAttribute("array", "true");
        visitsNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategy");
        visitsNode.setAttribute("name", "PSStrategy");

        int numOfClasses = sn.nclasses;
        Station istStation = sn.stations.get((int) sn.nodeToStation.get(ind));
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            visitsNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            switch (sn.sched.get(istStation)) {
                case PS:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.EPSStrategy");
                    subParameterNode.setAttribute("name", "EPSStrategy");
                    break;
                case DPS:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.DPSStrategy");
                    subParameterNode.setAttribute("name", "DPSStrategy");
                    break;
                case GPS:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.GPSStrategy");
                    subParameterNode.setAttribute("name", "GPSStrategy");
                    break;
                case LPS:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.EPSStrategy");
                    subParameterNode.setAttribute("name", "EPSStrategy");
                    break;
                case PSPRIO:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.EPSStrategyPriority");
                    subParameterNode.setAttribute("name", "EPSStrategyPriority");
                    break;
                case DPSPRIO:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.DPSStrategyPriority");
                    subParameterNode.setAttribute("name", "DPSStrategyPriority");
                    break;
                case GPSPRIO:
                    subParameterNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategies.GPSStrategyPriority");
                    subParameterNode.setAttribute("name", "GPSStrategyPriority");
                    break;
            }
            visitsNode.appendChild(subParameterNode);
            section.appendChild(visitsNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePreemptiveWeights(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEPREEMPTIVEWEIGHTS(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element visitsNode = simDoc.createElement("parameter");
        visitsNode.setAttribute("array", "true");
        visitsNode.setAttribute("classPath", "java.lang.Double");
        visitsNode.setAttribute("name", "serviceWeights");

        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        boolean[] exportClasses = getExportableClasses();
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            visitsNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.Double");
            subParameterNode.setAttribute("name", "serviceWeight");

            Element valueNode2 = simDoc.createElement("value");
            valueNode2.appendChild(simDoc.createTextNode(String.valueOf(sn.schedparam.get(i, r))));

            subParameterNode.appendChild(valueNode2);
            visitsNode.appendChild(subParameterNode);
            section.appendChild(visitsNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePutStrategies(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element queuePutStrategyNode = simDoc.createElement("parameter");
        queuePutStrategyNode.setAttribute("array", "true");
        queuePutStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategy");
        queuePutStrategyNode.setAttribute("name", "QueuePutStrategy");

        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        Station istStation = sn.stations.get(i);
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            queuePutStrategyNode.appendChild(refClassNode2);
            Element subParameterNode2 = simDoc.createElement("subParameter");
            switch (sn.sched.get(istStation)) {
                case SIRO:
                    subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy");
                    subParameterNode2.setAttribute("name", "RandStrategy");
                    break;
                case LCFS:
                    subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy");
                    subParameterNode2.setAttribute("name", "HeadStrategy");
                    break;
                default:
                    subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy");
                    subParameterNode2.setAttribute("name", "TailStrategy");
                    break;
            }
            queuePutStrategyNode.appendChild(subParameterNode2);
            section.appendChild(queuePutStrategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePutStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEPUTSTRATEGY(SIMDOC, SECTION, CURRENTNODE)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element queuePutStrategyNode = simDoc.createElement("parameter");
        queuePutStrategyNode.setAttribute("array", "true");
        queuePutStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategy");
        queuePutStrategyNode.setAttribute("name", "QueuePutStrategy");

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();
        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            queuePutStrategyNode.appendChild(refClassNode2);

            Element subParameterNode2 = simDoc.createElement("subParameter");
            // if not a station treat as FCFS
            if (sn.isstation.get(ind, 0) != 1) {
                subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy");
                subParameterNode2.setAttribute("name", "TailStrategy");
            } else { // if a station
                switch (sn.sched.get(sn.stations.get((int) sn.nodeToStation.get(ind)))) {
                    case SIRO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy");
                        subParameterNode2.setAttribute("name", "RandStrategy");
                        break;
                    case LJF:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LJFStrategy");
                        subParameterNode2.setAttribute("name", "LJFStrategy");
                        break;
                    case SJF:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.SJFStrategy");
                        subParameterNode2.setAttribute("name", "SJFStrategy");
                        break;
                    case LEPT:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LEPTStrategy");
                        subParameterNode2.setAttribute("name", "LEPTStrategy");
                        break;
                    case SEPT:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.SEPTStrategy");
                        subParameterNode2.setAttribute("name", "SEPTStrategy");
                        break;
                    case SRPT:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.SRPTStrategy");
                        subParameterNode2.setAttribute("name", "SRPTStrategy");
                        break;
                    case SRPTPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.SRPTStrategyPriority");
                        subParameterNode2.setAttribute("name", "SRPTStrategyPriority");
                        break;
                    case LCFS:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy");
                        subParameterNode2.setAttribute("name", "HeadStrategy");
                        break;
                    case LCFSPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategyPriority");
                        subParameterNode2.setAttribute("name", "HeadStrategyPriority");
                        break;
                    case LCFSPR:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategy");
                        subParameterNode2.setAttribute("name", "LCFSPRStrategy");
                        break;
                    case LCFSPI:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LCFSPIStrategy");
                        subParameterNode2.setAttribute("name", "LCFSPIStrategy");
                        break;
                    case FCFSPR:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.FCFSPRStrategy");
                        subParameterNode2.setAttribute("name", "FCFSPRStrategy");
                        break;
                    case FCFSPI:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.FCFSPIStrategy");
                        subParameterNode2.setAttribute("name", "FCFSPIStrategy");
                        break;
                    case LCFSPRPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategyPriority");
                        subParameterNode2.setAttribute("name", "LCFSPRStrategyPriority");
                        break;
                    case LCFSPIPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LCFSPIStrategyPriority");
                        subParameterNode2.setAttribute("name", "LCFSPIStrategyPriority");
                        break;
                    case FCFSPRPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.FCFSPRStrategyPriority");
                        subParameterNode2.setAttribute("name", "FCFSPRStrategyPriority");
                        break;
                    case FCFSPIPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.FCFSPIStrategyPriority");
                        subParameterNode2.setAttribute("name", "FCFSPIStrategyPriority");
                        break;
                    case HOL:
                    case FCFSPRIO:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority");
                        subParameterNode2.setAttribute("name", "TailStrategyPriority");
                        break;
                    case EDD:
                        // Note: JMT does not natively support EDD yet. This generates XML for future compatibility.
                        // For now, this will likely fall back to FCFS behavior in JMT execution.
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.EDDStrategy");
                        subParameterNode2.setAttribute("name", "EDDStrategy");
                        break;
                    case EDF:
                        // Note: JMT does not natively support EDF yet. This generates XML for future compatibility.
                        // For now, this will likely fall back to FCFS behavior in JMT execution.
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.EDFStrategy");
                        subParameterNode2.setAttribute("name", "EDFStrategy");
                        break;
                    case LPS:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LPSStrategy");
                        subParameterNode2.setAttribute("name", "LPSStrategy");
                        break;
                    default: // treat as FCFS -this is required for PS
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy");
                        subParameterNode2.setAttribute("name", "TailStrategy");
                        break;
                }
            }
            queuePutStrategyNode.appendChild(subParameterNode2);
            section.appendChild(queuePutStrategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveImpatience(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEIMPATIENCE(SIMDOC, SECTION, IND)
        // Generates XML for impatience (reneging) distributions for JMT Queue sections
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element impatienceNode = simDoc.createElement("parameter");
        impatienceNode.setAttribute("array", "true");
        impatienceNode.setAttribute("classPath", "jmt.engine.NetStrategies.ImpatienceStrategies.Impatience");
        impatienceNode.setAttribute("name", "Impatience");

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();
        int ist = (int) sn.nodeToStation.get(ind);

        // Check if this node has a valid station mapping
        // Cache nodes and other special nodes may not have stations (ist = -1)
        if (ist < 0) {
            // For non-station nodes, don't generate impatience parameters at all
            // JMT expects the array to have entries for all classes, so an empty array causes errors
            return new DocumentSectionPair(simDoc, section);
        }

        Station istStation = sn.stations.get(ist);

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            impatienceNode.appendChild(refClassNode);

            JobClass rstJobClass = sn.jobclasses.get(r);

            // Check if impatience is defined for this station-class pair
            boolean hasImpatience = false;
            ProcessType procType = null;
            if (sn.impatienceType != null && sn.impatienceType.containsKey(istStation)) {
                Map<JobClass, ProcessType> classMap = sn.impatienceType.get(istStation);
                if (classMap != null && classMap.containsKey(rstJobClass)) {
                    hasImpatience = true;
                    procType = classMap.get(rstJobClass);
                }
            }

            Element impatienceStrategyNode = simDoc.createElement("subParameter");
            if (!hasImpatience) {
                // No impatience defined - use null with Impatience interface
                impatienceStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ImpatienceStrategies.Impatience");
                impatienceStrategyNode.setAttribute("name", "Impatience");
                Element subParValue = simDoc.createElement("value");
                subParValue.appendChild(simDoc.createTextNode("null"));
                impatienceStrategyNode.appendChild(subParValue);
            } else {
                // Reneging impatience - use Reneging class
                impatienceStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ImpatienceStrategies.Reneging");
                impatienceStrategyNode.setAttribute("name", "Reneging");
                // Impatience is defined - generate distribution XML
                Element distributionNode = simDoc.createElement("subParameter");
                String javaClass = "";
                String javaParClass = "";

                int impatiencePhaseCount = sn.impatiencePhases.get(istStation).get(rstJobClass);
                // jmt.engine.random.HyperExpPar carries exactly (p, lambda1,
                // lambda2), so it can only express a 2-phase hyper-exponential.
                // An n>2 HyperExp is exported as the acyclic phase-type built
                // from its OWN PH representation (alpha = p, T = -diag(lambda)),
                // i.e. the same distribution rather than a moment-matched refit.
                // This mirrors the service path in saveServiceStrategy.
                boolean emitPhaseType = procType == ProcessType.PH
                        || procType == ProcessType.APH
                        || procType == ProcessType.COXIAN
                        || (procType == ProcessType.HYPEREXP && impatiencePhaseCount > 2);

                if (emitPhaseType) {
                    javaClass = "jmt.engine.random.PhaseTypeDistr";
                    javaParClass = "jmt.engine.random.PhaseTypePar";
                } else {
                    switch (procType) {
                        case DET:
                            javaClass = "jmt.engine.random.DeterministicDistr";
                            javaParClass = "jmt.engine.random.DeterministicDistrPar";
                            break;
                        case ERLANG:
                            javaClass = "jmt.engine.random.Erlang";
                            javaParClass = "jmt.engine.random.ErlangPar";
                            break;
                        case EXP:
                            javaClass = "jmt.engine.random.Exponential";
                            javaParClass = "jmt.engine.random.ExponentialPar";
                            break;
                        case GAMMA:
                            javaClass = "jmt.engine.random.GammaDistr";
                            javaParClass = "jmt.engine.random.GammaDistrPar";
                            break;
                        case HYPEREXP:
                            javaClass = "jmt.engine.random.HyperExp";
                            javaParClass = "jmt.engine.random.HyperExpPar";
                            break;
                        case PARETO:
                            javaClass = "jmt.engine.random.Pareto";
                            javaParClass = "jmt.engine.random.ParetoPar";
                            break;
                        case WEIBULL:
                            javaClass = "jmt.engine.random.Weibull";
                            javaParClass = "jmt.engine.random.WeibullPar";
                            break;
                        case LOGNORMAL:
                            javaClass = "jmt.engine.random.Lognormal";
                            javaParClass = "jmt.engine.random.LognormalPar";
                            break;
                        case UNIFORM:
                            javaClass = "jmt.engine.random.Uniform";
                            javaParClass = "jmt.engine.random.UniformPar";
                            break;
                        default:
                            throw new RuntimeException("Unsupported impatience distribution type: " + procType);
                    }
                }

                distributionNode.setAttribute("classPath", javaClass);
                if (emitPhaseType) {
                    distributionNode.setAttribute("name", "Phase-Type");
                } else {
                    switch (procType) {
                        case EXP:
                            distributionNode.setAttribute("name", "Exponential");
                            break;
                        case HYPEREXP:
                            distributionNode.setAttribute("name", "Hyperexponential");
                            break;
                        default:
                            distributionNode.setAttribute("name", procType.toString());
                            break;
                    }
                }
                impatienceStrategyNode.appendChild(distributionNode);

                // Create distribution parameters
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");

                // Get impatience parameters
                Matrix impatienceMu = sn.impatienceMu.get(istStation).get(rstJobClass);
                Matrix impatiencePhi = sn.impatiencePhi.get(istStation).get(rstJobClass);

                if (emitPhaseType) {
                    MatrixCell impatienceProc = sn.impatienceProc.get(istStation).get(rstJobClass);
                    Matrix impatiencePie = sn.impatiencePie.get(istStation).get(rstJobClass);
                    int phases = impatiencePhaseCount;
                    Matrix PH = impatienceProc.get(0);

                    // Alpha vector
                    Element subParNodeAlpha = simDoc.createElement("subParameter");
                    subParNodeAlpha.setAttribute("array", "true");
                    subParNodeAlpha.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlpha.setAttribute("name", "alpha");

                    Element subParNodeAlphaVec = simDoc.createElement("subParameter");
                    subParNodeAlphaVec.setAttribute("array", "true");
                    subParNodeAlphaVec.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlphaVec.setAttribute("name", "vector");

                    for (int k = 0; k < phases; k++) {
                        Element subParNodeAlphaElem = simDoc.createElement("subParameter");
                        subParNodeAlphaElem.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlphaElem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(alphaEntry(impatiencePie, k)))));
                        subParNodeAlphaElem.appendChild(subParValue);
                        subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
                    }
                    subParNodeAlpha.appendChild(subParNodeAlphaVec);

                    // T matrix
                    Element subParNodeT = simDoc.createElement("subParameter");
                    subParNodeT.setAttribute("array", "true");
                    subParNodeT.setAttribute("classPath", "java.lang.Object");
                    subParNodeT.setAttribute("name", "T");

                    for (int k = 0; k < phases; k++) {
                        Element subParNodeTvec = simDoc.createElement("subParameter");
                        subParNodeTvec.setAttribute("array", "true");
                        subParNodeTvec.setAttribute("classPath", "java.lang.Object");
                        subParNodeTvec.setAttribute("name", "vector");

                        for (int j = 0; j < phases; j++) {
                            Element subParNodeTElem = simDoc.createElement("subParameter");
                            subParNodeTElem.setAttribute("classPath", "java.lang.Double");
                            subParNodeTElem.setAttribute("name", "entry");
                            Element subParValue = simDoc.createElement("value");
                            if (k == j) {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -Math.abs(PH.get(k, j)))));
                            } else {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(PH.get(k, j)))));
                            }
                            subParNodeTElem.appendChild(subParValue);
                            subParNodeTvec.appendChild(subParNodeTElem);
                        }
                        subParNodeT.appendChild(subParNodeTvec);
                    }

                    distrParNode.appendChild(subParNodeAlpha);
                    distrParNode.appendChild(subParNodeT);
                } else {
                switch (procType) {
                    case DET: {
                        Element subParNodeT = simDoc.createElement("subParameter");
                        subParNodeT.setAttribute("classPath", "java.lang.Double");
                        subParNodeT.setAttribute("name", "t");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / impatienceMu.get(0, 0))));
                        subParNodeT.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeT);
                        break;
                    }
                    case EXP: {
                        Element subParNodeLambda = simDoc.createElement("subParameter");
                        subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda.setAttribute("name", "lambda");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", impatienceMu.get(0, 0))));
                        subParNodeLambda.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeLambda);
                        break;
                    }
                    case ERLANG: {
                        int phases = sn.impatiencePhases.get(istStation).get(rstJobClass);
                        Element subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", impatienceMu.get(0, 0) * phases)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);

                        Element subParNodeR = simDoc.createElement("subParameter");
                        subParNodeR.setAttribute("classPath", "java.lang.Long");
                        subParNodeR.setAttribute("name", "r");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.valueOf(phases)));
                        subParNodeR.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeR);
                        break;
                    }
                    case HYPEREXP: {
                        MatrixCell impatienceProc = sn.impatienceProc.get(istStation).get(rstJobClass);
                        Matrix impatiencePie = sn.impatiencePie.get(istStation).get(rstJobClass);

                        Element subParNodeP = simDoc.createElement("subParameter");
                        subParNodeP.setAttribute("classPath", "java.lang.Double");
                        subParNodeP.setAttribute("name", "p");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", impatiencePie.get(0, 0))));
                        subParNodeP.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeP);

                        Element subParNodeLambda1 = simDoc.createElement("subParameter");
                        subParNodeLambda1.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda1.setAttribute("name", "lambda1");
                        Element subParValue1 = simDoc.createElement("value");
                        subParValue1.appendChild(simDoc.createTextNode(String.format("%.12f", -impatienceProc.get(0).get(0, 0))));
                        subParNodeLambda1.appendChild(subParValue1);
                        distrParNode.appendChild(subParNodeLambda1);

                        Element subParNodeLambda2 = simDoc.createElement("subParameter");
                        subParNodeLambda2.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda2.setAttribute("name", "lambda2");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", -impatienceProc.get(0).get(1, 1))));
                        subParNodeLambda2.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeLambda2);
                        break;
                    }
                    case GAMMA: {
                        double scv = impatiencePhi.get(0, 0);
                        Element subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / scv)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);

                        Element subParNodeBeta = simDoc.createElement("subParameter");
                        subParNodeBeta.setAttribute("classPath", "java.lang.Double");
                        subParNodeBeta.setAttribute("name", "beta");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", scv / impatienceMu.get(0, 0))));
                        subParNodeBeta.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeBeta);
                        break;
                    }
                    case PARETO: {
                        double scv = impatiencePhi.get(0, 0);
                        double shape = Math.sqrt(1.0 + 1.0 / scv) + 1.0;
                        double scale = 1.0 / impatienceMu.get(0, 0) * (shape - 1.0) / shape;

                        Element subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", shape)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);

                        Element subParNodeK = simDoc.createElement("subParameter");
                        subParNodeK.setAttribute("classPath", "java.lang.Double");
                        subParNodeK.setAttribute("name", "k");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", scale)));
                        subParNodeK.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeK);
                        break;
                    }
                    case WEIBULL: {
                        double scv = impatiencePhi.get(0, 0);
                        double c = Math.sqrt(scv);
                        double rval = Math.pow(c, -1.086); // Justus approximation (1976)
                        double alpha = 1.0 / impatienceMu.get(0, 0) / gamma(1.0 + 1.0 / rval);

                        Element subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);

                        Element subParNodeR = simDoc.createElement("subParameter");
                        subParNodeR.setAttribute("classPath", "java.lang.Double");
                        subParNodeR.setAttribute("name", "r");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", rval)));
                        subParNodeR.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeR);
                        break;
                    }
                    case LOGNORMAL: {
                        double scv = impatiencePhi.get(0, 0);
                        double c = Math.sqrt(scv);
                        double mu = Math.log(1.0 / impatienceMu.get(0, 0) / Math.sqrt(c * c + 1.0));
                        double sigma = Math.sqrt(Math.log(c * c + 1.0));

                        Element subParNodeMu = simDoc.createElement("subParameter");
                        subParNodeMu.setAttribute("classPath", "java.lang.Double");
                        subParNodeMu.setAttribute("name", "mu");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", mu)));
                        subParNodeMu.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeMu);

                        Element subParNodeSigma = simDoc.createElement("subParameter");
                        subParNodeSigma.setAttribute("classPath", "java.lang.Double");
                        subParNodeSigma.setAttribute("name", "sigma");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", sigma)));
                        subParNodeSigma.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeSigma);
                        break;
                    }
                    case UNIFORM: {
                        double mean = 1.0 / impatienceMu.get(0, 0);
                        double b = 2.0 * mean; // Uniform [0, b] with mean = b/2

                        Element subParNodeMin = simDoc.createElement("subParameter");
                        subParNodeMin.setAttribute("classPath", "java.lang.Double");
                        subParNodeMin.setAttribute("name", "min");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode("0.0"));
                        subParNodeMin.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeMin);

                        Element subParNodeMax = simDoc.createElement("subParameter");
                        subParNodeMax.setAttribute("classPath", "java.lang.Double");
                        subParNodeMax.setAttribute("name", "max");
                        Element subParValue2 = simDoc.createElement("value");
                        subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", b)));
                        subParNodeMax.appendChild(subParValue2);
                        distrParNode.appendChild(subParNodeMax);
                        break;
                    }
                }
                }

                impatienceStrategyNode.appendChild(distrParNode);
            }

            impatienceNode.appendChild(impatienceStrategyNode);
        }

        section.appendChild(impatienceNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveRetrialDistributions(DocumentSectionPair documentSectionPair, int ind) {
        // Generates XML for retrial delay distributions for JMT Queue sections
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element retrialNode = simDoc.createElement("parameter");
        retrialNode.setAttribute("array", "true");
        retrialNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        retrialNode.setAttribute("name", "retrialDistributions");

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();
        int ist = (int) sn.nodeToStation.get(ind);

        if (ist < 0) {
            return new DocumentSectionPair(simDoc, section);
        }

        Station istStation = sn.stations.get(ist);

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            retrialNode.appendChild(refClassNode);

            JobClass rstJobClass = sn.jobclasses.get(r);

            // Check if retrial is defined for this station-class pair
            boolean hasRetrial = false;
            ProcessType procType = null;
            if (sn.retrialType != null && sn.retrialType.containsKey(istStation)) {
                Map<JobClass, ProcessType> classMap = sn.retrialType.get(istStation);
                if (classMap != null && classMap.containsKey(rstJobClass)) {
                    procType = classMap.get(rstJobClass);
                    if (procType != null && procType != ProcessType.DISABLED) {
                        hasRetrial = true;
                    }
                }
            }

            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            if (!hasRetrial) {
                // No retrial - write null ServiceTimeStrategy
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");
                Element nullDistrNode = simDoc.createElement("subParameter");
                nullDistrNode.setAttribute("classPath", "jmt.engine.random.Exponential");
                nullDistrNode.setAttribute("name", "Exponential");
                serviceTimeStrategyNode.appendChild(nullDistrNode);
                Element nullDistrParNode = simDoc.createElement("subParameter");
                nullDistrParNode.setAttribute("classPath", "jmt.engine.random.ExponentialPar");
                nullDistrParNode.setAttribute("name", "distrPar");
                Element subParNodeLambda = simDoc.createElement("subParameter");
                subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                subParNodeLambda.setAttribute("name", "lambda");
                Element subParValue = simDoc.createElement("value");
                subParValue.appendChild(simDoc.createTextNode("1.000000000000"));
                subParNodeLambda.appendChild(subParValue);
                nullDistrParNode.appendChild(subParNodeLambda);
                serviceTimeStrategyNode.appendChild(nullDistrParNode);
            } else {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

                Element distributionNode = simDoc.createElement("subParameter");
                String javaClass = "";
                String javaParClass = "";

                // The (D0,D1) representation of the retrial delay, when the delay
                // is Markovian. It carries the phase count and, for a phase-type
                // renewal process, the entry vector alpha.
                MatrixCell retrialProc = null;
                if (sn.retrialProc != null && sn.retrialProc.containsKey(istStation)) {
                    retrialProc = sn.retrialProc.get(istStation).get(rstJobClass);
                }
                int retrialPhaseCount = 1;
                if (retrialProc != null && retrialProc.get(0) != null) {
                    retrialPhaseCount = retrialProc.get(0).getNumCols();
                }

                // jmt.engine.random.HyperExpPar carries exactly (p, lambda1,
                // lambda2), so it can only express a 2-phase hyper-exponential.
                // An n>2 HyperExp, and any other phase-type retrial delay, is
                // exported as the acyclic phase-type built from its OWN (D0,D1),
                // i.e. the same distribution rather than a moment-matched refit.
                // This mirrors the service path in saveServiceStrategy.
                boolean emitPhaseType = procType == ProcessType.PH
                        || procType == ProcessType.APH
                        || procType == ProcessType.COXIAN
                        || (procType == ProcessType.HYPEREXP && retrialPhaseCount > 2);
                if (emitPhaseType && (retrialProc == null || retrialProc.get(0) == null || retrialProc.get(1) == null)) {
                    throw new RuntimeException("Retrial delay of type " + procType + " at station "
                            + istStation.getName() + " has no (D0,D1) representation, so it cannot be"
                            + " exported to JMT as a phase-type.");
                }

                if (emitPhaseType) {
                    javaClass = "jmt.engine.random.PhaseTypeDistr";
                    javaParClass = "jmt.engine.random.PhaseTypePar";
                } else {
                    switch (procType) {
                        case DET:
                            javaClass = "jmt.engine.random.DeterministicDistr";
                            javaParClass = "jmt.engine.random.DeterministicDistrPar";
                            break;
                        case ERLANG:
                            javaClass = "jmt.engine.random.Erlang";
                            javaParClass = "jmt.engine.random.ErlangPar";
                            break;
                        case EXP:
                            javaClass = "jmt.engine.random.Exponential";
                            javaParClass = "jmt.engine.random.ExponentialPar";
                            break;
                        case GAMMA:
                            javaClass = "jmt.engine.random.GammaDistr";
                            javaParClass = "jmt.engine.random.GammaDistrPar";
                            break;
                        case HYPEREXP:
                            javaClass = "jmt.engine.random.HyperExp";
                            javaParClass = "jmt.engine.random.HyperExpPar";
                            break;
                        case PARETO:
                            javaClass = "jmt.engine.random.Pareto";
                            javaParClass = "jmt.engine.random.ParetoPar";
                            break;
                        case WEIBULL:
                            javaClass = "jmt.engine.random.Weibull";
                            javaParClass = "jmt.engine.random.WeibullPar";
                            break;
                        case LOGNORMAL:
                            javaClass = "jmt.engine.random.Lognormal";
                            javaParClass = "jmt.engine.random.LognormalPar";
                            break;
                        case UNIFORM:
                            javaClass = "jmt.engine.random.Uniform";
                            javaParClass = "jmt.engine.random.UniformPar";
                            break;
                        default:
                            javaClass = "jmt.engine.random.Exponential";
                            javaParClass = "jmt.engine.random.ExponentialPar";
                            break;
                    }
                }

                distributionNode.setAttribute("classPath", javaClass);
                if (emitPhaseType) {
                    distributionNode.setAttribute("name", "Phase-Type");
                } else {
                    switch (procType) {
                        case EXP:
                            distributionNode.setAttribute("name", "Exponential");
                            break;
                        case HYPEREXP:
                            distributionNode.setAttribute("name", "Hyperexponential");
                            break;
                        default:
                            distributionNode.setAttribute("name", procType.toString());
                            break;
                    }
                }
                serviceTimeStrategyNode.appendChild(distributionNode);

                // Create distribution parameters
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");

                Matrix retrialMu = sn.retrialMu.get(istStation).get(rstJobClass);

                if (emitPhaseType) {
                    Matrix PH = retrialProc.get(0);
                    Matrix pie = Map_pie.map_pie(retrialProc.get(0), retrialProc.get(1));
                    int phases = retrialPhaseCount;

                    // Alpha vector
                    Element subParNodeAlpha = simDoc.createElement("subParameter");
                    subParNodeAlpha.setAttribute("array", "true");
                    subParNodeAlpha.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlpha.setAttribute("name", "alpha");

                    Element subParNodeAlphaVec = simDoc.createElement("subParameter");
                    subParNodeAlphaVec.setAttribute("array", "true");
                    subParNodeAlphaVec.setAttribute("classPath", "java.lang.Object");
                    subParNodeAlphaVec.setAttribute("name", "vector");

                    for (int k = 0; k < phases; k++) {
                        Element subParNodeAlphaElem = simDoc.createElement("subParameter");
                        subParNodeAlphaElem.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlphaElem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(alphaEntry(pie, k)))));
                        subParNodeAlphaElem.appendChild(subParValue);
                        subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
                    }
                    subParNodeAlpha.appendChild(subParNodeAlphaVec);

                    // T matrix
                    Element subParNodeT = simDoc.createElement("subParameter");
                    subParNodeT.setAttribute("array", "true");
                    subParNodeT.setAttribute("classPath", "java.lang.Object");
                    subParNodeT.setAttribute("name", "T");

                    for (int k = 0; k < phases; k++) {
                        Element subParNodeTvec = simDoc.createElement("subParameter");
                        subParNodeTvec.setAttribute("array", "true");
                        subParNodeTvec.setAttribute("classPath", "java.lang.Object");
                        subParNodeTvec.setAttribute("name", "vector");

                        for (int j = 0; j < phases; j++) {
                            Element subParNodeTElem = simDoc.createElement("subParameter");
                            subParNodeTElem.setAttribute("classPath", "java.lang.Double");
                            subParNodeTElem.setAttribute("name", "entry");
                            Element subParValue = simDoc.createElement("value");
                            if (k == j) {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -Math.abs(PH.get(k, j)))));
                            } else {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(PH.get(k, j)))));
                            }
                            subParNodeTElem.appendChild(subParValue);
                            subParNodeTvec.appendChild(subParNodeTElem);
                        }
                        subParNodeT.appendChild(subParNodeTvec);
                    }

                    distrParNode.appendChild(subParNodeAlpha);
                    distrParNode.appendChild(subParNodeT);
                } else if (procType == ProcessType.HYPEREXP) {
                    // 2-phase only: HyperExpPar carries exactly (p, lambda1,
                    // lambda2). Emitting a bare `lambda` here, as the generic
                    // fallback below did, silently degraded the retrial delay
                    // to an exponential of the same mean.
                    Matrix PH = retrialProc.get(0);
                    Matrix pie = Map_pie.map_pie(retrialProc.get(0), retrialProc.get(1));

                    Element subParNodeP = simDoc.createElement("subParameter");
                    subParNodeP.setAttribute("classPath", "java.lang.Double");
                    subParNodeP.setAttribute("name", "p");
                    Element subParValueP = simDoc.createElement("value");
                    subParValueP.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(alphaEntry(pie, 0)))));
                    subParNodeP.appendChild(subParValueP);
                    distrParNode.appendChild(subParNodeP);

                    Element subParNodeLambda1 = simDoc.createElement("subParameter");
                    subParNodeLambda1.setAttribute("classPath", "java.lang.Double");
                    subParNodeLambda1.setAttribute("name", "lambda1");
                    Element subParValue1 = simDoc.createElement("value");
                    subParValue1.appendChild(simDoc.createTextNode(String.format("%.12f", -PH.get(0, 0))));
                    subParNodeLambda1.appendChild(subParValue1);
                    distrParNode.appendChild(subParNodeLambda1);

                    Element subParNodeLambda2 = simDoc.createElement("subParameter");
                    subParNodeLambda2.setAttribute("classPath", "java.lang.Double");
                    subParNodeLambda2.setAttribute("name", "lambda2");
                    Element subParValue2 = simDoc.createElement("value");
                    subParValue2.appendChild(simDoc.createTextNode(String.format("%.12f", -PH.get(1, 1))));
                    subParNodeLambda2.appendChild(subParValue2);
                    distrParNode.appendChild(subParNodeLambda2);
                } else {
                switch (procType) {
                    case DET: {
                        Element subParNodeT = simDoc.createElement("subParameter");
                        subParNodeT.setAttribute("classPath", "java.lang.Double");
                        subParNodeT.setAttribute("name", "t");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / retrialMu.get(0, 0))));
                        subParNodeT.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeT);
                        break;
                    }
                    case EXP: {
                        Element subParNodeLambda = simDoc.createElement("subParameter");
                        subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda.setAttribute("name", "lambda");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", retrialMu.get(0, 0))));
                        subParNodeLambda.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeLambda);
                        break;
                    }
                    default: {
                        // For other distributions, use rate as exponential fallback
                        Element subParNodeLambda = simDoc.createElement("subParameter");
                        subParNodeLambda.setAttribute("classPath", "java.lang.Double");
                        subParNodeLambda.setAttribute("name", "lambda");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", retrialMu.get(0, 0))));
                        subParNodeLambda.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeLambda);
                        break;
                    }
                }
                }

                serviceTimeStrategyNode.appendChild(distrParNode);
            }

            retrialNode.appendChild(serviceTimeStrategyNode);
        }

        section.appendChild(retrialNode);
        return new DocumentSectionPair(simDoc, section);
    }

    /**
     * Reads the k-th entry of a phase-type entry-probability vector, which the
     * refresh routines store as a row vector (Map_pie returns 1-by-n) while the
     * hand-built fallbacks store a 1-by-1. Indexing it as a column silently read
     * out of range for k > 0, so the orientation is resolved explicitly here.
     *
     * @param pie the entry-probability vector, as a row or column Matrix
     * @param k   the zero-based phase index
     * @return the entry probability of phase k
     */
    private static double alphaEntry(Matrix pie, int k) {
        if (pie.getNumRows() == 1) {
            return pie.get(0, k);
        }
        return pie.get(k, 0);
    }

    private double gamma(double x) {
        // Simplified gamma function approximation using Lanczos approximation
        return org.apache.commons.math3.special.Gamma.gamma(x);
    }

    /**
     * Per-class station capacities that must be exported as JMT blocking-region
     * class constraints. Entry (i,r) is finite only when sn.classcap(i,r) is a
     * genuine buffer limit that the rest of the exported model does not already
     * imply; it is Inf otherwise, so that models without a per-class capacity
     * emit byte-identical XML.
     *
     * refreshCapacity derives classcap(i,r) as
     * min(chain population of r, station.classCap(r), station.cap), so a value
     * equal to the chain population or to the station total is not a per-class
     * buffer at all: the closed population, respectively the Queue "size"
     * written by saveBufferCapacity, already enforces it. Only a strictly
     * tighter value carries information that JMT would otherwise lose.
     *
     * @return an (nstations,nclasses) Matrix, Inf where no constraint is needed
     */
    private Matrix classCapConstraints() {
        Matrix classCapCon = new Matrix(sn.nstations, sn.nclasses);
        classCapCon.fill(Inf);
        if (sn.classcap == null || sn.classcap.isEmpty()) {
            return classCapCon;
        }

        // Population bound implied by the closed classes of each chain (Inf when open)
        Matrix chainpop = new Matrix(1, sn.nclasses);
        chainpop.fill(Inf);
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            double chainCap = 0;
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                chainCap += sn.njobs.get(0, (int) inchain_c.get(0, idx));
            }
            if (chainCap >= Integer.MAX_VALUE) {
                chainCap = Inf;
            }
            for (int idx = 0; idx < inchain_c.length(); idx++) {
                chainpop.set(0, (int) inchain_c.get(0, idx), chainCap);
            }
        }

        for (int i = 0; i < sn.nstations; i++) {
            int ind = (int) sn.stationToNode.get(i);
            // A Source has no buffer, and a Place carries its per-class capacities
            // in the JMT Storage section, which does have a capacities array.
            if (sn.nodetype.get(ind) == NodeType.Source || sn.nodetype.get(ind) == NodeType.Place) {
                continue;
            }
            for (int r = 0; r < sn.nclasses; r++) {
                // A class disabled at the station gets classcap 0 from
                // refreshCapacity; it is not routed there, and a zero constraint
                // is not a buffer limit.
                if (Double.isNaN(sn.rates.get(i, r))) {
                    continue;
                }
                // Utils.isInf also reads Integer.MAX_VALUE as infinite, which is
                // the "unset capacity" sentinel refreshCapacity leaves behind.
                double ccap = sn.classcap.get(i, r);
                if (!Utils.isInf(ccap) && ccap < FastMath.min(sn.cap.get(i, 0), chainpop.get(0, r))) {
                    classCapCon.set(i, r, ccap);
                }
            }
        }
        return classCapCon;
    }

    /**
     * Asserts that LINE's sn.classcap semantics for class r at station i is
     * expressible as a JMT blocking-region constraint, throwing a descriptive
     * exception when it is not. LINE enforces classcap in
     * State.afterEventStation by not enabling the arrival; sn.droprule is read
     * there only to select BAS/BBS/RSRD blocking, so WaitingQueue and Drop
     * behave identically for a capacity limit. The resulting reference
     * semantics, verified against SolverCTMC, is
     * <ul>
     * <li>open class: the arrival is lost and the upstream is unaffected</li>
     * <li>closed class: the job stays at the upstream station and the
     * population is conserved</li>
     * </ul>
     * which is the same predicate the CTMC analyzer applies in canDropClass.
     * Only the open case has a faithful counterpart, dropThisClass="true": JMT
     * discards the arrival at the region input station, exactly as LINE does.
     * Every constraint this class emits is therefore a loss constraint.
     *
     * @param i the station index
     * @param r the class index
     */
    private void assertClassCapExpressible(int i, int r) {
        // The closed case is rejected rather than mapped to dropThisClass=false.
        // JMT does not hold a blocked job at the upstream station: it parks it in
        // the region's synthetic input station, where it belongs to no station
        // and leaves the upstream server free. Against SolverCTMC on a two-chain
        // closed model that costs the population count (sum of queue lengths 4.79
        // instead of 5) and moves throughput by ~29%, so the region does not
        // express this constraint at all.
        if (!Utils.isInf(sn.njobs.get(0, r))) {
            throw new RuntimeException(String.format("Station %s carries a finite capacity %d for closed class %s. "
                    + "LINE holds a blocked closed job at its upstream station, whereas JMT can only express a "
                    + "per-class capacity as a blocking region, which parks the job in the region input station "
                    + "instead, freeing the upstream server and losing it from the population count. SolverJMT "
                    + "therefore cannot reproduce this model; use SolverCTMC, SolverSSA or SolverLDES, or express "
                    + "the limit as the station capacity.", sn.nodenames.get((int) sn.stationToNode.get(i)),
                    (int) sn.classcap.get(i, r), sn.classnames.get(r)));
        }
        // The blocking strategies below have no blocking-region counterpart
        // either: a region drops or defers the arrival, it cannot hold a job in
        // the upstream server.
        DropStrategy dropRule = sn.droprule.get(sn.stations.get(i)).get(sn.jobclasses.get(r));
        if (dropRule == DropStrategy.BlockingAfterService || dropRule == DropStrategy.BlockingBeforeService
                || dropRule == DropStrategy.ReServiceOnRejection || dropRule == DropStrategy.Retrial
                || dropRule == DropStrategy.RetrialWithLimit) {
            throw new RuntimeException(String.format("Station %s applies drop strategy \"%s\" to class %s and also "
                    + "carries a finite capacity %d for it. JMT exports a per-class capacity as a blocking region, "
                    + "which can only drop or defer an arrival and cannot reproduce that strategy. Remove the "
                    + "per-class capacity or use the station capacity instead, which is exported with its drop "
                    + "strategy.", sn.nodenames.get((int) sn.stationToNode.get(i)), DropStrategy.toText(dropRule),
                    sn.classnames.get(r), (int) sn.classcap.get(i, r)));
        }
    }

    public ElementDocumentPair saveRegions(ElementDocumentPair elementDocumentPair) {
        Document simDoc = elementDocumentPair.simDoc;
        Element simElem = elementDocumentPair.simElem;

        // Per-class station capacities (sn.classcap) have no JMT counterpart on
        // the Queue section: all jmt.engine.NodeSections.Queue constructors take
        // a scalar "size" and only the Storage (Place) section carries a
        // per-class capacities array. The per-class constraint is therefore
        // exported as a synthetic single-station blocking region, since
        // jmt.engine.QueueNet.BlockingRegion is the only JMT object with a
        // maxCapacityPerClass vector.
        Matrix classCapCon = classCapConstraints();
        boolean[] covered = new boolean[sn.nstations]; // stations already inside a region emitted below

        for (int r = 0; r < simModel.getRegions().size(); r++) {
            Region region = simModel.getRegions().get(r);
            Element blockingRegion = simDoc.createElement("blockingRegion");
            blockingRegion.setAttribute("name", "FCRegion" + (r + 1));
            blockingRegion.setAttribute("type", "default");

            for (Node node : region.getNodes()) {
                Element regionNode = simDoc.createElement("regionNode");
                regionNode.setAttribute("nodeName", node.getName());
                blockingRegion.appendChild(regionNode);
            }

            // Stations of this user-defined region, and the per-class station
            // capacities that have to be merged into it. A JMT node belongs to at
            // most one blocking region (jmt.engine.QueueNet.NetNode holds a single
            // region reference), so an overlapping synthetic region is not an
            // option. Merging is faithful only when the region spans exactly the
            // constrained station: a region constrains the jobs inside the WHOLE
            // node set, so for a multi-node region no choice of maxJobsPerClass
            // reproduces a per-station limit. The drop rule is shared by every
            // constraint of a class within a region, so the region rule and the
            // capacity rule must also agree.
            List<Integer> regionStations = new ArrayList<Integer>();
            for (Node node : region.getNodes()) {
                int istr = (int) sn.nodeToStation.get(node.getNodeIndex());
                if (istr >= 0) {
                    regionStations.add(istr);
                    covered[istr] = true;
                }
            }
            Matrix regionClassCap = new Matrix(1, sn.nclasses);
            regionClassCap.fill(Inf);
            for (int idx = 0; idx < regionStations.size(); idx++) {
                int istr = regionStations.get(idx);
                for (int c = 0; c < sn.nclasses; c++) {
                    if (Utils.isInf(classCapCon.get(istr, c))) {
                        continue;
                    }
                    if (regionStations.size() > 1) {
                        throw new RuntimeException(String.format("Station %s carries a finite capacity %d for class %s "
                                + "and also belongs to the multi-station finite capacity region FCRegion%d. JMT "
                                + "constrains a blocking region as a whole and allows a node to belong to only one "
                                + "region, so a per-station class capacity cannot be expressed alongside it. Remove "
                                + "the per-class capacity or shrink the region to that single station.",
                                sn.nodenames.get((int) sn.stationToNode.get(istr)), (int) classCapCon.get(istr, c),
                                sn.classnames.get(c), r + 1));
                    }
                    assertClassCapExpressible(istr, c);
                    boolean regionDrops = region.getDropStrategy(sn.jobclasses.get(c)) == DropStrategy.Drop;
                    if (!regionDrops) {
                        throw new RuntimeException(String.format("Station %s carries a finite capacity %d for class %s "
                                + "and also belongs to region FCRegion%d, whose drop rule for that class is \"%s\". A "
                                + "JMT blocking region admits a single drop rule per class, shared by all of its "
                                + "constraints, so the two cannot be merged: the per-class capacity of an open class is "
                                + "a loss constraint. Set the region drop rule for that class to DROP.",
                                sn.nodenames.get((int) sn.stationToNode.get(istr)), (int) classCapCon.get(istr, c),
                                sn.classnames.get(c), r + 1,
                                DropStrategy.toText(region.getDropStrategy(sn.jobclasses.get(c)))));
                    }
                    regionClassCap.set(0, c, classCapCon.get(istr, c));
                }
            }

            Element globalConstraint = simDoc.createElement("globalConstraint");
            globalConstraint.setAttribute("maxJobs", String.valueOf(region.getGlobalMaxJobs()));
            blockingRegion.appendChild(globalConstraint);

            Element globalMemoryConstraint = simDoc.createElement("globalMemoryConstraint");
            double regionGlobalMaxMem = region.getGlobalMaxMemory();
            // keep integral budgets formatted as integers so the XML stays stable
            globalMemoryConstraint.setAttribute("maxMemory", regionGlobalMaxMem == Math.rint(regionGlobalMaxMem)
                    ? String.valueOf((long) regionGlobalMaxMem) : String.valueOf(regionGlobalMaxMem));
            blockingRegion.appendChild(globalMemoryConstraint);

            for (int c = 0; c < sn.nclasses; c++) {
                JobClass cstJobClass = sn.jobclasses.get(c);
                // The two constraints apply to the same node set here, so the
                // tighter one subsumes the other and min() is exact.
                double classMaxJobs = region.getClassMaxJobs(cstJobClass) == Region.UNBOUNDED
                        ? regionClassCap.get(0, c)
                        : FastMath.min(region.getClassMaxJobs(cstJobClass), regionClassCap.get(0, c));
                if (!Utils.isInf(classMaxJobs) && classMaxJobs != Region.UNBOUNDED) {
                    Element classConstraint = simDoc.createElement("classConstraint");
                    classConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    classConstraint.setAttribute("maxJobsPerClass", String.valueOf((int) classMaxJobs));
                    blockingRegion.appendChild(classConstraint);
                }

                if (region.getClassMaxMemory(cstJobClass) != Region.UNBOUNDED) {
                    Element classMemoryConstraint = simDoc.createElement("classMemoryConstraint");
                    classMemoryConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    classMemoryConstraint.setAttribute("maxMemoryPerClass", String.valueOf(region.getClassMaxMemory(cstJobClass)));
                    blockingRegion.appendChild(classMemoryConstraint);
                }

                // Always write dropRules element - JMT defaults to drop when not specified
                Element dropRuleElem = simDoc.createElement("dropRules");
                dropRuleElem.setAttribute("jobClass", region.classes.get(c).getName());
                DropStrategy dropStrategy = region.getDropStrategy(cstJobClass);
                if (dropStrategy == DropStrategy.Drop) {
                    dropRuleElem.setAttribute("dropThisClass", "true");
                } else {
                    dropRuleElem.setAttribute("dropThisClass", "false");
                }
                blockingRegion.appendChild(dropRuleElem);

                double regionClassSize = region.getClassSize(cstJobClass);
                if (regionClassSize != 1.0) {
                    Element classMemoryConstraint = simDoc.createElement("classSize");
                    classMemoryConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    // keep integral sizes formatted as integers so the XML stays stable
                    classMemoryConstraint.setAttribute("size", regionClassSize == Math.rint(regionClassSize)
                            ? String.valueOf((long) regionClassSize) : String.valueOf(regionClassSize));
                    blockingRegion.appendChild(classMemoryConstraint);
                }

                // Per-class weight toward the region capacity constraint. Drives the
                // JMT "FCR Capacity" measure (weighted occupation / Total Weight).
                double regionClassWeight = region.getClassWeight(cstJobClass);
                if (regionClassWeight != 1.0) {
                    Element classWeightElem = simDoc.createElement("classWeight");
                    classWeightElem.setAttribute("jobClass", region.classes.get(c).getName());
                    classWeightElem.setAttribute("weight", regionClassWeight == Math.rint(regionClassWeight)
                            ? String.valueOf((long) regionClassWeight) : String.valueOf(regionClassWeight));
                    blockingRegion.appendChild(classWeightElem);
                }
            }
            simElem.appendChild(blockingRegion);
        }

        // Lastly, export the per-class station capacities that no region above
        // already carries, each as a synthetic single-station blocking region.
        // The region node set is the constrained station alone, so the region
        // occupation of class c is by construction the number of class-c jobs at
        // that station, which is exactly what sn.classcap bounds in
        // State.afterEventStation. The global constraint is left unbounded: the
        // station total sn.cap is already exported as the Queue "size" (see
        // saveBufferCapacity), and refreshCapacity guarantees
        // classcap(i,r) <= cap(i), so the region can never contradict the total.
        int classCapRegionIdx = 0;
        for (int i = 0; i < sn.nstations; i++) {
            if (covered[i]) {
                continue;
            }
            boolean anyConstraint = false;
            for (int c = 0; c < sn.nclasses; c++) {
                if (!Utils.isInf(classCapCon.get(i, c))) {
                    anyConstraint = true;
                    break;
                }
            }
            if (!anyConstraint) {
                continue;
            }
            // Validate every constrained class before emitting anything, so that
            // an inexpressible capacity throws instead of half-writing a region.
            for (int c = 0; c < sn.nclasses; c++) {
                if (!Utils.isInf(classCapCon.get(i, c))) {
                    assertClassCapExpressible(i, c);
                }
            }
            classCapRegionIdx++;

            Element blockingRegion = simDoc.createElement("blockingRegion");
            blockingRegion.setAttribute("name", "ClassCapRegion" + classCapRegionIdx);
            blockingRegion.setAttribute("type", "default");

            Element regionNode = simDoc.createElement("regionNode");
            regionNode.setAttribute("nodeName", sn.nodenames.get((int) sn.stationToNode.get(i)));
            blockingRegion.appendChild(regionNode);

            Element globalConstraint = simDoc.createElement("globalConstraint");
            globalConstraint.setAttribute("maxJobs", "-1");
            blockingRegion.appendChild(globalConstraint);

            Element globalMemoryConstraint = simDoc.createElement("globalMemoryConstraint");
            globalMemoryConstraint.setAttribute("maxMemory", "-1");
            blockingRegion.appendChild(globalMemoryConstraint);

            for (int c = 0; c < sn.nclasses; c++) {
                if (!Utils.isInf(classCapCon.get(i, c))) {
                    Element classConstraint = simDoc.createElement("classConstraint");
                    classConstraint.setAttribute("jobClass", sn.classnames.get(c));
                    classConstraint.setAttribute("maxJobsPerClass", String.valueOf((int) classCapCon.get(i, c)));
                    blockingRegion.appendChild(classConstraint);
                }
            }

            // Only the constrained classes need a drop rule: an unconstrained
            // class has neither a class constraint nor a global one, so
            // BlockingRegion.isBlocked is always false for it and its flag would
            // never be read. Every constraint emitted here is a loss constraint
            // (assertClassCapExpressible rejects the rest).
            for (int c = 0; c < sn.nclasses; c++) {
                if (!Utils.isInf(classCapCon.get(i, c))) {
                    Element dropRuleElem = simDoc.createElement("dropRules");
                    dropRuleElem.setAttribute("jobClass", sn.classnames.get(c));
                    dropRuleElem.setAttribute("dropThisClass", "true");
                    blockingRegion.appendChild(dropRuleElem);
                }
            }

            simElem.appendChild(blockingRegion);
        }

        return new ElementDocumentPair(simElem, simDoc);
    }

    public DocumentSectionPair saveRoutingStrategy(DocumentSectionPair documentSectionPair, int ind) {
//        [SIMDOC, SECTION] = SAVEROUTINGSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategy");
        strategyNode.setAttribute("name", "RoutingStrategy");

        int K = sn.nclasses;
        Node indNode = simModel.getNodes().get(ind);
        boolean[] exportClasses = getExportableClasses();
        // since the class switch node always outputs to a single node, it is faster to translate it to RAND. Also some problems with sn.rt value otherwise.
        if (sn.nodetype.get(ind) == NodeType.ClassSwitch) {
            for (JobClass jobClass : sn.routing.get(indNode).keySet()) {
                sn.routing.get(indNode).put(jobClass, RoutingStrategy.RAND);
            }
        }
        for (int r = 0; r < K; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode);

            Element concStratNode = simDoc.createElement("subParameter");
            Element concStratNode2;
            Matrix conn_i;
            Matrix conn_i_find;

            RoutingStrategy routingStrategy = sn.routing.get(indNode).get(sn.jobclasses.get(r));
            switch (routingStrategy) {
                case RAND:
                    // MATLAB uses RandomStrategy for RAND routing (see saveRoutingStrategy.m lines 31-34)
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy");
                    concStratNode.setAttribute("name", "Random");
                    break;
                case RROBIN:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy");
                    concStratNode.setAttribute("name", "Round Robin");
                    break;
                case JSQ:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.ShortestQueueLengthRoutingStrategy");
                    concStratNode.setAttribute("name", "Join the Shortest Queue (JSQ)");
                    break;
                case SQ: {
                    // Power of K routing strategy implementation
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy");
                    concStratNode.setAttribute("name", "Power of k");

                    // Read d from NodeParam (populated in Network.refreshLocalVars).
                    // withMemory is always false: LINE implements SQ(d) only, while JMT's
                    // flag selects Anselmi & Dufour SQ(d,N), a different policy.
                    int kVal = 2;
                    boolean memoryVal = false;
                    if (sn.nodeparam != null) {
                        jline.lang.NodeParam np = sn.nodeparam.get(indNode);
                        if (np != null) {
                            jline.lang.JobClass jc = sn.jobclasses.get(r);
                            if (np.d != null && np.d.containsKey(jc)) {
                                kVal = np.d.get(jc);
                            }
                        }
                    }

                    concStratNode2 = simDoc.createElement("subParameter");
                    concStratNode2.setAttribute("classPath", "java.lang.Integer");
                    concStratNode2.setAttribute("name", "k");
                    Element kValueNode = simDoc.createElement("value");
                    kValueNode.appendChild(simDoc.createTextNode(Integer.toString(kVal)));
                    concStratNode2.appendChild(kValueNode);
                    concStratNode.appendChild(concStratNode2);

                    Element memNode = simDoc.createElement("subParameter");
                    memNode.setAttribute("classPath", "java.lang.Boolean");
                    memNode.setAttribute("name", "withMemory");
                    Element memValueNode = simDoc.createElement("value");
                    memValueNode.appendChild(simDoc.createTextNode(memoryVal ? "true" : "false"));
                    memNode.appendChild(memValueNode);
                    concStratNode.appendChild(memNode);
                    break;
                }
//                    Element concStratNode2ValueNode = simDoc.createElement("value");
//                    concStratNode2ValueNode.appendChild(simDoc.createTextNode(String.format("%d", sn.nodeparam.get(indNode).d.get(sn.jobclasses.get(r)))));
//                    concStratNode2.appendChild(concStratNode2ValueNode);
//                    concStratNode3 = simDoc.createElement("subParameter");
//                    concStratNode3.setAttribute("classPath", "java.lang.Boolean");
//                    concStratNode3.setAttribute("name", "withMemory");
//                    Element concStratNode3ValueNode = simDoc.createElement("value");
//                    if (!sn.nodeparam.get(indNode).m.get(sn.jobclasses.get(r)).isEmpty()) {
//                        concStratNode3ValueNode.appendChild(simDoc.createTextNode("true"));
//                    } else {
//                        concStratNode3ValueNode.appendChild(simDoc.createTextNode("false"));
//                    }
//                    concStratNode3.appendChild(concStratNode3ValueNode);
//                    concStratNode.appendChild(concStratNode2);
//                    concStratNode.appendChild(concStratNode3);
                case WRROBIN:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.WeightedRoundRobinStrategy");
                    concStratNode.setAttribute("name", "Weighted Round Robin");
                    concStratNode2 = simDoc.createElement("subParameter");
                    concStratNode2.setAttribute("array", "true");
                    concStratNode2.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.WeightEntry");
                    concStratNode2.setAttribute("name", "WeightEntryArray");

                    // linked stations
                    conn_i = new Matrix(0, 0);
                    Matrix.extractRows(simModel.getStruct().connmatrix, ind, ind + 1, conn_i);
                    conn_i_find = conn_i.find();
                    for (int idx = 0; idx < conn_i_find.length(); idx++) {
                        int j = (int) conn_i_find.get(idx);
                        double weight = sn.nodeparam.get(indNode).weights.get(sn.jobclasses.get(r)).get(j);

                        Element concStratNode3 = simDoc.createElement("subParameter");
                        concStratNode3.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.WeightEntry");
                        concStratNode3.setAttribute("name", "WeightEntry");

                        Element concStratNode4Station = simDoc.createElement("subParameter");
                        concStratNode4Station.setAttribute("classPath", "java.lang.String");
                        concStratNode4Station.setAttribute("name", "stationName");

                        Element concStratNode4StationValueNode = simDoc.createElement("value");
                        concStratNode4StationValueNode.appendChild(simDoc.createTextNode(sn.nodenames.get(j)));
                        concStratNode4Station.appendChild(concStratNode4StationValueNode);
                        concStratNode3.appendChild(concStratNode4Station);
                        Element concStratNode4Weight = simDoc.createElement("subParameter");
                        concStratNode4Weight.setAttribute("classPath", "java.lang.Integer");
                        concStratNode4Weight.setAttribute("name", "weight");
                        Element concStratNode4WeightValueNode = simDoc.createElement("value");
                        concStratNode4WeightValueNode.appendChild(simDoc.createTextNode(String.format("%d", (int) weight)));
                        concStratNode4Weight.appendChild(concStratNode4WeightValueNode);
                        concStratNode3.appendChild(concStratNode4Station);
                        concStratNode3.appendChild(concStratNode4Weight);
                        concStratNode2.appendChild(concStratNode3);
                    }
                    concStratNode.appendChild(concStratNode2);
                    break;
                case PROB:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy");
                    concStratNode.setAttribute("name", "Probabilities");
                    concStratNode2 = simDoc.createElement("subParameter");
                    concStratNode2.setAttribute("array", "true");
                    concStratNode2.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                    concStratNode2.setAttribute("name", "EmpiricalEntryArray");

                    // linked stations
                    conn_i = new Matrix(0, 0);
                    Matrix.extractRows(sn.connmatrix, ind, ind + 1, conn_i);
                    conn_i_find = conn_i.find();
                    for (int idx = 0; idx < conn_i_find.length(); idx++) {
                        int j = (int) conn_i_find.get(idx);
                        double probRouting = sn.rtnodes.get(ind * K + r, j * K + r);
                        if (probRouting > 0) {
                            Element concStratNode3 = simDoc.createElement("subParameter");
                            concStratNode3.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                            concStratNode3.setAttribute("name", "EmpiricalEntry");
                            Element concStratNode4Station = simDoc.createElement("subParameter");
                            concStratNode4Station.setAttribute("classPath", "java.lang.String");
                            concStratNode4Station.setAttribute("name", "stationName");
                            Element concStratNode4StationValueNode = simDoc.createElement("value");
                            concStratNode4StationValueNode.appendChild(simDoc.createTextNode(sn.nodenames.get(j)));
                            concStratNode4Station.appendChild(concStratNode4StationValueNode);
                            concStratNode3.appendChild(concStratNode4Station);
                            Element concStratNode4Probability = simDoc.createElement("subParameter");
                            concStratNode4Probability.setAttribute("classPath", "java.lang.Double");
                            concStratNode4Probability.setAttribute("name", "probability");
                            Element concStratNode4ProbabilityValueNode = simDoc.createElement("value");
                            concStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode(String.format("%12.12f", probRouting)));
                            concStratNode4Probability.appendChild(concStratNode4ProbabilityValueNode);
                            concStratNode3.appendChild(concStratNode4Station);
                            concStratNode3.appendChild(concStratNode4Probability);
                            concStratNode2.appendChild(concStratNode3);
                        }
                    }
                    concStratNode.appendChild(concStratNode2);
                    break;
                case DISABLED:
                default:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy");
                    concStratNode.setAttribute("name", "Disabled");
            }
            strategyNode.appendChild(concStratNode);
            section.appendChild(strategyNode);
        }
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveServerVisits(DocumentSectionPair documentSectionPair) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element visitsNode = simDoc.createElement("parameter");
        visitsNode.setAttribute("array", "true");
        visitsNode.setAttribute("classPath", "java.lang.Integer");
        visitsNode.setAttribute("name", "numberOfVisits");

        int numOfClasses = sn.nclasses;
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            visitsNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.Integer");
            subParameterNode.setAttribute("name", "numberOfVisits");

            Element valueNode2 = simDoc.createElement("value");
            valueNode2.appendChild(simDoc.createTextNode("1"));

            subParameterNode.appendChild(valueNode2);
            visitsNode.appendChild(subParameterNode);
            section.appendChild(visitsNode);
        }

        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveServiceStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVESERVICESTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        strategyNode.setAttribute("name", "ServiceStrategy");

        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        Station istStation = sn.stations.get(i);
        boolean[] exportClasses = getExportableClasses();

        for (int r = 0; r < numOfClasses; r++) {
            if (!exportClasses[r]) continue;
            JobClass rstJobClass = sn.jobclasses.get(r);
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode2);
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            ProcessType type = sn.procid.get(istStation).get(rstJobClass);
            if (type == ProcessType.DISABLED) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "DisabledServiceTimeStrategy");
            } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.IMMEDIATE) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
            } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.PH
                    || sn.procid.get(istStation).get(rstJobClass) == ProcessType.APH
                    || sn.procid.get(istStation).get(rstJobClass) == ProcessType.COX2
                    || sn.procid.get(istStation).get(rstJobClass) == ProcessType.COXIAN
                    || (sn.procid.get(istStation).get(rstJobClass) == ProcessType.HYPEREXP
                        && sn.phases.get(i, r) > 2)) {
                // jmt.engine.random.HyperExpPar is 2-phase only, so an n-phase
                // hyper-exponential is emitted as the phase-type it already is
                // (same alpha and T, not a moment-matched refit). This mirrors
                // the n-phase HyperExp handling in saveTimingStrategies.
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.PhaseTypeDistr");
                distributionNode.setAttribute("name", "Phase-Type");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.PhaseTypePar");
                distrParNode.setAttribute("name", "distrPar");

                Element subParNodeAlpha = simDoc.createElement("subParameter");
                subParNodeAlpha.setAttribute("array", "true");
                subParNodeAlpha.setAttribute("classPath", "java.lang.Object");
                subParNodeAlpha.setAttribute("name", "alpha");
                Element subParNodeAlphaVec = simDoc.createElement("subParameter");
                subParNodeAlphaVec.setAttribute("array", "true");
                subParNodeAlphaVec.setAttribute("classPath", "java.lang.Object");
                subParNodeAlphaVec.setAttribute("name", "vector");

                MatrixCell PH = sn.proc.get(istStation).get(rstJobClass);
                Matrix alpha = sn.pie.get(istStation).get(rstJobClass);
                alpha.absEq();
                for (int k = 0; k < sn.phases.get(i, r); k++) {
                    Element subParNodeAlphaElem = simDoc.createElement("subParameter");
                    subParNodeAlphaElem.setAttribute("classPath", "java.lang.Double");
                    subParNodeAlphaElem.setAttribute("name", "entry");
                    Element subParValue = simDoc.createElement("value");
                    subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha.get(k))));
                    subParNodeAlphaElem.appendChild(subParValue);
                    subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
                }
                Element subParNodeT = simDoc.createElement("subParameter");
                subParNodeT.setAttribute("array", "true");
                subParNodeT.setAttribute("classPath", "java.lang.Object");
                subParNodeT.setAttribute("name", "T");
                Matrix T = PH.get(0);

                for (int k = 0; k < sn.phases.get(i, r); k++) {
                    Element subParNodeTvec = simDoc.createElement("subParameter");
                    subParNodeTvec.setAttribute("array", "true");
                    subParNodeTvec.setAttribute("classPath", "java.lang.Object");
                    subParNodeTvec.setAttribute("name", "vector");

                    for (int j = 0; j < sn.phases.get(i, r); j++) {
                        Element subParNodeTElem = simDoc.createElement("subParameter");
                        subParNodeTElem.setAttribute("classPath", "java.lang.Double");
                        subParNodeTElem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        if (k == j) {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * FastMath.abs(T.get(k, j)))));
                        } else {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", FastMath.abs(T.get(k, j)))));
                        }
                        subParNodeTElem.appendChild(subParValue);
                        subParNodeTvec.appendChild(subParNodeTElem);
                    }
                    subParNodeT.appendChild(subParNodeTvec);
                }
                subParNodeAlpha.appendChild(subParNodeAlphaVec);
                distrParNode.appendChild(subParNodeAlpha);
                distrParNode.appendChild(subParNodeT);
                serviceTimeStrategyNode.appendChild(distributionNode);
                serviceTimeStrategyNode.appendChild(distrParNode);
            } else if (sn.procid.get(istStation).get(rstJobClass) == ProcessType.MAP) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                distributionNode.setAttribute("name", "Burst (MAP)");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                distrParNode.setAttribute("name", "distrPar");

                MatrixCell MAP = sn.proc.get(istStation).get(rstJobClass);

                Element subParNodeD0 = simDoc.createElement("subParameter");
                subParNodeD0.setAttribute("array", "true");
                subParNodeD0.setAttribute("classPath", "java.lang.Object");
                subParNodeD0.setAttribute("name", "D0");
                Matrix D0 = MAP.get(0);

                for (int k = 0; k < sn.phases.get(i, r); k++) {
                    Element subParNodeD0vec = simDoc.createElement("subParameter");
                    subParNodeD0vec.setAttribute("array", "true");
                    subParNodeD0vec.setAttribute("classPath", "java.lang.Object");
                    subParNodeD0vec.setAttribute("name", "vector");

                    for (int j = 0; j < sn.phases.get(i, r); j++) {
                        Element subParNodeD0Elem = simDoc.createElement("subParameter");
                        subParNodeD0Elem.setAttribute("classPath", "java.lang.Double");
                        subParNodeD0Elem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D0.get(k, j))));
                        subParNodeD0Elem.appendChild(subParValue);
                        subParNodeD0vec.appendChild(subParNodeD0Elem);
                    }
                    subParNodeD0.appendChild(subParNodeD0vec);
                }
                distrParNode.appendChild(subParNodeD0);

                Element subParNodeD1 = simDoc.createElement("subParameter");
                subParNodeD1.setAttribute("array", "true");
                subParNodeD1.setAttribute("classPath", "java.lang.Object");
                subParNodeD1.setAttribute("name", "D1");
                Matrix D1 = MAP.get(1);

                for (int k = 0; k < sn.phases.get(i, r); k++) {
                    Element subParNodeD1vec = simDoc.createElement("subParameter");
                    subParNodeD1vec.setAttribute("array", "true");
                    subParNodeD1vec.setAttribute("classPath", "java.lang.Object");
                    subParNodeD1vec.setAttribute("name", "vector");

                    for (int j = 0; j < sn.phases.get(i, r); j++) {
                        Element subParNodeD1Elem = simDoc.createElement("subParameter");
                        subParNodeD1Elem.setAttribute("classPath", "java.lang.Double");
                        subParNodeD1Elem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D1.get(k, j))));
                        subParNodeD1Elem.appendChild(subParValue);
                        subParNodeD1vec.appendChild(subParNodeD1Elem);
                    }
                    subParNodeD1.appendChild(subParNodeD1vec);
                }
                distrParNode.appendChild(subParNodeD1);
                serviceTimeStrategyNode.appendChild(distributionNode);
                serviceTimeStrategyNode.appendChild(distrParNode);
            } else {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

                Element distributionNode = simDoc.createElement("subParameter");
                String javaClass = "";
                String javaParClass = "";
                switch (sn.procid.get(istStation).get(rstJobClass)) {
                    case DET:
                        javaClass = "jmt.engine.random.DeterministicDistr";
                        javaParClass = "jmt.engine.random.DeterministicDistrPar";
                        break;
                    case COX2:
                    case COXIAN:
                        javaClass = "jmt.engine.random.CoxianDistr";
                        javaParClass = "jmt.engine.random.CoxianPar";
                        break;
                    case ERLANG:
                        javaClass = "jmt.engine.random.Erlang";
                        javaParClass = "jmt.engine.random.ErlangPar";
                        break;
                    case EXP:
                        javaClass = "jmt.engine.random.Exponential";
                        javaParClass = "jmt.engine.random.ExponentialPar";
                        break;
                    case GAMMA:
                        javaClass = "jmt.engine.random.GammaDistr";
                        javaParClass = "jmt.engine.random.GammaDistrPar";
                        break;
                    case HYPEREXP:
                        javaClass = "jmt.engine.random.HyperExp";
                        javaParClass = "jmt.engine.random.HyperExpPar";
                        break;
                    case PARETO:
                        javaClass = "jmt.engine.random.Pareto";
                        javaParClass = "jmt.engine.random.ParetoPar";
                        break;
                    case WEIBULL:
                        javaClass = "jmt.engine.random.Weibull";
                        javaParClass = "jmt.engine.random.WeibullPar";
                        break;
                    case LOGNORMAL:
                        javaClass = "jmt.engine.random.Lognormal";
                        javaParClass = "jmt.engine.random.LognormalPar";
                        break;
                    case UNIFORM:
                        javaClass = "jmt.engine.random.Uniform";
                        javaParClass = "jmt.engine.random.UniformPar";
                        break;
                    case MMPP2:
                        javaClass = "jmt.engine.random.MMPP2Distr";
                        javaParClass = "jmt.engine.random.MMPP2Par";
                        break;
                    case REPLAYER:
                    case TRACE:
                        javaClass = "jmt.engine.random.Replayer";
                        javaParClass = "jmt.engine.random.ReplayerPar";
                        break;
                }
                distributionNode.setAttribute("classPath", javaClass);
                switch (sn.procid.get(istStation).get(rstJobClass)) {
                    case REPLAYER:
                    case TRACE:
                        distributionNode.setAttribute("name", "Replayer");
                        break;
                    case EXP:
                        distributionNode.setAttribute("name", "Exponential");
                        break;
                    case HYPEREXP:
                        distributionNode.setAttribute("name", "Hyperexponential");
                        break;
                    default:
                        distributionNode.setAttribute("name", ProcessType.toText(sn.procid.get(istStation).get(rstJobClass)));
                        break;
                }
                serviceTimeStrategyNode.appendChild(distributionNode);

                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");

                Element subParNodeAlpha = simDoc.createElement("subParameter");
                Element subParValue = simDoc.createElement("value");
                double c;
                switch (sn.procid.get(istStation).get(rstJobClass)) {
                    case DET:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "t");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / sn.rates.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case EXP:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case HYPEREXP:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "p");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.pie.get(istStation).get(rstJobClass).get(0))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda2");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).get(1, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case ERLANG:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(i, r) * sn.phases.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Long");
                        subParNodeAlpha.setAttribute("name", "r");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.valueOf((int) sn.phases.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case MMPP2:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda0");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.proc.get(istStation).get(rstJobClass).get(1).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.proc.get(istStation).get(rstJobClass).get(1).get(1, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma0");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.proc.get(istStation).get(rstJobClass).get(0).get(0, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.proc.get(istStation).get(rstJobClass).get(0).get(1, 0))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case GAMMA:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1.0 / sn.scv.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "beta");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.scv.get(i, r) / sn.rates.get(i, r))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case PARETO:
                        double shape = FastMath.sqrt(1 + 1.0 / sn.scv.get(i, r)) + 1;
                        double scale = 1.0 / sn.rates.get(i, r) * (shape - 1) / shape;
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", shape)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "k");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", scale)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case WEIBULL:
                        c = FastMath.sqrt(sn.scv.get(i, r));
                        double rval = FastMath.pow(c, -1.086); //Justus approximation (1976)
                        double alpha = 1.0 / sn.rates.get(i, r) / Maths.gammaFunction(1 + 1.0 / rval);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "r");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", rval)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case LOGNORMAL:
                        c = FastMath.sqrt(sn.scv.get(i, r));
                        double mu = FastMath.log(1 / sn.rates.get(i, r) / FastMath.sqrt(c * c + 1));
                        double sigma = FastMath.sqrt(Math.log(c * c + 1));
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "mu");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", mu)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sigma)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case UNIFORM:
                        double maxVal = (Math.sqrt(12 * sn.scv.get(i, r) / FastMath.pow(sn.rates.get(i, r), 2)) + 2 / sn.rates.get(i, r)) / 2;
                        double minVal = 2 / sn.rates.get(i, r) - maxVal;
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "min");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", minVal)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "max");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", maxVal)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case REPLAYER:
                    case TRACE:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.String");
                        subParNodeAlpha.setAttribute("name", "fileName");
                        subParValue.appendChild(simDoc.createTextNode(((ServiceNodeParam) sn.nodeparam.get(simModel.getNodes().get(i))).fileName.get(r)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                }
                serviceTimeStrategyNode.appendChild(distrParNode);
            }
            strategyNode.appendChild(serviceTimeStrategyNode);
            section.appendChild(strategyNode);
        }

        return new DocumentSectionPair(simDoc, section);
    }

    // Save the timing strategy for the given node
    public DocumentSectionPair saveTimingStrategies(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        strategyNode.setAttribute("name", "timingStrategies");

        Node istNode = simModel.getNodes().get(ind);
        int numOfModes = ((TransitionNodeParam) sn.nodeparam.get(istNode)).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            int nphases = (int) ((TransitionNodeParam) sn.nodeparam.get(simModel.getNodes().get(ind))).firingphases.get(m);
            Element timimgStrategyNode = simDoc.createElement("subParameter");
            Mode mstMode = ((Transition) istNode).getModes().get(m);

            if (((TransitionNodeParam) sn.nodeparam.get(istNode)).timing.get(m) == TimingStrategy.IMMEDIATE) {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
            } else if (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode) == ProcessType.APH || ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode) == ProcessType.COX2 || ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode) == ProcessType.COXIAN || (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingphases.get(m) > 2 && ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode) == ProcessType.HYPEREXP)) {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "timingStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.PhaseTypeDistr");
                distributionNode.setAttribute("name", "Phase-Type");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.PhaseTypePar");
                distrParNode.setAttribute("name", "distrPar");

                Element subParNodeAlpha = simDoc.createElement("subParameter");
                subParNodeAlpha.setAttribute("array", "true");
                subParNodeAlpha.setAttribute("classPath", "java.lang.Object");
                subParNodeAlpha.setAttribute("name", "alpha");
                Element subParNodeAlphaVec = simDoc.createElement("subParameter");
                subParNodeAlphaVec.setAttribute("array", "true");
                subParNodeAlphaVec.setAttribute("classPath", "java.lang.Object");
                subParNodeAlphaVec.setAttribute("name", "vector");
                MatrixCell PH = ((TransitionNodeParam) sn.nodeparam.get(simModel.getNodes().get(ind))).firingproc.get(mstMode);
                Matrix alpha = Map_pie.map_pie(PH.get(0), PH.get(1));
                alpha.absEq();
                for (int k = 0; k < nphases; k++) {
                    Element subParNodeAlphaElem = simDoc.createElement("subParameter");
                    subParNodeAlphaElem.setAttribute("classPath", "java.lang.Double");
                    subParNodeAlphaElem.setAttribute("name", "entry");
                    Element subParValue = simDoc.createElement("value");
                    subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", alpha.get(k))));
                    subParNodeAlphaElem.appendChild(subParValue);
                    subParNodeAlphaVec.appendChild(subParNodeAlphaElem);
                }

                Element subParNodeT = simDoc.createElement("subParameter");
                subParNodeT.setAttribute("array", "true");
                subParNodeT.setAttribute("classPath", "java.lang.Object");
                subParNodeT.setAttribute("name", "T");

                Matrix T = PH.get(0);

                for (int k = 0; k < nphases; k++) {
                    Element subParNodeTvec = simDoc.createElement("subParameter");
                    subParNodeTvec.setAttribute("array", "true");
                    subParNodeTvec.setAttribute("classPath", "java.lang.Object");
                    subParNodeTvec.setAttribute("name", "vector");
                    for (int j = 0; j < nphases; j++) {
                        Element subParNodeTElem = simDoc.createElement("subParameter");
                        subParNodeTElem.setAttribute("classPath", "java.lang.Double");
                        subParNodeTElem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        if (k == j) {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * FastMath.abs(T.get(k, j)))));
                        } else {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", FastMath.abs(T.get(k, j)))));
                        }
                        subParNodeTElem.appendChild(subParValue);
                        subParNodeTvec.appendChild(subParNodeTElem);
                    }
                    subParNodeT.appendChild(subParNodeTvec);
                }
                subParNodeAlpha.appendChild(subParNodeAlphaVec);
                distrParNode.appendChild(subParNodeAlpha);
                distrParNode.appendChild(subParNodeT);
                timimgStrategyNode.appendChild(distributionNode);
                timimgStrategyNode.appendChild(distrParNode);
            } else if (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode) == ProcessType.MAP) {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "timingStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                distributionNode.setAttribute("name", "Burst (MAP)");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                distrParNode.setAttribute("name", "distrPar");

                MatrixCell MAP = ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingproc.get(mstMode);

                Element subParNodeD0 = simDoc.createElement("subParameter");
                subParNodeD0.setAttribute("array", "true");
                subParNodeD0.setAttribute("classPath", "java.lang.Object");
                subParNodeD0.setAttribute("name", "D0");
                Matrix D0 = MAP.get(0);
                for (int k = 0; k < nphases; k++) {
                    Element subParNodeD0vec = simDoc.createElement("subParameter");
                    subParNodeD0vec.setAttribute("array", "true");
                    subParNodeD0vec.setAttribute("classPath", "java.lang.Object");
                    subParNodeD0vec.setAttribute("name", "vector");
                    for (int j = 0; j < nphases; j++) {
                        Element subParNodeD0Elem = simDoc.createElement("subParameter");
                        subParNodeD0Elem.setAttribute("classPath", "java.lang.Double");
                        subParNodeD0Elem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D0.get(k, j))));
                        subParNodeD0Elem.appendChild(subParValue);
                        subParNodeD0vec.appendChild(subParNodeD0Elem);
                    }
                    subParNodeD0.appendChild(subParNodeD0vec);
                }
                distrParNode.appendChild(subParNodeD0);

                Element subParNodeD1 = simDoc.createElement("subParameter");
                subParNodeD1.setAttribute("array", "true");
                subParNodeD1.setAttribute("classPath", "java.lang.Object");
                subParNodeD1.setAttribute("name", "D1");
                Matrix D1 = MAP.get(1);
                for (int k = 0; k < nphases; k++) {
                    Element subParNodeD1vec = simDoc.createElement("subParameter");
                    subParNodeD1vec.setAttribute("array", "true");
                    subParNodeD1vec.setAttribute("classPath", "java.lang.Object");
                    subParNodeD1vec.setAttribute("name", "vector");
                    for (int j = 0; j < nphases; j++) {
                        Element subParNodeD1Elem = simDoc.createElement("subParameter");
                        subParNodeD1Elem.setAttribute("classPath", "java.lang.Double");
                        subParNodeD1Elem.setAttribute("name", "entry");
                        Element subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", D1.get(k, j))));
                        subParNodeD1Elem.appendChild(subParValue);
                        subParNodeD1vec.appendChild(subParNodeD1Elem);
                    }
                    subParNodeD1.appendChild(subParNodeD1vec);
                }
                distrParNode.appendChild(subParNodeD1);
                timimgStrategyNode.appendChild(distributionNode);
                timimgStrategyNode.appendChild(distrParNode);
            } else {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "timingStrategy");

                Element distributionNode = simDoc.createElement("subParameter");
                String javaClass = "";
                String javaParClass = "";
                switch (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode)) {
                    case DET:
                        javaClass = "jmt.engine.random.DeterministicDistr";
                        javaParClass = "jmt.engine.random.DeterministicDistrPar";
                        break;
                    case COX2:
                    case COXIAN:
                        javaClass = "jmt.engine.random.CoxianDistr";
                        javaParClass = "jmt.engine.random.CoxianPar";
                        break;
                    case ERLANG:
                        javaClass = "jmt.engine.random.Erlang";
                        javaParClass = "jmt.engine.random.ErlangPar";
                        break;
                    case EXP:
                        javaClass = "jmt.engine.random.Exponential";
                        javaParClass = "jmt.engine.random.ExponentialPar";
                        break;
                    case GAMMA:
                        javaClass = "jmt.engine.random.GammaDistr";
                        javaParClass = "jmt.engine.random.GammaDistrPar";
                        break;
                    case HYPEREXP:
                        javaClass = "jmt.engine.random.HyperExp";
                        javaParClass = "jmt.engine.random.HyperExpPar";
                        break;
                    case PARETO:
                        javaClass = "jmt.engine.random.Pareto";
                        javaParClass = "jmt.engine.random.ParetoPar";
                        break;
                    case WEIBULL:
                        javaClass = "jmt.engine.random.Weibull";
                        javaParClass = "jmt.engine.random.WeibullPar";
                        break;
                    case LOGNORMAL:
                        javaClass = "jmt.engine.random.Lognormal";
                        javaParClass = "jmt.engine.random.LognormalPar";
                        break;
                    case UNIFORM:
                        javaClass = "jmt.engine.random.Uniform";
                        javaParClass = "jmt.engine.random.UniformPar";
                        break;
                    case MMPP2:
                        javaClass = "jmt.engine.random.MMPP2Distr";
                        javaParClass = "jmt.engine.random.MMPP2Par";
                        break;
                    case REPLAYER:
                    case TRACE:
                        javaClass = "jmt.engine.random.Replayer";
                        javaParClass = "jmt.engine.random.ReplayerPar";
                        break;
                }
                distributionNode.setAttribute("classPath", javaClass);
                switch (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode)) {
                    case REPLAYER:
                    case TRACE:
                        distributionNode.setAttribute("name", "Replayer");
                        break;
                    case EXP:
                        distributionNode.setAttribute("name", "Exponential");
                        break;
                    case HYPEREXP:
                        distributionNode.setAttribute("name", "Hyperexponential");
                        break;
                    default:
                        distributionNode.setAttribute("name", ProcessType.toText(((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode)));
                }
                timimgStrategyNode.appendChild(distributionNode);

                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");
                Element subParNodeAlpha = simDoc.createElement("subParameter");
                Element subParValue = simDoc.createElement("value");

                MatrixCell matrixMap = ((TransitionNodeParam) sn.nodeparam.get(istNode)).firingproc.get(mstMode);

                switch (((TransitionNodeParam) sn.nodeparam.get(istNode)).firingprocid.get(mstMode)) {
                    case DET:
                        // Network.refreshPetriNetNodes stores firingproc as
                        // (singleton(mean), singleton(scv)) for non-Markovian
                        // non-Pareto distributions. Read mean directly.
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "t");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(0).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case EXP:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Map_lambda.map_lambda(matrixMap.get(0), matrixMap.get(1)))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case HYPEREXP:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "p");
                        subParValue = simDoc.createElement("value");
                        Matrix pie = Map_pie.map_pie(matrixMap.get(0), matrixMap.get(1));
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", pie.get(0))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * matrixMap.get(0).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda2");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * matrixMap.get(0).get(1, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case ERLANG:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue = simDoc.createElement("value");
                        double timingrate = Map_lambda.map_lambda(matrixMap.get(0), matrixMap.get(1));
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", timingrate * nphases)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Long");
                        subParNodeAlpha.setAttribute("name", "r");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%d", nphases)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case MMPP2:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda0");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(1).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(1).get(1, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma0");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(0).get(0, 1))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(0).get(1, 0))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case PARETO:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(0).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "k");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(1).value())));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case GAMMA:
                        // firingproc storage is (singleton(mean), singleton(scv));
                        // see Network.refreshPetriNetNodes for non-Markovian
                        // non-Pareto distributions.
                        double gammaMean = matrixMap.get(0).value();
                        double gammaScv = matrixMap.get(1).value();
                        double gammaAlpha = 1.0 / gammaScv;  // shape
                        double gammaBeta = gammaScv * gammaMean;  // scale = scv * mean
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", gammaAlpha)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "beta");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", gammaBeta)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case WEIBULL:
                        // Weibull shape from SCV using Justus approximation (1976).
                        // firingproc storage is (singleton(mean), singleton(scv)).
                        double weibullMean = matrixMap.get(0).value();
                        double weibullScv = matrixMap.get(1).value();
                        double weibullR = FastMath.pow(FastMath.sqrt(weibullScv), -1.086);
                        double weibullAlpha = weibullMean / Maths.gammaFunction(1 + 1.0 / weibullR);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "alpha");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", weibullAlpha)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "r");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", weibullR)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case LOGNORMAL:
                        // firingproc storage is (singleton(mean), singleton(scv));
                        // see Network.refreshPetriNetNodes.
                        double lognormalMean = matrixMap.get(0).value();
                        double lognormalScv = matrixMap.get(1).value();
                        double lognormalMu = FastMath.log(lognormalMean / FastMath.sqrt(lognormalScv + 1));
                        double lognormalSigma = FastMath.sqrt(FastMath.log(lognormalScv + 1));
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "mu");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", lognormalMu)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "sigma");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", lognormalSigma)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case UNIFORM:
                        // firingproc storage is (singleton(mean), singleton(scv));
                        // see Network.refreshPetriNetNodes.
                        double uniformMean = matrixMap.get(0).value();
                        double uniformScv = matrixMap.get(1).value();
                        double uniformMax = (FastMath.sqrt(12 * uniformScv * uniformMean * uniformMean) + 2 * uniformMean) / 2;
                        double uniformMin = 2 * uniformMean - uniformMax;
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "min");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", uniformMin)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "max");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", uniformMax)));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case REPLAYER:
                    case TRACE:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.String");
                        subParNodeAlpha.setAttribute("name", "fileName");
                        subParValue = simDoc.createElement("value");
                        // Use fileName from ServiceNodeParam parent class if available
                        TransitionNodeParam transParam = (TransitionNodeParam) sn.nodeparam.get(istNode);
                        String traceFileName = "";
                        if (transParam != null && transParam.fileName != null && m < transParam.fileName.size()) {
                            traceFileName = transParam.fileName.get(m);
                        }
                        subParValue.appendChild(simDoc.createTextNode(traceFileName));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    default:
                        throw new RuntimeException(String.format("Unsupported firing distribution type for mode %d", m));
                }
                timimgStrategyNode.appendChild(distrParNode);
            }
            strategyNode.appendChild(timimgStrategyNode);
        }
        section.appendChild(strategyNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveTotalCapacity(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "totalCapacity");
        Element valueNode = simDoc.createElement("value");
        if (sn.isstation.get(ind) == 0 || Utils.isInf(sn.cap.get((int) sn.nodeToStation.get(ind)))) {
            valueNode.appendChild(simDoc.createTextNode("-1"));
        } else {
            if (sn.cap.get((int) sn.nodeToStation.get(ind)) == Integer.MAX_VALUE) {
                valueNode.appendChild(simDoc.createTextNode("-1"));
            } else if (Utils.isInf(sn.cap.get((int) sn.nodeToStation.get(ind)))) {
                valueNode.appendChild(simDoc.createTextNode("-1"));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) sn.cap.get((int) sn.nodeToStation.get(ind)))));
            }
        }
        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public ElementDocumentPair saveXMLHeader(String logPath) throws ParserConfigurationException {
        String xmlnsXsi = "http://www.w3.org/2001/XMLSchema-instance";
        String fname = simFileName + ".jsimg";
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        Document simDoc = docBuilder.newDocument();
        Element simElem = simDoc.createElement("sim");
        simElem.setAttribute("xmlns:xsi", xmlnsXsi);
        simElem.setAttribute("name", fname);
        // simElem.setAttribute("timestamp", ""Tue Jan 1 00:00:01 GMT+00:00 2000"");
        simElem.setAttribute("xsi:noNamespaceSchemaLocation", "SIMmodeldefinition.xsd");
        simElem.setAttribute("disableStatisticStop", "true");
        simElem.setAttribute("logDecimalSeparator", ".");
        simElem.setAttribute("logDelimiter", ";");
        simElem.setAttribute("logPath", logPath);
        simElem.setAttribute("logReplaceMode", "0");
        simElem.setAttribute("maxSamples", String.valueOf(maxSamples));
        simElem.setAttribute("maxEvents", String.valueOf(maxEvents));

        if (!Utils.isInf(maxSimulatedTime)) {
            simElem.setAttribute("maxSimulated", String.format("%.3f", maxSimulatedTime));
        }
        simElem.setAttribute("polling", "1.0");
        simElem.setAttribute("seed", String.valueOf(seed));
        simDoc.appendChild(simElem);
        return new ElementDocumentPair(simElem, simDoc);
    }
}
