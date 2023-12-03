package jline.solvers.jmt;

import jline.api.SN;
import jline.io.SysUtils;
import jline.util.Maths;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.sections.Section;
import jline.lang.state.State;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverHandles;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.Matrix;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static java.lang.Double.*;
import static jline.lib.KPCToolbox.map_lambda;
import static jline.lib.KPCToolbox.map_pie;
import static jline.api.SN.snGetDemandsChain;
import static jline.io.SysUtils.*;

public class SolverJMT extends NetworkSolver {
    private String jmtPath;
    private String filePath;
    private String fileName;
    private double maxSimulatedTime;
    private int maxSamples;
    private int maxEvents;
    private int seed;
    private double simConfInt;
    private double simMaxRelErr;

    private SolverJMTResult jmtResult;

    public static final String XSI_NO_NAMESPACE_SCHEMA_LOCATION = "Archive.xsd";
    public static final String FILE_FORMAT = "jsimg";
    public static final String JSIMG_PATH = "";

    public SolverJMT(Network model) {
        this(model, SolverJMT.defaultOptions());
    }

    public SolverJMT(Network model, SolverOptions options) {
        super(model, "SolverJMT", options);
        this.simConfInt = 0.99;
        this.simMaxRelErr = 0.03;
        this.maxEvents = -1;
        this.jmtPath = jmtGetPath();
    }

    public SolverJMT(Network model, SolverOptions options, String jmtPath) {
        super(model, "SolverJMT", options);
        this.simConfInt = 0.99;
        this.simMaxRelErr = 0.03;
        this.maxEvents = -1;
        this.jmtPath = jmtGetPath(jmtPath);
    }

    public SolverJMT(Network model, String jmtPath) {
        this(model, SolverJMT.defaultOptions(), jmtPath);
    }

    public NetworkStruct getStruct() {
        if (this.sn == null)
            this.sn = this.model.getStruct(true);
        return this.sn;
    }

    public DocumentSectionPair saveArrivalStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        strategyNode.setAttribute("name", "ServiceStrategy");
        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);

        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode2);
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
            serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");

            if (Double.isFinite(sn.njobs.get(r))) {   // Closed
                Element subParValue = simDoc.createElement("value");
                subParValue.appendChild(simDoc.createTextNode("null"));
                serviceTimeStrategyNode.appendChild(subParValue);
            } else { // Open
                Station istStation = sn.stations.get(i);
                JobClass rstJobClass = sn.jobclasses.get(r);
                if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.DISABLED) {
                    Element subParValue = simDoc.createElement("value");
                    subParValue.appendChild(simDoc.createTextNode("null"));
                    serviceTimeStrategyNode.appendChild(subParValue);
                } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.IMMEDIATE) {
                    serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                    serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
                } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.PH
                        || sn.proctype.get(istStation).get(rstJobClass) == ProcessType.APH
                        || sn.proctype.get(istStation).get(rstJobClass) == ProcessType.COXIAN
                        || sn.proctype.get(istStation).get(rstJobClass) == ProcessType.HYPEREXP) {
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
                    Map<Integer, Matrix> PH = sn.proc.get(istStation).get(rstJobClass);
                    Matrix alpha = sn.pie.get(istStation).get(rstJobClass);
                    alpha.abs();
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
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * Math.abs(T.get(k, j)))));
                            } else {
                                subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(T.get(k, j)))));
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
                } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.MAP
                        || sn.proctype.get(istStation).get(rstJobClass) == ProcessType.MMPP2) {
                    Element distributionNode = simDoc.createElement("subParameter");
                    distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                    distributionNode.setAttribute("name", "Burst (MAP)");
                    Element distrParNode = simDoc.createElement("subParameter");
                    distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                    distrParNode.setAttribute("name", "distrPar");

                    Map<Integer, Matrix> MAP = sn.proc.get(istStation).get(rstJobClass);

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
                    Element distributionNode = simDoc.createElement("subParameter");
                    String javaClass = "";
                    String javaParClass = "";
                    switch (sn.proctype.get(istStation).get(rstJobClass)) {
                        case DET:
                            javaClass = "jmt.engine.random.DeterministicDistr";
                            javaParClass = "jmt.engine.random.DeterministicDistrPar";
                            break;
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
                    switch (sn.proctype.get(istStation).get(rstJobClass)) {
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
                            distributionNode.setAttribute("name", ProcessType.toText(sn.proctype.get(istStation).get(rstJobClass)));
                            break;
                    }
                    serviceTimeStrategyNode.appendChild(distributionNode);

                    Element distrParNode = simDoc.createElement("subParameter");
                    distrParNode.setAttribute("classPath", javaParClass);
                    distrParNode.setAttribute("name", "distrPar");

                    Element subParNodeAlpha = simDoc.createElement("subParameter");
                    Element subParValue = simDoc.createElement("value");
                    double c = Math.sqrt(sn.scv.get(i, r));

                    switch (sn.proctype.get(istStation).get(rstJobClass)) {
                        case DET:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "t");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(i, r))));
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
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).get(0, 0))));
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
                            subParValue.appendChild(simDoc.createTextNode(String.format("%d", (int) sn.phases.get(i, r))));
                            subParNodeAlpha.appendChild(subParValue);
                            distrParNode.appendChild(subParNodeAlpha);
                            break;
                        case GAMMA:
                            subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                            subParNodeAlpha.setAttribute("name", "alpha");
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1 / sn.scv.get(i, r))));
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
                            double shape = Math.sqrt(1 + 1 / sn.scv.get(i, r)) + 1;
                            double scale = 1 / sn.rates.get(i, r) * (shape - 1) / shape;
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
                            double rval = Math.pow(c, -1.086); //Justus approximation (1976)
                            double alpha = 1 / sn.rates.get(i, r) / Maths.gammaFunction(1 + 1 / rval);
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
                            double mu = Math.log(1 / sn.rates.get(i, r) / Math.sqrt(c * c + 1));
                            double sigma = Math.sqrt(Math.log(c * c + 1));
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
                            double maxVal = (Math.sqrt(12 * sn.scv.get(i, r) / Math.pow(sn.rates.get(i, r), 2)) + 2 / sn.rates.get(i, r)) / 2;
                            double minVal = 2 / sn.rates.get(i, r) - maxVal;
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
                            String istFileName = sn.nodeparam.get(sn.nodes.get(ind)).fileName;
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

    public DocumentSectionPair saveBufferCapacity(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = getStruct();
        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "size");
        Element valueNode = simDoc.createElement("value");
        int istStation = (int) sn.nodeToStation.get(ind);
        double cap;
        if (istStation >= 0) {
            cap = sn.cap.get(istStation);
        } else {
            cap = POSITIVE_INFINITY;
        }
        if (sn.isstation.get(ind) == 0 || Double.isInfinite(cap)) {
            valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
        } else {
            if (cap == Arrays.stream(sn.njobs.nz_values).sum()) {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) cap)));
            }
        }
        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);
        return new DocumentSectionPair(simDoc, section);
    }


    public DocumentSectionPair saveDropStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;

        Element schedStrategyNode = simDoc.createElement("parameter");
        schedStrategyNode.setAttribute("array", "true");
        schedStrategyNode.setAttribute("classPath", "java.lang.String");
        schedStrategyNode.setAttribute("name", "dropStrategies");
        double i = sn.nodeToStation.get(ind);
        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            schedStrategyNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.String");
            subParameterNode.setAttribute("name", "dropStrategy");

            Element valueNode2 = simDoc.createElement("value");
            if (Double.isNaN(i) || i < 0 || sn.droprule.get(sn.stations.get((int) i)).get(sn.jobclasses.get(r)) == DropStrategy.Drop) {
                valueNode2.appendChild(simDoc.createTextNode("drop"));
            } else {
                valueNode2.appendChild(simDoc.createTextNode(DropStrategy.toText(sn.droprule.get(sn.stations.get((int) i)).get(sn.jobclasses.get(r)))));
            }
            subParameterNode.appendChild(valueNode2);
            schedStrategyNode.appendChild(subParameterNode);
            section.appendChild(schedStrategyNode);
        }
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

    public DocumentSectionPair savePutStrategies(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element queuePutStrategyNode = simDoc.createElement("parameter");
        queuePutStrategyNode.setAttribute("array", "true");
        queuePutStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategy");
        queuePutStrategyNode.setAttribute("name", "QueuePutStrategy");

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        Station istStation = sn.stations.get(i);

        for (int r = 0; r < numOfClasses; r++) {
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

    public DocumentSectionPair saveNumberOfServers(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVENUMBEROFSERVERS(SIMDOC, SECTION, CURRENTNODE)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "maxJobs");

        NetworkStruct sn = getStruct();
        Element valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) sn.nservers.get((int) sn.nodeToStation.get(ind)))));

        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);

        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePreemptiveStrategy(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element visitsNode = simDoc.createElement("parameter");
        visitsNode.setAttribute("array", "true");
        visitsNode.setAttribute("classPath", "jmt.engine.NetStrategies.PSStrategy");
        visitsNode.setAttribute("name", "PSStrategy");

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        Station istStation = sn.stations.get((int) sn.nodeToStation.get(ind));

        for (int r = 0; r < numOfClasses; r++) {
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

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            visitsNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.Double");
            subParameterNode.setAttribute("name", "serviceWeight");

            Element valueNode2 = simDoc.createElement("value");
            valueNode2.appendChild(simDoc.createTextNode(String.valueOf((int) sn.schedparam.get(i, r))));

            subParameterNode.appendChild(valueNode2);
            visitsNode.appendChild(subParameterNode);
            section.appendChild(visitsNode);
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

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        for (int r = 0; r < numOfClasses; r++) {
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
                    case LCFS:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy");
                        subParameterNode2.setAttribute("name", "HeadStrategy");
                        break;
                    case LCFSPR:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategy");
                        subParameterNode2.setAttribute("name", "LCFSPRStrategy");
                        break;
                    case HOL:
                        subParameterNode2.setAttribute("classPath", "jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority");
                        subParameterNode2.setAttribute("name", "TailStrategyPriority");
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

    public DocumentSectionPair saveRoutingStrategy(DocumentSectionPair documentSectionPair, int ind) {
//        [SIMDOC, SECTION] = SAVEROUTINGSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategy");
        strategyNode.setAttribute("name", "RoutingStrategy");

        NetworkStruct sn = getStruct();
        int K = sn.nclasses;
        Node istNode = sn.nodes.get(ind);
        // since the class switch node always outputs to a single node, it is faster to translate it to RAND. Also some problems with sn.rt value otherwise.
        if (sn.nodetypes.get(ind) == NodeType.ClassSwitch) {
            for (JobClass jobClass : sn.routing.get(istNode).keySet()) {
                sn.routing.get(istNode).put(jobClass, RoutingStrategy.RAND);
            }
        }
        for (int r = 0; r < K; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode);

            Element concStratNode = simDoc.createElement("subParameter");
            Matrix conn_i;
            Matrix conn_i_find;

            RoutingStrategy routingStrategy = sn.routing.get(istNode).get(sn.jobclasses.get(r));
            switch (routingStrategy) {
                case RAND:
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
                case KCHOICES:
                    // TODO: not implemented
                    throw new RuntimeException("KCHOICES has not yet been implemented in JLNE.");
//                    concStratNode = simDoc.createElement("subParameter");
//                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy");
//                    concStratNode.setAttribute("name", "Power of k");
//                    concStratNode2 = simDoc.createElement("subParameter");
//                    concStratNode2.setAttribute("classPath", "java.lang.Integer");
//                    concStratNode2.setAttribute("name", "k");
//                    Element concStratNode2ValueNode = simDoc.createElement("value");
//                    concStratNode2ValueNode.appendChild(simDoc.createTextNode(String.format("%d", sn.nodeparam.get(istNode).k.get(sn.jobclasses.get(r)))));
//                    concStratNode2.appendChild(concStratNode2ValueNode);
//                    concStratNode3 = simDoc.createElement("subParameter");
//                    concStratNode3.setAttribute("classPath", "java.lang.Boolean");
//                    concStratNode3.setAttribute("name", "withMemory");
//                    Element concStratNode3ValueNode = simDoc.createElement("value");
//                    if (!sn.nodeparam.get(istNode).withMemory.get(sn.jobclasses.get(r)).isEmpty()) {
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
                    Element concStratNode2 = simDoc.createElement("subParameter");
                    concStratNode2.setAttribute("array", "true");
                    concStratNode2.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.WeightEntry");
                    concStratNode2.setAttribute("name", "WeightEntryArray");

                    // linked stations
                    conn_i = new Matrix(0, 0);
                    Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
                    conn_i_find = conn_i.find();
                    for (int idx = 0; idx < conn_i_find.length(); idx++) {
                        int j = (int) conn_i_find.get(idx);
                        double weight = sn.nodeparam.get(istNode).weights.get(sn.jobclasses.get(r)).get(j);

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
                    Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
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
                default:
                    concStratNode = simDoc.createElement("subParameter");
                    concStratNode.setAttribute("classPath", "jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy");
                    concStratNode.setAttribute("name", "Random");
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

        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;

        for (int r = 0; r < numOfClasses; r++) {
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

        NetworkStruct sn = this.getStruct();
        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);
        Station istStation = sn.stations.get(i);

        for (int r = 0; r < numOfClasses; r++) {
            JobClass rstJobClass = sn.jobclasses.get(r);
            Element refClassNode2 = simDoc.createElement("refClass");
            refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            strategyNode.appendChild(refClassNode2);
            Element serviceTimeStrategyNode = simDoc.createElement("subParameter");
            ProcessType type = sn.proctype.get(istStation).get(rstJobClass);
            if (type == ProcessType.DISABLED) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "DisabledServiceTimeStrategy");
            } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.IMMEDIATE) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
            } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.PH || sn.proctype.get(istStation).get(rstJobClass) == ProcessType.APH) {
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

                Map<Integer, Matrix> PH = sn.proc.get(istStation).get(rstJobClass);
                Matrix alpha = sn.pie.get(istStation).get(rstJobClass);
                alpha.abs();
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
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * Math.abs(T.get(k, j)))));
                        } else {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(T.get(k, j)))));
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
            } else if (sn.proctype.get(istStation).get(rstJobClass) == ProcessType.MAP) {
                serviceTimeStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                serviceTimeStrategyNode.setAttribute("name", "ServiceTimeStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                distributionNode.setAttribute("name", "Burst (MAP)");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                distrParNode.setAttribute("name", "distrPar");

                Map<Integer, Matrix> MAP = sn.proc.get(istStation).get(rstJobClass);

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
                switch (sn.proctype.get(istStation).get(rstJobClass)) {
                    case DET:
                        javaClass = "jmt.engine.random.DeterministicDistr";
                        javaParClass = "jmt.engine.random.DeterministicDistrPar";
                        break;
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
                switch (sn.proctype.get(istStation).get(rstJobClass)) {
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
                        distributionNode.setAttribute("name", ProcessType.toText(sn.proctype.get(istStation).get(rstJobClass)));
                        break;
                }
                serviceTimeStrategyNode.appendChild(distributionNode);

                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");

                Element subParNodeAlpha = simDoc.createElement("subParameter");
                Element subParValue = simDoc.createElement("value");
                double c;
                switch (sn.proctype.get(istStation).get(rstJobClass)) {
                    case DET:
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "t");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.rates.get(i, r))));
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
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * sn.proc.get(istStation).get(rstJobClass).get(0).get(0, 0))));
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
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", sn.proc.get(istStation).get(rstJobClass).get(1).get(0, 0))));
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
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", 1 / sn.scv.get(i, r))));
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
                        double shape = Math.sqrt(1 + 1 / sn.scv.get(i, r)) + 1;
                        double scale = 1 / sn.rates.get(i, r) * (shape - 1) / shape;
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
                        c = Math.sqrt(sn.scv.get(i, r));
                        double rval = Math.pow(c, -1.086); //Justus approximation (1976)
                        double alpha = 1 / sn.rates.get(i, r) / Maths.gammaFunction(1 + 1 / rval);
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
                        c = Math.sqrt(sn.scv.get(i, r));
                        double mu = Math.log(1 / sn.rates.get((int) i, r) / Math.sqrt(c * c + 1));
                        double sigma = Math.sqrt(Math.log(c * c + 1));
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
                        double maxVal = (Math.sqrt(12 * sn.scv.get(i, r) / Math.pow(sn.rates.get(i, r), 2)) + 2 / sn.rates.get(i, r)) / 2;
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
                        subParValue.appendChild(simDoc.createTextNode(sn.nodeparam.get(sn.nodes.get(i)).fileName));
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

    public DocumentSectionPair saveClassSwitchStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVECLASSSWITCHSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element paramNode = simDoc.createElement("parameter");
        paramNode.setAttribute("array", "true");
        paramNode.setAttribute("classPath", "java.lang.Object");
        paramNode.setAttribute("name", "matrix");

        NetworkStruct sn = getStruct();
        int K = sn.nclasses;
        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
        Matrix conn_i_find = conn_i.find();
        int j = (int) conn_i_find.get(0);
        for (int r = 0; r < K; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            paramNode.appendChild(refClassNode);

            Element subParNodeRow = simDoc.createElement("subParameter");
            subParNodeRow.setAttribute("array", "true");
            subParNodeRow.setAttribute("classPath", "java.lang.Float");
            subParNodeRow.setAttribute("name", "row");
            for (int s = 0; s < K; s++) {
                refClassNode = simDoc.createElement("refClass");
                refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(s)));
                subParNodeRow.appendChild(refClassNode);

                Element subParNodeCell = simDoc.createElement("subParameter");
                subParNodeCell.setAttribute("classPath", "java.lang.Float");
                subParNodeCell.setAttribute("name", "cell");
                Element valNode = simDoc.createElement("value");
                valNode.appendChild(simDoc.createTextNode(String.format("%12.12f", sn.rtnodes.get(ind * K + r, j * K + s))));
                subParNodeCell.appendChild(valNode);
                subParNodeRow.appendChild(subParNodeCell);
            }
            paramNode.appendChild(subParNodeRow);
        }
        section.appendChild(paramNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveLogTunnel(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVELOGTUNNEL(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        List<String> loggerNodesCP = new ArrayList<>();
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
        if (!sn.nodeparam.get(istNode).filePath.endsWith(File.separator)) {
            sn.nodeparam.get(istNode).filePath = sn.nodeparam.get(istNode).filePath + File.separator;
        }

        List<String> loggerNodesValues = new ArrayList<>();
        loggerNodesValues.add(sn.nodeparam.get(istNode).fileName);
        loggerNodesValues.add(sn.nodeparam.get(istNode).filePath);
        loggerNodesValues.add(sn.nodeparam.get(istNode).startTime);
        loggerNodesValues.add(sn.nodeparam.get(istNode).loggerName);
        loggerNodesValues.add(sn.nodeparam.get(istNode).timestamp);
        loggerNodesValues.add(sn.nodeparam.get(istNode).jobID);
        loggerNodesValues.add(sn.nodeparam.get(istNode).jobClass);
        loggerNodesValues.add(sn.nodeparam.get(istNode).timeSameClass);
        loggerNodesValues.add(sn.nodeparam.get(istNode).timeAnyClass);
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

    public DocumentSectionPair saveForkStrategy(DocumentSectionPair documentSectionPair, int ind) {
        // [SIMDOC, SECTION] = SAVEFORKSTRATEGY(SIMDOC, SECTION, NODEIDX)
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);

        Element jplNode = simDoc.createElement("parameter");
        jplNode.setAttribute("classPath", "java.lang.Integer");
        jplNode.setAttribute("name", "jobsPerLink");
        Element valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) sn.nodeparam.get(istNode).fanOut)));
        jplNode.appendChild(valueNode);
        section.appendChild(jplNode);

        Element blockNode = simDoc.createElement("parameter");
        blockNode.setAttribute("classPath", "java.lang.Integer");
        blockNode.setAttribute("name", "block");
        valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode("-1"));
        blockNode.appendChild(valueNode);
        section.appendChild(blockNode);

        Element issimplNode = simDoc.createElement("parameter");
        issimplNode.setAttribute("classPath", "java.lang.Boolean");
        issimplNode.setAttribute("name", "isSimplifiedFork");
        valueNode = simDoc.createElement("value");
        valueNode.appendChild(simDoc.createTextNode("true"));
        issimplNode.appendChild(valueNode);
        section.appendChild(issimplNode);

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategy");
        strategyNode.setAttribute("name", "ForkStrategy");

        int numOfClasses = sn.nclasses;
        for (int r = 0; r < numOfClasses; r++) {
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
                    Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
                    Matrix conn_i_find = conn_i.find();

                    Element classStratNode3 = simDoc.createElement("subParameter");
                    Element classStratNode4 = simDoc.createElement("subParameter");
                    Element classStratNode4Station = simDoc.createElement("subParameter");
                    Element classStratNode4StationValueNode = simDoc.createElement("value");

                    for (int idx = 0; idx < conn_i_find.length(); idx++) {
                        int k = (int) conn_i_find.get(idx);
                        classStratNode3 = simDoc.createElement("subParameter");
                        classStratNode3.setAttribute("classPath", "jmt.engine.NetStrategies.ForkStrategies.OutPath");
                        classStratNode3.setAttribute("name", "OutPathEntry");
                        classStratNode4 = simDoc.createElement("subParameter");
                        classStratNode4.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                        classStratNode4.setAttribute("name", "outUnitProbability");
                        classStratNode4Station = simDoc.createElement("subParameter");
                        classStratNode4Station.setAttribute("classPath", "java.lang.String");
                        classStratNode4Station.setAttribute("name", "stationName");
                        classStratNode4StationValueNode = simDoc.createElement("value");
                        classStratNode4StationValueNode.appendChild(simDoc.createTextNode(String.format("%s", sn.nodenames.get(k))));
                    }
                    classStratNode4Station.appendChild(classStratNode4StationValueNode);
                    classStratNode3.appendChild(classStratNode4Station);
                    Element classStratNode4Probability = simDoc.createElement("subParameter");
                    classStratNode4Probability.setAttribute("classPath", "java.lang.Double");
                    classStratNode4Probability.setAttribute("name", "probability");
                    Element classStratNode4ProbabilityValueNode = simDoc.createElement("value");
                    classStratNode4ProbabilityValueNode.appendChild(simDoc.createTextNode("1.0"));
                    classStratNode4Probability.appendChild(classStratNode4ProbabilityValueNode);

                    Element classStratNode4b = simDoc.createElement("subParameter");
                    classStratNode4b.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                    classStratNode4b.setAttribute("array", "true");
                    classStratNode4b.setAttribute("name", "JobsPerLinkDis");
                    Element classStratNode5b = simDoc.createElement("subParameter");
                    classStratNode5b.setAttribute("classPath", "jmt.engine.random.EmpiricalEntry");
                    classStratNode5b.setAttribute("name", "EmpiricalEntry");
                    Element classStratNode5bStation = simDoc.createElement("subParameter");
                    classStratNode5bStation.setAttribute("classPath", "java.lang.String");
                    classStratNode5bStation.setAttribute("name", "numbers");
                    Element classStratNode5bStationValueNode = simDoc.createElement("value");
                    classStratNode5bStationValueNode.appendChild(simDoc.createTextNode(String.valueOf((int) sn.nodeparam.get(istNode).fanOut)));
                    classStratNode5bStation.appendChild(classStratNode5bStationValueNode);
                    classStratNode4b.appendChild(classStratNode5bStation);
                    Element classStratNode5bProbability = simDoc.createElement("subParameter");
                    classStratNode5bProbability.setAttribute("classPath", "java.lang.Double");
                    classStratNode5bProbability.setAttribute("name", "probability");
                    Element classStratNode5bProbabilityValueNode = simDoc.createElement("value");
                    classStratNode5bProbabilityValueNode.appendChild(simDoc.createTextNode("1.0"));
                    classStratNode5bProbability.appendChild(classStratNode5bProbabilityValueNode);

                    classStratNode4.appendChild(classStratNode4Station);
                    classStratNode4.appendChild(classStratNode4Probability);
                    classStratNode3.appendChild(classStratNode4);
                    classStratNode5b.appendChild(classStratNode5bStation);
                    classStratNode5b.appendChild(classStratNode5bProbability);
                    classStratNode4b.appendChild(classStratNode5b);
                    classStratNode3.appendChild(classStratNode4b);
                    classStratNode2.appendChild(classStratNode3);
            }
            classStratNode.appendChild(classStratNode2);
            strategyNode.appendChild(classStratNode);
        }
        section.appendChild(strategyNode);
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

        NetworkStruct sn = this.getStruct();
        int numOfClasses = sn.nclasses;

        Element refClassNode2 = simDoc.createElement("refClass");
        for (int r = 0; r < numOfClasses; r++) {
            Node istNode = sn.nodes.get(ind);
            JobClass rstJobClass = sn.jobclasses.get(r);
            switch (sn.nodeparam.get(istNode).joinStrategy.get(rstJobClass)) {
                case STD:
                    refClassNode2 = simDoc.createElement("refClass");
                    refClassNode2.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    strategyNode.appendChild(refClassNode2);

                    Element joinStrategyNode = simDoc.createElement("subParameter");
                    joinStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.JoinStrategies.NormalJoin");
                    joinStrategyNode.setAttribute("name", "Standard Join");
                    Element reqNode = simDoc.createElement("subParameter");
                    reqNode.setAttribute("classPath", "java.lang.Integer");
                    reqNode.setAttribute("name", "numRequired");
                    Element valueNode = simDoc.createElement("value");
                    valueNode.appendChild(simDoc.createTextNode(String.valueOf(sn.nodeparam.get(istNode).fanIn.get(rstJobClass).intValue())));
                    reqNode.appendChild(valueNode);
                    joinStrategyNode.appendChild(reqNode);
                    strategyNode.appendChild(joinStrategyNode);
                    section.appendChild(strategyNode);
                    break;
                case Quorum:
                case Guard:
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
                    valueNode.appendChild(simDoc.createTextNode(String.valueOf(sn.nodeparam.get(istNode).joinRequired.get(rstJobClass))));
                    reqNode.appendChild(valueNode);
                    joinStrategyNode.appendChild(reqNode);
                    strategyNode.appendChild(joinStrategyNode);
                    section.appendChild(strategyNode);
                    break;

            }
        }
        return new DocumentSectionPair(simDoc, section);
    }


    public ElementDocumentPair saveClasses(ElementDocumentPair elementDocumentPair) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;

        // [SIMELEM, SIMDOC] = SAVECLASSES(SIMELEM, SIMDOC)
        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;

        for (int r = 0; r < numOfClasses; r++) {
            Element userClass = simDoc.createElement("userClass");
            userClass.setAttribute("name", sn.classnames.get(r));
            if (Double.isInfinite(sn.njobs.get(r))) {
                userClass.setAttribute("type", "open");
            } else {
                userClass.setAttribute("type", "closed");
            }
            userClass.setAttribute("priority", String.valueOf((int) sn.classprio.get(r)));

            int refStatIndex = (int) sn.refstat.get(r);
            Map<Integer, Matrix> integerMatrixMap = sn.proc.get(sn.stations.get(refStatIndex)).get(sn.jobclasses.get(r));

            if (!integerMatrixMap.isEmpty()) {
                if (Double.isFinite(sn.njobs.get(r))) {
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

    public ElementDocumentPair saveLinks(ElementDocumentPair elementDocumentPair) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;

        NetworkStruct sn = getStruct();

        for (int j = 0; j < sn.connmatrix.numCols; j++) {
            for (int i = 0; i < sn.connmatrix.numRows; i++) {
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

    public ElementDocumentPair saveRegions(ElementDocumentPair elementDocumentPair) {
        Document simDoc = elementDocumentPair.simDoc;
        Element simElem = elementDocumentPair.simElem;

        NetworkStruct sn = getStruct();

        for (int r = 0; r < model.getRegions().size(); r++) {
            FiniteCapacityRegion region = model.getRegions().get(r);
            Element blockingRegion = simDoc.createElement("blockingRegion");
            blockingRegion.setAttribute("name", "FCRegion" + (r + 1));
            blockingRegion.setAttribute("type", "default");

            for (Node node : region.getNodes()) {
                Element regionNode = simDoc.createElement("regionNode");
                regionNode.setAttribute("nodeName", node.getName());
                blockingRegion.appendChild(regionNode);
            }

            Element globalConstraint = simDoc.createElement("globalConstraint");
            globalConstraint.setAttribute("maxJobs", String.valueOf(region.getGlobalMaxJobs()));
            blockingRegion.appendChild(globalConstraint);

            Element globalMemoryConstraint = simDoc.createElement("globalMemoryConstraint");
            globalMemoryConstraint.setAttribute("maxMemory", String.valueOf(region.getGlobalMaxMemory()));
            blockingRegion.appendChild(globalMemoryConstraint);

            for (int c = 0; c < sn.nclasses; c++) {
                JobClass cstJobClass = sn.jobclasses.get(c);
                if (region.classMaxJobs.get(cstJobClass) != FiniteCapacityRegion.UNBOUNDED) {
                    Element classConstraint = simDoc.createElement("classConstraint");
                    classConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    classConstraint.setAttribute("maxJobsPerClass", String.valueOf(region.classMaxJobs.get(cstJobClass)));
                    blockingRegion.appendChild(classConstraint);
                }

                if (region.classMaxMemory.get(cstJobClass) != FiniteCapacityRegion.UNBOUNDED) {
                    Element classMemoryConstraint = simDoc.createElement("classMemoryConstraint");
                    classMemoryConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    classMemoryConstraint.setAttribute("maxMemoryPerClass", String.valueOf(region.classMaxMemory.get(cstJobClass)));
                    blockingRegion.appendChild(classMemoryConstraint);
                }

                if (region.dropRule.get(cstJobClass)) {
                    Element dropRule = simDoc.createElement("dropRules");
                    dropRule.setAttribute("jobClass", region.classes.get(c).getName());
                    dropRule.setAttribute("dropThisClass", "true");
                    blockingRegion.appendChild(dropRule);
                }

                if (region.classSize.get(cstJobClass) != 1) {
                    Element classMemoryConstraint = simDoc.createElement("classSize");
                    classMemoryConstraint.setAttribute("jobClass", region.classes.get(c).getName());
                    classMemoryConstraint.setAttribute("size", String.valueOf(region.classSize.get(cstJobClass)));
                    blockingRegion.appendChild(classMemoryConstraint);
                }
            }
            simElem.appendChild(blockingRegion);
        }

        return new ElementDocumentPair(simElem, simDoc);
    }

    public ElementDocumentPair saveMetric(ElementDocumentPair elementDocumentPair, Map<Station, Map<JobClass, SolverHandles.Metric>> handles) {
        Element simElem = elementDocumentPair.simElem;
        Document simDoc = elementDocumentPair.simDoc;

        for (int i = 0; i < handles.keySet().size(); i++) {
            Station istStation = sn.stations.get(i);
            for (int r = 0; r < handles.get(istStation).keySet().size(); r++) {
                JobClass rstJobClass = sn.jobclasses.get(r);
                SolverHandles.Metric currentPerformanceIndex = handles.get(istStation).get(rstJobClass);
                if (!currentPerformanceIndex.isDisabled) {
                    Element performanceNode = simDoc.createElement("measure");
                    performanceNode.setAttribute("alpha", String.format("%.2f", 1 - this.simConfInt));
                    performanceNode.setAttribute("name", "Performance_" + (i + 1));
                    performanceNode.setAttribute("nodeType", "station");
                    performanceNode.setAttribute("precision", String.format("%.2f", this.simMaxRelErr));

                    if (currentPerformanceIndex.station == null) {
                        performanceNode.setAttribute("referenceNode", "");
                    } else {
                        performanceNode.setAttribute("referenceNode", currentPerformanceIndex.station.getName());
                    }
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
        if (this.handles == null) {
            this.getAvgHandles();
        }
        SolverHandles handles = this.handles;

        ElementDocumentPair res = saveMetric(elementDocumentPair, handles.Q);
        res = saveMetric(res, handles.U);
        res = saveMetric(res, handles.R);
        res = saveMetric(res, handles.T);
        res = saveMetric(res, handles.A);

        // JMT ResidT is inconsistently defined with LINE"s on some
        // difficult class switching cases, hence we recompute it at the
        // level of the NetworkSolver class to preserve consistency.
        return res;
    }

    public ElementDocumentPair saveXMLHeader(String logPath) throws ParserConfigurationException {
        String xmlnsXsi = "http://www.w3.org/2001/XMLSchema-instance";
        String fname = getFileName() + ".jsimg";
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
        simElem.setAttribute("maxSamples", String.valueOf(this.maxSamples));
        simElem.setAttribute("maxEvents", String.valueOf(this.maxEvents));

        if (!Double.isInfinite(this.maxSimulatedTime)) {
            simElem.setAttribute("maxSimulated", String.format("%.3f", this.maxSimulatedTime));
        }
        simElem.setAttribute("polling", "1.0");
        simElem.setAttribute("seed", String.valueOf(this.options.seed));
        simDoc.appendChild(simElem);
        return new ElementDocumentPair(simElem, simDoc);
    }

    public DocumentSectionPair saveTotalCapacity(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = getStruct();
        Element sizeNode = simDoc.createElement("parameter");
        sizeNode.setAttribute("classPath", "java.lang.Integer");
        sizeNode.setAttribute("name", "totalCapacity");
        Element valueNode = simDoc.createElement("value");
        if (sn.isstation.get(ind) == 0 || Double.isInfinite(sn.cap.get((int) sn.nodeToStation.get(ind)))) {
            valueNode.appendChild(simDoc.createTextNode("-1"));
        } else {
            if (Double.isInfinite(sn.cap.get((int) sn.nodeToStation.get(ind)))) {
                valueNode.appendChild(simDoc.createTextNode("-1"));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf((int) sn.cap.get((int) sn.nodeToStation.get(ind)))));
            }
        }
        sizeNode.appendChild(valueNode);
        section.appendChild(sizeNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair savePlaceCapacities(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element placeCapacityNode = simDoc.createElement("parameter");
        placeCapacityNode.setAttribute("array", "true");
        placeCapacityNode.setAttribute("classPath", "java.lang.Integer");
        placeCapacityNode.setAttribute("name", "capacities");
        NetworkStruct sn = getStruct();
        int numOfClasses = sn.nclasses;
        double i = sn.nodeToStation.get(ind);
        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            placeCapacityNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.Integer");
            subParameterNode.setAttribute("name", "capacity");

            Element valueNode2 = simDoc.createElement("value");
            if (Double.isInfinite(sn.classcap.get((int) i, r))) {
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

    public DocumentSectionPair saveDropRule(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        NetworkStruct sn = this.getStruct();
        Element schedStrategyNode = simDoc.createElement("parameter");
        schedStrategyNode.setAttribute("array", "true");
        schedStrategyNode.setAttribute("classPath", "java.lang.String");
        schedStrategyNode.setAttribute("name", "droprules");

        int numOfClasses = sn.nclasses;
        int i = (int) sn.nodeToStation.get(ind);

        for (int r = 0; r < numOfClasses; r++) {
            Element refClassNode = simDoc.createElement("refClass");
            refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
            schedStrategyNode.appendChild(refClassNode);

            Element subParameterNode = simDoc.createElement("subParameter");
            subParameterNode.setAttribute("classPath", "java.lang.String");
            subParameterNode.setAttribute("name", "droprule");

            Element valueNode2 = simDoc.createElement("value");

            valueNode2.appendChild(simDoc.createTextNode(DropStrategy.toText(sn.droprule.get(sn.stations.get(i)).get(sn.jobclasses.get(r)))));
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

        NetworkStruct sn = getStruct();
        int numOfNodes = sn.nnodes;
        int numOfClasses = sn.nclasses;
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subEnablingConditionNode = simDoc.createElement("subParameter");
            subEnablingConditionNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
            subEnablingConditionNode.setAttribute("name", "enablingCondition");

            Element subEnablingVectorsNode = simDoc.createElement("subParameter");
            subEnablingVectorsNode.setAttribute("array", "true");
            subEnablingVectorsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionVector");
            subEnablingVectorsNode.setAttribute("name", "enablingVectors");

            for (int k = 0; k < numOfNodes; k++) {
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

                boolean exists = false;

                for (int r = 0; r < numOfClasses; r++) {
                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(sn.classnames.get(r)));
                    subEnablingEntriesNode.appendChild(refClassNode);

                    Element subParameterNode = simDoc.createElement("subParameter");
                    subParameterNode.setAttribute("classPath", "java.lang.Integer");
                    subParameterNode.setAttribute("name", "enablingEntry");

                    Element valueNode2 = simDoc.createElement("value");

                    if (Double.isInfinite(sn.nodeparam.get(istNode).enabling.get(m).get(k, r))) {
                        valueNode2.appendChild(simDoc.createTextNode("-1"));
                        exists = true;
                    } else if (sn.nodeparam.get(istNode).enabling.get(m).get(k, r) > 0) {
                        valueNode2.appendChild(simDoc.createTextNode(String.valueOf(sn.nodeparam.get(istNode).enabling.get(m).get(k, r))));
                        exists = true;
                    } else if (Double.isFinite(sn.nodeparam.get(istNode).inhibiting.get(m).get(k, r)) && sn.nodeparam.get(istNode).inhibiting.get(m).get(k, r) > 0) {
                        valueNode2.appendChild(simDoc.createTextNode("0"));
                        exists = true;
                    }

                    subParameterNode.appendChild(valueNode2);
                    subEnablingEntriesNode.appendChild(subParameterNode);
                    subEnablingVectorNode.appendChild(subEnablingEntriesNode);

                }
                if (exists) {
                    subEnablingVectorsNode.appendChild(subEnablingVectorNode);
                }
            }
            subEnablingConditionNode.appendChild(subEnablingVectorsNode);
            enablingNode.appendChild(subEnablingConditionNode);
        }
        section.appendChild(enablingNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveInhibitingConditions(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element inhibitingConditionsNode = simDoc.createElement("parameter");
        inhibitingConditionsNode.setAttribute("array", "true");
        inhibitingConditionsNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
        inhibitingConditionsNode.setAttribute("name", "inhibitingConditions");

        NetworkStruct sn = getStruct();

        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractColumn(this.sn.connmatrix, ind, conn_i);
        Matrix inputs = conn_i.find();
        List<String> connections = new ArrayList<>();
        for (int idx = 0; idx < inputs.length(); idx++) {
            int i = (int) inputs.get(idx);
            connections.add(sn.nodenames.get((int) inputs.get(i)));
        }
        int numOfInputs = connections.size();

        int numOfClasses = sn.nclasses;
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;

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

                    if (Double.isInfinite(sn.nodeparam.get(istNode).inhibiting.get(m).get(k, r))) {
                        valueNode2.appendChild(simDoc.createTextNode("-1"));
                    } else {
                        valueNode2.appendChild(simDoc.createTextNode(String.valueOf(sn.nodeparam.get(istNode).inhibiting.get(m).get((int) inputs.get(k), r))));
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

    public DocumentSectionPair saveFiringOutcomes(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element firingOutcomesNode = simDoc.createElement("parameter");
        firingOutcomesNode.setAttribute("array", "true");
        firingOutcomesNode.setAttribute("classPath", "jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix");
        firingOutcomesNode.setAttribute("name", "firingOutcomes");

        NetworkStruct sn = getStruct();

        Matrix conn_i = new Matrix(0, 0);
        Matrix.extractRows(this.sn.connmatrix, ind, ind + 1, conn_i);
        Matrix outputs = conn_i.find();
        List<String> connections = new ArrayList<>();
        for (int idx = 0; idx < outputs.length(); idx++) {
            int i = (int) outputs.get(idx);
            connections.add(sn.nodenames.get((int) outputs.get(i)));
        }
        int numOfOutput = connections.size();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;
        int numOfClasses = this.model.getClasses().size();

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
                    JobClass currentClass = this.model.getClasses().get(j);

                    Element refClassNode = simDoc.createElement("refClass");
                    refClassNode.appendChild(simDoc.createTextNode(currentClass.getName()));
                    subFiringEntriesNode.appendChild(refClassNode);

                    Element subParameterNode = simDoc.createElement("subParameter");
                    subParameterNode.setAttribute("classPath", "java.lang.Integer");
                    subParameterNode.setAttribute("name", "firingEntry");

                    Element valueNode2 = simDoc.createElement("value");
                    valueNode2.appendChild(simDoc.createTextNode(String.valueOf(sn.nodeparam.get(istNode).firing.get(m).get((int) outputs.get(k), j))));

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

    public DocumentSectionPair saveModeNames(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element modeNamesNode = simDoc.createElement("parameter");
        modeNamesNode.setAttribute("classPath", "java.lang.String");
        modeNamesNode.setAttribute("name", "modeNames");
        modeNamesNode.setAttribute("array", "true");

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;
        for (int m = 0; m < numOfModes; m++) {
            Element subModeNameNode = simDoc.createElement("subParameter");
            subModeNameNode.setAttribute("classPath", "java.lang.String");
            subModeNameNode.setAttribute("name", "modeName");

            Element valueNode = simDoc.createElement("value");
            valueNode.appendChild(simDoc.createTextNode(sn.nodeparam.get(istNode).modenames.get(m)));

            subModeNameNode.appendChild(valueNode);
            modeNamesNode.appendChild(subModeNameNode);
        }
        section.appendChild(modeNamesNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveNumbersOfServers(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element numbersOfServersNode = simDoc.createElement("parameter");
        numbersOfServersNode.setAttribute("classPath", "java.lang.Integer");
        numbersOfServersNode.setAttribute("name", "numbersOfServers");
        numbersOfServersNode.setAttribute("array", "true");

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;
        for (int m = 0; m < numOfModes; m++) {
            Element subNumberOfServersNode = simDoc.createElement("subParameter");
            subNumberOfServersNode.setAttribute("classPath", "java.lang.Integer");
            subNumberOfServersNode.setAttribute("name", "numberOfServers");

            Element valueNode = simDoc.createElement("value");
            double nmodeservers = sn.nodeparam.get(istNode).nmodeservers.get(m);

            if (Double.isInfinite(nmodeservers)) {
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

    public DocumentSectionPair saveTimingStrategies(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element strategyNode = simDoc.createElement("parameter");
        strategyNode.setAttribute("array", "true");
        strategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategy");
        strategyNode.setAttribute("name", "timingStrategies");

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            int nphases = sn.nodeparam.get(sn.nodes.get(ind)).firingphases.get(m);
            Element timimgStrategyNode = simDoc.createElement("subParameter");
            JobClass mstJobClass = sn.jobclasses.get(m);

            if (sn.nodeparam.get(istNode).firingid.get(mstJobClass) == TimingStrategy.IMMEDIATE) {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "ZeroServiceTimeStrategy");
            } else if (sn.nodeparam.get(istNode).firingprocid.get(mstJobClass) == ProcessType.APH || sn.nodeparam.get(istNode).firingprocid.get(mstJobClass) == ProcessType.COXIAN || (sn.nodeparam.get(istNode).firingphases.get(m) > 2 && sn.nodeparam.get(istNode).firingprocid.get(mstJobClass) == ProcessType.HYPEREXP)) {
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
                Map<Integer, Matrix> PH = sn.nodeparam.get(sn.nodes.get(ind)).firingproc.get(mstJobClass);
                Matrix alpha = map_pie(PH.get(0), PH.get(1));
                alpha.abs();
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
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * Math.abs(T.get(k, j)))));
                        } else {
                            subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", Math.abs(T.get(k, j)))));
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
            } else if (sn.nodeparam.get(istNode).firingprocid.get(mstJobClass) == ProcessType.MAP) {
                timimgStrategyNode.setAttribute("classPath", "jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy");
                timimgStrategyNode.setAttribute("name", "timingStrategy");
                Element distributionNode = simDoc.createElement("subParameter");
                distributionNode.setAttribute("classPath", "jmt.engine.random.MAPDistr");
                distributionNode.setAttribute("name", "Burst (MAP)");
                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", "jmt.engine.random.MAPPar");
                distrParNode.setAttribute("name", "distrPar");

                Map<Integer, Matrix> MAP = sn.nodeparam.get(istNode).firingproc.get(mstJobClass);

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
                switch (sn.nodeparam.get(istNode).firingprocid.get(mstJobClass)) {
                    case DET:
                        javaClass = "jmt.engine.random.DeterministicDistr";
                        javaParClass = "jmt.engine.random.DeterministicDistrPar";
                        break;
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
                switch (sn.nodeparam.get(istNode).firingprocid.get(mstJobClass)) {
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
                        distributionNode.setAttribute("name", ProcessType.toText(sn.nodeparam.get(istNode).firingprocid.get(mstJobClass)));
                }
                timimgStrategyNode.appendChild(distributionNode);

                Element distrParNode = simDoc.createElement("subParameter");
                distrParNode.setAttribute("classPath", javaParClass);
                distrParNode.setAttribute("name", "distrPar");
                Element subParNodeAlpha = simDoc.createElement("subParameter");
                Element subParValue = simDoc.createElement("value");

                Map<Integer, Matrix> matrixMap = sn.nodeparam.get(istNode).firingproc.get(mstJobClass);

                switch (sn.nodeparam.get(istNode).firingprocid.get(mstJobClass)) {
                    case DET:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "t");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", map_lambda(matrixMap.get(0), matrixMap.get(1)))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case EXP:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", map_lambda(matrixMap.get(0), matrixMap.get(1)))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        break;
                    case HYPEREXP:
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "p");
                        subParValue = simDoc.createElement("value");
                        Matrix pie = map_pie(matrixMap.get(0), matrixMap.get(1));
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", pie.get(0))));
                        subParNodeAlpha.appendChild(subParValue);
                        distrParNode.appendChild(subParNodeAlpha);
                        subParNodeAlpha = simDoc.createElement("subParameter");
                        subParNodeAlpha.setAttribute("classPath", "java.lang.Double");
                        subParNodeAlpha.setAttribute("name", "lambda1");
                        subParValue = simDoc.createElement("value");
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", -1 * matrixMap.get(0).get(0, 0))));
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
                        double timingrate = map_lambda(matrixMap.get(0), matrixMap.get(1));
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
                        subParValue.appendChild(simDoc.createTextNode(String.format("%.12f", matrixMap.get(1).get(0, 0))));
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
                    case GAMMA:
                    case PARETO:
                    case WEIBULL:
                    case LOGNORMAL:
                    case UNIFORM:
                    case REPLAYER:
                    case TRACE:
                    default:
                        // TODO: not implemented
                        throw new RuntimeException(String.format("Unsupported firing distribution for mode %d", m));
                }
                timimgStrategyNode.appendChild(distrParNode);
            }
            strategyNode.appendChild(timimgStrategyNode);
        }
        section.appendChild(strategyNode);
        return new DocumentSectionPair(simDoc, section);
    }

    public DocumentSectionPair saveFiringPriorities(DocumentSectionPair documentSectionPair, int ind) {
        Document simDoc = documentSectionPair.simDoc;
        Element section = documentSectionPair.section;

        Element firingPrioritiesNode = simDoc.createElement("parameter");
        firingPrioritiesNode.setAttribute("classPath", "java.lang.Integer");
        firingPrioritiesNode.setAttribute("name", "firingPriorities");
        firingPrioritiesNode.setAttribute("array", "true");

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subFiringPriorityNode = simDoc.createElement("subParameter");
            subFiringPriorityNode.setAttribute("classPath", "java.lang.Integer");
            subFiringPriorityNode.setAttribute("name", "firingPriority");

            Element valueNode = simDoc.createElement("value");
            double firingPrio = sn.nodeparam.get(istNode).firingprio.get(m);
            if (Double.isInfinite(firingPrio)) {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(-1)));
            } else {
                valueNode.appendChild(simDoc.createTextNode(String.valueOf(firingPrio)));
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

        NetworkStruct sn = getStruct();
        Node istNode = sn.nodes.get(ind);
        int numOfModes = sn.nodeparam.get(istNode).nmodes;

        for (int m = 0; m < numOfModes; m++) {
            Element subFiringWeightNode = simDoc.createElement("subParameter");
            subFiringWeightNode.setAttribute("classPath", "java.lang.Double");
            subFiringWeightNode.setAttribute("name", "firingWeight");

            Element valueNode = simDoc.createElement("value");
            double firingWeights = sn.nodeparam.get(istNode).fireweight.get(m);

            if (Double.isInfinite(firingWeights)) {
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

    public String getFileName() {
        return fileName;
    }

    public void setJmtJarPath(String path) {
        this.jmtPath = path;
    }

    public String getJmtJarPath() {
        return jmtPath;
    }

    public String getFilePath() {
        return filePath;
    }

    public enum ViewMode {
        JSIMW,
        JSIMG
    }

    public static void viewModel(String filename, ViewMode viewMode) {
        viewModel(jmtGetPath(), filename, viewMode);
    }

    public static void viewModel(String jmtPath, String filename, ViewMode viewMode) {
        Path path = Paths.get(filename).getParent();
        if (path == null) {
            filename = Paths.get(java.lang.System.getProperty("user.dir"), filename).toString();
        }

        String redirectOutput = " > /dev/null";
        if (java.lang.System.getProperty("os.name").startsWith("Windows")) {
            redirectOutput = " > nul 2>&1";
        }

        String cmd = String.format(
                "java -cp %s jmt.commandline.Jmt %s %s %s",
                jmtPath, viewMode.toString().toLowerCase(), filename, redirectOutput
        );

        java.lang.System.out.println("JMT view model command: " + cmd);
        String output = SysUtils.system(cmd);
        java.lang.System.out.println("JMT view model command output: " + output);
    }

    public void jsimwView(String jmtPath) {
        try {
            jsimwView(jmtPath, SolverJMT.defaultOptions());
        } catch (ParserConfigurationException e) {
            throw new RuntimeException(e);
        }
    }

    public void jsimwView(String jmtPath, SolverOptions options) throws ParserConfigurationException {
        if (this.enableChecks && !supports(this.model)) {
            throw new RuntimeException("This model contains features not supported by the solver.");
        }

        if (options.samples < 5000) {
            java.lang.System.err.println("JMT requires at least 5000 samples for each metric. Setting the samples to 5000.");
            options.samples = 5000;
        }

        this.seed = options.seed;
        this.maxSamples = options.samples;

        NetworkStruct sn = getStruct();
        this.writeJSIM(sn);

        String fileName = this.getFilePath() + File.separator + this.getFileName() + ".jsim";
        java.lang.System.out.println("JMT Model: " + fileName);

        viewModel(jmtPath, fileName, ViewMode.JSIMW);
    }


    public void jsimgView() {
        jsimgView(jmtGetPath(), this.options);
    }

    public void jsimgView(SolverOptions options) {
        jsimgView(jmtGetPath(), options);
    }

    public void jsimgView(String jmtPath, SolverOptions options) {
//        if (this.enableChecks && !supports(this.model)) {
//            throw new RuntimeException("This model contains features not supported by the solver.");
//        }

        if (options == null) {
            options = defaultOptions();
        }

        if (options.samples == 0) {
            options.samples = 10000;
        } else if (options.samples < 5000) {
            // line_warning
            java.lang.System.err.println("JMT requires at least 5000 samples for each metric. Setting the samples to 5000.");
            options.samples = 5000;
        }

        // set seed and maxSamples
        this.seed = options.seed;
        this.maxSamples = options.samples;

        NetworkStruct sn = getStruct();
        try {
            this.writeJSIM(sn);
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        }

        String fileName = this.getFilePath() + File.separator + this.getFileName() + ".jsim";
        java.lang.System.out.println("JMT Model: " + fileName);

        viewModel(jmtPath, fileName, ViewMode.JSIMG);
    }

    public void jsimwView() throws ParserConfigurationException {
        jsimwView(jmtGetPath(), this.options);
    }

    public void jsimgView(String jmtPath) {
        jsimgView(jmtPath, SolverJMT.defaultOptions());
    }

    public String writeJSIM(NetworkStruct sn, String outputFileName) throws ParserConfigurationException {
        ElementDocumentPair xml = saveXMLHeader(this.model.getLogPath());
        xml = saveClasses(xml);
        int numOfClasses = sn.nclasses;
        int numOfNodes = sn.nnodes;

        for (int i = 0; i < numOfNodes; i++) {
            Node currentNode = this.model.getNodes().get(i);
            Element node = xml.simDoc.createElement("node");
            node.setAttribute("name", currentNode.getName());
            List<Section> nodeSections = Arrays.asList(currentNode.getInput(), currentNode.getServer(), currentNode.getOutput());
            for (int j = 0; j < nodeSections.size(); j++) {
                Element xml_section = xml.simDoc.createElement("section");
                Section currentSection = nodeSections.get(j);
                if (currentSection != null && !currentSection.getClassName().startsWith("Generic ")) {
                    xml_section.setAttribute("className", currentSection.getClassName());
                    DocumentSectionPair simXML = new DocumentSectionPair(xml.simDoc, xml_section);
                    switch (currentSection.getClassName()) {
                        case "Buffer":
                            simXML.section.setAttribute("className", "Queue"); // overwrite with JMT class name
                            simXML = saveBufferCapacity(simXML, i);
                            simXML = saveDropStrategy(simXML, i); // unfinished
                            simXML = saveGetStrategy(simXML);
                            simXML = savePutStrategy(simXML, i);
                            break;
                        case "Server":
                        case "jline.Server":
                        case "PreemptiveServer":
                            simXML = saveNumberOfServers(simXML, i);
                            simXML = saveServerVisits(simXML);
                            simXML = saveServiceStrategy(simXML, i);
                            break;
                        case "SharedServer":
                            simXML.section.setAttribute("className", "PSServer"); // overwrite with JMT class name
                            simXML = saveNumberOfServers(simXML, i);
                            simXML = saveServerVisits(simXML);
                            simXML = saveServiceStrategy(simXML, i);
                            simXML = savePreemptiveStrategy(simXML, i);
                            simXML = savePreemptiveWeights(simXML, i);
                            break;
                        case "InfiniteServer":
                        case "jline.InfiniteServer":
                            simXML.section.setAttribute("className", "Delay"); // overwrite with JMT class name
                            simXML = saveServiceStrategy(simXML, i);
                            break;
                        case "RandomSource":
                        case "jline.RandomSource":
                            simXML = saveArrivalStrategy(simXML, i);
                            break;
                        case "Dispatcher":
                        case "jline.Dispatcher":
                        case "ClassSwitchDispatcher":
                            simXML.section.setAttribute("className", "Router"); // overwrite with JMT class name
                            simXML = saveRoutingStrategy(simXML, i);
                            break;
                        case "StatelessClassSwitcher":
                        case "jline.StatelessClassSwitcher":
                            simXML.section.setAttribute("className", "ClassSwitch"); // overwrite with JMT class name
                            simXML = saveClassSwitchStrategy(simXML, i);
                            break;
                        case "LogTunnel":
                            simXML = saveLogTunnel(simXML, i);
                            break;
                        case "Joiner":
                            simXML.section.setAttribute("className", "Join"); // overwrite with JMT class name
                            simXML = saveJoinStrategy(simXML, i);
                            break;
                        case "Forker":
                            simXML.section.setAttribute("className", "Fork"); // overwrite with JMT class name
                            simXML = saveForkStrategy(simXML, i);
                            break;
                        case "Storage":
                            simXML.section.setAttribute("className", "Storage"); // overwrite with JMT class name
                            simXML = saveTotalCapacity(simXML, i);
                            simXML = savePlaceCapacities(simXML, i);
                            simXML = saveDropRule(simXML, i);
                            simXML = saveGetStrategy(simXML);
                            simXML = savePutStrategies(simXML, i);
                            break;
                        case "Enabling":
                            simXML.section.setAttribute("className", "Enabling"); // overwrite with JMT class name
                            simXML = saveEnablingConditions(simXML, i);
                            simXML = saveInhibitingConditions(simXML, i);
                            break;
                        case "Firing":
                            simXML.section.setAttribute("className", "Firing"); // overwrite with JMT class name
                            simXML = saveFiringOutcomes(simXML, i);
                            break;
                        case "Timing":
                            simXML.section.setAttribute("className", "Timing"); // overwrite with JMT class name
                            simXML = saveModeNames(simXML, i);
                            simXML = saveNumbersOfServers(simXML, i);
                            simXML = saveTimingStrategies(simXML, i);
                            simXML = saveFiringPriorities(simXML, i);
                            simXML = saveFiringWeights(simXML, i);
                            break;
                    }
                    node.appendChild(simXML.section);
                }
            }
            xml.simElem.appendChild(node);
        }
        xml = saveMetrics(xml);
        xml = saveLinks(xml);
        xml = saveRegions(xml);

        boolean hasReferenceNodes = false;
        Element preloadNode = xml.simDoc.createElement("preload");
        Map<StatefulNode, Matrix> s0 = sn.state;
        int numOfStations = sn.nstations;

        for (int i = 0; i < numOfStations; i++) {
            boolean isReferenceNode = false;
            int nodeIndex = (int) sn.stationToNode.get(i);
            int isf = (int) sn.stationToStateful.get(i);
            Element stationPopulationsNode = null;

            if (sn.nodetypes.get(nodeIndex) != NodeType.Source && sn.nodetypes.get(nodeIndex) != NodeType.Join) {
                State.StateMarginalStatistics sms = State.toMarginal(sn, nodeIndex, s0.get(sn.stations.get(isf)), null, null, null, null, null);
                stationPopulationsNode = xml.simDoc.createElement("stationPopulations");
                stationPopulationsNode.setAttribute("stationName", sn.nodenames.get(nodeIndex));

                for (int r = 0; r < numOfClasses; r++) {
                    Element classPopulationNode = xml.simDoc.createElement("classPopulation");

                    if (Double.isInfinite(sn.njobs.get(r))) {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(r))));
                        classPopulationNode.setAttribute("refClass", sn.classnames.get(r));
                        stationPopulationsNode.appendChild(classPopulationNode);
                    } else {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(r))));
                        classPopulationNode.setAttribute("refClass", sn.classnames.get(r));
                        stationPopulationsNode.appendChild(classPopulationNode);
                    }
                }
            }

            if (isReferenceNode) {
                preloadNode.appendChild(stationPopulationsNode);
            }
            hasReferenceNodes = hasReferenceNodes || isReferenceNode;
        }

        if (hasReferenceNodes) {
            xml.simElem.appendChild(preloadNode);
        }
        try {
            SysUtils.writeXML(outputFileName, xml.simDoc);
        } catch (Exception e) {
            e.printStackTrace();
            try {
                SysUtils.writeXML(outputFileName, xml.simDoc);
            } catch (Exception retryException) {
                retryException.printStackTrace();
            }
        }
        return outputFileName;
    }

    public String writeJSIM(NetworkStruct sn) throws ParserConfigurationException {
        String outputFileName = getJSIMTempPath();
        return writeJSIM(sn, outputFileName);
    }

    public SolverResult getResults() {
        SolverOptions options = this.options;
        SolverResult solverResult = new SolverResult();
        switch (options.method) {
            case "jsim":
            case "default":
                this.jmtResult = getResultsJSIM();
                break;
            default:
                this.jmtResult = getResultsJMVA();
                break;
        }

        NetworkStruct sn = getStruct();
        solverResult.QN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.UN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.RN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.TN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.AN = new Matrix(sn.nstations, sn.nclasses);
        solverResult.WN = new Matrix(sn.nstations, sn.nclasses);

        for (int m = 0; m < this.jmtResult.metrics.size(); m++) {
            Metric metric = this.jmtResult.metrics.get(m);
            String stationName = metric.getStation();
            List<Integer> indList = findInd(stationName, sn.nodenames);
            List<Integer> istStations = new ArrayList<>();
            for (int ind : indList) {
                istStations.add((int) sn.nodeToStation.get(ind));
            }
            List<Integer> rList = new ArrayList<>();
            String metricClass = metric.getClassName();
            for (int i = 0; i < sn.classnames.size(); i++) {
                String className = sn.classnames.get(i);
                if (metricClass.equalsIgnoreCase(className)) {
                    rList.add(i);
                }
            }
            boolean open = true;
            for (int r : rList) {
                if (Double.isFinite(sn.njobs.get(r))) {
                    open = false;
                }
            }
            int sumJobs = 0;
            for (int r : rList) {
                sumJobs += sn.njobs.get(r);
            }
            switch (metric.getMetricType()) {
                case QLen:
                    if (open) {
                        solverResult.QN = setValues(solverResult.QN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.QN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.QN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.QN, istStations, rList, 0);
                    }
                    break;
                case Util:
                    if (open) {
                        solverResult.UN = setValues(solverResult.UN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.UN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.UN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.UN, istStations, rList, 0);
                    }
                    break;
                case RespT:
                    if (open) {
                        solverResult.RN = setValues(solverResult.RN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.RN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.RN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.RN, istStations, rList, 0);
                    }
                    break;
                case ResidT:
                    if (open) {
                        solverResult.WN = setValues(solverResult.WN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.WN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.WN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.WN, istStations, rList, 0);
                    }
                    break;
                case ArvR:
                    if (open) {
                        solverResult.AN = setValues(solverResult.AN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.AN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.AN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.AN, istStations, rList, 0);
                    }
                    break;
                case Tput:
                    if (open) {
                        solverResult.TN = setValues(solverResult.TN, istStations, rList, metric.getMeanValue());
                    } else {    // closed
                        solverResult.TN = metric.getAnalyzedSamples() > sumJobs ? setValues(solverResult.TN, istStations, rList, metric.getMeanValue()) : setValues(solverResult.TN, istStations, rList, 0);
                    }
                    break;
            }
        }
        this.result = solverResult;
        return solverResult;
    }

    private List<Integer> findInd(String term, List<String> list) {
        List<Integer> res = new ArrayList<>();
        for (int i = 0; i < list.size(); i++) {
            String nodeName = list.get(i);
            if (nodeName.equalsIgnoreCase(term)) {
                res.add(i);
            }
        }
        return res;
    }

    private Matrix setValues(Matrix matrix, List<Integer> istStations, List<Integer> rList, double value) {
        for (int i : istStations) {
            for (int r : rList) {
                matrix.set(i, r, value);
            }
        }
        return matrix;
    }

    public SolverJMTResult getResultsJSIM() {
        SolverJMTResult result = new SolverJMTResult();
        String fileName = this.getFileName() + ".jsim-result.jsim";
        Path filePath = Paths.get(getFilePath(), fileName);
        File file = new File(filePath.toUri());
        if (file.exists() && !file.isDirectory()) {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            try {
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(file);
                NodeList measure = doc.getElementsByTagName("measure");
                for (int i = 0; i < measure.getLength(); i++) {
                    NamedNodeMap attributes = measure.item(i).getAttributes();
                    Metric metric = new Metric(attributes);
                    result.metrics.add(metric);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            java.lang.System.err.println("JMT did not output a result file, the simulation has likely failed.");
        }
        return result;
    }

    public SolverJMTResult getResultsJMVA() {
        SolverJMTResult result = new SolverJMTResult();
        String fileName = this.getFileName() + ".jmva-result.jmva";
        Path filePath = Paths.get(getFilePath(), fileName);
        File file = new File(filePath.toUri());
        if (file.exists() && !file.isDirectory()) {
            DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
            try {
                DocumentBuilder builder = factory.newDocumentBuilder();
                Document doc = builder.parse(file);
                NodeList measure = doc.getElementsByTagName("measure");
                for (int i = 0; i < measure.getLength(); i++) {
                    NamedNodeMap attributes = measure.item(i).getAttributes();
                    Metric metric = new Metric(attributes);
                    result.metrics.add(metric);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else {
            java.lang.System.err.println("JMT did not output a result file, the analysis has likely failed.");
        }
        return result;
    }

    public double getProbNormConstAggr() throws ParserConfigurationException {
        Double lNormConst;
        if (GlobalConstants.DummyMode) {
            lNormConst = NaN;
            return lNormConst;
        }
        switch (options.method) {
            case "jmva":
            case "jmva.recal":
            case "jmva.comom":
            case "jmva.ls":
                this.runAnalyzer();
                lNormConst = this.jmtResult.logNormConstAggr;
                break;
            default:
                lNormConst = NaN;
                java.lang.System.err.println("Selected solver method does not compute normalizing constants. Choose either jmva.recal, jmva.comom, or jmva.ls.");
        }
        return lNormConst;
    }

    public double getProbAggr(Node node, Matrix state_a) {
        double Pr = NaN;
        if (GlobalConstants.DummyMode) {
            return Pr;
        }

        NetworkStruct sn = getStruct();
        double stationStateAggr = this.sampleAggr(node);
        java.lang.System.err.println("getProbAggr() has not yet been implemented in JLINE.");
        return Pr;
    }

    public double getProbAggr(Node node) {
        double Pr = NaN;
        if (GlobalConstants.DummyMode) {
            return Pr;
        }
        Matrix state_a = sn.state.get(sn.stations.get((int) sn.stationToStateful.get((int) sn.nodeToStation.get(node.getNodeIdx()))));
        return getProbAggr(node, state_a);
    }

    // TODO: add sampling method
    public double sampleAggr(Node node, int numEvents, boolean markActivePassive) throws IOException {
        if (GlobalConstants.DummyMode) {
            return NaN;
        }
        // TODO: not implemented
        throw new RuntimeException("sampleAggr() has not yet been implemented in JLINE.");
//        NetworkStruct sn = getStruct();
//        Map<Station, Map<JobClass, SolverHandles.Metric>> Q = getAvgQLenHandles();
//        Network modelCopy = this.model;
//        modelCopy.resetNetwork();
//
//        int numberOfNodes = modelCopy.getNumberOfNodes();
//        int numberOfClasses = modelCopy.getNumberOfClasses();
//        boolean[][] isNodeClassLogged = new boolean[numberOfNodes][numberOfClasses];
//        int ind = this.model.getNodeIndex(node);
//
//        for (int i = 0; i < numberOfNodes; i++) {
//            for (int j = 0; j < numberOfClasses; j++) {
//                isNodeClassLogged[i][j] = i == ind;
//            }
//        }
//        Map<JobClass, Map<JobClass, Matrix>> Plinked = sn.rtorig;
//        String logpath = JMT.lineTempName();
//        modelCopy.linkAndLog();
    }

    // TODO: add sampling method
    public double sampleAggr(Node node, int numEvents) throws IOException {
        return sampleAggr(node, numEvents, false);
    }

    // TODO: add sampling method
    public double sampleAggr(Node node) {
        // TODO: not implemented
        throw new RuntimeException("JMT does not allow to fix the number of events for individual nodes. The number of returned events may be inaccurate.");
    }

    public double[] getTranProbAggr(Node node) {
//        double Pi_t;
//        double SSnode_a;
//        if (GlobalConstants.DummyMode){
//            Pi_t = NaN;
//            SSnode_a = NaN;
//            return new double[] {Pi_t, SSnode_a};
//        }
//        System.err.println("Method not implemented yet.");
//        SolverOptions options = this.options;
//        int initSeed = this.options.seed;
//        if (options.timespan != null && Double.isFinite(options.timespan[1])){
//            List<Double> tu = new ArrayList<>();
//             options.iter_max
//            double [] TranSysStateAggr = new SampleSysAggr[options.get("iter_max")];
//        }
        // TODO: implementation
        throw new RuntimeException("getTranProbAggr() has not yet been implemented in JLINE.");
    }

    public NetworkStruct sampleSysAggr(int numEvents, boolean markActivePassive) {
        if (GlobalConstants.DummyMode) {
            return null;
        }
        // TODO: not implemented
        throw new RuntimeException("sampleSysAggr() has not yet been implemented in JLINE.");
//        NetworkStruct sn = getStruct();
//        numEvents -= 1;
//        Map<Station, Map<JobClass, SolverHandles.Metric>> Q = getAvgQLenHandles();
//        Matrix statStateAggr = new Matrix(sn.nstations, 1);
//        Network modelCopy = this.model;
//        modelCopy.resetNetwork();
    }

    public NetworkStruct sampleSysAggr(int numEvents) {
        boolean markActivePassive = false;
        return sampleSysAggr(numEvents, markActivePassive);
    }

    public NetworkStruct sampleSysAggr() {
        int numEvents = this.options.samples;
        return sampleSysAggr(numEvents);
    }

    public double probSysStateAggr() {
        if (GlobalConstants.DummyMode) {
            return NaN;
        }
        // TODO: not implemented
        throw new RuntimeException("probSysStateAggr() has not yet been implemented in JLINE.");
    }

    public Matrix getCdfRespT() {
        Map<Station, Map<JobClass, SolverHandles.Metric>> R = getAvgRespTHandles();
        return getCdfRespT(R);
    }

    public Matrix getCdfRespT(Map<Station, Map<JobClass, SolverHandles.Metric>> R) {
        // TODO: not implemented
        throw new RuntimeException("getCdfRespT() has not yet been implemented in JLINE.");
    }

    public Matrix getTranCdfRespT() {
        Map<Station, Map<JobClass, SolverHandles.Metric>> R = getAvgRespTHandles();
        return getCdfRespT(R);
    }

    public Matrix getTranCdfRespT(Map<Station, Map<JobClass, SolverHandles.Metric>> R) {
        // TODO: not implemented
        throw new RuntimeException("getTranCdfRespT() has not yet been implemented in JLINE.");
    }

    public Matrix getTranCdfPassT() {
        NetworkStruct sn = getStruct();
        Matrix RD = new Matrix(sn.nstations, sn.nclasses);
        if (GlobalConstants.DummyMode) {
            return RD;
        }
        Map<Station, Map<JobClass, SolverHandles.Metric>> R = getAvgRespTHandles();
        return getTranCdfPassT(R);
    }

    public Matrix getTranCdfPassT(Map<Station, Map<JobClass, SolverHandles.Metric>> R) {
        // TODO: not implemented
        throw new RuntimeException("getTranCdfPassT() has not yet been implemented in JLINE.");
//        NetworkStruct sn = getStruct();
//        Matrix RD = new Matrix(sn.nstations,sn.nclasses);
//        if (GlobalConstants.DummyMode) {
//            return RD;
//        }
//        Network cdfmodel = this.model;
//        cdfmodel.resetNetwork();
//        int numberOfNodes = cdfmodel.getNumberOfNodes();
//        int numberOfClasses = cdfmodel.getNumberOfClasses();
//        boolean[][] isNodeClassLogged = new boolean[numberOfNodes][numberOfClasses];
//        for (int i = 0; i < numberOfNodes; i++) {
//            for (int r = 0; r < numberOfClasses; r++) {
//                isNodeClassLogged[i][r] = false;
//            }
//        }
//        for (int i = 0; i < numberOfNodes; i++) {
//            for (int r = 0; r < numberOfClasses; r++) {
//                if (!R.get(sn.stations.get((int) sn.nodeToStation.get(i))).get(sn.jobclasses.get(r)).isDisabled){
//                    int ni = this.model.getNodeIndex(cdfmodel.getStationFromIndex(i));
//                    isNodeClassLogged[ni][r] = true;
//                }
//            }
//        }
    }

    public static String listValidMethods() {
        return "default,jsim,jmva,jmva.amva,jmva.mva,jmva.recal,jmva.comom,jmva.chow,jmva.bs,jmva.aql,jmva.lin,jmva.dmlin,jmva.ls";
    }

    public static FeatureSet getFeatureSet() {
        FeatureSet s = new FeatureSet();
        // TODO: update with the features supported.
        String[] features = {"Sink", "Source", "Router",
                "ClassSwitch", "Delay", "Queue", "Fork", "Join", "Forker", "Joiner", "Logger",
                "Coxian", "APH", "Erlang", "Exp", "HyperExp", "Det", "Gamma", "Lognormal", "MAP", "MMPP2", "Normal", "PH",
                "Pareto", "Weibull", "Replayer", "Uniform",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
                "Server", "Sink", "RandomSource", "ServiceTunnel",
                "CacheClassSwitcher", "Cache", "LogTunnel", "Linkage", "Enabling", "Timing", "Firing", "Storage", "Place",
                "Transition",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS", "SchedStrategy_GPS", "SchedStrategy_SIRO", "SchedStrategy_HOL",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_SEPT", "SchedStrategy_LEPT", "SchedStrategy_SJF", "SchedStrategy_LJF",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND", "RoutingStrategy_RROBIN", "RoutingStrategy_WRROBIN",
                "RoutingStrategy_KCHOICES", "SchedStrategy_EXT",
                "ClosedClass", "OpenClass",};
        s.setTrue(features);
        return s;
    }

    public static boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverJMT.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    public static String writeJMVA(NetworkStruct sn, String outputFileName, SolverOptions options) {
        try {
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document mvaDoc = dBuilder.newDocument();

            Element mvaElem = mvaDoc.createElement("model");
            mvaDoc.appendChild(mvaElem);
            mvaElem.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
            mvaElem.setAttribute("xsi:noNamespaceSchemaLocation", "JMTmodel.xsd");

            Element algTypeElement = mvaDoc.createElement("algType");

            switch (options.method) {
                case "jmva.recal":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "RECAL");
                    break;
                case "jmva.comom":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "CoMoM");
                    break;
                case "jmva.chow":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Chow");
                    break;
                case "jmva.bs":
                case "jmva.amva":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Bard-Schweitzer");
                    break;
                case "jmva.aql":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "AQL");
                    break;
                case "jmva.lin":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "Linearizer");
                    break;
                case "jmva.dmlin":
                    setAlgTypeName(algTypeElement, sn.nservers, options.method, "De Souza-Muntz Linearizer");
                    break;
                case "jmva.ls":
                    algTypeElement.setAttribute("name", "Logistic Sampling");
                    break;
                default:
                    algTypeElement.setAttribute("name", "MVA");
                    break;
            }

            algTypeElement.setAttribute("tolerance", "1.0E-7");
            algTypeElement.setAttribute("maxSamples", String.valueOf(options.samples));

            int M = sn.nstations;    //number of stations
            Matrix NK = sn.njobs.transpose();  // initial population per class
            int C = sn.nchains;
            Matrix SCV = sn.scv;
            Matrix ST = new Matrix(sn.rates.numRows, sn.rates.numCols);

            for (int i = 0; i < sn.rates.numRows; i++) {
                for (int j = 0; j < sn.rates.numCols; j++) {

                    double currentRate = sn.rates.get(i, j);
                    if (Double.isNaN(currentRate)) {
                        ST.set(i, j, 0);
                    } else {
                        ST.set(i, j, 1.0 / currentRate);
                    }

                    if (Double.isNaN(SCV.get(i, j))) {
                        SCV.set(i, j, 1);
                    }
                }
            }

            SN.snGetDemandsChainReturn snGetDemandsChainReturn = snGetDemandsChain(sn);

            Element parametersElem = mvaDoc.createElement("parameters");
            Element classesElem = mvaDoc.createElement("classes");
            classesElem.setAttribute("number", String.valueOf(sn.nchains));
            Element stationsElem = mvaDoc.createElement("stations");
            int numberOfStations = sn.nstations - countNodesWithType(sn.nodetypes, NodeType.Source);
            stationsElem.setAttribute("number", String.valueOf(numberOfStations));
            Element refStationsElem = mvaDoc.createElement("ReferenceStation");
            refStationsElem.setAttribute("number", String.valueOf(sn.nchains));
            Element algParamsElem = mvaDoc.createElement("algParams");

            boolean[] sourceid = new boolean[sn.nodetypes.size()];
            for (int i = 0; i < sn.nodetypes.size(); i++) {
                sourceid[i] = sn.nodetypes.get(i) == NodeType.Source;
            }

            for (int c = 0; c < sn.nchains; c++) {
                Element classElem;
                double sumOfNJobs = 0.0;
                for (int i = 0; i < sn.chains.numCols; i++) {
                    sumOfNJobs += sn.njobs.get((int) sn.chains.get(c, i));
                }
                if (Double.isFinite(sumOfNJobs)) {
                    classElem = mvaDoc.createElement("closedclass");
                    classElem.setAttribute("population", String.valueOf(snGetDemandsChainReturn.Nchain.get(c)));
                    classElem.setAttribute("name", String.format("Chain%02d", c + 1));
                } else {
                    double rateSum = 0.0;
                    for (int i = 0; i < sourceid.length; i++) {
                        if (sourceid[i]) {
                            for (int j = 0; j < sn.chains.numCols; j++) {
                                rateSum += sn.rates.get(i, (int) sn.chains.get(c, j));
                            }
                        }
                    }
                    classElem = mvaDoc.createElement("openclass");
                    classElem.setAttribute("rate", String.valueOf(rateSum));
                    classElem.setAttribute("name", String.format("Chain%02d", c + 1));
                }
                classesElem.appendChild(classElem);
            }

            boolean[] isLoadDep = new boolean[sn.nstations];
            for (int i = 0; i < sn.nstations; i++) {
                Element statElem = null;
                NodeType currentNodeType = sn.nodetypes.get((int) sn.stationToNode.get(i));
                switch (currentNodeType) {
                    case Delay:
                        statElem = mvaDoc.createElement("delaystation");
                        statElem.setAttribute("name", sn.nodenames.get((int) sn.stationToNode.get(i)));
                        break;
                    case Queue:
                        if (sn.nservers.get(i) == 1) {
                            isLoadDep[i] = false;
                            statElem = mvaDoc.createElement("listation");
                        } else {
                            isLoadDep[i] = true;
                            statElem = mvaDoc.createElement("ldstation");
                        }
                        statElem.setAttribute("name", sn.nodenames.get((int) sn.stationToNode.get(i)));
                        statElem.setAttribute("servers", "1");
                        break;
                    default:
                        continue;
                }

                Element srvTimesElem = mvaDoc.createElement("servicetimes");
                for (int c = 0; c < sn.nchains; c++) {
                    if (isLoadDep[i]) {
                        Element statSrvTimeElem = mvaDoc.createElement("servicetimes");
                        statSrvTimeElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));
                        String ldSrvString = String.valueOf(snGetDemandsChainReturn.STchain.get(i, c));
                        if (anyElementIsInfinity(NK.toArray1D())) {
                            throw new RuntimeException("JMVA does not support open classes in load-dependent models;");
                        }

                        for (int n = 1; n <= Arrays.stream(NK.toArray1D()).sum(); n++) {
                            ldSrvString = String.format("%s;%s", ldSrvString, snGetDemandsChainReturn.STchain.get(i, c) / Math.min(n, sn.nservers.get(i)));
                        }

                        statSrvTimeElem.appendChild(mvaDoc.createTextNode(ldSrvString));
                        srvTimesElem.appendChild(statSrvTimeElem);
                    } else {
                        Element statSrvTimeElem = mvaDoc.createElement("servicetime");
                        statSrvTimeElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));
                        statSrvTimeElem.appendChild(mvaDoc.createTextNode(String.valueOf(snGetDemandsChainReturn.STchain.get(i, c))));
                        srvTimesElem.appendChild(statSrvTimeElem);
                    }
                }
                statElem.appendChild(srvTimesElem);
                Element visitsElem = mvaDoc.createElement("visits");
                for (int c = 0; c < sn.nchains; c++) {
                    Element statVisitElem = mvaDoc.createElement("visit");
                    statVisitElem.setAttribute("customerclass", String.format("Chain%02d", c + 1));

                    double val;
                    if (snGetDemandsChainReturn.STchain.get(i, c) > 0) {
                        val = snGetDemandsChainReturn.Lchain.get(i, c) / snGetDemandsChainReturn.STchain.get(i, c);
                    } else {
                        val = 0;
                    }

                    statVisitElem.appendChild(mvaDoc.createTextNode(String.valueOf(val)));
                    visitsElem.appendChild(statVisitElem);
                }
                statElem.appendChild(visitsElem);
                stationsElem.appendChild(statElem);
            }

            int[] refstatchain = new int[C];
            for (int c = 0; c < sn.nchains; c++) {
                Matrix inchain = sn.inchain.get(c);
                refstatchain[c] = (int) sn.refstat.get((int) inchain.get(0));
            }
            for (int c = 0; c < sn.nchains; c++) {
                Element classRefElem = mvaDoc.createElement("Class");
                classRefElem.setAttribute("name", String.format("Chain%d", c + 1));
                classRefElem.setAttribute("refStation", sn.nodenames.get((int) sn.stationToNode.get(refstatchain[c])));
                refStationsElem.appendChild(classRefElem);
            }

            Element compareAlgsElem = mvaDoc.createElement("compareAlgs");
            compareAlgsElem.setAttribute("value", "false");
            algParamsElem.appendChild(algTypeElement);
            algParamsElem.appendChild(compareAlgsElem);

            parametersElem.appendChild(classesElem);
            parametersElem.appendChild(stationsElem);
            parametersElem.appendChild(refStationsElem);
            mvaElem.appendChild(parametersElem);
            mvaElem.appendChild(algParamsElem);
            writeXML(outputFileName, mvaDoc);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return outputFileName;
    }

    private static void setAlgTypeName(Element algTypeElement, Matrix nservers, String method, String name) {
        double maxFiniteValue = Double.NEGATIVE_INFINITY;  // initial value set to negative infinity

        for (int i = 0; i < nservers.numRows; i++) {
            if (Double.isFinite(nservers.get(i, 0))) {
                double currentValue = nservers.get(i, 0);
                if (currentValue > maxFiniteValue) {
                    maxFiniteValue = currentValue;
                }
            }
        }

        if (maxFiniteValue > 1) {
            throw new RuntimeException(method + " does not support multi-server stations.");
        } else {
            algTypeElement.setAttribute("name", name);
        }
    }

    private static int countNodesWithType(List<NodeType> nodetypes, NodeType type) {
        int count = 0;
        for (NodeType nodeType : nodetypes) {
            if (nodeType == type) count++;
        }
        return count;
    }

    private static boolean anyElementIsInfinity(double[] array) {
        for (double val : array) {
            if (Double.isInfinite(val)) {
                return true;
            }
        }
        return false;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }


    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public double getMaxSimulatedTime() {
        return maxSimulatedTime;
    }

    public void setMaxSimulatedTime(double maxSimulatedTime) {
        this.maxSimulatedTime = maxSimulatedTime;
    }

    public int getMaxSamples() {
        return maxSamples;
    }

    public void setMaxSamples(int maxSamples) {
        this.maxSamples = maxSamples;
    }

    public int getMaxEvents() {
        return maxEvents;
    }

    public void setMaxEvents(int maxEvents) {
        this.maxEvents = maxEvents;
    }

    public int getSeed() {
        return seed;
    }

    public void setSeed(int seed) {
        this.seed = seed;
    }

    public double getSimConfInt() {
        return simConfInt;
    }

    public void setSimConfInt(double simConfInt) {
        this.simConfInt = simConfInt;
    }

    public double getSimMaxRelErr() {
        return simMaxRelErr;
    }

    public void setSimMaxRelErr(double simMaxRelErr) {
        this.simMaxRelErr = simMaxRelErr;
    }

    private String getJSIMTempPath() {
        if (this.filePath == null || this.fileName == null) {
            try {
                this.filePath = SysUtils.lineTempName("jsim");
            } catch (IOException ioe) {
                ioe.printStackTrace();
                throw new RuntimeException("Unable to get JSIM temp path");
            }
            this.fileName = "model";
        }
        String fname = this.fileName + ".jsim";
        return filePath + File.separator + fname;
    }

    private String getJMVATempPath() {
        if (this.filePath == null || this.fileName == null) {
            try {
                this.filePath = SysUtils.lineTempName("jmva");
            } catch (IOException ioe) {
                ioe.printStackTrace();
                throw new RuntimeException("Unable to get JMVA temp path");
            }
            this.fileName = "model";
        }
        String fname = this.fileName + ".jmva";
        return filePath + File.separator + fname;
    }

    protected boolean hasAvgResults() {
        return hasResults();
    }

    @Override
    protected void runAnalyzer() throws ParserConfigurationException {
        long startTime = java.lang.System.currentTimeMillis();
//        if (this.enableChecks && !SolverJMT.supports(this.model)) {
//            throw new RuntimeException("This model contains features not supported by the solver.");
//        }

        if (this.model == null)
            throw new RuntimeException("Model is not provided");
        if (this.options == null)
            this.options = new SolverOptions(SolverType.JMT);
        if (options.verbose == null) {
            options.verbose = VerboseLevel.values()[0];
        }
        if (options.samples == 0) {
            options.samples = 10000;
        } else if (options.samples < 5000) {
            if (!options.method.equalsIgnoreCase("jmva.ls")) {
                java.lang.System.err.printf(String.format("JMT requires at least 5000 samples for each metric, the current value is %d. Starting the simulation with 5000 samples.%n", options.samples));
            }
            options.samples = 5000;
        }
        if (options.seed == 0) {
            options.seed = Math.toIntExact(Math.round((Math.random() * (1e6 - 1)) + 1));
        }
        this.seed = options.seed;
        if (options.timespan == null) {
            options.timespan = new double[]{0.0, Double.POSITIVE_INFINITY};
        } else {
            this.maxSimulatedTime = options.timespan[1];
        }

        if (!this.model.hasInitState()) {
            this.model.initDefault();
        }
        this.maxSamples = options.samples;
        NetworkStruct sn = this.getStruct();
        String fname = "";
        String cmd = "";
        String cmdOutput = "";
        long runTime = 0;
        SolverResult solverResult = null;


        switch (options.method) {
            case "jsim":
            case "default":
                fname = this.writeJSIM(sn);
                cmd = String.format(
                        "java -cp %s jmt.commandline.Jmt sim %s -seed %s --illegal-access=permit",
                        this.jmtPath, fname, options.seed
                );

                if (options.verbose != VerboseLevel.SILENT) {
                    java.lang.System.out.println("JMT Model: " + fname);
                }
                if (options.verbose == VerboseLevel.DEBUG) {
                    java.lang.System.out.println("JMT Command: " + cmd);
                }
                cmdOutput = SysUtils.system(cmd);
                if (!cmdOutput.isEmpty()) {
                    java.lang.System.out.println("JMT Command output: " + cmdOutput);
                }
                runTime = java.lang.System.currentTimeMillis() - startTime;

                solverResult = getResults();
                if (!options.keep) {
//                    try {
//                        //JMT.removeDirectory(Paths.get(this.getFilePath()));
//                    } catch (IOException ioe) {
//                        ioe.printStackTrace();
//                    }
                }
                solverResult.runtime = runTime / 1000.0;
                this.result = solverResult;
                if (options.verbose != VerboseLevel.values()[0]) {
                    java.lang.System.out.println("JMT Analysis completed. Runtime: " + result.runtime + " seconds");
                }
                break;
            case "closing":
                // TODO: implementation
                throw new RuntimeException("method: closing has not yet been implemented in JLINE.");

//                SolverOptions options = this.options;
//                int initSeed = this.options.seed;
//                if (options.timespan.length > 1){
//                    double[] initTimeSpan = this.options.timespan;
//                    this.options.timespan[0] = this.options.timespan[1];
//                    if (Double.isFinite(this.options.timespan[1])){
//                        double [] tu;
//                        for (int it = 0; it < options.iter_max; it++) {
//                            this.options.seed = initSeed + it - 1;
//                            TranSysStateAggr{it} = sampleSysAggr(self);
//                            if isempty(tu)
//                            tu = TranSysStateAggr{it}.t;
//                        else
//                            % we need to limit the time series at the minimum
//                            % as otherwise the predictor of the state cannot
//                            % take into account constraints that exist on the
//                            % state space
//                            tumax = min(max(tu),max(TranSysStateAggr{it}.t));
//                            tu = union(tu, TranSysStateAggr{it}.t);
//                            tu = tu(tu<=tumax);
//                            end
//                        }
//                        Matrix QNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//                        Matrix UNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//                        Matrix TNt = new Matrix(sn.nstations, sn.nclasses, tu.length, 2);
//
//                    }
//                }

            case "jmva":
            case "jmva.amva":
            case "jmva.mva":
            case "jmva.recal":
            case "jmva.comom":
            case "jmva.chow":
            case "jmva.bs":
            case "jmva.aql":
            case "jmva.lin":
            case "jmva.dmlin":
            case "jmva.ls":
            case "jmt.jmva":
            case "jmt.jmva.mva":
            case "jmt.jmva.amva":
            case "jmt.jmva.recal":
            case "jmt.jmva.comom":
            case "jmt.jmva.chow":
            case "jmt.jmva.bs":
            case "jmt.jmva.aql":
            case "jmt.jmva.lin":
            case "jmt.jmva.dmlin":
            case "jmt.jmva.ls":
                fname = writeJMVA(sn, getJMVATempPath(), this.options);
                cmd = String.format(
                        "java -cp %s jmt.commandline.Jmt mva %s -seed %s --illegal-access=permit",
                        this.jmtPath, fname, this.options.seed
                );

                if (this.options.verbose != VerboseLevel.SILENT) {
                    java.lang.System.out.println("JMT Model: " + fname);
                    java.lang.System.out.println("JMT Command: " + cmd);
                }
                cmdOutput = SysUtils.system(cmd);
                java.lang.System.out.println("JMT Command output: " + cmdOutput);
                runTime = java.lang.System.currentTimeMillis() - startTime;
                solverResult.runtime = runTime / 1000.0;
                this.result = solverResult;
                java.lang.System.out.println("JMT Analysis completed. Runtime: " + runTime / 1000.0 + " seconds");

                getResults();
                if (!this.options.keep) {
//                    try {
//                        //JMT.removeDirectory(Paths.get(this.getFilePath()));
//                    } catch (IOException ioe) {
//                        ioe.printStackTrace();
//                    }
                }
                break;
            default:
                java.lang.System.err.println("Warning: This solver does not support the specified method. Setting to default.");
                this.options.method = "default";
                runAnalyzer();
        }
    }


    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.JMT);
    }
}
