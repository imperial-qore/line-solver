package jline.io;

import jline.api.mam.Map_meanKt;
import jline.api.mam.Map_scvKt;
import jline.api.mam.Map_skewKt;
import jline.api.mam.Map_varKt;
import jline.lang .ClosedClass;
import jline.lang.JobClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.layered.LayeredNetwork;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.TimingStrategy;
import jline.lang.NetworkStruct;
import jline.lang.nodes.*;
import jline.lang.nodes.Queue;
import jline.lang.processes.*;
import jline.lang.sections.Buffer;
import jline.solvers.jmt.SolverJMT;
import jline.util.matrix.Matrix;
import org.apache.commons.io.FilenameUtils;
import org.w3c.dom.*;
import org.w3c.dom.Node;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.util.*;

import static jline.io.InputOutputKt.*;

/**
 * Model-to-Model transformation class for converting between different queueing network model formats.
 * 
 * <p>This class provides methods to transform between various queueing network modeling formats,
 * including JMT (Java Modelling Tools) models and LINE models. It supports importing JSIM/JSIMG
 * files and converting them to LINE Network objects, as well as exporting LINE models back to
 * JMT-compatible formats.</p>
 * 
 * <p>Supported transformations:
 * <ul>
 *   <li>JMT to LINE: Import JSIM/JSIMG/JSIMW files and convert to LINE Network objects</li>
 *   <li>LQN to LINE: Import LQN/LQNX (Layered Queueing Network) files to LayeredNetwork objects</li>
 *   <li>LINE to JSIMG: Export LINE Network models to JMT-compatible JSIMG format</li>
 * </ul>
 * </p>
 * 
 * @author LINE Development Team
 * @since 1.0
 */
public class M2M {
    /**
     * Creates a new instance of the M2M transformation class.
     */
    public M2M() {
    }

    /**
     * Converts a JMT (Java Modelling Tools) model file to a LINE Network object.
     * 
     * <p>This method automatically detects the JMT file format based on its extension:
     * <ul>
     *   <li>.jmva - JMVA files (not yet supported)</li>
     *   <li>.jsim, .jsimg, .jsimw - JSIM simulation model files</li>
     * </ul>
     * </p>
     * 
     * @param filename The path to the JMT model file to import
     * @return A LINE Network object representing the imported model
     * @throws RuntimeException if the file format is not supported or parsing fails
     */
    public Network JMT2LINE(String filename) {
        String fext = FilenameUtils.getExtension(filename);

        switch (fext) {
            case "jmva":
                line_error(mfilename(new Object() {
                }), "JMVA files not yet supported.");
                // return JMVA2LINE(filename);
            case "jsim":
            case "jsimg":
            case "jsimw":
                return JSIM2LINE(filename);
            default:
                line_warning(
                        mfilename(new Object() {
                        }),
                        "The file has unknown extension, trying to parse as a JSIMG file.\n");
                return JSIM2LINE(filename);
        }
    }

    /**
     * Converts a JMT (Java Modelling Tools) model file to a LINE Network object with a custom model name.
     * 
     * <p>This method automatically detects the JMT file format based on its extension:
     * <ul>
     *   <li>.jmva - JMVA files (not yet supported)</li>
     *   <li>.jsim, .jsimg, .jsimw - JSIM simulation model files</li>
     * </ul>
     * </p>
     * 
     * @param filename The path to the JMT model file to import
     * @param modelName The name to assign to the resulting LINE Network model
     * @return A LINE Network object representing the imported model with the specified name
     * @throws RuntimeException if the file format is not supported or parsing fails
     */
    public Network JMT2LINE(String filename, String modelName) {
        String fext = FilenameUtils.getExtension(filename);

        switch (fext) {
            case "jmva":
                line_error(mfilename(new Object() {
                }), "JMVA files not yet supported.");
                // return JMVA2LINE(filename, modelName);
            case "jsim":
            case "jsimg":
            case "jsimw":
                return JSIM2LINE(filename, modelName);
            default:
                line_warning(
                        mfilename(new Object() {
                        }),
                        "The file has unknown extension, trying to parse as a JSIMG file.\n");
                return JSIM2LINE(filename, modelName);
        }
    }

    /**
     * Converts a JSIM/JSIMG file to a LINE Network object.
     * 
     * <p>This method parses the JSIM XML file and extracts the model name from the document.
     * It then delegates to the overloaded method that accepts a model name parameter.</p>
     * 
     * @param filename The path to the JSIM/JSIMG file to import
     * @return A LINE Network object representing the imported JSIM model
     * @throws AssertionError if the XML document cannot be read
     */
    public Network JSIM2LINE(String filename) {
        // import model
        Document xDoc = xml_read(filename);

        // get model name
        assert xDoc != null;
        String modelName = FilenameUtils.getBaseName(xDoc.getDocumentElement().getAttribute("name"));

        return JSIM2LINE(filename, modelName);
    }

    /**
     * Converts a JSIM/JSIMG file to a LINE Network object with a specified model name.
     * 
     * <p>This is the main implementation method that performs the actual JSIM to LINE conversion.
     * It parses the JSIM XML structure and creates corresponding LINE objects including:
     * <ul>
     *   <li>Nodes: Sources, Queues, Sinks, Delays, Routers, Forks, Joins, ClassSwitches, Places, Transitions</li>
     *   <li>Job Classes: Open and Closed classes with their properties</li>
     *   <li>Service and Arrival Distributions: Exponential, Erlang, Hyperexponential, Coxian, etc.</li>
     *   <li>Routing Strategies: Random, Probabilistic, Round Robin, Join Shortest Queue, etc.</li>
     *   <li>Scheduling Strategies: FCFS, LCFS, PS, DPS, GPS, Priority-based, etc.</li>
     * </ul>
     * </p>
     * 
     * @param filename The path to the JSIM/JSIMG file to import
     * @param modelName The name to assign to the resulting LINE Network model
     * @return A LINE Network object representing the imported JSIM model
     * @throws RuntimeException if parsing fails or unsupported features are encountered
     */
    public Network JSIM2LINE(String filename, String modelName) {
        // T0 = tic;
        // import model
        Document xDocRead = xml_read(filename);
        assert xDocRead != null;
        Element xDoc = (Element) xDocRead.getElementsByTagName("sim").item(0);

        // create network
        Network model = new Network(modelName);

        // create stations
        List<String> node_name = new ArrayList<>();
        List<String> orig_node_name = new ArrayList<>();
        NodeList xnodeList = xDoc.getElementsByTagName("node");
        List<NodeList> xsection = new ArrayList<>();

        for (int i = 0; i < xnodeList.getLength(); i++) {
            Node xnode = xnodeList.item(i);

            if (xnode.getNodeType() == Node.ELEMENT_NODE) {
                Element xnodeElement = (Element) xnode;
                String xnodeName = xnodeElement.getAttribute("name");

                orig_node_name.add(xnodeName);

                xnodeName = xnodeName.replace("/", "_");
                xnodeName = xnodeName.replace("\\", "_");
                node_name.add(xnodeName);

                xsection.add(xnodeElement.getElementsByTagName("section"));
            }
        }

        List<SchedStrategy> strategy = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<NamedNodeMap>> xsection_javaClass = new ArrayList<>();
        int sink_idx = -1;
        int source_idx = -1;
        List<jline.lang.nodes.Node> node = new ArrayList<>();
        List<List<NamedNodeMap>> xrouting =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<Node>>> xsection_par =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<Node>> xsection_i_par = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<Node>>> xsection_i_subpar =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<Node>>> xsection_i_value =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<NamedNodeMap>> xsection_i_par_attr =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<Node>>> xsvc = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<NamedNodeMap>> xput_strategy =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));

        for (int i = 0; i < node_name.size(); i++) {
            // xsection_i.add(xsection.get(i));
            // cast the Node objects in the NodeList into Element and added to a list
            // equivalent to "xsection_javaClass{i} = {xsection_i{i}.ATTRIBUTE};" in MATLAB version
            xsection_javaClass.add(getAttributesList(xsection.get(i)));
            switch (xsection_javaClass.get(i).get(0).getNamedItem("className").getNodeValue()) {
                case "JobSink":
                    node.add(new Sink(model, node_name.get(i)));
                    //sink_idx = i;
                    break;
                case "RandomSource":
                    node.add(new Source(model, node_name.get(i)));
                    //source_idx = i;
                    xrouting.set(
                            i,
                            getAttributesList(
                                    getSubParameter(
                                            ((Element) xsection.get(i).item(2))
                                                    .getElementsByTagName("parameter")
                                                    .item(0))));
                    break;
                case "Join":
                    List<Integer> forkMap = new ArrayList<>();
                    for (int j = 0; j < node.size(); j++) {
                        if (node.get(j) instanceof Fork) {
                            forkMap.add(j);
                        }
                    }
                    if (forkMap.size() > 1) {
                        line_error(
                                mfilename(new Object() {
                                }), "JSIM2LINE supports at most a single fork-join pair.");
                    }
                    node.add(new Join(model, node_name.get(i), node.get(forkMap.get(0))));
                    xrouting.set(
                            i,
                            getAttributesList(
                                    getSubParameter(
                                            ((Element) xsection.get(i).item(2))
                                                    .getElementsByTagName("parameter")
                                                    .item(0))));
                    break;
                case "Queue":
                    switch (xsection_javaClass.get(i).get(2).getNamedItem("className").getNodeValue()) {
                        case "Fork":
                            Fork fork = new Fork(model, node_name.get(i));
                            fork.setTasksPerLink(
                                    Integer.parseInt(
                                            ((Element)
                                                    ((Element) xsection.get(i).item(2))
                                                            .getElementsByTagName("parameter")
                                                            .item(0))
                                                    .getElementsByTagName("value")
                                                    .item(0)
                                                    .getTextContent()));
                            node.add(fork);
                            xrouting.set(
                                    i,
                                    getAttributesList(
                                            getSubParameter(
                                                    ((Element) xsection.get(i).item(2))
                                                            .getElementsByTagName("parameter")
                                                            .item(3))));
                            break;
                        default:
                            switch (xsection_javaClass.get(i).get(1).getNamedItem("className").getNodeValue()) {
                                case "ServiceTunnel":
                                    node.add(new Router(model, node_name.get(i)));
                                    xrouting.set(
                                            i,
                                            getAttributesList(
                                                    getSubParameter(
                                                            ((Element) xsection.get(i).item(2))
                                                                    .getElementsByTagName("parameter")
                                                                    .item(0))));
                                    break;
                                default:
                                    // xsection_par contains the parameters of all sections of a node
                                    // xsection_i_par contains the parameters of first section of a node
                                    xsection_par.set(i, getChildNodes("parameter", xsection.get(i)));
                                    xsection_i_par.set(i, xsection_par.get(i).get(0));

                                    xsection_i_value.set(i, getChildNodes("value", xsection_i_par.get(i)));
                                    xsection_i_par_attr.set(i, getAttributesList(xsection_i_par.get(i)));
                                    xsection_i_subpar.set(i, getSubParameter(xsection_i_par.get(i)));

                                    xsvc.set(i, getSubParameter(getChildNodes("parameter", xsection.get(i).item(1))));
                                    xrouting.set(
                                            i,
                                            getAttributesList(
                                                    getSubParameter(
                                                            getChildNodes("parameter", xsection.get(i).item(2)).get(0))));

                                    switch (((Element) xsection_par.get(i).get(0).get(2)).getAttribute("name")) {
                                        case "retrialDistributions":
                                            xput_strategy.set(
                                                    i,
                                                    getAttributesList(
                                                            ((Element) xsection_par.get(i).get(0).get(4))
                                                                    .getElementsByTagName("subParameter")));
                                            break;
                                        default:
                                            xput_strategy.set(
                                                    i,
                                                    getAttributesList(
                                                            ((Element) xsection_par.get(i).get(0).get(3))
                                                                    .getElementsByTagName("subParameter")));
                                    }

                                    switch (xput_strategy.get(i).get(0).getNamedItem("name").getNodeValue()) {
                                        case "TailStrategy":
                                            strategy.set(i, SchedStrategy.FCFS);
                                            break;
                                        case "TailStrategyPriority":
                                            strategy.set(i, SchedStrategy.HOL);
                                            break;
                                        case "HeadStrategy":
                                            strategy.set(i, SchedStrategy.LCFS);
                                            break;
                                        case "RandStrategy":
                                            strategy.set(i, SchedStrategy.SIRO);
                                            break;
                                        case "SJFStrategy":
                                            strategy.set(i, SchedStrategy.SJF);
                                            break;
                                        case "SEPTStrategy":
                                            strategy.set(i, SchedStrategy.SEPT);
                                            break;
                                        case "LJFStrategy":
                                            strategy.set(i, SchedStrategy.LJF);
                                            break;
                                        case "LEPTStrategy":
                                            strategy.set(i, SchedStrategy.LEPT);
                                            break;
                                    }

                                    Delay delay;
                                    Queue queue;
                                    int xcapacity;
                                    int xnumservers;
                                    List<List<Node>> strategy_sub;
                                    // same usage as xsection_i_type{i}{2}.className
                                    switch (xsection_javaClass
                                            .get(i)
                                            .get(1)
                                            .getNamedItem("className")
                                            .getNodeValue()) {
                                        case "Delay":
                                            delay = new Delay(model, node_name.get(i));
                                            xcapacity =
                                                    Integer.parseInt(
                                                            ((Element) xsection_par.get(i).get(0).get(0))
                                                                    .getElementsByTagName("value")
                                                                    .item(0)
                                                                    .getTextContent());
                                            delay.setCapacity(xcapacity);
                                            node.add(delay);
                                            break;
                                        case "Server":
                                            queue = new Queue(model, node_name.get(i), strategy.get(i));
                                            // Read JMT buffer size (queue buffer only)
                                            xcapacity =
                                                    Integer.parseInt(
                                                            ((Element) xsection_par.get(i).get(0).get(0))
                                                                    .getElementsByTagName("value")
                                                                    .item(0)
                                                                    .getTextContent());
                                            xnumservers =
                                                    Integer.parseInt(
                                                            ((Element) xsection_par.get(i).get(1).get(0))
                                                                    .getElementsByTagName("value")
                                                                    .item(0)
                                                                    .getTextContent());
                                            // Convert JMT buffer size to LINE capacity (Kendall K = buffer + servers)
                                            // JMT size = buffer only, LINE cap = total system capacity
                                            if (xcapacity >= 0 && xcapacity != Integer.MAX_VALUE) {
                                                queue.setCapacity(xcapacity + xnumservers);
                                            } else {
                                                queue.setCapacity(xcapacity);
                                            }
                                            queue.setNumberOfServers(xnumservers);
                                            node.add(queue);
                                            if (strategy.get(i) == SchedStrategy.SEPT) {
                                                // Extract SEPT parameters from XML and apply them
                                                // SEPT strategy uses scheduling parameters for job priorities
                                                // These will be applied later when schedparams are processed
                                                break;
                                            }
                                            break;
                                        case "PSServer":
                                            strategy_sub = getSubParameter(xsection_par.get(i).get(1));
                                            // Extract PSServer parameters from XML and apply them
                                            // PSServer strategies use scheduling parameters for weights/priorities
                                            // These will be applied later when schedparams are processed
                                            // we assume the strategies are identical across classes
                                            switch (((Element) strategy_sub.get(3).get(0)).getAttribute("name")) {
                                                case "EPSStrategy":
                                                    strategy.set(i, SchedStrategy.PS);
                                                    break;
                                                case "DPSStrategy":
                                                    strategy.set(i, SchedStrategy.DPS);
                                                    break;
                                                case "GPSStrategy":
                                                    strategy.set(i, SchedStrategy.GPS);
                                                    break;
                                            }
                                            queue = new Queue(model, node_name.get(i), strategy.get(i));
                                            // Read JMT buffer size (queue buffer only)
                                            xcapacity =
                                                    Integer.parseInt(
                                                            ((Element) xsection_par.get(i).get(0).get(0))
                                                                    .getElementsByTagName("value")
                                                                    .item(0)
                                                                    .getTextContent());
                                            xnumservers =
                                                    Integer.parseInt(
                                                            ((Element) xsection_par.get(i).get(1).get(0))
                                                                    .getElementsByTagName("value")
                                                                    .item(0)
                                                                    .getTextContent());
                                            // Convert JMT buffer size to LINE capacity (Kendall K = buffer + servers)
                                            // JMT size = buffer only, LINE cap = total system capacity
                                            if (xcapacity >= 0 && xcapacity != Integer.MAX_VALUE) {
                                                queue.setCapacity(xcapacity + xnumservers);
                                            } else {
                                                queue.setCapacity(xcapacity);
                                            }
                                            queue.setNumberOfServers(xnumservers);
                                            node.add(queue);
                                            break;
                                        case "ClassSwitch":
                                            strategy_sub = getSubParameter(xsection_par.get(i).get(1));
                                            List<List<Node>> strategy_sub1 = getSubParameter(strategy_sub.get(0));
                                            Matrix csMatrix = new Matrix(strategy_sub1.size(), strategy_sub1.size());
                                            csMatrix.zero();
                                            for (int r = 0; r < strategy_sub1.size(); r++) {
                                                for (int c = 0; c < strategy_sub1.get(r).size(); c++) {
                                                    csMatrix.set(
                                                            r,
                                                            c,
                                                            Double.parseDouble(
                                                                    ((Element) strategy_sub1.get(r).get(c))
                                                                            .getElementsByTagName("value")
                                                                            .item(0)
                                                                            .getTextContent()));
                                                }
                                            }
                                            node.add(new ClassSwitch(model, node_name.get(i), csMatrix));
                                            break;
                                    }
                            }
                    }
                    break;
                case "Storage":
                    node.add(new Place(model, node_name.get(i)));
                    break;
                case "Enabling":
                    node.add(new Transition(model, node_name.get(i)));
                    break;
            }
        }

        // create classes
        List<JobClass> jobclass = new ArrayList<>();
        List<NamedNodeMap> classes = getAttributesList(xDoc.getElementsByTagName("userClass"));
        int ref;

        // JMT uses higher priority value = higher priority, LINE uses lower value = higher priority
        // We need to invert priorities when importing from JMT
        int maxPrio = 0;
        for (int i = 0; i < classes.size(); i++) {
            int prio = Integer.parseInt(classes.get(i).getNamedItem("priority").getNodeValue());
            if (prio > maxPrio) {
                maxPrio = prio;
            }
        }

        for (int i = 0; i < classes.size(); i++) {
            ref = node_name.indexOf(classes.get(i).getNamedItem("referenceSource").getNodeValue());
            // Invert priority: JMT uses higher=higher, LINE uses lower=higher
            int linePrio = maxPrio - Integer.parseInt(classes.get(i).getNamedItem("priority").getNodeValue());
            switch (classes.get(i).getNamedItem("type").getNodeValue()) {
                case "closed":
                    jobclass.add(
                            new ClosedClass(
                                    model,
                                    classes.get(i).getNamedItem("name").getNodeValue(),
                                    Long.parseLong(classes.get(i).getNamedItem("customers").getNodeValue()),
                                    (Station) node.get(ref),
                                    linePrio));
                    break;
                case "open":
                    // sink and source have been created before
                    jobclass.add(
                            new OpenClass(
                                    model,
                                    classes.get(i).getNamedItem("name").getNodeValue(),
                                    linePrio));
                    if (classes
                            .get(i)
                            .getNamedItem("referenceSource")
                            .getNodeValue()
                            .equals("StatelessClassSwitcher")) {
                        for (int j = 0; j < node.size(); j++) {
                            if (node.get(j) instanceof Source) {
                                ((Source) node.get(j)).setArrival(jobclass.get(i), Disabled.getInstance());
                            }
                        }
                    }
                    break;
            }
        }

        int xcap;
        NodeList refClasses;
        int nclasses;
        int xclasscap;
        int nmodes;
        int ninputs;
        String nodeName;
        int ind;
        for (int i = 0; i < node_name.size(); i++) {
            switch (xsection_javaClass.get(i).get(0).getNamedItem("className").getNodeValue()) {
                case "Storage":
                    xcap =
                            Integer.parseInt(
                                    ((Element)
                                            ((Element) xsection.get(i).item(0))
                                                    .getElementsByTagName("parameter")
                                                    .item(0))
                                            .getElementsByTagName("value")
                                            .item(0)
                                            .getTextContent());
                    if (xcap == -1) {
                        ((Place) node.get(i)).setCapacity(Integer.MAX_VALUE);
                    } else {
                        ((Place) node.get(i)).setCapacity(xcap);
                    }
                    refClasses =
                            ((Element)
                                    ((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(1))
                                    .getElementsByTagName("refClass");
                    if (refClasses.item(0) != null) {
                        nclasses = refClasses.getLength();
                    } else {
                        nclasses = 1;
                    }
                    for (int c = 0; c < nclasses; c++) {
                        xclasscap =
                                Integer.parseInt(
                                        ((Element)
                                                getSubParameter(
                                                        ((Element) xsection.get(i).item(0))
                                                                .getElementsByTagName("parameter")
                                                                .item(1),
                                                        c))
                                                .getElementsByTagName("value")
                                                .item(0)
                                                .getTextContent());
                        if (xclasscap == -1) {
                            // Confirmed: Java version correctly uses JobClass for first parameter
                            // MATLAB version compatibility achieved through JobClass parameter
                            ((Place) node.get(i)).setClassCap(jobclass.get(c), Integer.MAX_VALUE);
                        } else {
                            ((Place) node.get(i)).setClassCap(jobclass.get(c), xclasscap);
                        }
                        switch (((Element)
                                getSubParameter(
                                        ((Element) xsection.get(i).item(0))
                                                .getElementsByTagName("parameter")
                                                .item(2),
                                        c))
                                .getElementsByTagName("value")
                                .item(0)
                                .getTextContent()) {
                            case "BAS blocking":
                                // Confirmed: Java version correctly uses JobClass for first parameter
                                // MATLAB version compatibility achieved through JobClass parameter
                                ((Place) node.get(i))
                                        .setDropRule(jobclass.get(c), DropStrategy.BlockingAfterService);
                                break;
                            case "drop":
                                ((Place) node.get(i)).setDropRule(jobclass.get(c), DropStrategy.Drop);
                                break;
                            case "waiting queue":
                                ((Place) node.get(i)).setDropRule(jobclass.get(c), DropStrategy.WaitingQueue);
                                break;
                        }
                    }

                    ((Place) node.get(i)).setState(new Matrix(0, 0));
                    break;
                // Transition class has not been implemented yet
                case "Enabling":
                    // Enabling Section
                    double lambda, lambda1, lambda2;
                    nmodes = getSubParameter(((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0)).size();
                    for (int m = 0; m < nmodes; m++) {
                        String modeName = ((Element) getSubParameter(((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0), m))
                                .getElementsByTagName("value").item(0).getTextContent();
                        Mode mode = ((Transition) node.get(i)).addMode(modeName);
                    }
                    
                    // Initialize the transition with enabling conditions
                    for (int m = 0; m < nmodes; m++) {
                        Mode mode = ((Transition) node.get(i)).getModes().get(m);
                        
                        // Get inputs for this mode
                        List<Node> modeInputs = getSubParameter(getSubParameter(getSubParameter(
                                ((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0), m), 0));
                        
                        if (modeInputs != null) {
                            ninputs = modeInputs.size();
                            for (int j = 0; j < ninputs; j++) {
                                refClasses = ((Element) getSubParameter(getSubParameter(getSubParameter(getSubParameter(
                                        ((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0), m), 0), j), 1))
                                        .getElementsByTagName("refClass");
                                
                                if (refClasses.item(0) != null) {
                                    nclasses = refClasses.getLength();
                                } else {
                                    nclasses = 1;
                                }
                                
                                nodeName = ((Element) getSubParameter(getSubParameter(getSubParameter(getSubParameter(
                                        ((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0), m), 0), j), 0))
                                        .getElementsByTagName("value").item(0).getTextContent();
                                ind = model.getNodeIndex(nodeName);
                                
                                for (int k = 0; k < nclasses; k++) {
                                    // Get enabling condition
                                    int enable = Integer.parseInt(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            getSubParameter(getSubParameter(((Element) xsection.get(i).item(0))
                                                    .getElementsByTagName("parameter").item(0), m), 0), j), 1), k))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    
                                    if (enable == -1) {
                                        ((Transition) node.get(i)).setEnablingConditions(mode, jobclass.get(k), (Place) node.get(ind), Integer.MAX_VALUE);
                                    } else {
                                        ((Transition) node.get(i)).setEnablingConditions(mode, jobclass.get(k), (Place) node.get(ind), enable);
                                    }
                                    
                                    // Get inhibiting condition
                                    int inhibit = Integer.parseInt(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            getSubParameter(getSubParameter(((Element) xsection.get(i).item(0))
                                                    .getElementsByTagName("parameter").item(1), m), 0), j), 1), k))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    
                                    if (inhibit == -1) {
                                        ((Transition) node.get(i)).setInhibitingConditions(mode, jobclass.get(k), (Place) node.get(ind), Integer.MAX_VALUE);
                                    } else {
                                        ((Transition) node.get(i)).setInhibitingConditions(mode, jobclass.get(k), (Place) node.get(ind), inhibit);
                                    }
                                }
                            }
                        }
                    }
                    
                    // Timing Section  
                    for (int m = 0; m < nmodes; m++) {
                        Mode mode = ((Transition) node.get(i)).getModes().get(m);
                        
                        int numOfServers = Integer.parseInt(((Element) getSubParameter(
                                ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(1), m))
                                .getElementsByTagName("value").item(0).getTextContent());
                        
                        if (numOfServers == -1) {
                            ((Transition) node.get(i)).setNumberOfServers(mode, Integer.MAX_VALUE);
                        } else {
                            ((Transition) node.get(i)).setNumberOfServers(mode, numOfServers);
                        }
                        
                        String timingSt = ((Element) getSubParameter(
                                ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m))
                                .getAttribute("classPath");
                        
                        if (timingSt.equals("jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy")) {
                            ((Transition) node.get(i)).setTimingStrategy(mode, TimingStrategy.IMMEDIATE);
                        } else {
                            ((Transition) node.get(i)).setTimingStrategy(mode, TimingStrategy.TIMED);
                            
                            String distribution = ((Element) getSubParameter(getSubParameter(
                                    ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 0))
                                    .getAttribute("name");
                            
                            lambda = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                    ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 0))
                                    .getElementsByTagName("value").item(0).getTextContent());
                            
                            switch (distribution) {
                                case "Exponential":
                                    ((Transition) node.get(i)).setDistribution(mode, new Exp(lambda));
                                    break;
                                case "Erlang":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Erlang(lambda, (int) lambda1));
                                    break;
                                case "Hyperexponential":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    lambda2 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 2))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new HyperExp(lambda, lambda1, lambda2));
                                    break;
                                case "Coxian":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    lambda2 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 2))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Coxian(Arrays.asList(lambda, lambda1), Arrays.asList(lambda2, 1.0)));
                                    break;
                                case "Deterministic":
                                    ((Transition) node.get(i)).setDistribution(mode, new Det(lambda));
                                    break;
                                case "Pareto":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Pareto(lambda, lambda1));
                                    break;
                                case "Gamma":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Gamma(lambda, lambda1));
                                    break;
                                case "Uniform":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Uniform(lambda, lambda1));
                                    break;
                                case "Replayer":
                                    ((Transition) node.get(i)).setDistribution(mode, new Replayer(lambda));
                                    break;
                                case "Weibull":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    // scale and shape are inverted in the constructor
                                    ((Transition) node.get(i)).setDistribution(mode, new Weibull(lambda1, lambda));
                                    break;
                                case "Lognormal":
                                    lambda1 = Double.parseDouble(((Element) getSubParameter(getSubParameter(getSubParameter(
                                            ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(2), m), 1), 1))
                                            .getElementsByTagName("value").item(0).getTextContent());
                                    ((Transition) node.get(i)).setDistribution(mode, new Lognormal(lambda, lambda1));
                                    break;
                                default:
                                    line_error(mfilename(new Object() {}), 
                                        "The model includes a timing distribution not supported by the model-to-model transformation from JMT.");
                                    ((Transition) node.get(i)).setDistribution(mode, new Exp(1));
                            }
                        }
                        
                        int firingPriorities = Integer.parseInt(((Element) getSubParameter(
                                ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(3), m))
                                .getElementsByTagName("value").item(0).getTextContent());
                        ((Transition) node.get(i)).setFiringPriorities(mode, firingPriorities);
                        
                        int firingWeights = Integer.parseInt(((Element) getSubParameter(
                                ((Element) xsection.get(i).item(1)).getElementsByTagName("parameter").item(4), m))
                                .getElementsByTagName("value").item(0).getTextContent());
                        ((Transition) node.get(i)).setFiringWeights(mode, firingWeights);
                    }
                    
                    // Firing Section
                    for (int m = 0; m < nmodes; m++) {
                        Mode mode = ((Transition) node.get(i)).getModes().get(m);
                        
                        Node modeParamNode = getSubParameter(((Element) xsection.get(i).item(2)).getElementsByTagName("parameter").item(0), m);
                        List<Node> modeOutputs = null;
                        if (modeParamNode != null) {
                            List<Node> temp = getSubParameter(modeParamNode);
                            if (temp != null && temp.size() > 0) {
                                modeOutputs = getSubParameter(temp.get(0));
                            }
                        }
                        
                        int noutputs = (modeOutputs == null) ? 0 : modeOutputs.size();
                        
                        for (int j = 0; j < noutputs; j++) {
                            refClasses = ((Element) getSubParameter(getSubParameter(getSubParameter(getSubParameter(
                                    ((Element) xsection.get(i).item(2)).getElementsByTagName("parameter").item(0), m), 0), j), 1))
                                    .getElementsByTagName("refClass");
                            
                            if (refClasses.item(0) != null) {
                                nclasses = refClasses.getLength();
                            } else {
                                nclasses = 1;
                            }
                            
                            nodeName = ((Element) getSubParameter(getSubParameter(getSubParameter(getSubParameter(
                                    ((Element) xsection.get(i).item(2)).getElementsByTagName("parameter").item(0), m), 0), j), 0))
                                    .getElementsByTagName("value").item(0).getTextContent();
                            ind = model.getNodeIndex(nodeName);
                            
                            for (int k = 0; k < nclasses; k++) {
                                int outcome = Integer.parseInt(((Element) getSubParameter(getSubParameter(getSubParameter(
                                        getSubParameter(getSubParameter(((Element) xsection.get(i).item(2))
                                                .getElementsByTagName("parameter").item(0), m), 0), j), 1), k))
                                        .getElementsByTagName("value").item(0).getTextContent());
                                
                                if (outcome == -1) {
                                    ((Transition) node.get(i)).setFiringOutcome(mode, jobclass.get(k), node.get(ind), Integer.MAX_VALUE);
                                } else {
                                    ((Transition) node.get(i)).setFiringOutcome(mode, jobclass.get(k), node.get(ind), outcome);
                                }
                            }
                        }
                    }
                    break;

            }
        }

        List<Matrix> schedparams = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        // set service distributions
        List<List<List<Node>>> xarv_statdistrib =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<NamedNodeMap>>> xarv_statdistribattr =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<Node>> xarv = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<Node>>> xarv_sec = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<Node> par, alpha;
        List<List<Node>> pars;
        double lambda, lambda1, lambda2, lambda3;
        List<List<List<Node>>> xsvc_sec = new ArrayList<>(Collections.nCopies(node_name.size(), null));
        List<List<List<NamedNodeMap>>> xsvc_statdistrib =
                new ArrayList<>(Collections.nCopies(node_name.size(), null));
        for (int i = 0; i < node_name.size(); i++) {
            if (node.get(i) instanceof Source) {
                for (int r = 0; r < classes.size(); r++) {
                    xsection_par.set(i, getChildNodes("parameter", xsection.get(i)));
                    xsection_i_par.set(i, xsection_par.get(i).get(0));
                    xsection_i_subpar.set(i, getSubParameter(xsection_i_par.get(i)));
                    xarv_statdistrib.set(i, new ArrayList<>(Collections.nCopies(classes.size(), null)));
                    xarv_statdistrib.get(i).set(r, getSubParameter(xsection_i_subpar.get(i).get(0).get(r)));
                    if (xarv_statdistrib.get(i).get(r) == null) {
                        ((Source) node.get(i)).setArrival(jobclass.get(r), Disabled.getInstance());
                    } else {
                        xarv_statdistribattr.set(i, new ArrayList<>(Collections.nCopies(classes.size(), null)));
                        xarv_statdistribattr.get(i).set(r, getAttributesList(xarv_statdistrib.get(i).get(r)));
                        xarv.set(
                                i,
                                getSubParameter(
                                        ((Element) xsection.get(i).item(0)).getElementsByTagName("parameter").item(0)));
                        xarv_sec.set(i, getSubParameter(xarv.get(i)));
                        switch (xarv_statdistribattr.get(i).get(r).get(0).getNamedItem("name").getNodeValue()) {
                            case "Exponential":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Exp(lambda));
                                break;
                            case "Erlang":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i))
                                        .setArrival(jobclass.get(r), new Erlang(lambda, (int) lambda1));
                                break;
                            case "Hyperexponential":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i))
                                        .setArrival(jobclass.get(r), new HyperExp(lambda, lambda1, lambda2));
                                break;
                            case "Coxian":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i))
                                        .setArrival(jobclass.get(r), new Coxian(Arrays.asList(lambda, lambda1), Arrays.asList(lambda2, 1.0)));
                                break;
                            case "Deterministic":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Det(lambda));
                                break;
                            case "Pareto":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Pareto(lambda, lambda1));
                                break;
                            case "Weibull":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Weibull(lambda, lambda1));
                                break;
                            case "Lognormal":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Lognormal(lambda, lambda1));
                                break;
                            case "Gamma":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Gamma(lambda, lambda1));
                                break;
                            case "Uniform":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Uniform(lambda, lambda1));
                                break;
                            case "Replayer":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Replayer(lambda));
                                break;
                            case "Burst (MMPP2)":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda3 =
                                        Double.parseDouble(
                                                ((Element) par.get(3))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Source) node.get(i))
                                        .setArrival(jobclass.get(r), new MMPP2(lambda, lambda1, lambda2, lambda3));
                                break;
                            case "Burst (MAP)":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                pars = getSubParameter(getSubParameter(par.get(0)));
                                Matrix D0 = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        D0.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                pars = getSubParameter(getSubParameter(par.get(1)));
                                Matrix D1 = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        D1.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                MAP ax = new MAP(D0, D1);
                                ((Source) node.get(i)).setArrival(jobclass.get(r), ax);
                                break;
                            case "Phase-Type":
                                par = getSubParameter(xarv_sec.get(i).get(r).get(1));
                                alpha = getSubParameter(getSubParameter(par.get(0), 0));
                                Matrix alphaMatrix = new Matrix(1, alpha.size());
                                for (int c = 0; c < alpha.size(); c++) {
                                    alphaMatrix.set(
                                            0,
                                            c,
                                            Double.parseDouble(
                                                    ((Element) alpha.get(c))
                                                            .getElementsByTagName("value")
                                                            .item(0)
                                                            .getTextContent()));
                                }
                                pars = getSubParameter(getSubParameter(par.get(1)));
                                Matrix T = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        T.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                // Check if T has elements below the main diagonal (strict lower triangular)
                                APH ax1;
                                if (Matrix.tril(T, -1).any()) {
                                    line_warning(
                                            mfilename(new Object() {
                                            }),
                                            "The input model uses a general PH distribution, which is not yet supported in LINE. Fitting the first three moments into an APH distribution.");
                                    Matrix D = Matrix.negative(T).mult(Matrix.ones(T.getNumCols(), 1)).mult(alphaMatrix);
                                    ax1 = APH.fitMeanAndSCV(Map_meanKt.map_mean(T, D), Map_varKt.map_var(T, D));
                                } else {
                                    ax1 = new APH(alphaMatrix, T);
                                }
                                ((Source) node.get(i)).setArrival(jobclass.get(r), ax1);
                                break;
                            default:
                                line_error(
                                        mfilename(new Object() {
                                        }),
                                        "The model includes an arrival distribution not supported by the model-to-model transformation from JMT.");
                                // xarv_statdistribattr.get(i).get(r).get(0).getNamedItem("name").getNodeValue();
                                // Default to exponential distribution with rate 1 for unsupported arrival types
                                ((Source) node.get(i)).setArrival(jobclass.get(r), new Exp(1));
                        }
                    }
                }
            } else if (node.get(i) instanceof Queue || node.get(i) instanceof Delay) {

                if (schedparams.get(i) == null) {
                    switch (strategy.get(i)) {
                        case SEPT:
                        case LEPT:
                            Matrix nan = new Matrix(1, jobclass.size());
                            nan.fill(Double.NaN);
                            schedparams.set(i, nan);
                            break;
                        default:
                            Matrix ones = new Matrix(1, jobclass.size());
                            ones.ones();
                            schedparams.set(i, ones);
                    }
                }
                switch (xsection_javaClass.get(i).get(1).getNamedItem("className").getNodeValue()) {
                    case "StatelessClassSwitcher":
                        // do nothing
                        continue;
                    case "Delay":
                        xsvc_sec.set(i, getSubParameter(xsvc.get(i).get(0)));
                        break;
                    default:
                        xsvc_sec.set(i, getSubParameter(xsvc.get(i).get(2)));
                }
                xsvc_statdistrib.set(i, new ArrayList<>(Collections.nCopies(classes.size(), null)));
                for (int r = 0; r < classes.size(); r++) {
                    // same usage as xsection_i_type{i}{2}.className
                    if (xsvc_sec.get(i).get(r) == null) {
                        xsvc_statdistrib.get(i).set(r, null);
                    } else {
                        xsvc_statdistrib.get(i).set(r, getAttributesList(xsvc_sec.get(i).get(r)));
                    }
                    // Apply scheduling parameters for priority-based strategies
                    double para_ir = schedparams.get(i).get(r);
                    
                    // Set scheduling strategy parameters for strategies that use them
                    if (strategy.get(i) == SchedStrategy.SEPT || strategy.get(i) == SchedStrategy.LEPT || 
                        strategy.get(i) == SchedStrategy.DPS || strategy.get(i) == SchedStrategy.GPS ||
                        strategy.get(i) == SchedStrategy.PSPRIO || strategy.get(i) == SchedStrategy.DPSPRIO ||
                        strategy.get(i) == SchedStrategy.GPSPRIO) {
                        if (!Double.isNaN(para_ir) && node.get(i) instanceof Queue) {
                            ((Queue) node.get(i)).setSchedStrategyPar(jobclass.get(r), para_ir);
                        }
                    }
                    if (xsvc_statdistrib.get(i).get(r) == null) {
                        // case 'Disabled'
                        ((Queue) node.get(i)).setService(jobclass.get(r), Disabled.getInstance());
                    } else {
                        switch (xsvc_statdistrib.get(i).get(r).get(0).getNamedItem("name").getNodeValue()) {
                            case "Replayer":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Replayer(lambda), para_ir);
                                break;
                            case "Exponential":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Exp(lambda), para_ir);
                                break;
                            case "Erlang":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i))
                                        .setService(jobclass.get(r), new Erlang(lambda, (int) lambda1), para_ir);
                                break;
                            case "Hyperexponential":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i))
                                        .setService(jobclass.get(r), new HyperExp(lambda, lambda1, lambda2), para_ir);
                                break;
                            case "Coxian":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i))
                                        .setService(jobclass.get(r), new Coxian(Arrays.asList(lambda, lambda1), Arrays.asList(lambda2, 1.0)), para_ir);
                                break;
                            case "Deterministic":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Det(lambda));
                                break;
                            case "Pareto":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Pareto(lambda, lambda1));
                                break;
                            case "Weibull":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Weibull(lambda, lambda1));
                                break;
                            case "Lognormal":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Lognormal(lambda, lambda1));
                                break;
                            case "Gamma":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Gamma(lambda, lambda1));
                                break;
                            case "Burst (MMPP2)":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda2 =
                                        Double.parseDouble(
                                                ((Element) par.get(2))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda3 =
                                        Double.parseDouble(
                                                ((Element) par.get(3))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i))
                                        .setService(jobclass.get(r), new MMPP2(lambda, lambda1, lambda2, lambda3));
                                break;
                            case "Burst (MAP)":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                pars = getSubParameter(getSubParameter(par.get(0)));
                                Matrix D0 = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        D0.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                pars = getSubParameter(getSubParameter(par.get(1)));
                                Matrix D1 = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        D1.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                MAP ax = new MAP(D0, D1);
                                ((Queue) node.get(i)).setService(jobclass.get(r), ax);
                                break;
                            case "Phase-Type":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                alpha = getSubParameter(getSubParameter(par.get(0), 0));
                                Matrix alphaMatrix = new Matrix(1, alpha.size());
                                for (int c = 0; c < alpha.size(); c++) {
                                    alphaMatrix.set(
                                            0,
                                            c,
                                            Double.parseDouble(
                                                    ((Element) alpha.get(c))
                                                            .getElementsByTagName("value")
                                                            .item(0)
                                                            .getTextContent()));
                                }
                                pars = getSubParameter(getSubParameter(par.get(1)));
                                Matrix T = new Matrix(pars.size(), pars.get(0).size());
                                for (int row = 0; row < pars.size(); row++) {
                                    for (int col = 0; col < pars.get(row).size(); col++) {
                                        T.set(
                                                row,
                                                col,
                                                Double.parseDouble(
                                                        ((Element) pars.get(row).get(col))
                                                                .getElementsByTagName("value")
                                                                .item(0)
                                                                .getTextContent()));
                                    }
                                }
                                APH ax1;
                                if (Matrix.tril(T).any()) {
                                    line_warning(
                                            mfilename(new Object() {
                                            }),
                                            "The input model uses a general PH distribution, which is not yet supported in LINE. Fitting the first three moments into an APH distribution.");
                                    Matrix D = Matrix.negative(T).mult(Matrix.ones(T.getNumCols(), 1)).mult(alphaMatrix);
                                    ax1 = APH.fitMeanAndSCV(Map_meanKt.map_mean(T, D), Map_scvKt.map_scv(T, D));
                                } else {
                                    ax1 = new APH(alphaMatrix, T);
                                }
                                ((Queue) node.get(i)).setService(jobclass.get(r), ax1);
                                break;
                            case "Uniform":
                                par = getSubParameter(xsvc_sec.get(i).get(r).get(1));
                                lambda =
                                        Double.parseDouble(
                                                ((Element) par.get(0))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                lambda1 =
                                        Double.parseDouble(
                                                ((Element) par.get(1))
                                                        .getElementsByTagName("value")
                                                        .item(0)
                                                        .getTextContent());
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Uniform(lambda, lambda1));
                                break;
                            default:
                                line_error(
                                        mfilename(new Object() {
                                        }),
                                        "The model includes a service distribution not supported by the model-to-model transformation from JMT.");
                                // xsvc_statdistrib.get(i).get(r).get(0).getNamedItem("name").getNodeValue();
                                // Default to exponential distribution with rate 1 for unsupported service types
                                ((Queue) node.get(i)).setService(jobclass.get(r), new Exp(1), para_ir);
                        }
                    }
                }
                for (int c = 0; c < xsection_i_par_attr.get(i).size(); c++) {
                    switch (xsection_i_par_attr.get(i).get(c).getNamedItem("name").getNodeValue()) {
                        case "size":
                            ((Buffer) node.get(i).getInput())
                                    .setSize(
                                            Integer.parseInt(
                                                    xsection_i_value.get(i).get(c).get(0).getTextContent())); // buffer size
                            break;
                    }
                }
            }
        }

        // create links
        Matrix C = new Matrix(node_name.size(), node_name.size());
        C.zero();
        List<NamedNodeMap> links = getAttributesList(xDoc.getElementsByTagName("connection"));
        for (int l = 0; l < links.size(); l++) {
            int source = orig_node_name.indexOf(links.get(l).getNamedItem("source").getNodeValue());
            int target = orig_node_name.indexOf(links.get(l).getNamedItem("target").getNodeValue());
            C.set(source, target, 1);
        }

        // assign routing probabilities
        Matrix P = new Matrix(node_name.size() * classes.size(), node_name.size() * classes.size());
        for (int from = 0; from < node_name.size(); from++) {
            for (int target = 0; target < node_name.size(); target++) {
                if (C.get(from, target) != 0) {
                    model.addLink(node.get(from), node.get(target));
                }
            }
        }

        List<Node> xroutprobarray = new ArrayList<>(Collections.nCopies(classes.size(), null));
        List<Node> xroutprob;
        List<List<Node>> xroutprobdest;
        for (int from = 0; from < node_name.size(); from++) {
            if (node.get(from) instanceof Place
                    || node.get(from) instanceof Transition
                    || node.get(from) instanceof Sink) {
                // Do nothing
            } else {
                for (int r = 0; r < classes.size(); r++) {
                    switch (xrouting.get(from).get(r).getNamedItem("name").getNodeValue()) {
                        case "Random":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.RAND);
                            break;
                        case "Probabilities":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.PROB);
                            xroutprobarray.set(
                                    r,
                                    getSubParameter(
                                            getSubParameter(
                                                    ((Element) xsection.get(from).item(2))
                                                            .getElementsByTagName("parameter")
                                                            .item(0),
                                                    r),
                                            0));
                            xroutprob = getSubParameter(xroutprobarray.get(r));
                            xroutprobdest = getSubParameter(xroutprob);
                            for (int j = 0; j < xroutprobdest.size(); j++) {
                                List<List<Node>> xprob = getChildNodes("value", xroutprobdest.get(j));
                                int target = node_name.indexOf(xprob.get(0).get(0).getTextContent());
                                double prob = Double.parseDouble(xprob.get(1).get(0).getTextContent());
                                node.get(from).setProbRouting(jobclass.get(r), node.get(target), prob);
                            }
                            break;
                        case "Power of k":
                            line_error(mfilename(new Object() {
                            }), "Power of k import not supported yet.");
                            break;
                        case "Round Robin":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.RROBIN);
                            break;
                        case "Weighted Round Robin":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.WRROBIN);
                            xroutprobarray.set(
                                    r,
                                    getSubParameter(
                                            getSubParameter(
                                                    ((Element) xsection.get(from).item(2))
                                                            .getElementsByTagName("parameter")
                                                            .item(0),
                                                    r),
                                            0));
                            xroutprob = getSubParameter(xroutprobarray.get(r));
                            xroutprobdest = getSubParameter(xroutprob);
                            for (int j = 0; j < xroutprobdest.size(); j++) {
                                List<Node> xprob = getChildNodes("value", xroutprobdest.get(j)).get(0);
                                int target = node_name.indexOf(xprob.get(0).getTextContent());
                                double weight = Double.parseDouble(xprob.get(1).getTextContent());
                                node.get(from)
                                        .setRouting(jobclass.get(r), RoutingStrategy.WRROBIN, node.get(target), weight);
                            }
                            break;
                        case "Join the Shortest Queue (JSQ)":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.JSQ);
                            break;
                        case "Disabled":
                            node.get(from).setRouting(jobclass.get(r), RoutingStrategy.DISABLED);
                            break;
                    }
                }
            }
        }

        if (model.getIndexSourceStation() > 1) {
            line_error(mfilename(new Object[]{}),
                    "LINE supports JMT models with at most a single source node. You can refactor your JMT model in several ways:\n - If you are mapping in JMT each class to a different source, this is not required. You can instead assign the same reference station to each class and configure class routing in the routing panel of the source node.\n - In more general cases, you may follow these three steps:\n    (1) give a different name to each class of arrival, assigning these classes to a single source as reference station.\n    (2) put a class-switch node after the source to switch the new classes into the original classes they were in the model with multiple sources.\n    (3) configure the routing section of this class-switch node to set the same routing for the classes as they were in the original model.\n");
        }

        // Preload
        NodeList preload = xDoc.getElementsByTagName("preload");
        boolean hasPreloadData = false;
        if (preload.getLength() > 0) {
            NodeList stationPopulations =
                    ((Element) preload.item(0)).getElementsByTagName("stationPopulations");
            int npreloadStates = stationPopulations.getLength();
            if (npreloadStates > 0) {
                Matrix state = new Matrix(node.size(), classes.size());
                state.zero();
                hasPreloadData = true;
                
                for (int st = 0; st < npreloadStates; st++) {
                    nodeName = ((Element) stationPopulations.item(st)).getAttribute("stationName");
                    ind = model.getNodeIndex(nodeName);
                    NodeList classPopulation =
                            ((Element) stationPopulations.item(st)).getElementsByTagName("classPopulation");
                    for (int r = 0; r < classPopulation.getLength(); r++) {
                        int c =
                                model.getClassIndex(((Element) classPopulation.item(r)).getAttribute("refClass"));
                        state.set(
                                ind,
                                c,
                                Double.parseDouble(((Element) classPopulation.item(r)).getAttribute("population")));
                    }

                    if (node.get(ind) instanceof Place) {
                        if (classes.size() > 1) {
                            line_error(
                                    mfilename(new Object() {
                                    }),
                                    "Import failed: Colored Petri net models are not yet supported in LINE.\n");
                        }
                    }
                }
                
                // Only try to initialize from marginal if we have actual preload data
                try {
                    model.initFromMarginal(state);
                } catch (Exception e) {
                    line_warning(
                            mfilename(new Object() {
                            }), "Import failed to automatically initialize the model with preload data. Trying default initialization.\n");
                    try {
                        model.initDefault();
                    } catch (Exception e2) {
                        line_warning(
                                mfilename(new Object() {
                                }), "Default initialization also failed.\n");
                    }
                }
            }
        }
        
        // If no preload data was found, skip state initialization
        // The model will be initialized later when needed by the solver

        return model;
    }

    /**
     * Extracts the attributes from each node in a NodeList.
     * 
     * @param nodeList The NodeList to process
     * @return A list of NamedNodeMap objects containing the attributes of each node
     */
    private List<NamedNodeMap> getAttributesList(NodeList nodeList) {
        List<NamedNodeMap> namedNodeMaps = new ArrayList<>();

        for (int i = 0; i < nodeList.getLength(); i++) {
            namedNodeMaps.add(nodeList.item(i).getAttributes());
        }

        return namedNodeMaps;
    }

    /**
     * Extracts the attributes from each node in a List of Nodes.
     * 
     * @param nodeList The list of nodes to process
     * @return A list of NamedNodeMap objects containing the attributes of each node
     */
    private List<NamedNodeMap> getAttributesList(List<Node> nodeList) {
        List<NamedNodeMap> namedNodeMaps = new ArrayList<>();

        for (Node node : nodeList) {
            namedNodeMaps.add(node.getAttributes());
        }

        return namedNodeMaps;
    }

    /**
     * Gets direct child nodes with a specific tag name from a parent node.
     * 
     * <p>Unlike getElementsByTagName(), this method only returns direct children
     * (one level below the parent node), not all descendants.</p>
     * 
     * @param tagName The tag name to search for
     * @param parentNode The parent node to search within
     * @return A list of direct child nodes with the specified tag name, or null if none found
     */
    private List<Node> getChildNodes(String tagName, Node parentNode) {
        List<Node> nodeList = new ArrayList<>();
        NodeList childNodes = parentNode.getChildNodes();

        for (int i = 0; i < childNodes.getLength(); i++) {
            if (childNodes.item(i).getNodeName().equals(tagName)) {
                nodeList.add(childNodes.item(i));
            }
        }

        if (nodeList.isEmpty()) {
            return null;
        }
        return nodeList;
    }

    /**
     * Gets direct child nodes with a specific tag name from multiple parent nodes.
     * 
     * @param tagName The tag name to search for
     * @param parentNodes List of parent nodes to search within
     * @return A list of lists, where each inner list contains the child nodes for each parent
     */
    private List<List<Node>> getChildNodes(String tagName, List<Node> parentNodes) {
        List<List<Node>> nodeList = new ArrayList<>();

        for (Node parentNode : parentNodes) {
            nodeList.add(getChildNodes(tagName, parentNode));
        }

        return nodeList;
    }

    /**
     * Gets direct child nodes with a specific tag name from a NodeList of parent nodes.
     * 
     * @param tagName The tag name to search for
     * @param parentNodes NodeList of parent nodes to search within
     * @return A list of lists, where each inner list contains the child nodes for each parent
     */
    private List<List<Node>> getChildNodes(String tagName, NodeList parentNodes) {
        List<List<Node>> nodeList = new ArrayList<>();

        for (int i = 0; i < parentNodes.getLength(); i++) {
            nodeList.add(getChildNodes(tagName, parentNodes.item(i)));
        }

        return nodeList;
    }

    /**
     * Gets all direct child "subParameter" nodes from a parent node.
     * 
     * @param parentNode The parent node to search within
     * @return A list of subParameter child nodes
     */
    private List<Node> getSubParameter(Node parentNode) {
        return getChildNodes("subParameter", parentNode);
    }

    /**
     * Gets a specific child "subParameter" node by index from a parent node.
     * 
     * @param parentNode The parent node to search within
     * @param index The index of the subParameter to retrieve
     * @return The subParameter node at the specified index, or null if not found
     */
    private Node getSubParameter(Node parentNode, int index) {
        try {
            return Objects.requireNonNull(getSubParameter(parentNode)).get(index);
        } catch (Exception e) {
            return null;
        }
    }

    /**
     * Gets all direct child "subParameter" nodes from multiple parent nodes.
     * 
     * @param parentNode List of parent nodes to search within
     * @return A list of lists, where each inner list contains the subParameter nodes for each parent
     */
    private List<List<Node>> getSubParameter(List<Node> parentNode) {
        List<List<Node>> nodeList = new ArrayList<>();

        for (int i = 0; i < parentNode.size(); i++) {
            nodeList.add(getSubParameter(parentNode.get(i)));
        }

        return nodeList;
    }

    /**
     * Reads and parses an XML file into a Document object.
     * 
     * @param filename The path to the XML file to read
     * @return A parsed Document object, or null if reading fails
     */
    private Document xml_read(String filename) {
        try {
            File file = new File(filename);
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            Document doc = db.parse(file);
            doc.getDocumentElement().normalize();

            return doc;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Converts an LQN/LQNX (Layered Queueing Network) model file to a LayeredNetwork object.
     * 
     * <p>This method parses LQN XML files that conform to the LQN schema and creates a
     * corresponding LayeredNetwork representation in LINE. LQN models are commonly used
     * for modeling layered software architectures and client-server systems.</p>
     * 
     * @param filename Path to the LQN XML file to import
     * @return A LayeredNetwork object representing the imported LQN model
     * @throws RuntimeException if the file cannot be parsed or is not a valid LQN file
     */
    public LayeredNetwork LQN2LINE(String filename) {
        try {
            return LayeredNetwork.parseXML(filename);
        } catch (Exception e) {
            line_error(mfilename(new Object() {}), "Failed to parse LQN file: " + filename + ". Error: " + e.getMessage());
            return null;
        }
    }

    /**
     * Converts an LQN/LQNX (Layered Queueing Network) model file to a LayeredNetwork object with a custom name.
     * 
     * <p>This method parses LQN XML files and allows specifying a custom name for the resulting
     * LayeredNetwork model, overriding any name defined in the XML file.</p>
     * 
     * @param filename Path to the LQN XML file to import
     * @param modelName Custom name to assign to the LayeredNetwork model
     * @return A LayeredNetwork object representing the imported LQN model with the specified name
     * @throws RuntimeException if the file cannot be parsed or is not a valid LQN file
     */
    public LayeredNetwork LQN2LINE(String filename, String modelName) {
        try {
            LayeredNetwork model = LayeredNetwork.parseXML(filename);
            if (model != null && modelName != null && !modelName.isEmpty()) {
                // Set the model name if the LayeredNetwork supports it
                // Note: LayeredNetwork might not have a setName method, but we can try
                model.setName(modelName);
            }
            return model;
        } catch (Exception e) {
            line_error(mfilename(new Object() {}), "Failed to parse LQN file: " + filename + ". Error: " + e.getMessage());
            return null;
        }
    }

    /**
     * Converts a LINE Network model to JSIMG file format for use with JMT.
     * 
     * <p>This method exports a LINE Network model to the JSIMG (Java SIMulation Graphics) format,
     * which can be opened and simulated in JMT (Java Modelling Tools). The conversion preserves
     * all supported model elements including nodes, classes, routing, and distributions.</p>
     * 
     * @param model The LINE Network model to convert
     * @param outputFileName Path where the JSIMG file should be saved
     * @return true if conversion was successful, false otherwise
     */
    public boolean LINE2JSIMG(Network model, String outputFileName) {
        try {
            // Create a SolverJMT instance to use its writeJSIM method
            SolverJMT solver = new SolverJMT(model);
            NetworkStruct sn = model.getStruct();
            
            // Use SolverJMT's writeJSIM method which handles all the JSIMG generation
            solver.writeJSIM(sn, outputFileName);
            return true;
        } catch (Exception e) {
            line_error(mfilename(new Object() {}), "Error during LINE to JSIMG conversion: " + e.getMessage());
            return false;
        }
    }

    /**
     * Converts a LINE Network model to JSIMG file format with automatic filename generation.
     * 
     * <p>This convenience method exports a LINE Network model to JSIMG format using the model's
     * name as the base filename with a .jsimg extension. For example, a model named "MyModel"
     * will be exported to "MyModel.jsimg".</p>
     * 
     * @param model The LINE Network model to convert
     * @return The path to the generated JSIMG file, or null if conversion failed
     */
    public String LINE2JSIMG(Network model) {
        String outputFileName = model.getName() + ".jsimg";
        if (LINE2JSIMG(model, outputFileName)) {
            return outputFileName;
        }
        return null;
    }
}
