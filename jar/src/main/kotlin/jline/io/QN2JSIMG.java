/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.state.State;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Logger;
import jline.lang.nodes.Node;
import jline.lang.sections.Section;
import jline.solvers.jmt.handlers.SaveHandlers;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static jline.lang.state.ToMarginal.toMarginal;

/**
 * Writes a LINE Network model to JMT JSIMG format.
 *
 * <p>This class provides static utility methods for exporting queueing network
 * models to JMT's JSIM graphical format (.jsimg files). The JSIMG format is
 * compatible with JMT's simulation and analysis tools.
 *
 * <p>Example usage:
 * <pre>
 * Network model = new Network("example");
 * // ... define model ...
 * String filename = QN2JSIMG.writeJSIM(model);
 * </pre>
 *
 * @see jline.lang.Network
 * @see SaveHandlers
 */
public class QN2JSIMG {

    /**
     * Writes a Network model to JSIMG format at a temporary location.
     *
     * @param model The Network model to export
     * @return Path to the generated JSIMG file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public static String writeJSIM(Network model) throws ParserConfigurationException {
        return writeJSIM(model, model.getStruct(), null);
    }

    /**
     * Writes a Network model to JSIMG format at the specified location.
     *
     * @param model The Network model to export
     * @param outputFileName The output file path (null for temp file)
     * @return Path to the generated JSIMG file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public static String writeJSIM(Network model, String outputFileName) throws ParserConfigurationException {
        return writeJSIM(model, model.getStruct(), outputFileName);
    }

    /**
     * Writes a Network model with a specific structure to JSIMG format.
     *
     * @param model The Network model to export
     * @param sn The network structure to use
     * @param outputFileName The output file path (null for temp file)
     * @return Path to the generated JSIMG file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public static String writeJSIM(Network model, NetworkStruct sn, String outputFileName) throws ParserConfigurationException {
        SaveHandlers saveHandlers = new SaveHandlers(model);
        return writeJSIM(model, sn, outputFileName, saveHandlers);
    }

    /**
     * Core implementation: Writes a Network model to JSIMG format using custom SaveHandlers.
     *
     * @param model The Network model to export
     * @param sn The network structure to use
     * @param outputFileName The output file path (null for temp file)
     * @param saveHandlers The SaveHandlers instance to use for XML generation
     * @return Path to the generated JSIMG file
     * @throws ParserConfigurationException if XML parsing fails
     */
    public static String writeJSIM(Network model, NetworkStruct sn, String outputFileName, SaveHandlers saveHandlers)
            throws ParserConfigurationException {

        // Generate default temp path if not provided
        if (outputFileName == null || outputFileName.isEmpty()) {
            try {
                String filePath = SysUtilsKt.lineTempName("jsim");
                outputFileName = filePath + File.separator + "model.jsimg";
            } catch (IOException ioe) {
                throw new RuntimeException("Unable to create temp path for JSIM file", ioe);
            }
        }

        // Update SaveHandlers with the current network structure
        saveHandlers.updateNetworkStruct(sn);

        ElementDocumentPair xml = saveHandlers.saveXMLHeader(model.getLogPath());
        xml = saveHandlers.saveClasses(xml);
        int numOfClasses = sn.nclasses;
        int numOfNodes = sn.nnodes;

        for (int i = 0; i < numOfNodes; i++) {
            Node currentNode = model.getNodes().get(i);
            org.w3c.dom.Element node = xml.simDoc.createElement("node");
            node.setAttribute("name", currentNode.getName());
            List<Section> nodeSections = Arrays.asList(currentNode.getInput(), currentNode.getServer(), currentNode.getOutput());
            for (int j = 0; j < nodeSections.size(); j++) {
                org.w3c.dom.Element xml_section = xml.simDoc.createElement("section");
                Section currentSection = nodeSections.get(j);
                // For Logger nodes, we need to include Generic sections too
                boolean isLogger = currentNode instanceof Logger;
                if (currentSection != null && (isLogger || !currentSection.getClassName().startsWith("Generic "))) {
                    xml_section.setAttribute("className", currentSection.getClassName());

                    // Override className for SRPT strategies - JMT requires PreemptiveServer
                    if (currentSection.getClassName().equals("Server") && currentNode instanceof jline.lang.nodes.Queue) {
                        jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) currentNode;
                        SchedStrategy sched = queue.getSchedStrategy();
                        if (sched == SchedStrategy.SRPT || sched == SchedStrategy.SRPTPRIO) {
                            xml_section.setAttribute("className", "PreemptiveServer");
                        }
                    }

                    DocumentSectionPair simXML = new DocumentSectionPair(xml.simDoc, xml_section);
                    switch (currentSection.getClassName()) {
                        case "Buffer":
                            simXML.section.setAttribute("className", "Queue"); // overwrite with JMT class name
                            simXML = saveHandlers.saveBufferCapacity(simXML, i);
                            simXML = saveHandlers.saveDropStrategy(simXML, i); // unfinished
                            simXML = saveHandlers.saveGetStrategy(simXML, i);
                            simXML = saveHandlers.savePutStrategy(simXML, i);
                            simXML = saveHandlers.saveImpatience(simXML, i);
                            break;
                        case "Server":
                        case "jline.Server":
                        case "PreemptiveServer":
                            simXML = saveHandlers.saveNumberOfServers(simXML, i);
                            simXML = saveHandlers.saveServerVisits(simXML);
                            simXML = saveHandlers.saveServiceStrategy(simXML, i);
                            simXML = saveHandlers.saveDelayOffStrategy(simXML, i);
                            break;
                        case "PollingServer":
                            simXML = saveHandlers.setPollingServerClassName(simXML, i);
                            simXML = saveHandlers.saveNumberOfServers(simXML, i);
                            simXML = saveHandlers.saveServerVisits(simXML);
                            simXML = saveHandlers.saveServiceStrategy(simXML, i);
                            simXML = saveHandlers.saveSwitchoverStrategy(simXML, i);
                            break;
                        case "SharedServer":
                            simXML.section.setAttribute("className", "PSServer"); // overwrite with JMT class name
                            simXML = saveHandlers.saveNumberOfServers(simXML, i);
                            simXML = saveHandlers.saveServerVisits(simXML);
                            simXML = saveHandlers.saveServiceStrategy(simXML, i);
                            simXML = saveHandlers.saveDelayOffStrategy(simXML, i);
                            simXML = saveHandlers.savePreemptiveStrategy(simXML, i);
                            simXML = saveHandlers.savePreemptiveWeights(simXML, i);
                            break;
                        case "InfiniteServer":
                        case "jline.InfiniteServer":
                            simXML.section.setAttribute("className", "Delay"); // overwrite with JMT class name
                            simXML = saveHandlers.saveServiceStrategy(simXML, i);
                            break;
                        case "RandomSource":
                        case "jline.RandomSource":
                            simXML = saveHandlers.saveArrivalStrategy(simXML, i);
                            break;
                        case "Dispatcher":
                        case "jline.Dispatcher":
                        case "ClassSwitchDispatcher":
                            simXML.section.setAttribute("className", "Router"); // overwrite with JMT class name
                            simXML = saveHandlers.saveRoutingStrategy(simXML, i);
                            break;
                        case "StatelessClassSwitcher":
                        case "jline.StatelessClassSwitcher":
                            simXML.section.setAttribute("className", "ClassSwitch"); // overwrite with JMT class name
                            simXML = saveHandlers.saveClassSwitchStrategy(simXML, i);
                            break;
                        case "Cache":
                            simXML.section.setAttribute("className", "Cache"); // overwrite with JMT class name
                            simXML = saveHandlers.saveCacheStrategy(simXML, i);
                            break;
                        case "LogTunnel":
                            simXML = saveHandlers.saveLogTunnel(simXML, i);
                            break;
                        case "Joiner":
                            simXML.section.setAttribute("className", "Join"); // overwrite with JMT class name
                            simXML = saveHandlers.saveJoinStrategy(simXML, i);
                            break;
                        case "Forker":
                            simXML.section.setAttribute("className", "Fork"); // overwrite with JMT class name
                            simXML = saveHandlers.saveForkStrategy(simXML, i);
                            break;
                        case "Storage":
                            simXML.section.setAttribute("className", "Storage"); // overwrite with JMT class name
                            simXML = saveHandlers.saveTotalCapacity(simXML, i);
                            simXML = saveHandlers.savePlaceCapacities(simXML, i);
                            simXML = saveHandlers.saveDropRule(simXML, i);
                            simXML = saveHandlers.saveGetStrategy(simXML);
                            simXML = saveHandlers.savePutStrategies(simXML, i);
                            break;
                        case "Enabling":
                            simXML.section.setAttribute("className", "Enabling"); // overwrite with JMT class name
                            simXML = saveHandlers.saveEnablingConditions(simXML, i);
                            simXML = saveHandlers.saveInhibitingConditions(simXML, i);
                            break;
                        case "Firing":
                            simXML.section.setAttribute("className", "Firing"); // overwrite with JMT class name
                            simXML = saveHandlers.saveFiringOutcomes(simXML, i);
                            break;
                        case "Timing":
                            simXML.section.setAttribute("className", "Timing"); // overwrite with JMT class name
                            simXML = saveHandlers.saveModeNames(simXML, i);
                            simXML = saveHandlers.saveNumbersOfServers(simXML, i);
                            simXML = saveHandlers.saveTimingStrategies(simXML, i);
                            simXML = saveHandlers.saveFiringPriorities(simXML, i);
                            simXML = saveHandlers.saveFiringWeights(simXML, i);
                            break;
                        case "Generic Input":
                            // For Logger nodes, create a Queue section for the input
                            if (currentNode instanceof Logger) {
                                simXML.section.setAttribute("className", "Queue");
                                simXML = saveHandlers.saveBufferCapacity(simXML, i);
                                simXML = saveHandlers.saveDropStrategy(simXML, i);
                                simXML = saveHandlers.saveGetStrategy(simXML, i);
                                simXML = saveHandlers.savePutStrategy(simXML, i);
                                simXML = saveHandlers.saveImpatience(simXML, i);
                            }
                            break;
                        case "Generic Output":
                            // For Logger nodes, output section is handled by Router section below
                            if (currentNode instanceof Logger) {
                                continue; // Skip adding this section, Router will be added separately
                            }
                            break;
                    }
                    node.appendChild(simXML.section);
                }
            }

            // Special handling for Logger nodes - they need Router sections even with Generic Output
            if (currentNode instanceof Logger) {
                org.w3c.dom.Element routerSection = xml.simDoc.createElement("section");
                routerSection.setAttribute("className", "Router");
                DocumentSectionPair routerXML = new DocumentSectionPair(xml.simDoc, routerSection);
                routerXML = saveHandlers.saveRoutingStrategy(routerXML, i);
                node.appendChild(routerXML.section);
            }

            xml.simElem.appendChild(node);
        }
        xml = saveHandlers.saveMetrics(xml);
        xml = saveHandlers.saveLinks(xml);
        xml = saveHandlers.saveRegions(xml);

        boolean hasReferenceNodes = false;
        org.w3c.dom.Element preloadNode = xml.simDoc.createElement("preload");
        Map<jline.lang.nodes.StatefulNode, Matrix> s0 = sn.state;
        int numOfStations = sn.nstations;

        for (int i = 0; i < numOfStations; i++) {
            boolean isReferenceNode = false;
            int nodeIndex = (int) sn.stationToNode.get(i);
            int isf = (int) sn.stationToStateful.get(i);
            org.w3c.dom.Element stationPopulationsNode = null;

            if (sn.nodetype.get(nodeIndex) != NodeType.Source && sn.nodetype.get(nodeIndex) != NodeType.Join) {
                State.StateMarginalStatistics sms = toMarginal(sn, nodeIndex, s0.get(model.getStatefulNodes().get(isf)), null, null, null, null, null);
                stationPopulationsNode = xml.simDoc.createElement("stationPopulations");
                stationPopulationsNode.setAttribute("stationName", sn.nodenames.get(nodeIndex));

                // TODO: here we assume that the current state is the one on top of the state data structure
                // however it could be that stateprior places the mass on another (or multiple other) states
                for (int r = 0; r < numOfClasses; r++) {
                    org.w3c.dom.Element classPopulationNode = xml.simDoc.createElement("classPopulation");

                    if (Double.isInfinite(sn.njobs.get(r)) || sn.njobs.get(r) == Integer.MAX_VALUE) {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(0, r))));
                        classPopulationNode.setAttribute("refClass", sn.classnames.get(r));
                        stationPopulationsNode.appendChild(classPopulationNode);
                    } else {
                        isReferenceNode = true;
                        classPopulationNode.setAttribute("population", String.valueOf(Math.round(sms.nir.get(0, r))));
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
            InputOutputKt.writeXML(outputFileName, xml.simDoc);
        } catch (Exception e) {
            e.printStackTrace();
            try {
                InputOutputKt.writeXML(outputFileName, xml.simDoc);
            } catch (Exception retryException) {
                retryException.printStackTrace();
            }
        }
        return outputFileName;
    }
}
