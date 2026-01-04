package jline.solvers.jmt;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;

/**
 * Test fixtures and helper methods for SolverJMT tests.
 *
 * Contains utility methods shared across test methods including:
 * - XML parsing helpers for JSIM files
 * - Parameter extraction from generated XML
 */
public class SolverJMTTestFixtures {

    /**
     * Helper method to parse generated JSIM XML file
     */
    protected Document parseJsimXml(String jsimPath) throws ParserConfigurationException, IOException, SAXException {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        return builder.parse(new File(jsimPath));
    }

    /**
     * Helper method to find impatience parameter node for a specific queue and class
     */
    protected Element findImpatienceParameter(Document doc, String queueName, String className) {
        NodeList nodes = doc.getElementsByTagName("node");
        for (int i = 0; i < nodes.getLength(); i++) {
            Element node = (Element) nodes.item(i);
            if (node.getAttribute("name").equals(queueName)) {
                NodeList sections = node.getElementsByTagName("section");
                for (int j = 0; j < sections.getLength(); j++) {
                    Element section = (Element) sections.item(j);
                    if (section.getAttribute("className").equals("Queue")) {
                        NodeList parameters = section.getElementsByTagName("parameter");
                        for (int k = 0; k < parameters.getLength(); k++) {
                            Element param = (Element) parameters.item(k);
                            if (param.getAttribute("name").equals("Impatience")) {
                                // Find the refClass for this class
                                NodeList refClasses = param.getElementsByTagName("refClass");
                                for (int m = 0; m < refClasses.getLength(); m++) {
                                    Element refClass = (Element) refClasses.item(m);
                                    if (refClass.getTextContent().equals(className)) {
                                        // Found the right class, get the next subParameter
                                        NodeList subParams = param.getElementsByTagName("subParameter");
                                        if (subParams.getLength() > m) {
                                            return (Element) subParams.item(m);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return null;
    }
}
