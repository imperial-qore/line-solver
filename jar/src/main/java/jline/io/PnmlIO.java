/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodes.Place;
import jline.lang.nodes.Station;
import jline.lang.nodes.Transition;
import jline.lang.processes.Det;
import jline.lang.processes.Disabled;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Gamma;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Lognormal;
import jline.lang.processes.Pareto;
import jline.lang.processes.Uniform;
import jline.lang.processes.Weibull;
import jline.util.NamedParam;
import jline.util.matrix.Matrix;

/**
 * PNML (ISO/IEC 15909-2) place/transition nets, read and written.
 *
 * <p>Java twin of matlab/src/io/pnml_save.m and pnml_load.m. The grammar written is
 * http://www.pnml.org/version-2009/grammar/ptnet, so that a LINE net can be read by the
 * tools built around that corpus (GreatSPN, TINA, the Model Checking Contest harnesses)
 * and a net from that corpus can be analysed here.</p>
 *
 * <p>THE P/T GRAMMAR IS UNCOLOURED, so what it can carry is narrower than what LINE can
 * express, and the difference is REFUSED rather than approximated: more than one job
 * class, an open class or a Source/Sink, a queueing place, a firing-rate dependence, and
 * any distribution outside the scalar-parameter families listed in {@link #distToXml}.</p>
 *
 * <p>TIMING RIDES IN A TOOLSPECIFIC BLOCK, which is where the grammar puts what it does
 * not define. Each LINE MODE becomes one PNML transition, so that the arcs of a mode are
 * the arcs of a transition as the grammar requires; the block records which LINE
 * transition and mode the PNML transition came from, so the reader regroups the modes the
 * writer split. A reader that ignores the block still sees a correct untimed P/T net, and
 * a P/T net with no such block is read with every transition TIMED and EXPONENTIAL AT
 * RATE 1, the convention of the stochastic Petri net literature.</p>
 *
 * @since LINE 3.0
 */
public class PnmlIO {

    private PnmlIO() {
    }

    /**
     * Write the Petri net held by a Network to a PNML place/transition file.
     *
     * @param model    network holding only places and transitions, with one closed class
     * @param filename output path
     * @throws IOException if the file cannot be written
     */
    public static void save(Network model, String filename) throws IOException {
        if (model == null) {
            throw new IllegalArgumentException("pnml_save expects a Network holding a Petri net.");
        }
        List<JobClass> classes = model.getClasses();
        if (classes.size() != 1) {
            throw new IllegalArgumentException("The PNML place/transition grammar is UNCOLOURED, so it cannot carry a net with "
                    + classes.size() + " job classes: its tokens are indistinguishable. Export a single-class net, or use "
                    + "LineModelIO for the full model.");
        }
        if (classes.get(0) instanceof OpenClass) {
            throw new IllegalArgumentException("The PNML place/transition grammar has no unbounded token source, so an open "
                    + "class cannot be represented. Close the class, or use LineModelIO.");
        }

        List<Place> places = new ArrayList<Place>();
        List<Transition> transitions = new ArrayList<Transition>();
        List<jline.lang.nodes.Node> nodes = model.getNodes();
        for (int i = 0; i < nodes.size(); i++) {
            jline.lang.nodes.Node nd = nodes.get(i);
            if (nd instanceof Place) {
                Place pl = (Place) nd;
                if (pl.isQueueing()) {
                    throw new IllegalArgumentException("Place " + pl.getName() + " is a QUEUEING place, i.e. a station with an "
                            + "embedded queue, which the place/transition grammar cannot represent.");
                }
                places.add(pl);
            } else if (nd instanceof Transition) {
                transitions.add((Transition) nd);
            } else {
                throw new IllegalArgumentException("Node " + nd.getName() + " is a " + nd.getClass().getSimpleName()
                        + ". A PNML place/transition net holds only places and transitions; a Source, a Sink or a queueing "
                        + "station has no counterpart in the grammar.");
            }
        }
        if (places.isEmpty()) {
            throw new IllegalArgumentException("The model holds no Place, so there is no Petri net to write.");
        }

        long[] marking = initialMarking(model, places, classes.get(0));

        String netName = model.getName();
        if (netName == null || netName.isEmpty()) {
            netName = "net";
        }

        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">\n");
        sb.append("  <net id=\"").append(escape(netName)).append("\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">\n");
        sb.append("    <name><text>").append(escape(netName)).append("</text></name>\n");
        sb.append("    <page id=\"page0\">\n");

        for (int p = 0; p < places.size(); p++) {
            sb.append("      <place id=\"").append(escape(places.get(p).getName())).append("\">\n");
            sb.append("        <name><text>").append(escape(places.get(p).getName())).append("</text></name>\n");
            sb.append("        <initialMarking><text>").append(marking[p]).append("</text></initialMarking>\n");
            sb.append("      </place>\n");
        }

        List<String> arcs = new ArrayList<String>();
        int arcId = 0;
        for (int t = 0; t < transitions.size(); t++) {
            Transition tr = transitions.get(t);
            List<Mode> modes = tr.getModes();
            if (modes.isEmpty()) {
                throw new IllegalArgumentException("Transition " + tr.getName()
                        + " declares no mode, so it has no firing behaviour to write.");
            }
            for (int m = 0; m < modes.size(); m++) {
                Mode mode = modes.get(m);
                if (tr.firingRateDependence.get(mode) != null) {
                    throw new IllegalArgumentException("Transition " + tr.getName() + " mode " + (m + 1)
                            + " declares a marking-dependent firing rate, which no PNML element can carry.");
                }
                String tid = modes.size() == 1 ? tr.getName() : tr.getName() + "." + mode.getName();
                TimingStrategy timing = tr.timingStrategies.get(mode);
                if (timing == null) {
                    timing = TimingStrategy.TIMED;
                }
                sb.append("      <transition id=\"").append(escape(tid)).append("\">\n");
                sb.append("        <name><text>").append(escape(tid)).append("</text></name>\n");
                sb.append("        <toolspecific tool=\"LINE\" version=\"3.0\">\n");
                sb.append("          <mode transition=\"").append(escape(tr.getName()))
                        .append("\" name=\"").append(escape(mode.getName()))
                        .append("\" timing=\"").append(timing == TimingStrategy.IMMEDIATE ? "immediate" : "timed")
                        .append("\" servers=\"").append(num(tr.getNumberOfModeServers(mode)))
                        .append("\" priority=\"").append(num(tr.firingPriorities.get(0, m)))
                        .append("\" weight=\"").append(num(tr.firingWeights.get(0, m)))
                        .append("\">\n");
                if (timing != TimingStrategy.IMMEDIATE) {
                    sb.append(distToXml(tr.distributions.get(mode), tr.getName(), m + 1));
                }
                sb.append("          </mode>\n");
                sb.append("        </toolspecific>\n");
                sb.append("      </transition>\n");

                Matrix enabling = tr.enablingConditions.get(mode);
                Matrix inhibiting = tr.inhibitingConditions.get(mode);
                Matrix firing = tr.firingOutcomes.get(mode);
                for (int p = 0; p < places.size(); p++) {
                    int pidx = model.getNodeIndex(places.get(p));
                    double w = enabling == null ? 0.0 : enabling.get(pidx, 0);
                    if (w > 0) {
                        arcId++;
                        arcs.add(arc(arcId, places.get(p).getName(), tid, (long) Math.round(w), false));
                    }
                    double inh = inhibiting == null ? Double.POSITIVE_INFINITY : inhibiting.get(pidx, 0);
                    if (!Double.isInfinite(inh) && !Double.isNaN(inh)) {
                        arcId++;
                        arcs.add(arc(arcId, places.get(p).getName(), tid, (long) Math.round(inh), true));
                    }
                    double f = firing == null ? 0.0 : firing.get(pidx, 0);
                    if (f > 0) {
                        arcId++;
                        arcs.add(arc(arcId, tid, places.get(p).getName(), (long) Math.round(f), false));
                    }
                }
            }
        }

        for (int a = 0; a < arcs.size(); a++) {
            sb.append(arcs.get(a));
        }
        sb.append("    </page>\n");
        sb.append("  </net>\n");
        sb.append("</pnml>\n");

        BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename),
                Charset.forName("UTF-8")));
        try {
            out.write(sb.toString());
        } finally {
            out.close();
        }
    }

    /**
     * Read the first net of a PNML place/transition document.
     *
     * @param filename input path
     * @return the equivalent LINE network
     */
    public static Network load(String filename) {
        return load(filename, null);
    }

    /**
     * Read one net of a PNML place/transition document.
     *
     * @param filename input path
     * @param netId    id of the net to read, or null for the first
     * @return the equivalent LINE network
     */
    public static Network load(String filename, String netId) {
        File f = new File(filename);
        if (!f.exists()) {
            throw new IllegalArgumentException("File " + filename + " cannot be found.");
        }
        Document doc;
        try {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            doc = db.parse(f);
        } catch (Exception e) {
            throw new RuntimeException("Cannot parse " + filename + ": " + e.getMessage(), e);
        }

        List<Element> nets = children(doc.getDocumentElement(), "net");
        if (nets.isEmpty()) {
            throw new IllegalArgumentException(filename + " holds no <net> element.");
        }
        Element net = null;
        for (int i = 0; i < nets.size(); i++) {
            if (netId == null || netId.isEmpty() || netId.equals(nets.get(i).getAttribute("id"))) {
                net = nets.get(i);
                break;
            }
        }
        if (net == null) {
            throw new IllegalArgumentException(filename + " holds no net with id \"" + netId + "\".");
        }
        String nettype = net.getAttribute("type");
        if (nettype != null && !nettype.isEmpty() && nettype.indexOf("ptnet") < 0) {
            throw new IllegalArgumentException("Net \"" + net.getAttribute("id") + "\" declares type " + nettype
                    + ". Only the place/transition grammar (http://www.pnml.org/version-2009/grammar/ptnet) is read: a "
                    + "coloured or a symmetric net carries token colours that a single-class LINE net cannot hold.");
        }

        // Places, transitions and arcs may sit directly under <net> or under any
        // <page>; the grammar allows both and tools differ, so the whole subtree is
        // searched rather than one level of it.
        List<Element> placeElems = descendants(net, "place");
        List<Element> transElems = descendants(net, "transition");
        List<Element> arcElems = descendants(net, "arc");
        if (placeElems.isEmpty()) {
            throw new IllegalArgumentException("Net \"" + net.getAttribute("id") + "\" holds no place.");
        }

        List<String> placeNames = new ArrayList<String>();
        List<Long> placeMarking = new ArrayList<Long>();
        long total = 0;
        for (int i = 0; i < placeElems.size(); i++) {
            placeNames.add(elementId(placeElems.get(i)));
            long mk = (long) textNumber(placeElems.get(i), "initialMarking", 0.0);
            placeMarking.add(Long.valueOf(mk));
            total += mk;
        }
        if (total == 0) {
            throw new IllegalArgumentException("The initial marking of this net is empty. A LINE closed class needs tokens "
                    + "to hold, and a net with none has no reachable behaviour to analyse.");
        }

        int nt = transElems.size();
        List<String> transIds = new ArrayList<String>();
        List<String> modeOwner = new ArrayList<String>();
        List<String> modeName = new ArrayList<String>();
        List<TimingStrategy> modeTiming = new ArrayList<TimingStrategy>();
        List<Integer> modeServers = new ArrayList<Integer>();
        List<Integer> modePriority = new ArrayList<Integer>();
        List<Double> modeWeight = new ArrayList<Double>();
        List<Distribution> modeDist = new ArrayList<Distribution>();
        for (int i = 0; i < nt; i++) {
            Element te = transElems.get(i);
            String id = elementId(te);
            transIds.add(id);
            String owner = id;
            String name = "Mode1";
            TimingStrategy timing = TimingStrategy.TIMED;
            int servers = 1;
            int priority = 1;
            double weight = 1.0;
            Distribution dist = new Exp(1);
            Element spec = lineToolspecific(te);
            if (spec != null) {
                String ownerAttr = spec.getAttribute("transition");
                if (ownerAttr != null && !ownerAttr.isEmpty()) {
                    owner = ownerAttr;
                }
                String nameAttr = spec.getAttribute("name");
                if (nameAttr != null && !nameAttr.isEmpty()) {
                    name = nameAttr;
                }
                if ("immediate".equalsIgnoreCase(spec.getAttribute("timing"))) {
                    timing = TimingStrategy.IMMEDIATE;
                }
                double sv = attrNumber(spec, "servers", 1.0);
                servers = Double.isInfinite(sv) ? Integer.MAX_VALUE : (int) Math.round(sv);
                priority = (int) Math.round(attrNumber(spec, "priority", 1.0));
                weight = attrNumber(spec, "weight", 1.0);
                if (timing != TimingStrategy.IMMEDIATE) {
                    dist = distFromXml(spec);
                }
            }
            modeOwner.add(owner);
            modeName.add(name);
            modeTiming.add(timing);
            modeServers.add(Integer.valueOf(servers));
            modePriority.add(Integer.valueOf(priority));
            modeWeight.add(Double.valueOf(weight));
            modeDist.add(dist);
        }

        List<String> ownerNames = new ArrayList<String>();
        int[] ownerOf = new int[nt];
        for (int i = 0; i < nt; i++) {
            int k = ownerNames.indexOf(modeOwner.get(i));
            if (k < 0) {
                ownerNames.add(modeOwner.get(i));
                k = ownerNames.size() - 1;
            }
            ownerOf[i] = k;
        }

        String name = netId;
        if (name == null || name.isEmpty()) {
            name = net.getAttribute("id");
        }
        if (name == null || name.isEmpty()) {
            name = "pnml";
        }

        Network model = new Network(name);
        List<Place> placeObj = new ArrayList<Place>();
        for (int i = 0; i < placeNames.size(); i++) {
            placeObj.add(new Place(model, placeNames.get(i)));
        }
        List<Transition> transObj = new ArrayList<Transition>();
        for (int i = 0; i < ownerNames.size(); i++) {
            transObj.add(new Transition(model, ownerNames.get(i)));
        }

        // The reference station is the first place holding tokens, so the class starts
        // where the marking says it does.
        int refIdx = 0;
        for (int i = 0; i < placeMarking.size(); i++) {
            if (placeMarking.get(i).longValue() > 0) {
                refIdx = i;
                break;
            }
        }
        ClosedClass jobclass = new ClosedClass(model, "Class1", (int) total, (Station) placeObj.get(refIdx), 0);

        List<Mode> modeObj = new ArrayList<Mode>();
        for (int i = 0; i < nt; i++) {
            Transition tr = transObj.get(ownerOf[i]);
            Mode mode = tr.addMode(modeName.get(i));
            modeObj.add(mode);
            tr.setTimingStrategy(mode, modeTiming.get(i));
            tr.setNumberOfServers(mode, modeServers.get(i));
            tr.setFiringPriorities(mode, modePriority.get(i).intValue());
            tr.setFiringWeights(mode, modeWeight.get(i).doubleValue());
            if (modeTiming.get(i) != TimingStrategy.IMMEDIATE) {
                tr.setDistribution(mode, modeDist.get(i));
            }
        }

        RoutingMatrix routing = model.initRoutingMatrix();
        for (int a = 0; a < arcElems.size(); a++) {
            Element ae = arcElems.get(a);
            String src = ae.getAttribute("source");
            String tgt = ae.getAttribute("target");
            double w = textNumber(ae, "inscription", 1.0);
            if (w <= 0) {
                throw new IllegalArgumentException("Arc " + src + " -> " + tgt + " carries a non-positive inscription " + w + ".");
            }
            int ip = placeNames.indexOf(src);
            int it = transIds.indexOf(tgt);
            if (ip >= 0 && it >= 0) {
                Transition tr = transObj.get(ownerOf[it]);
                if (isInhibitor(ae)) {
                    tr.setInhibitingConditions(modeObj.get(it), jobclass, placeObj.get(ip), (int) Math.round(w));
                } else {
                    tr.setEnablingConditions(modeObj.get(it), jobclass, placeObj.get(ip), (int) Math.round(w));
                }
                routing.set(jobclass, jobclass, placeObj.get(ip), tr, 1.0);
                continue;
            }
            it = transIds.indexOf(src);
            ip = placeNames.indexOf(tgt);
            if (it >= 0 && ip >= 0) {
                Transition tr = transObj.get(ownerOf[it]);
                tr.setFiringOutcome(modeObj.get(it), jobclass, placeObj.get(ip), (int) Math.round(w));
                routing.set(jobclass, jobclass, tr, placeObj.get(ip), 1.0);
                continue;
            }
            throw new IllegalArgumentException("Arc " + src + " -> " + tgt + " connects two places or two transitions, "
                    + "which the place/transition grammar does not allow.");
        }
        model.link(routing);

        for (int i = 0; i < placeObj.size(); i++) {
            placeObj.get(i).setMarking((int) placeMarking.get(i).longValue());
        }
        return model;
    }

    /**
     * Token count of each place in the initial marking. The state set on the place is
     * authoritative; a place with no state holds the class population when it is the
     * reference station of the class, which is the default LINE itself applies.
     */
    private static long[] initialMarking(Network model, List<Place> places, JobClass jobclass) {
        long[] marking = new long[places.size()];
        Station refstat = null;
        if (jobclass instanceof ClosedClass) {
            refstat = ((ClosedClass) jobclass).getReferenceStation();
        }
        for (int p = 0; p < places.size(); p++) {
            Matrix st = places.get(p).getState();
            if (st != null && !st.isEmpty()) {
                marking[p] = Math.round(st.get(0, 0));
            } else if (refstat != null && refstat.getName().equals(places.get(p).getName())) {
                marking[p] = Math.round(((ClosedClass) jobclass).getPopulation());
            }
        }
        return marking;
    }

    private static String arc(int id, String source, String target, long weight, boolean inhibitor) {
        StringBuilder sb = new StringBuilder();
        sb.append("      <arc id=\"a").append(id).append("\" source=\"").append(escape(source))
                .append("\" target=\"").append(escape(target)).append("\">\n");
        if (inhibitor) {
            // An inhibitor arc is not in the P/T grammar itself; <type value="inhibitor"/>
            // is the extension GreatSPN, TINA and PIPE all read, so it is the one written.
            sb.append("        <type value=\"inhibitor\"/>\n");
        }
        sb.append("        <inscription><text>").append(weight).append("</text></inscription>\n");
        sb.append("      </arc>\n");
        return sb.toString();
    }

    /**
     * Timing of a timed mode, as a distribution name and its scalar parameters. The
     * parameter NAMES are LINE's own, so the block is self-describing and the reader
     * reconstructs the object by name rather than by position. A distribution whose
     * parameters are not scalars -- a phase-type or a MAP given by matrices -- is REFUSED
     * here: writing its mean rate instead would produce a file that reads back as a
     * different model.
     */
    private static String distToXml(Distribution dist, String trname, int m) {
        if (dist == null) {
            throw new IllegalArgumentException("Transition " + trname + " mode " + m
                    + " is timed but carries no distribution.");
        }
        if (dist instanceof Immediate) {
            return "            <distribution name=\"Immediate\"/>\n";
        }
        if (dist instanceof Disabled) {
            return "            <distribution name=\"Disabled\"/>\n";
        }
        int nparams = dist.getNumParams(0);
        if (nparams <= 0) {
            throw new IllegalArgumentException("Transition " + trname + " mode " + m + " holds a "
                    + dist.getClass().getSimpleName() + ", which declares no scalar parameter and therefore cannot be "
                    + "written to PNML.");
        }
        StringBuilder sb = new StringBuilder();
        sb.append("            <distribution name=\"").append(escape(dist.getClass().getSimpleName())).append("\">\n");
        for (int i = 1; i <= nparams; i++) {
            NamedParam np = dist.getParam(i);
            Object v = np == null ? null : np.getValue();
            if (!(v instanceof Number)) {
                throw new IllegalArgumentException("Transition " + trname + " mode " + m + " holds a "
                        + dist.getClass().getSimpleName() + " whose parameter \"" + (np == null ? "?" : np.getName())
                        + "\" is not a scalar. The PNML timing block carries scalar parameters only; a "
                        + "matrix-parameterized law (PH, APH, MAP, MMPP2, ME, RAP) has no PNML representation and "
                        + "writing its mean rate instead would read back as a different model.");
            }
            sb.append("              <parameter name=\"").append(escape(np.getName()))
                    .append("\" value=\"").append(num(((Number) v).doubleValue())).append("\"/>\n");
        }
        sb.append("            </distribution>\n");
        return sb.toString();
    }

    /**
     * Rebuild the distribution of a timed mode from the toolspecific block. The
     * constructors are named ONE BY ONE rather than applied positionally: LINE stores
     * Weibull as (alpha=scale, r=shape) while its constructor takes (shape, scale), so a
     * positional rebuild would silently transpose the two.
     */
    private static Distribution distFromXml(Element modeElem) {
        List<Element> dists = children(modeElem, "distribution");
        if (dists.isEmpty()) {
            return new Exp(1);
        }
        Element d = dists.get(0);
        String name = d.getAttribute("name");
        Map<String, Double> p = new LinkedHashMap<String, Double>();
        List<Element> params = children(d, "parameter");
        for (int i = 0; i < params.size(); i++) {
            p.put(params.get(i).getAttribute("name"), Double.valueOf(attrNumber(params.get(i), "value", Double.NaN)));
        }
        if ("Immediate".equals(name)) {
            return new Immediate();
        }
        if ("Disabled".equals(name)) {
            return new Disabled();
        }
        if ("Exp".equals(name)) {
            return new Exp(need(p, "lambda", name));
        }
        if ("Det".equals(name)) {
            return new Det(need(p, "t", name));
        }
        if ("Erlang".equals(name)) {
            return new Erlang(need(p, "alpha", name), (int) Math.round(need(p, "r", name)));
        }
        if ("HyperExp".equals(name)) {
            return new HyperExp(need(p, "p", name), need(p, "lambda1", name), need(p, "lambda2", name));
        }
        if ("Uniform".equals(name)) {
            return new Uniform(need(p, "min", name), need(p, "max", name));
        }
        if ("Gamma".equals(name)) {
            return new Gamma(need(p, "alpha", name), need(p, "beta", name));
        }
        if ("Pareto".equals(name)) {
            return new Pareto(need(p, "alpha", name), need(p, "k", name));
        }
        if ("Weibull".equals(name)) {
            return new Weibull(need(p, "r", name), need(p, "alpha", name));
        }
        if ("Lognormal".equals(name)) {
            return new Lognormal(need(p, "mu", name), need(p, "sigma", name));
        }
        throw new IllegalArgumentException("The timing block names distribution \"" + name + "\", which the PNML reader "
                + "does not construct. The families it reads are Exp, Det, Erlang, HyperExp, Uniform, Gamma, Pareto, "
                + "Weibull, Lognormal, Immediate and Disabled, which are the ones the writer writes.");
    }

    private static double need(Map<String, Double> p, String key, String distName) {
        Double v = p.get(key);
        if (v == null || v.isNaN()) {
            throw new IllegalArgumentException("Distribution " + distName + " is missing the parameter \"" + key + "\".");
        }
        return v.doubleValue();
    }

    private static Element lineToolspecific(Element transElem) {
        List<Element> blocks = children(transElem, "toolspecific");
        for (int i = 0; i < blocks.size(); i++) {
            if (!"LINE".equalsIgnoreCase(blocks.get(i).getAttribute("tool"))) {
                continue;
            }
            List<Element> modes = children(blocks.get(i), "mode");
            if (!modes.isEmpty()) {
                return modes.get(0);
            }
        }
        return null;
    }

    private static boolean isInhibitor(Element arcElem) {
        List<Element> types = children(arcElem, "type");
        for (int i = 0; i < types.size(); i++) {
            if ("inhibitor".equalsIgnoreCase(types.get(i).getAttribute("value"))) {
                return true;
            }
        }
        return "inhibitor".equalsIgnoreCase(arcElem.getAttribute("type"));
    }

    private static String elementId(Element elem) {
        String id = elem.getAttribute("id");
        if (id == null || id.isEmpty()) {
            List<Element> names = children(elem, "name");
            if (!names.isEmpty()) {
                id = textOf(names.get(0));
            }
        }
        if (id == null || id.isEmpty()) {
            throw new IllegalArgumentException("A place or transition carries neither an id nor a name.");
        }
        return id;
    }

    private static String textOf(Element elem) {
        List<Element> texts = children(elem, "text");
        if (!texts.isEmpty()) {
            return texts.get(0).getTextContent().trim();
        }
        return elem.getTextContent().trim();
    }

    private static double textNumber(Element elem, String labelName, double defaultValue) {
        List<Element> labels = children(elem, labelName);
        if (labels.isEmpty()) {
            return defaultValue;
        }
        String txt = textOf(labels.get(0));
        if (txt.isEmpty()) {
            return defaultValue;
        }
        // Some tools write the marking as "3" and some as "1`3" (a coloured multiset of
        // one colour); the plain integer is the one the grammar defines.
        int tick = txt.lastIndexOf('`');
        if (tick >= 0) {
            txt = txt.substring(tick + 1);
        }
        try {
            return Double.parseDouble(txt.trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Label <" + labelName + "> holds \"" + txt + "\", which is not a number.");
        }
    }

    private static double attrNumber(Element elem, String key, double defaultValue) {
        String txt = elem.getAttribute(key);
        if (txt == null || txt.trim().isEmpty()) {
            return defaultValue;
        }
        txt = txt.trim();
        if ("Inf".equalsIgnoreCase(txt)) {
            return Double.POSITIVE_INFINITY;
        }
        if ("-Inf".equalsIgnoreCase(txt)) {
            return Double.NEGATIVE_INFINITY;
        }
        try {
            return Double.parseDouble(txt);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Attribute " + key + " holds \"" + txt + "\", which is not a number.");
        }
    }

    /** Direct children with the given local name, ignoring any namespace prefix. */
    private static List<Element> children(Element elem, String localName) {
        List<Element> out = new ArrayList<Element>();
        NodeList kids = elem.getChildNodes();
        for (int i = 0; i < kids.getLength(); i++) {
            Node nd = kids.item(i);
            if (nd.getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            if (localName.equals(local(nd.getNodeName()))) {
                out.add((Element) nd);
            }
        }
        return out;
    }

    /** Descendants with the given local name, in document order. */
    private static List<Element> descendants(Element elem, String localName) {
        List<Element> out = new ArrayList<Element>();
        NodeList kids = elem.getChildNodes();
        for (int i = 0; i < kids.getLength(); i++) {
            Node nd = kids.item(i);
            if (nd.getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            if (localName.equals(local(nd.getNodeName()))) {
                out.add((Element) nd);
            } else {
                out.addAll(descendants((Element) nd, localName));
            }
        }
        return out;
    }

    private static String local(String name) {
        int k = name.lastIndexOf(':');
        return k < 0 ? name : name.substring(k + 1);
    }

    /**
     * Numbers are written in the shortest form that reads back exactly, so an integral
     * count does not acquire a decimal point and an infinite server count keeps the
     * spelling the reader expects.
     */
    private static String num(double v) {
        if (Double.isInfinite(v)) {
            return v > 0 ? "Inf" : "-Inf";
        }
        if (v == Math.rint(v) && Math.abs(v) < 9.007199254740992E15) {
            return Long.toString((long) v);
        }
        return Double.toString(v);
    }

    private static String num(int v) {
        if (v == Integer.MAX_VALUE) {
            return "Inf";
        }
        return Integer.toString(v);
    }

    private static String escape(String s) {
        if (s == null) {
            return "";
        }
        String out = s.replace("&", "&amp;");
        out = out.replace("<", "&lt;");
        out = out.replace(">", "&gt;");
        out = out.replace("\"", "&quot;");
        return out;
    }
}
