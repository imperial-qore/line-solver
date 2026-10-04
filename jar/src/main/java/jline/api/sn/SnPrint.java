package jline.api.sn;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import jline.GlobalConstants;
import jline.lang.Event;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.JobClass;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SnPrint {
    private SnPrint() {}

    public static void snPrint(NetworkStruct sn) {
        System.out.println("nstations: " + sn.nstations);
        System.out.println("nstateful: " + sn.nstateful);
        System.out.println("nnodes: " + sn.nnodes);
        System.out.println("nclasses: " + sn.nclasses);
        System.out.println("nclosedjobs: " + sn.nclosedjobs);
        System.out.println("nchains: " + sn.nchains);

        printMatrix("refstat", sn.refstat);
        printMatrix("njobs", sn.njobs);
        printMatrix("nservers", sn.nservers);
        printMatrix("connmatrix", sn.connmatrix);
        printMatrix("scv", sn.scv);
        printMatrix("isstation", sn.isstation);
        printMatrix("isstateful", sn.isstateful);
        printMatrix("isstatedep", sn.isstatedep);
        printMatrix("nodeToStateful", sn.nodeToStateful);
        printMatrix("nodeToStation", sn.nodeToStation);
        printMatrix("stationToNode", sn.stationToNode);
        printMatrix("stationToStateful", sn.stationToStateful);
        printMatrix("statefulToStation", sn.statefulToStation);
        printMatrix("statefulToNode", sn.statefulToNode);
        printMatrix("rates", sn.rates);
        printMatrix("classprio", sn.classprio);
        printMatrix("phases", sn.phases);
        printMatrix("phasessz", sn.phasessz);
        printMatrix("phaseshift", sn.phaseshift);
        printMatrix("schedparam", sn.schedparam);
        printMatrix("chains", sn.chains);
        printMatrix("rt", sn.rt);
        printMatrix("nvars", sn.nvars);
        printMatrix("rtnodes", sn.rtnodes);
        printMatrix("csmask", sn.csmask);
        printMatrix("isslc", sn.isslc);
        printMatrix("cap", sn.cap);
        printMatrix("classcap", sn.classcap);
        printMatrix("refclass", sn.refclass);
        printMatrix("lldscaling", sn.lldscaling);
        printMatrix("fj", sn.fj);

        printList("nodetype", sn.nodetype);
        printList("classnames", sn.classnames);
        printList("nodenames", sn.nodenames);

        printMapContents("rtorig", sn.rtorig);
        printMapContents("lst", sn.lst);
        printMapContents("state", sn.state);
        printMapContents("stateprior", sn.stateprior);
        printMapContents("space", sn.space);
        printMapContents("routing", sn.routing);
        printMapContents("procid", sn.procid);
        printMapContents("mu", sn.mu);
        printMapContents("phi", sn.phi);
        printMapContents("proc", sn.proc);
        printMapContents("pie", sn.pie);
        printMapContents("sched", sn.sched);
        printMapContents("inchain", sn.inchain);
        printMapContents("visits", sn.visits);
        printMapContents("nodevisits", sn.nodevisits);
        printMapContents("droprule", sn.droprule);
        printMapContents("nodeparam", sn.nodeparam);
        printMapContents("sync", sn.sync);
        printMapContents("gsync", sn.gsync);
        printMapContents("cdscaling", sn.cdscaling);
        printMapContents("jdscaling", sn.jdscaling);

        if (sn.stations != null) {
            List<String> stationNames = new ArrayList<String>();
            for (Station s : sn.stations) stationNames.add(s.getName());
            printList("stations", stationNames);
        } else {
            printList("stations", null);
        }
        if (sn.stateful != null) {
            List<String> statefulNames = new ArrayList<String>();
            for (StatefulNode s : sn.stateful) statefulNames.add(s.getName());
            printList("stateful", statefulNames);
        } else {
            printList("stateful", null);
        }
        if (sn.jobclasses != null) {
            List<String> jobClassNames = new ArrayList<String>();
            for (JobClass jc : sn.jobclasses) jobClassNames.add(jc.getName());
            printList("jobclasses", jobClassNames);
        } else {
            printList("jobclasses", null);
        }
        if (sn.nodes != null) {
            List<String> nodeNames = new ArrayList<String>();
            for (Node n : sn.nodes) nodeNames.add(n.getName());
            printList("nodes", nodeNames);
        } else {
            printList("nodes", null);
        }
    }

    private static String printMatrixCompact(Matrix matrix) {
        if (matrix == null) return "null";
        if (matrix.isEmpty()) return "[]";
        StringBuilder sb = new StringBuilder("[");
        int numRows = matrix.getNumRows();
        int numCols = matrix.getNumCols();
        for (int i = 0; i < numRows; i++) {
            if (i > 0) sb.append("; ");
            for (int j = 0; j < numCols; j++) {
                if (j > 0) sb.append(" ");
                double value = matrix.get(i, j);
                if (Double.isNaN(value)) sb.append("NaN");
                else if (value == (double) GlobalConstants.MaxInt) sb.append("Inf");
                else if (Double.isInfinite(value)) sb.append("Inf");
                else if (matrix.isInteger()) sb.append(Math.round(value));
                else if (value == Math.floor(value)) sb.append((int) value);
                else sb.append(value);
            }
        }
        sb.append("]");
        return sb.toString();
    }

    private static void printList(String name, List<?> list) {
        System.out.print(name + ": ");
        if (list == null || list.isEmpty()) { System.out.println("[]"); return; }
        System.out.print("[");
        for (int i = 0; i < list.size(); i++) {
            if (i > 0) System.out.print(", ");
            Object item = list.get(i);
            if (item instanceof String) System.out.print("\"" + item + "\"");
            else System.out.print(item);
        }
        System.out.println("]");
    }

    private static void printMatrix(String name, Matrix matrix) {
        if (matrix == null) System.out.println(name + ": null");
        else if (matrix.isEmpty()) System.out.println(name + ": []");
        else { System.out.print(name + ": "); System.out.print(printMatrixCompact(matrix)); System.out.println(); }
    }

    private static String keyName(Object key) {
        if (key == null) return "null";
        if (key instanceof String) return (String) key;
        if (key instanceof Integer) return key.toString();
        try {
            Method m = key.getClass().getMethod("getName");
            return m.invoke(key).toString();
        } catch (Exception e) {
            return key.getClass().getSimpleName();
        }
    }

    private static void printMapContents(String name, Map<?, ?> map) {
        if (map == null) { System.out.println(name + ": null"); return; }
        if (map.isEmpty()) { System.out.println(name + ": {}"); return; }
        System.out.print(name + ": {");
        boolean first = true;
        TreeMap<String, Map.Entry<?, ?>> sorted = new TreeMap<String, Map.Entry<?, ?>>();
        for (Map.Entry<?, ?> e : map.entrySet()) sorted.put(keyName(e.getKey()), e);
        for (Map.Entry<String, Map.Entry<?, ?>> entry : sorted.entrySet()) {
            if (!first) System.out.print(", ");
            first = false;
            Object key = entry.getValue().getKey();
            if (key == null) System.out.print("null");
            else if (key instanceof String) System.out.print("\"" + key + "\"");
            else if (key instanceof Integer) System.out.print(key);
            else {
                try {
                    Method m = key.getClass().getMethod("getName");
                    System.out.print("\"" + m.invoke(key) + "\"");
                } catch (Exception ex) {
                    System.out.print("\"" + key.getClass().getSimpleName() + "\"");
                }
            }
            System.out.print(": ");
            printValue(entry.getValue().getValue());
        }
        System.out.println("}");
    }

    private static void printValue(Object value) {
        if (value == null) { System.out.print("null"); return; }
        if (value instanceof Map<?, ?>) {
            Map<?, ?> m = (Map<?, ?>) value;
            if (m.isEmpty()) { System.out.print("{}"); return; }
            System.out.print("{");
            boolean innerFirst = true;
            TreeMap<String, Map.Entry<?, ?>> sorted = new TreeMap<String, Map.Entry<?, ?>>();
            for (Map.Entry<?, ?> e : m.entrySet()) sorted.put(keyName(e.getKey()), e);
            for (Map.Entry<String, Map.Entry<?, ?>> entry : sorted.entrySet()) {
                if (!innerFirst) System.out.print(", ");
                innerFirst = false;
                Object k = entry.getValue().getKey();
                if (k == null) System.out.print("null");
                else if (k instanceof String) System.out.print("\"" + k + "\"");
                else if (k instanceof Integer || k instanceof Long) System.out.print(k);
                else {
                    try {
                        Method mm = k.getClass().getMethod("getName");
                        System.out.print("\"" + mm.invoke(k) + "\"");
                    } catch (Exception ex) {
                        if (k instanceof Number) System.out.print(k);
                        else System.out.print(k.getClass().getSimpleName());
                    }
                }
                System.out.print(": ");
                printValue(entry.getValue().getValue());
            }
            System.out.print("}");
        } else if (value instanceof Matrix) {
            System.out.print(printMatrixCompact((Matrix) value));
        } else if (value instanceof String) {
            System.out.print("\"" + value + "\"");
        } else if (value instanceof Boolean) {
            System.out.print(value);
        } else if (value instanceof Integer) {
            int iv = (Integer) value;
            if (iv == GlobalConstants.MaxInt) System.out.print("Inf");
            else System.out.print(iv);
        } else if (value instanceof Double) {
            double dv = (Double) value;
            if (Double.isNaN(dv)) System.out.print("NaN");
            else if (dv == (double) GlobalConstants.MaxInt) System.out.print("Inf");
            else if (Double.isInfinite(dv)) System.out.print("Inf");
            else System.out.print(dv);
        } else if (value instanceof Float) {
            float fv = (Float) value;
            if (Float.isNaN(fv)) System.out.print("NaN");
            else if (fv == (float) GlobalConstants.MaxInt) System.out.print("Inf");
            else if (Float.isInfinite(fv)) System.out.print("Inf");
            else System.out.print(fv);
        } else if (value instanceof Long) {
            long lv = (Long) value;
            if (lv == (long) GlobalConstants.MaxInt) System.out.print("Inf");
            else System.out.print(lv);
        } else if (value instanceof Enum<?>) {
            System.out.print(value.toString());
        } else if (value instanceof List<?>) {
            List<?> list = (List<?>) value;
            System.out.print("[");
            for (int i = 0; i < list.size(); i++) {
                if (i > 0) System.out.print(", ");
                Object item = list.get(i);
                if (item instanceof Matrix) System.out.print(printMatrixCompact((Matrix) item));
                else printValue(item);
            }
            System.out.print("]");
        } else if (value instanceof MatrixCell) {
            MatrixCell mc = (MatrixCell) value;
            System.out.print("{");
            for (int i = 0; i < mc.size(); i++) {
                if (i > 0) System.out.print(", ");
                System.out.print("[" + i + "]: ");
                System.out.print(printMatrixCompact(mc.get(i)));
            }
            System.out.print("}");
        } else if (value instanceof Event) {
            Event ev = (Event) value;
            System.out.print("\"(" + ev.getEvent().name() + ": node: " + ev.getNode() + ", class: " + ev.getJobClass());
            if (!Double.isNaN(ev.getProb()) && ev.getProb() != 1.0) System.out.print(", prob: " + ev.getProb());
            if (!Double.isNaN(ev.getT())) System.out.print(", t: " + ev.getT());
            if (!Double.isNaN(ev.getJob())) System.out.print(", job: " + ev.getJob());
            System.out.print(")\"");
        } else if (value instanceof Sync) {
            Sync sync = (Sync) value;
            System.out.print("{");
            System.out.print("\"active\": ");
            printValue(sync.active);
            System.out.print(", \"passive\": ");
            printValue(sync.passive);
            System.out.print("}");
        } else {
            String str = value.toString();
            if (str.contains("Lambda") && str.contains("$")) {
                System.out.print("<Function>");
            } else if (str.contains("@") && str.contains(".")) {
                try {
                    Method m = value.getClass().getMethod("getName");
                    System.out.print(m.invoke(value));
                } catch (Exception e) {
                    System.out.print(value.getClass().getSimpleName());
                }
            } else {
                System.out.print(str);
            }
        }
    }

    /** Stochastic network Print algorithms. */
    public static final class SnprintAlgo {}
}
