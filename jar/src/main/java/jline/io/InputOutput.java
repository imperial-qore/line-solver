/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import java.io.File;
import java.lang.reflect.Method;
import java.util.logging.ConsoleHandler;
import java.util.logging.Formatter;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.apache.commons.io.FileUtils;
import org.w3c.dom.Document;

import jline.GlobalConstants;
import jline.VerboseLevel;

public final class InputOutput {
    private InputOutput() {}

    /**
     * Custom formatter that removes timestamps and formats messages on single lines
     */
    private static final class SimpleLineFormatter extends Formatter {
        @Override
        public String format(LogRecord record) {
            String level;
            if (record.getLevel() == Level.SEVERE) {
                level = "SEVERE";
            } else if (record.getLevel() == Level.WARNING) {
                level = "WARNING";
            } else if (record.getLevel() == Level.INFO) {
                level = "INFO";
            } else {
                level = record.getLevel().getName();
            }
            return level + ": " + record.getMessage() + "\n";
        }
    }

    private static String lastWarning = "";
    private static boolean suppressedWarnings = false;
    private static long suppressedWarningTic = System.currentTimeMillis();
    private static boolean loggerConfigured = false;

    /**
     * Check if running under Maven test
     */
    private static boolean isRunningUnderMavenTest() {
        if (System.getProperty("surefire.test.class.path") != null) {
            return true;
        }
        if (System.getProperty("maven.test.skip") != null) {
            return true;
        }
        if (System.getProperty("basedir") != null && System.getProperty("surefire.real.class.path") != null) {
            return true;
        }
        StackTraceElement[] stack = Thread.currentThread().getStackTrace();
        for (StackTraceElement el : stack) {
            String cn = el.getClassName();
            if (cn.contains("org.junit") || cn.contains("org.testng") || cn.contains("maven.surefire")) {
                return true;
            }
        }
        return false;
    }

    /**
     * Configure logger to use custom formatting without timestamps globally
     */
    private static void configureLogger() {
        Level logLevel;
        if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
            logLevel = Level.OFF;
        } else {
            logLevel = Level.WARNING;
        }

        if (!loggerConfigured) {
            LogManager logManager = LogManager.getLogManager();

            Logger rootLogger = Logger.getLogger("");

            Handler[] handlers = rootLogger.getHandlers();
            for (Handler handler : handlers) {
                rootLogger.removeHandler(handler);
            }

            ConsoleHandler rootConsoleHandler = new ConsoleHandler();
            rootConsoleHandler.setFormatter(new SimpleLineFormatter());
            rootConsoleHandler.setLevel(logLevel);
            rootLogger.addHandler(rootConsoleHandler);
            rootLogger.setLevel(logLevel);

            Logger logger = Logger.getLogger(FileUtils.class.getName());
            logger.setUseParentHandlers(true);
            logger.setLevel(logLevel);

            loggerConfigured = true;
        } else {
            Logger rootLogger = Logger.getLogger("");
            rootLogger.setLevel(logLevel);
            Handler[] handlers = rootLogger.getHandlers();
            for (Handler handler : handlers) {
                handler.setLevel(logLevel);
            }
            Logger.getLogger(FileUtils.class.getName()).setLevel(logLevel);
        }
    }

    private static void line_printf(String message) {
        configureLogger();
        Logger.getLogger(FileUtils.class.getName()).log(Level.INFO, message);
    }

    /**
     * Print debug message if verbose level is DEBUG
     */
    public static void line_debug(VerboseLevel verbose, String message) {
        if (verbose == VerboseLevel.DEBUG) {
            line_printf("[DEBUG] " + message);
        }
    }

    /**
     * Emit a warning without the repeat suppression applied by
     * {@link #line_warning}.
     *
     * <p>{@code line_warning} keeps only the last message and hides an
     * identical repeat for 60 seconds. That is right for configuration notices
     * cast once per model, but wrong for a warning that reports a correctness
     * limitation of the analysis: solving several models in one session would
     * then flag only the first one, and the user would read the silence on the
     * others as a clean bill of health. Warnings that say "these numbers are
     * not exact" must be raised for every model they apply to, so they go
     * through here instead. Verbosity gating is unchanged: SILENT still
     * silences everything.</p>
     */
    public static void line_warning_always(String caller, String msg, Object... args) {
        if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
            return;
        }
        line_printf(String.format("[%s] %s", caller, String.format(msg, args)), Level.WARNING);
    }

    public static void line_warning(String caller, String msg, Object... args) {
        Logger logger = Logger.getLogger(FileUtils.class.getName());

        if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
            return;
        }

        String errmsg = String.format(msg, args);
        String finalmsg = String.format("[%s] %s", caller, errmsg);

        try {
            long currentTime = System.currentTimeMillis();

            if (finalmsg.compareTo(lastWarning) != 0 || (currentTime - suppressedWarningTic) > 60000) {
                line_printf(finalmsg, Level.WARNING);
                lastWarning = finalmsg;
                suppressedWarnings = false;
                suppressedWarningTic = currentTime;
            } else {
                if (!suppressedWarnings) {
                    line_printf(String.format("[%s] %s",
                                    caller,
                                    "Message casted more than once, repetitions will not be printed on screen for 60 seconds."),
                            Level.WARNING);
                    suppressedWarnings = true;
                    suppressedWarningTic = currentTime;
                }
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Exception in line_warning", e);
        }
    }

    private static void line_printf(String message, Level level) {
        configureLogger();
        Logger.getLogger(FileUtils.class.getName()).log(level, message);
    }

    public static void line_error(String caller, String msg) {
        String finalmsg = String.format("[%s] %s", caller, msg);
        if (GlobalConstants.Verbose != VerboseLevel.SILENT) {
            line_printf(finalmsg, Level.SEVERE);
        }
        throw new RuntimeException(finalmsg);
    }

    /**
     * Registry of the external-tool acknowledgements already printed in this
     * session (one JVM = one session). Distinct from the collective
     * library-attribution flag of the native solvers: each external tool a
     * wrapper solver delegates to must be acknowledged on its own.
     */
    private static final java.util.Set<String> toolAckShown =
            java.util.Collections.synchronizedSet(new java.util.HashSet<String>());

    /**
     * Print, once per session, the acknowledgement of the external tool that a
     * wrapper solver delegates to, together with the pointer to its official
     * website and the canonical paper to cite. The acknowledgement is
     * pull-based, like the library attribution: nothing is printed at the
     * default verbosity, and only a caller that asks for it explicitly, by
     * running at {@link VerboseLevel#DEBUG}, gets the line, at most once per
     * tool per session. {@code solver.citations()} and
     * {@link #line_citation(String)} are the quiet ways to obtain the same
     * reference.
     *
     * <p>The strings were verified against each upstream project's own pages;
     * do not reword an author list or a URL from memory. Mirror any edit in
     * the MATLAB ({@code matlab/src/io/line_ack.m}) and Python
     * ({@code line_solver/api/io/logging.py}) tables. The table covers the
     * in-tree wrapper solvers only; out-of-tree solvers carry their own
     * text.</p>
     *
     * <p>The machine-readable form of the same reference is
     * {@link #line_citation(String)}, which returns the BibTeX entry.</p>
     *
     * @param verbose  the caller's verbosity level
     * @param toolName the external tool, e.g. "JMT", "LQNS", "QNS"
     */
    /**
     * Prints the LINE startup banner, as MATLAB's lineStart and python's
     * line_solver.lineStart do.
     *
     * <p>Using the JAR as a library stays silent otherwise, so this is the
     * explicit entry point for a session banner. It also names where to obtain
     * the third-party dependencies and the algorithm references: attribution in
     * LINE is pull-based, in the spirit of Sage, so nothing is printed during a
     * solve.
     *
     * @return the LINE version string
     */
    public static String lineStart() {
        String version = new jline.lang.Model("").getVersion();
        System.out.printf("Starting LINE version %s: StdOut=console, VerboseLevel=%s, "
                + "CoarseTol=%.1e, FineTol=%.1e, Zero=%.1e, MaxInt=%d%n",
                version, GlobalConstants.Verbose, GlobalConstants.CoarseTol,
                GlobalConstants.FineTol, GlobalConstants.Zero, GlobalConstants.MaxInt);
        System.out.println("Type solver.libraries() for third-party dependencies, "
                + "solver.citations() for references.");
        return version;
    }

    public static void line_ack(VerboseLevel verbose, String toolName) {
        if (verbose != VerboseLevel.DEBUG || GlobalConstants.Verbose == VerboseLevel.SILENT) {
            return;
        }
        String key = toolName.toUpperCase();
        String msg = ackText(key);
        if (msg == null) {
            return;
        }
        if (!toolAckShown.add(key)) {
            return;
        }
        // see _kb/12-interfaces-and-docs.md (line_ack prints on System.out, not through the logger)
        System.out.println(msg);
        String cite = ackCitation(key);
        if (cite != null) {
            System.out.println("  Cite: " + cite);
        }
    }

    /**
     * Return the BibTeX entry for the canonical paper of an external tool that
     * a wrapper solver delegates to, so that the acknowledgement printed by
     * {@link #line_ack(VerboseLevel, String)} can be turned into a citation
     * without retyping it.
     *
     * <p>The keys match {@code doc/latex/biblio.bib} and
     * {@code BIBLIOGRAPHY.md}. Mirror any edit in
     * {@code matlab/src/io/line_citation.m} and
     * {@code line_solver.api.io.logging.line_citation}.</p>
     *
     * @param toolName the external tool, "JMT", "LQNS" or "QNS" ("QNS" shares
     *                 the LQNS reference, qnsolver being part of that
     *                 distribution)
     * @return the BibTeX entry, or the empty string for an unknown tool
     */
    public static String line_citation(String toolName) {
        String key = toolName.toUpperCase();
        if ("JMT".equals(key)) {
            return "@INPROCEEDINGS{BerCS07,\n"
                    + "  author = {M. Bertoli and G. Casale and G. Serazzi},\n"
                    + "  title = {The {JMT} Simulator for Performance Evaluation of Non-Product-Form\n"
                    + "\tQueueing Networks},\n"
                    + "  booktitle = {Proc. of the 40th Annual Simulation Symposium (ANSS)},\n"
                    + "  year = {2007},\n"
                    + "  pages = {3--10}\n"
                    + "}";
        } else if ("LQNS".equals(key) || "QNS".equals(key)) {
            return "@ARTICLE{fran.ea09,\n"
                    + "  author = {G. Franks and T. Al-Omari and M. Woodside and O. Das and S. Derisavi},\n"
                    + "  title = {Enhanced Modeling and Solution of Layered Queueing Networks},\n"
                    + "  journal = {IEEE Trans. Software Engineering},\n"
                    + "  year = {2009},\n"
                    + "  volume = {35},\n"
                    + "  pages = {148-161},\n"
                    + "  number = {2}\n"
                    + "}";
        }
        return "";
    }

    /**
     * One-line reference to the canonical paper of each tool, printed under the
     * acknowledgement. It describes the same work as
     * {@link #line_citation(String)} (keys {@code BerCS07} and
     * {@code fran.ea09} in {@code doc/latex/biblio.bib}), which is where the
     * citation key belongs: this line is for the reader.
     */
    private static String ackCitation(String toolName) {
        if ("JMT".equals(toolName)) {
            return "M. Bertoli, G. Casale, G. Serazzi. \"The JMT Simulator for Performance "
                    + "Evaluation of Non-Product-Form Queueing Networks\". Proc. of the 40th "
                    + "Annual Simulation Symposium (ANSS), pp. 3-10, 2007.";
        } else if ("LQNS".equals(toolName) || "QNS".equals(toolName)) {
            return "G. Franks, T. Al-Omari, M. Woodside, O. Das, S. Derisavi. \"Enhanced "
                    + "Modeling and Solution of Layered Queueing Networks\". IEEE Trans. "
                    + "Software Engineering, 35(2):148-161, 2009.";
        }
        return null;
    }

    private static String ackText(String toolName) {
        if ("JMT".equals(toolName)) {
            return "SolverJMT delegates to Java Modelling Tools (JMT), by M. Bertoli, "
                    + "G. Casale, G. Serazzi (Politecnico di Milano, Imperial College London). "
                    + "Please acknowledge the JMT authors: http://jmt.sourceforge.net/";
        } else if ("LQNS".equals(toolName)) {
            return "SolverLQNS delegates to LQNS/LQSIM, by G. Franks, M. Woodside et al. "
                    + "(Real-Time and Distributed Systems Group, Carleton University). "
                    + "Please acknowledge the LQNS authors: http://www.layeredqueues.org/";
        } else if ("QNS".equals(toolName)) {
            return "SolverQNS delegates to qnsolver, part of the LQNS distribution by "
                    + "G. Franks, M. Woodside et al. (Real-Time and Distributed Systems Group, "
                    + "Carleton University). Please acknowledge the LQNS authors: "
                    + "http://www.layeredqueues.org/";
        }
        return null;
    }

    public static String mfilename(Object obj) {
        Method enclosingMethod = obj.getClass().getEnclosingMethod();
        if (enclosingMethod != null) {
            return enclosingMethod.getName();
        } else {
            return obj.getClass().getSimpleName();
        }
    }

    /**
     * Writes the given XML Document object to a specified output file.
     */
    public static void writeXML(String outputFileName, Document doc) throws TransformerException {
        Transformer transformer = TransformerFactory.newInstance().newTransformer();
        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
        transformer.setOutputProperty(OutputKeys.VERSION, "1.0");
        transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
        transformer.setOutputProperty(OutputKeys.STANDALONE, "no");
        StreamResult streamResult = new StreamResult(new File(outputFileName));
        DOMSource source = new DOMSource(doc);
        transformer.transform(source, streamResult);
    }
}
