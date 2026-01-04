/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io

import jline.GlobalConstants
import jline.VerboseLevel
import org.apache.commons.io.FileUtils
import org.w3c.dom.Document
import java.io.File
import java.util.logging.Level
import java.util.logging.Logger
import java.util.logging.Formatter
import java.util.logging.LogRecord
import java.util.logging.Handler
import java.util.logging.ConsoleHandler
import java.util.logging.LogManager
import javax.xml.transform.OutputKeys
import javax.xml.transform.TransformerException
import javax.xml.transform.TransformerFactory
import javax.xml.transform.dom.DOMSource
import javax.xml.transform.stream.StreamResult

/**
 * Custom formatter that removes timestamps and formats messages on single lines
 */
private class SimpleLineFormatter : Formatter() {
    override fun format(record: LogRecord): String {
        val level = when (record.level) {
            Level.SEVERE -> "SEVERE"
            Level.WARNING -> "WARNING"
            Level.INFO -> "INFO"
            else -> record.level.name
        }
        return "$level: ${record.message}\n"
    }
}

/**
 * Functions to print on screen
 */
private var lastWarning = ""
private var suppressedWarnings = false
private var suppressedWarningTic = System.currentTimeMillis()
private var loggerConfigured = false

/**
 * Check if running under Maven test
 */
private fun isRunningUnderMavenTest(): Boolean {
    // Check for common Maven Surefire/Failsafe properties
    return System.getProperty("surefire.test.class.path") != null ||
           System.getProperty("maven.test.skip") != null ||
           System.getProperty("basedir") != null && System.getProperty("surefire.real.class.path") != null ||
           // Also check if any test framework is on the classpath
           Thread.currentThread().stackTrace.any { 
               it.className.contains("org.junit") || 
               it.className.contains("org.testng") ||
               it.className.contains("maven.surefire")
           }
}

/**
 * Configure logger to use custom formatting without timestamps globally
 */
private fun configureLogger() {
    // Determine logging level based on GlobalConstants.Verbose and test environment
    val logLevel = when {
        GlobalConstants.Verbose == VerboseLevel.SILENT -> Level.OFF
        else -> Level.WARNING
    }

    if (!loggerConfigured) {
        // Configure all existing loggers and handlers
        val logManager = LogManager.getLogManager()
        val loggerNames = logManager.loggerNames

        // Configure root logger first
        val rootLogger = Logger.getLogger("")

        // Remove all existing handlers and add our custom one
        rootLogger.handlers.forEach { handler ->
            rootLogger.removeHandler(handler)
        }

        val rootConsoleHandler = ConsoleHandler()
        rootConsoleHandler.formatter = SimpleLineFormatter()
        rootConsoleHandler.level = logLevel
        rootLogger.addHandler(rootConsoleHandler)
        rootLogger.level = logLevel

        // Configure our specific logger to use the root logger
        val logger = Logger.getLogger(FileUtils::class.java.getName())
        logger.useParentHandlers = true
        logger.level = logLevel

        loggerConfigured = true
    } else {
        // If already configured, just update the levels based on current GlobalConstants.Verbose
        val rootLogger = Logger.getLogger("")
        rootLogger.level = logLevel
        rootLogger.handlers.forEach { handler ->
            handler.level = logLevel
        }
        Logger.getLogger(FileUtils::class.java.getName()).level = logLevel
    }
}

private fun line_printf(message: String?) {
    configureLogger()
    Logger.getLogger(FileUtils::class.java.getName()).log(Level.INFO, message)
}

/**
 * Print debug message if verbose level is DEBUG
 */
fun line_debug(verbose: VerboseLevel, message: String) {
    if (verbose == VerboseLevel.DEBUG) {
        line_printf("[DEBUG] $message")
    }
}

fun line_warning(caller: String?, msg: String, vararg args: Any?) {
    val logger = Logger.getLogger(FileUtils::class.java.getName())

    // Check if global verbosity level allows warnings
    if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
        return
    }

    val errmsg = String.format(msg, *args)
    val finalmsg = String.format("[%s] %s", caller, errmsg)

    try {
        val currentTime = System.currentTimeMillis()

        if (finalmsg.compareTo(lastWarning) != 0 || (currentTime - suppressedWarningTic) > 60000) {
            line_printf(finalmsg, Level.WARNING)
            lastWarning = finalmsg
            suppressedWarnings = false
            suppressedWarningTic = currentTime
        } else {
            if (!suppressedWarnings) {
                line_printf(String.format("[%s] %s",
                    caller,
                    "Message casted more than once, repetitions will not be printed on screen for 60 seconds."),
                    Level.WARNING)
                suppressedWarnings = true
                suppressedWarningTic = currentTime
            }
        }
    } catch (e: Exception) {
        logger.log(Level.SEVERE, "Exception in line_warning", e)
    }
}

private fun line_printf(message: String?, level: Level?) {
    configureLogger()
    Logger.getLogger(FileUtils::class.java.getName()).log(level, message)
}


fun line_error(caller: String?, msg: String?) {
    val logger = Logger.getLogger(FileUtils::class.java.getName())
    
    // Check if global verbosity level allows errors
    if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
        // Still throw the exception but don't print to console
        val finalmsg = String.format("[%s] %s", caller, msg)
        throw RuntimeException(finalmsg)
    }
    
    val finalmsg = String.format("[%s] %s", caller, msg) 
    line_printf(finalmsg, Level.SEVERE)
    //System.exit(1);
    throw RuntimeException(finalmsg)
}


fun mfilename(obj: Any): String {
    val enclosingMethod = obj.javaClass.getEnclosingMethod()
    return if (enclosingMethod != null) {
        enclosingMethod.name
    } else {
        obj.javaClass.simpleName
    }
}

/**
 * Writes the given XML Document object to a specified output file.
 *
 * @param outputFileName The name of the output file where the XML data will be written.
 * @param doc            The Document object containing the XML data.
 * @throws TransformerException If an unrecoverable error occurs during the transformation.
 */

@Throws(TransformerException::class)
fun writeXML(outputFileName: String, doc: Document?) {
    val transformer = TransformerFactory.newInstance().newTransformer()
    transformer.setOutputProperty(OutputKeys.INDENT, "yes")
    transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2")
    transformer.setOutputProperty(OutputKeys.VERSION, "1.0")
    transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8")
    transformer.setOutputProperty(OutputKeys.STANDALONE, "no")
    val streamResult = StreamResult(File(outputFileName))
    val source = DOMSource(doc)
    transformer.transform(source, streamResult)
}

