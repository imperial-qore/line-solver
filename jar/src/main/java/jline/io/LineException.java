/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.io;

import jline.GlobalConstants;
import jline.VerboseLevel;

/**
 * The exception every deliberate LINE diagnostic throws.
 *
 * <p>It carries NO stack trace of its own. A {@code line_error} is a message
 * LINE wrote on purpose and its text already names the caller, so the frames
 * between the throw and whatever the user called ({@code runAnalyzerChecks} ->
 * {@code runAnalyzer} -> {@code getAvg} -> {@code getAvgTable} -> ...) are
 * noise: printing them buries a one-sentence refusal such as "finite station
 * capacity at station 'Queue2' is not supported by SolverNC" under a screen of
 * internal frames the user neither chose nor can act on. The JVM's default
 * handler, {@code printStackTrace()} and MATLAB's {@code lang='java'} bridge
 * then all report the message alone. This mirrors the empty-stack throw of
 * {@code matlab/src/io/line_error.m}.</p>
 *
 * <p>The suppression is decided ONCE, when the exception is constructed, from
 * {@link GlobalConstants#Verbose}: at {@link VerboseLevel#DEBUG} the trace is
 * writable and the full chain comes back, which is what a developer wants.
 * Genuine JVM faults ({@code IndexOutOfBoundsException},
 * {@code NullPointerException}) never come through here and are untouched.</p>
 *
 * <p>A CAUSE KEEPS ITS OWN TRACE. Only this wrapper's frames are dropped, so
 * the "Caused by:" block of a re-reported failure still points at the line that
 * actually broke -- the whole reason {@code line_error(caller, msg, cause)}
 * exists.</p>
 *
 * <p>It extends {@code RuntimeException}, which is what {@code line_error} threw
 * before, so every existing {@code catch (RuntimeException)} still catches it.</p>
 */
public class LineException extends RuntimeException {

    private static final long serialVersionUID = 1L;

    /**
     * @param message the already-formatted diagnostic
     */
    public LineException(String message) {
        this(message, null);
    }

    /**
     * @param message the already-formatted diagnostic
     * @param cause   the failure being re-reported, or null
     */
    public LineException(String message, Throwable cause) {
        // (message, cause, enableSuppression, writableStackTrace). Suppression
        // stays on so try-with-resources behaves as usual; only the trace goes.
        super(message, cause, true, GlobalConstants.Verbose == VerboseLevel.DEBUG);
    }
}
