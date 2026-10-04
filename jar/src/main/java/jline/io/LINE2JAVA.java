/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;

import java.io.PrintStream;

import static jline.io.InputOutput.line_error;

/**
 * Writes the Java (JLINE) source of a model, dispatching on its type.
 *
 * <p>Port of {@code matlab/src/io/LINE2JAVA.m}: a {@link Network} goes to {@link QN2JAVA} and a
 * {@link LayeredNetwork} to {@link LQN2JAVA}, in both cases under the model's own name.</p>
 */
public final class LINE2JAVA {

    private LINE2JAVA() {
    }

    /** Writes the generated source to standard output. */
    public static void write(Object model) {
        write(model, System.out);
    }

    /** Writes the generated source to a stream, which is flushed but not closed. */
    public static void write(Object model, PrintStream out) {
        if (model instanceof Network) {
            QN2JAVA.write((Network) model, ((Network) model).getName(), out);
        } else if (model instanceof LayeredNetwork) {
            LQN2JAVA.write((LayeredNetwork) model, ((LayeredNetwork) model).getName(), out);
        } else {
            refuse(model);
        }
    }

    /** Writes the generated source to {@code filename}, overwriting it. */
    public static void write(Object model, String filename) {
        if (model instanceof Network) {
            QN2JAVA.write((Network) model, ((Network) model).getName(), filename);
        } else if (model instanceof LayeredNetwork) {
            LQN2JAVA.write((LayeredNetwork) model, ((LayeredNetwork) model).getName(), filename);
        } else {
            refuse(model);
        }
    }

    /** Returns the generated source. */
    public static String generate(Object model) {
        if (model instanceof Network) {
            return QN2JAVA.generate((Network) model);
        } else if (model instanceof LayeredNetwork) {
            return LQN2JAVA.generate((LayeredNetwork) model);
        }
        refuse(model);
        return null;
    }

    private static void refuse(Object model) {
        line_error("LINE2JAVA", "LINE2JAVA expects a Network or a LayeredNetwork, got "
                + (model == null ? "null" : model.getClass().getName()) + ".");
    }
}
