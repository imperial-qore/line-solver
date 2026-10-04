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
 * Writes the MATLAB script of a model, dispatching on its type.
 *
 * <p>Port of {@code matlab/src/io/LINE2MATLAB.m}: a {@link Network} goes to {@link QN2MATLAB} and a
 * {@link LayeredNetwork} to {@link LQN2MATLAB}, in both cases under the model's own name.</p>
 */
public final class LINE2MATLAB {

    private LINE2MATLAB() {
    }

    /** Writes the generated script to standard output. */
    public static void write(Object model) {
        write(model, System.out);
    }

    /** Writes the generated script to a stream, which is flushed but not closed. */
    public static void write(Object model, PrintStream out) {
        if (model instanceof Network) {
            QN2MATLAB.write((Network) model, ((Network) model).getName(), out);
        } else if (model instanceof LayeredNetwork) {
            LQN2MATLAB.write((LayeredNetwork) model, ((LayeredNetwork) model).getName(), out);
        } else {
            refuse(model);
        }
    }

    /** Writes the generated script to {@code filename}, overwriting it. */
    public static void write(Object model, String filename) {
        if (model instanceof Network) {
            QN2MATLAB.write((Network) model, ((Network) model).getName(), filename);
        } else if (model instanceof LayeredNetwork) {
            LQN2MATLAB.write((LayeredNetwork) model, ((LayeredNetwork) model).getName(), filename);
        } else {
            refuse(model);
        }
    }

    /** Returns the generated script. */
    public static String generate(Object model) {
        if (model instanceof Network) {
            return QN2MATLAB.generate((Network) model);
        } else if (model instanceof LayeredNetwork) {
            return LQN2MATLAB.generate((LayeredNetwork) model);
        }
        refuse(model);
        return null;
    }

    private static void refuse(Object model) {
        line_error("LINE2MATLAB", "LINE2MATLAB expects a Network or a LayeredNetwork, got "
                + (model == null ? "null" : model.getClass().getName()) + ".");
    }
}
