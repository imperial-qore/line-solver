/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

/**
 * Run one example by name, the way {@code line-examples} runs one C++ example.
 *
 * <p>Every other example class exposes a {@code main} that runs its whole
 * family in sequence, which is the right shape for a person reading one file and
 * the wrong shape for a harness that needs the output of a single example
 * attributed to that example. {@code doc/latex/harvest-java.py} uses this to
 * capture per-example evidence for the Java edition of the manual, so a page can
 * show what its own listing prints rather than what its family prints.
 *
 * <pre>
 *   java -cp common/jline.jar jline.examples.ExampleRunner &lt;class&gt; &lt;method&gt;
 * </pre>
 *
 * <p>The method must be public, static and take no arguments -- the shape every
 * example demo already has. The exit status is 0 when the example returned and 1
 * when it threw, so a sweep is usable from a shell without parsing the output.
 */
public class ExampleRunner {

    public static void main(String[] args) {
        if (args.length != 2) {
            System.err.println("usage: ExampleRunner <fully.qualified.Class> <method>");
            System.exit(2);
            return;
        }
        String className = args[0];
        String methodName = args[1];

        Class<?> owner;
        try {
            owner = Class.forName(className);
        } catch (ClassNotFoundException e) {
            System.err.println("no such class: " + className);
            System.exit(2);
            return;
        }

        Method target = null;
        Method[] declared = owner.getMethods();
        for (int i = 0; i < declared.length; i++) {
            Method m = declared[i];
            if (m.getName().equals(methodName)
                    && m.getParameterTypes().length == 0
                    && Modifier.isStatic(m.getModifiers())) {
                target = m;
                break;
            }
        }
        if (target == null) {
            System.err.println("no no-arg static method " + methodName + " on " + className);
            System.exit(2);
            return;
        }

        try {
            Object result = target.invoke(null);
            // A builder returns the model rather than printing it. Showing the
            // model's own toString is closer to the reference than showing
            // nothing, and keeps a builder distinguishable from a silent demo.
            if (result != null) {
                System.out.println(result);
            }
        } catch (IllegalAccessException e) {
            System.err.println("cannot invoke " + className + "." + methodName + ": " + e);
            System.exit(1);
        } catch (InvocationTargetException e) {
            Throwable cause = e.getCause() == null ? e : e.getCause();
            System.err.println("FAILED " + methodName + ": " + cause);
            System.exit(1);
        }
    }
}
