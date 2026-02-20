/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;


import jline.lang.Element;

/**
 * Element of a LayeredNetwork model
 */
public class LayeredNetworkElement extends Element {
    public static final int ACTIVITY = 3;
    public static final int CALL = 4;
    public static final int ENTRY = 2;
    public static final int HOST = 0;
    public static final int PROCESSOR = 0;
    public static final int TASK = 1;
    public LayeredNetwork model;

    public LayeredNetworkElement(String name) {
        super(name);
    }

}
