package jline.lang.layerednetworks;


import jline.lang.Element;

/**
 * Element of a LayeredNetwork model
 */
public class LayeredNetworkElement extends Element {
    public LayeredNetwork model;
    public static final int HOST = 0;
    public static final int PROCESSOR = 0;
    public static final int TASK = 1;
    public static final int ENTRY = 2;
    public static final int ACTIVITY = 3;
    public static final int CALL = 4;


    public LayeredNetworkElement(String name) {
        super(name);
    }

}
