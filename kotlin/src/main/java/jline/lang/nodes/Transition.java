package jline.lang.nodes;

import org.apache.commons.lang3.NotImplementedException;

/**
 * Transition as in a stochastic Petri net model
 */
public class Transition extends Node{
    public Transition(String nodeName) {
        super(nodeName);
        throw new NotImplementedException("Transition node not implemented");
    }
}
