package jline.lang.nodes;

import org.apache.commons.lang3.NotImplementedException;

/**
 * Place as in a stochastic Petri net model
 */
public class Place extends Station{
    public Place(String name) {
        super(name);
        throw new NotImplementedException("Place station not implemented");
    }
}
