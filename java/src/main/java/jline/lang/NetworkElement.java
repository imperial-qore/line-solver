package jline.lang;

import java.io.Serializable;

/**
 * Class representing an element within a Network object
 */
public class NetworkElement extends Element {
    protected String name;
    public NetworkElement(String neName) {
        super(neName);
        this.name = neName;
    }

    public String getName() { return this.name; }
}
