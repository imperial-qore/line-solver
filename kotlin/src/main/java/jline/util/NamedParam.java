package jline.util;

import java.io.Serializable;

/**
 * A class to store a named parameter.
 */
public class NamedParam implements Serializable {
    protected String name;
    protected Object value;
    public NamedParam(String name, Object value) {
        this.name = name;
        this.value = value;
    }

    public Object getValue() {
        return this.value;
    }
}
