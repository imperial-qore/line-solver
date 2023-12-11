package jline.lang.sections;

import java.io.Serializable;

import jline.lang.NetworkElement;

/**
 * A general class modeling a node section
 */
public abstract class Section extends NetworkElement implements Serializable {
    String className;
    public Section(String className) {
        super("Section");
        this.className = className;
    }

    public String getClassName() {
        return className;
    }
}
