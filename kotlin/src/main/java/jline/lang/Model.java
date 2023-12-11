package jline.lang;

import java.io.Serializable;

/**
 * Class representing a model supported by the library
 */
public class Model implements Serializable {
    private String name;
    private final String lineVersion;

    public Model(String modelName) {
        this.name = modelName;
        this.lineVersion = "2.0.31";
    }

    public String getName() {
        return this.name;
    }
    public void setName(String setName) {
        this.name = setName;
    }
    public String getVersion() { return this.lineVersion; }
}