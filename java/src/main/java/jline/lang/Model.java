package jline.lang;

import jline.lang.constant.GlobalConstants;

import java.io.Serializable;

/**
 * Class representing a model supported by the library
 */
public class Model implements Serializable {
    private String name;
    private String lineVersion;

    public Model(String modelName) {
        this.name = modelName;
        this.setVersion(GlobalConstants.Version);
    }

    public String getName() {
        return this.name;
    }
    public void setName(String setName) {
        this.name = setName;
    }
    public String getVersion() { return this.lineVersion; }
    public void setVersion(String version) { this.lineVersion = version; }
}