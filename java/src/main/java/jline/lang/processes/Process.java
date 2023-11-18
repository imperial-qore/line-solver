package jline.lang.processes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import jline.util.NamedParam;

/**
 * A superclass for point processes
 */
public class Process implements Serializable {
    protected String name;
    protected int numParam;
    protected List<NamedParam> params;
    public Process(String name, int numParam) {
        this.name = name;
        this.numParam = numParam;
        this.params = new ArrayList<NamedParam>();
        for (int i = 0; i < this.numParam; i++) {
            this.params.add(new NamedParam("NULL_PARAM", null));
        }
    }

    public int getNumParams() {
        return this.numParam;
    }

    public void setParam(int id, String name, Object value) {
        this.params.set(id-1, new NamedParam(name, value));
    }

    public NamedParam getParam(int id) {
        return this.params.get(id-1);
    }
}