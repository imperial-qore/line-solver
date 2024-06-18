package jline.lang.state;

import jline.util.Matrix;

// Dataholder class for output of afterEvent* functions
public class EventResult {

    public final Matrix outspace;
    public final Matrix outrate;
    public final Matrix outprob;


    public EventResult(Matrix outspace, Matrix outrate, Matrix outprob) {
        this.outspace = outspace;
        this.outrate = outrate;
        this.outprob = outprob;
    }

}
