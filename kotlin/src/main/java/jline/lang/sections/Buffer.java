package jline.lang.sections;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import jline.lang.*;
import jline.lang.constant.SchedStrategyType;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

/**
 * Input buffer of a station
 */
public class Buffer extends InputSection implements Serializable {
    protected int size;
    public Buffer(List<JobClass> classes) {
        super("Buffer");

        this.size = -1;
        this.schedPolicy = SchedStrategyType.NP;
    }

    public void setSize(int size) {
        this.size = size;
    }
}
