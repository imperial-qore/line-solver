package jline.lang.sections;

import java.io.Serializable;
import java.util.List;
import jline.lang.*;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;
import jline.lang.*;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.sections.*;

/**
 * Output section of a queue that applies class switching upon job departure
 */
public class ClassSwitchOutputSection extends Dispatcher implements Serializable {
    public ClassSwitchOutputSection(List<JobClass> customerClasses) {
        super(customerClasses);
        this.className = "ClassSwitchDispatcher";
        this.isClassSwitch = true;
    }
}
