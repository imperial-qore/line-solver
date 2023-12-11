package jline.lang.layerednetworks;

import jline.lang.layerednetworks.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TaskTest {
    LayeredNetwork lqn = new LayeredNetwork("LQN");
    Processor processor = new Processor(lqn, "P1", 1, SchedStrategy.INF);
    Task task = new Task(lqn,"T1",1,SchedStrategy.PS);

    @Test
    void on() {
        task.on(processor);
        assertEquals(task.parent.getName(),"P1");
    }

    @Test
    void setAsReferenceTask() {
        task.setAsReferenceTask();
        assertEquals(task.scheduling,SchedStrategy.REF);
    }

    @Test
    void removeActivity() {
        Activity activity1 = new Activity(lqn, "A1", new Exp(1.0));
        activity1.on(task);
        Activity activity2 = new Activity(lqn, "A2", new Exp(1.0));
        activity2.on(task);
        task.removeActivity(0);
        assertEquals(task.activities.get(0).getName(),"A2");
    }

    @Test
    void setThinkTime() {
        task.setThinkTime(0.02);
        assertEquals(task.thinkTimeSCV,1.0);
    }
}