package jline.lang.layered;

import jline.lang.constant.SchedStrategy;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HostTest {
    LayeredNetwork lqn = new LayeredNetwork("LQN");
    Host host = new Host(lqn, "H1", 1, SchedStrategy.INF);
    @Test
    void addTask() {
        Task task1 = new Task(lqn, "T1", 1, SchedStrategy.REF);
        Task task2 = new Task(lqn, "T2", 1, SchedStrategy.PS);
        host.addTask(task1);
        host.addTask(task2);
        assertEquals(host.tasks.get(0).getName(),"T1");
        assertEquals(host.tasks.get(1).getName(),"T2");
    }
}