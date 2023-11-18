package jline.lang.layerednetworks;

import jline.lang.constant.SchedStrategy;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import jline.lang.layerednetworks.*;

class EntryTest {

    @Test
    void on() {
        LayeredNetwork lqn = new LayeredNetwork("LQN");
        Task task1 = new Task(lqn, "T1", 1, SchedStrategy.REF);
        Entry entry = new Entry(lqn, "E1");
        entry.on(task1);
        assertEquals(task1.entries.get(0).getName(),"E1");
        assertEquals(entry.parent.getName(),"T1");
    }
}