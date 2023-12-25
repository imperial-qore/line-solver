package jline.lang.layered;

import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ActivityTest {
    LayeredNetwork lqn = new LayeredNetwork("LQN");
    Activity activity = new Activity(lqn, "A1", new Exp(1.0));
    Task task = new Task(lqn,"T1",1, SchedStrategy.PS);

    @Test
    void on() {
        activity.on(task);
        assertEquals(activity.parent.getName(),"T1");
    }

    @Test
    void setHostDemand() {
        activity.setHostDemand(0.01);
        assertEquals(activity.hostDemandSCV,1.0);
    }

    @Test
    void repliesTo() throws Exception {
        activity.on(task);
        Entry entry = new Entry(lqn,"E1");
        activity.repliesTo(entry);
        assertEquals(entry.replyActivity.get(0),activity.getName());
    }

    @Test
    void boundTo() {
        Entry entry = new Entry(lqn,"E1");
        activity.boundTo(entry);
        assertEquals(activity.boundToEntry,"E1");
    }

    @Test
    void synchCall() {
        Entry entry = new Entry(lqn,"E1");
        activity.synchCall(entry);
        assertEquals(activity.syncCallDests.get(0),"E1");
    }

    @Test
    void asynchCall() {
        Entry entry = new Entry(lqn,"E1");
        activity.asynchCall(entry);
        assertEquals(activity.asyncCallDests.get(0),"E1");
    }
}