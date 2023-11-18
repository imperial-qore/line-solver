package jline.solvers.ctmc;

import jline.solvers.ssa.events.Event;
import jline.solvers.ssa.events.SynchedEvent;
import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Pair;

import java.util.Arrays;

public class EventData {
    private final Event event;
    private final Pair<SynchedEvent, Double> outputEventData;
    private final SSAStateMatrix initialNetworkState;
    private final SSAStateMatrix afterEventNetworkState;

    public EventData(Event event, Pair<SynchedEvent, Double> outputEventData, SSAStateMatrix initialNetworkState, SSAStateMatrix afterEventNetworkState) {
        this.event = event;
        this.outputEventData = outputEventData;
        this.initialNetworkState = initialNetworkState;
        this.afterEventNetworkState = afterEventNetworkState;
    }

    public Event getValue0() {
        return this.event;
    }

    public Pair<SynchedEvent, Double> getValue1() {
        return this.outputEventData;
    }

    public SSAStateMatrix getValue2() {
        return this.initialNetworkState;
    }

    public SSAStateMatrix getValue3() {
        return this.afterEventNetworkState;
    }

    @Override
    public int hashCode() {
        Object[] arrayList = new Object[4];
        if(event != null) {
            arrayList[0] = event.getNode().getNodeIdx();
        }
        if(outputEventData != null && outputEventData.getLeft().getNode() != null) {
            arrayList[1] = outputEventData.getLeft().getNode().getNodeIdx();
        }
        if(initialNetworkState != null) {
            arrayList[2] = initialNetworkState;
        }
        if(afterEventNetworkState != null) {
            arrayList[3] = afterEventNetworkState;
        }
        return Arrays.deepHashCode(arrayList);
    }

    @Override
    public boolean equals(Object o){
        assert o instanceof EventData;
        EventData that = (EventData) o;

        // null checks needed because of artificial first event in eventSpace
        if(this.event != that.event) {
            return false;
        }
        if(that.outputEventData != null && !that.outputEventData.getLeft().isDummy()) {
            if(this.outputEventData.getLeft() != that.outputEventData.getLeft()) {
                return false;
            }
        }
        if(that.initialNetworkState != null) {
            if(!this.initialNetworkState.equals(that.initialNetworkState)){
                return false;
            }
        } else {
            return false;
        }
        if(that.afterEventNetworkState != null) {
            return this.afterEventNetworkState.equals(that.afterEventNetworkState);
        } else {
            return false;
        }
    }
}
