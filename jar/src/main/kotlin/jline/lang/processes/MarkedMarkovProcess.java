/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.lang.Event;
import jline.lang.constant.EventType;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.List;
import java.util.Map;

/**
 * A class for continuous time Markov chain where transitions are labeled
 */
public class MarkedMarkovProcess extends MarkovProcess {
    protected MatrixCell eventFilt;
    protected List<Map<String, Object>> eventList;

    /**
     * Creates a MarkedCTMC with the specified generator, event filters, and events
     * @param infGen the infinitesimal generator matrix
     * @param eventFilt the event filter matrices
     * @param evs the event list
     */
    public MarkedMarkovProcess(Matrix infGen, MatrixCell eventFilt, List<Map<String, Object>> evs) {
        this(infGen, eventFilt, evs, true, null);
    }

    /**
     * Creates a MarkedCTMC with the specified generator, event filters, events, and finite flag
     * @param infGen the infinitesimal generator matrix
     * @param eventFilt the event filter matrices
     * @param evs the event list
     * @param isFinite whether the CTMC is finite
     */
    public MarkedMarkovProcess(Matrix infGen, MatrixCell eventFilt, List<Map<String, Object>> evs, boolean isFinite) {
        this(infGen, eventFilt, evs, isFinite, null);
    }

    /**
     * Creates a MarkedCTMC with the specified generator, event filters, events, finite flag, and state space
     * @param infGen the infinitesimal generator matrix
     * @param eventFilt the event filter matrices
     * @param evs the event list
     * @param isFinite whether the CTMC is finite
     * @param stateSpace the state space representation
     */
    public MarkedMarkovProcess(Matrix infGen, MatrixCell eventFilt, List<Map<String, Object>> evs, boolean isFinite, Matrix stateSpace) {
        super(infGen, isFinite, stateSpace);
        this.name = "MarkedCTMC";
        this.eventFilt = eventFilt;
        this.eventList = evs;
    }

    /**
     * Convert to MAP for a specific event
     * @param ev the event to convert to MAP
     * @return the MAP representation
     */
    public MAP toMAP(Event ev) {
        for (int e = 0; e < eventList.size(); e++) {
            Map<String, Object> eventItem = eventList.get(e);
            @SuppressWarnings("unchecked")
            List<Event> activeEvents = (List<Event>) eventItem.get("active");
            
            if (activeEvents != null && !activeEvents.isEmpty()) {
                Event activeEvent = activeEvents.get(0);
                if (activeEvent.getNode() == ev.getNode() && 
                    activeEvent.getJobClass() == ev.getJobClass() &&
                    activeEvent.getEvent() == ev.getEvent()) {
                    
                    Matrix D0 = infGen.sub(eventFilt.get(e));
                    Matrix D1 = eventFilt.get(e);
                    return new MAP(D0, D1);
                }
            }
        }
        throw new IllegalArgumentException("Event not found in event list");
    }

    /**
     * Convert to MAP for a specific event type, node, and class
     * @param evtype the event type
     * @param node the node index
     * @param jobclass the job class index
     * @return the MAP representation
     */
    public MAP toMAP(EventType evtype, int node, int jobclass) {
        Event ev = new Event(evtype, node, jobclass);
        return toMAP(ev);
    }

    /**
     * Solve for embedded probabilities for specified event set
     * @param evset the set of event indices to solve for (null for all events)
     * @return the embedded probabilities for each event
     */
    public MatrixCell embeddedSolve(int[] evset) {
        if (evset == null) {
            evset = new int[eventFilt.size()];
            for (int i = 0; i < eventFilt.size(); i++) {
                evset[i] = i;
            }
        }
        
        MatrixCell pie = new MatrixCell(evset.length);
        
        for (int i = 0; i < evset.length; i++) {
            int ev = evset[i];
            Matrix D0 = infGen.sub(eventFilt.get(ev));
            Matrix D1 = eventFilt.get(ev);
            MAP eMAP = new MAP(D0, D1);
            Matrix embeddedProb = eMAP.getEmbeddedProb();
            pie.set(i, embeddedProb);
        }
        
        return pie;
    }

    /**
     * Solve for embedded probabilities for all events
     * @return the embedded probabilities for each event
     */
    public MatrixCell embeddedSolve() {
        return embeddedSolve(null);
    }

    /**
     * Get the event filter matrices
     * @return the event filter matrices
     */
    public MatrixCell getEventFilt() {
        return eventFilt;
    }

    /**
     * Get the event list
     * @return the event list
     */
    public List<Map<String, Object>> getEventList() {
        return eventList;
    }

    /**
     * Create MarkedCTMC from sample system aggregation
     * @param sa the sample aggregation
     * @return MarkedCTMC constructed from samples
     */
    public static MarkedMarkovProcess fromSampleSysAggr(Object sa) {
        throw new UnsupportedOperationException("fromSampleSysAggr not yet implemented");
    }
}