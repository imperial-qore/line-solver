/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.lang.Event;
import jline.lang.constant.EventType;
import jline.solvers.ssa.SampleSysState;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static jline.api.mc.Ctmc_makeinfgenKt.ctmc_makeinfgen;
import static jline.api.mc.Dtmc_makestochasticKt.dtmc_makestochastic;

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
     * Create MarkedCTMC from sample system aggregation.
     * Matches MATLAB MarkedMarkovProcess.fromSampleSysAggr().
     *
     * Constructs a MarkedMarkovProcess from sampled state trajectories by:
     * 1. Computing the Cartesian product of per-node states to get system states
     * 2. Building a DTMC from transition counts
     * 3. Computing holding times
     * 4. Creating the infinitesimal generator
     * 5. Building event filter matrices from the event list
     *
     * @param sa the sample aggregation (SampleSysState from SSA/CTMC solver)
     * @return MarkedCTMC constructed from samples
     */
    @SuppressWarnings("unchecked")
    public static MarkedMarkovProcess fromSampleSysAggr(Object sa) {
        if (!(sa instanceof SampleSysState)) {
            throw new IllegalArgumentException(
                    "fromSampleSysAggr requires a SampleSysState object from SSA or CTMC solver");
        }
        SampleSysState sampleData = (SampleSysState) sa;

        // Step 1: Compute Cartesian product of per-node states
        List<Matrix> stateList = sampleData.state;
        Matrix sampleState = stateList.get(0);
        for (int r = 1; r < stateList.size(); r++) {
            sampleState = cartesianProduct(sampleState, stateList.get(r));
        }

        // Step 2: Get unique states and hash indices
        int nSamples = sampleState.getNumRows();
        int nCols = sampleState.getNumCols();

        // Map each unique row to an index
        LinkedHashMap<String, Integer> stateMap = new LinkedHashMap<String, Integer>();
        int[] stateHash = new int[nSamples];
        List<double[]> uniqueStatesList = new ArrayList<double[]>();

        for (int i = 0; i < nSamples; i++) {
            double[] row = new double[nCols];
            for (int c = 0; c < nCols; c++) {
                row[c] = sampleState.get(i, c);
            }
            String key = Arrays.toString(row);
            Integer idx = stateMap.get(key);
            if (idx == null) {
                idx = stateMap.size();
                stateMap.put(key, idx);
                uniqueStatesList.add(row);
            }
            stateHash[i] = idx;
        }

        int nStates = stateMap.size();

        // Build state space matrix
        Matrix stateSpace = new Matrix(nStates, nCols);
        for (int i = 0; i < nStates; i++) {
            double[] row = uniqueStatesList.get(i);
            for (int c = 0; c < nCols; c++) {
                stateSpace.set(i, c, row[c]);
            }
        }

        // Step 3: Build transition count matrix and holding times
        Matrix dtmc = new Matrix(nStates, nStates);
        double[] holdTime = new double[nStates];

        // Identify unique active events for event filter
        List<Event> events = sampleData.event;
        List<String> uniqueEventKeys = new ArrayList<String>();
        Map<String, Integer> eventKeyMap = new HashMap<String, Integer>();
        int[] evType = new int[nSamples];

        if (events != null && !events.isEmpty()) {
            for (int i = 0; i < events.size(); i++) {
                Event ev = events.get(i);
                String key = ev.getEvent().ordinal() + ":" + ev.getNode() + ":" + ev.getJobClass();
                Integer idx = eventKeyMap.get(key);
                if (idx == null) {
                    idx = uniqueEventKeys.size();
                    uniqueEventKeys.add(key);
                    eventKeyMap.put(key, idx);
                }
                if (i < nSamples) {
                    evType[i] = idx;
                }
            }
        }

        int nEventTypes = Math.max(uniqueEventKeys.size(), 1);
        Matrix[] eventFiltArr = new Matrix[nEventTypes];
        for (int e = 0; e < nEventTypes; e++) {
            eventFiltArr[e] = new Matrix(nStates, nStates);
        }

        // Count transitions
        for (int i = 1; i < nSamples; i++) {
            int from = stateHash[i - 1];
            int to = stateHash[i];
            dtmc.set(from, to, dtmc.get(from, to) + 1);
            holdTime[from] += sampleData.t.get(i, 0) - sampleData.t.get(i - 1, 0);

            if (events != null && !events.isEmpty() && i - 1 < events.size()) {
                int e = evType[i - 1];
                eventFiltArr[e].set(from, to, eventFiltArr[e].get(from, to) + 1);
            }
        }

        // Step 4: Normalize
        // holdTime = holdTime / sum(dtmc, 2)  (row sums)
        for (int i = 0; i < nStates; i++) {
            double rowSum = 0;
            for (int j = 0; j < nStates; j++) {
                rowSum += dtmc.get(i, j);
            }
            if (rowSum > 0) {
                holdTime[i] /= rowSum;
            }
        }

        // Normalize event filters
        for (int e = 0; e < nEventTypes; e++) {
            for (int i = 0; i < nStates; i++) {
                double rowSum = 0;
                for (int j = 0; j < nStates; j++) {
                    rowSum += dtmc.get(i, j);
                }
                if (rowSum > 0) {
                    for (int j = 0; j < nStates; j++) {
                        eventFiltArr[e].set(i, j, eventFiltArr[e].get(i, j) / rowSum);
                    }
                }
            }
        }

        // Step 5: Create infinitesimal generator
        // infGen = ctmc_makeinfgen(dtmc_makestochastic(dtmc) ./ (holdTime * ones(1,n)))
        Matrix stochDtmc = dtmc_makestochastic(dtmc);
        Matrix scaledDtmc = new Matrix(nStates, nStates);
        for (int i = 0; i < nStates; i++) {
            double ht = holdTime[i];
            if (ht > 0) {
                for (int j = 0; j < nStates; j++) {
                    scaledDtmc.set(i, j, stochDtmc.get(i, j) / ht);
                }
            }
        }
        Matrix infGen = ctmc_makeinfgen(scaledDtmc);

        // Step 6: Scale event filters by departure rates
        // eventFilt{e} = eventFilt{e} .* (-diag(infGen) * ones(1,n))
        MatrixCell eventFilt = new MatrixCell(nEventTypes);
        for (int e = 0; e < nEventTypes; e++) {
            for (int i = 0; i < nStates; i++) {
                double depRate = -infGen.get(i, i);
                for (int j = 0; j < nStates; j++) {
                    eventFiltArr[e].set(i, j, eventFiltArr[e].get(i, j) * depRate);
                }
            }
            eventFilt.set(e, eventFiltArr[e]);
        }

        // Step 7: Build event list
        List<Map<String, Object>> evList = new ArrayList<Map<String, Object>>();
        if (events != null) {
            for (Event ev : events) {
                Map<String, Object> evMap = new HashMap<String, Object>();
                List<Event> activeList = new ArrayList<Event>();
                activeList.add(ev);
                evMap.put("active", activeList);
                evMap.put("passive", new ArrayList<Event>());
                evList.add(evMap);
            }
        }

        return new MarkedMarkovProcess(infGen, eventFilt, evList, true, stateSpace);
    }

    /**
     * Compute Cartesian product of two state matrices (row-wise).
     * Each row of the result is the concatenation of a row from A with a row from B,
     * paired by row index.
     */
    private static Matrix cartesianProduct(Matrix A, Matrix B) {
        int nRows = A.getNumRows();
        int nColsA = A.getNumCols();
        int nColsB = B.getNumCols();
        Matrix result = new Matrix(nRows, nColsA + nColsB);
        for (int i = 0; i < nRows; i++) {
            for (int c = 0; c < nColsA; c++) {
                result.set(i, c, A.get(i, c));
            }
            for (int c = 0; c < nColsB; c++) {
                result.set(i, nColsA + c, B.get(i, c));
            }
        }
        return result;
    }
}