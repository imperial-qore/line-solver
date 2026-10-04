/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.constant.CallType;
import jline.lang.layered.Activity;
import jline.lang.layered.CacheTask;
import jline.lang.layered.Entry;
import jline.lang.layered.Host;
import jline.lang.layered.ItemEntry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.layered.Task;

import java.util.ArrayList;
import java.util.List;

import static jline.io.InputOutput.line_error;

/**
 * Refuses, by name, the layered constructs that the struct-driven {@link LQN2JAVA} generator cannot write.
 *
 * <p>{@code LQN2JAVA.m} writes processors, tasks with multiplicity, replication and think time, entries,
 * activities with their host demand and think time, bindings, calls and replies, and activity precedences. Anything
 * else a model carries would be dropped silently, and the program would rebuild a different model; this
 * check turns each such construct into an error that names it.</p>
 */
final class LqnChecks {

    private LqnChecks() {
    }

    static void refuseUnsupported(LayeredNetwork model, LayeredNetworkStruct sn, String caller) {
        List<String> found = new ArrayList<String>();
        for (Host h : model.getHosts().values()) {
            if (h.getQuantum() != 0.001 || h.getSpeedFactor() != 1) {
                found.add("processor " + h.getName() + " has a quantum or speed factor");
            }
            if (h.hasLinearConstraints() || h.hasRateDependence()) {
                found.add("processor " + h.getName() + " has admission constraints, rate dependence or server pools");
            }
        }
        for (Task t : model.getTasks().values()) {
            if (t.getClass() != Task.class) {
                found.add("task " + t.getName() + " is a " + t.getClass().getSimpleName());
            }
            if (t instanceof CacheTask || t.hasSetupDelayoff()) {
                found.add("task " + t.getName() + " has cache or setup/delay-off behaviour");
            }
            if (t.getPriority() != 0) {
                found.add("task " + t.getName() + " has priority " + t.getPriority());
            }
            if (t.getFanInValue() > 0 || !t.getFanOutMap().isEmpty()) {
                found.add("task " + t.getName() + " has fan-in or fan-out");
            }
            if (t.hasLinearConstraints() || t.hasRateDependence()) {
                found.add("task " + t.getName() + " has admission constraints, rate dependence or server pools");
            }
        }
        for (Entry e : model.getEntries().values()) {
            if (e instanceof ItemEntry) {
                found.add("entry " + e.getName() + " is an ItemEntry");
            }
            if (e.getArrival() != null) {
                found.add("entry " + e.getName() + " has an open arrival");
            }
            if (!"PH1PH2".equals(e.getType())) {
                found.add("entry " + e.getName() + " has type " + e.getType());
            }
        }
        for (Activity a : model.getActivities().values()) {
            if (a.getPhase() != 1) {
                found.add("activity " + a.getName() + " is in phase " + a.getPhase());
            }
            if (!a.getSyncCallGroups().isEmpty()) {
                found.add("activity " + a.getName() + " has a call group");
            }
            if (!"STOCHASTIC".equals(a.getCallOrder())) {
                found.add("activity " + a.getName() + " has call order " + a.getCallOrder());
            }
        }
        if (sn.calltype != null) {
            for (CallType ct : sn.calltype.values()) {
                if (ct == CallType.FWD) {
                    found.add("an entry forwards its calls");
                    break;
                }
            }
        }
        if (!found.isEmpty()) {
            StringBuilder b = new StringBuilder(caller).append(" cannot write this model: ");
            for (int i = 0; i < found.size(); i++) {
                b.append(i > 0 ? "; " : "").append(found.get(i));
            }
            b.append('.');
            line_error(caller, b.toString());
        }
    }
}
