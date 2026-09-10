/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Class for auxiliary information stored in Network objects
 */
public class NetworkAttribute extends ModelAttribute implements Serializable {
    private final Map<Integer, Integer[]> hosts;
    private final Map<Integer, Integer[]> tasks;
    private final Map<Integer, Integer[]> entries;
    private final Map<Integer, Integer[]> calls;
    private final Map<Integer, Integer[]> activities;
    private int clientIdx;
    private int serverIdx;
    private int sourceIdx;
    // LQN server element -> its station index in this layer, and the station
    // indices of the layer's host and task servers. Under 'srvn' layering there
    // is exactly one server, under 'flat' layering there is one per element.
    private final Map<Integer, Integer> serverIdxOf;
    private final java.util.List<Integer> hostStations;
    private final java.util.List<Integer> taskStations;

    public NetworkAttribute() {
        this.hosts = new HashMap<>();
        this.tasks = new HashMap<>();
        this.entries = new HashMap<>();
        this.calls = new HashMap<>();
        this.activities = new HashMap<>();
        this.serverIdxOf = new HashMap<>();
        this.hostStations = new java.util.ArrayList<>();
        this.taskStations = new java.util.ArrayList<>();
    }

    public Map<Integer, Integer> getServerIdxOf() {
        return serverIdxOf;
    }

    public void putServerIdxOf(int elemIdx, int stationIdx) {
        this.serverIdxOf.put(elemIdx, stationIdx);
    }

    public java.util.List<Integer> getHostStations() {
        return hostStations;
    }

    public java.util.List<Integer> getTaskStations() {
        return taskStations;
    }

    public void addActivities(Integer[] newActivities) {
        this.activities.put(activities.size() + 1, newActivities);
    }

    public void addCalls(Integer[] newCalls) {
        this.calls.put(calls.size() + 1, newCalls);
    }

    public void addEntries(Integer[] newEntries) {
        this.entries.put(entries.size() + 1, newEntries);
    }

    public void addHosts(Integer[] newHosts) {
        this.hosts.put(hosts.size() + 1, newHosts);
    }

    public void addTasks(Integer[] newTasks) {
        this.tasks.put(tasks.size() + 1, newTasks);
    }

    public Map<Integer, Integer[]> getActivities() {
        return activities;
    }

    public Map<Integer, Integer[]> getCalls() {
        return calls;
    }

    public int getClientIdx() {
        return clientIdx;
    }

    public void setClientIdx(int clientIdx) {
        this.clientIdx = clientIdx;
    }

    public Map<Integer, Integer[]> getEntries() {
        return entries;
    }

    public Map<Integer, Integer[]> getHosts() {
        return hosts;
    }

    public int getServerIdx() {
        return serverIdx;
    }

    public void setServerIdx(int serverIdx) {
        this.serverIdx = serverIdx;
    }

    public int getSourceIdx() {
        return sourceIdx;
    }

    public void setSourceIdx(int sourceIdx) {
        this.sourceIdx = sourceIdx;
    }

    public Map<Integer, Integer[]> getTasks() {
        return tasks;
    }
}
