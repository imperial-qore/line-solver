package jline.lang;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Class for auxiliary information stored in Network objects
 */
public class NetworkAttribute implements Serializable {
    private int clientIdx;
    private int serverIdx;
    private int sourceIdx;

    private final Map<Integer,Integer[]> hosts;
    private final Map<Integer,Integer[]> tasks;
    private final Map<Integer,Integer[]> entries;
    private final Map<Integer,Integer[]> calls;
    private final Map<Integer,Integer[]> activities;

    public NetworkAttribute(){
        this.hosts = new HashMap<>();
        this.tasks =new HashMap<>();
        this.entries = new HashMap<>();
        this.calls =new HashMap<>();
        this.activities =new HashMap<>();
    }

    public int getClientIdx() {
        return clientIdx;
    }

    public int getServerIdx() {
        return serverIdx;
    }

    public int getSourceIdx() {
        return sourceIdx;
    }

    public Map<Integer, Integer[]> getActivities() {
        return activities;
    }

    public Map<Integer, Integer[]> getCalls() {
        return calls;
    }

    public Map<Integer, Integer[]> getEntries() {
        return entries;
    }

    public Map<Integer, Integer[]> getHosts() {
        return hosts;
    }

    public Map<Integer, Integer[]> getTasks() {
        return tasks;
    }

    public void setSourceIdx(int sourceIdx) {
        this.sourceIdx = sourceIdx;
    }

    public void setClientIdx(int clientIdx) {
        this.clientIdx = clientIdx;
    }

    public void setServerIdx(int serverIdx) {
        this.serverIdx = serverIdx;
    }

    public void addActivities(Integer[] newActivities) {
        this.activities.put(activities.size()+1,newActivities);
    }

    public void addCalls(Integer[] newCalls) {
        this.calls.put(calls.size()+1,newCalls);
    }

    public void addEntries(Integer[] newEntries) {
        this.entries.put(entries.size()+1,newEntries);
    }

    public void addHosts(Integer[] newHosts) {
        this.hosts.put(hosts.size()+1,newHosts);
    }

    public void addTasks(Integer[] newTasks) {
        this.tasks.put(tasks.size()+1,newTasks);
    }
}
