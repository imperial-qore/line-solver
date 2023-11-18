package jline.solvers.ssa.state;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class LCFSBuffer extends SSAStateCell {
    /*
        Last come-first serve without preemption
     */

	public Deque<Integer> deque;
    public Deque<Integer> serverQueue;

    protected int nServers;
    protected int[] inQueue;
    protected int totalInQueue;

    protected PhaseList phaseList;

    private void moveToService() {
        if (this.deque.isEmpty()) {
            return;
        }

        int classToAdd = this.deque.removeFirst();
        this.phaseList.addToService(classToAdd);
        this.serverQueue.addLast(classToAdd);
    }

    public LCFSBuffer(int nClasses, int nServers, PhaseList phaseList) {
        this.deque = new ArrayDeque<Integer>();
        this.serverQueue = new ArrayDeque<Integer>();

        this.nServers = nServers;
        this.inQueue = new int[nClasses];

        for (int i = 0; i < nClasses; i++) {
            this.inQueue[i] = 0;
        }
        this.totalInQueue = 0;

        this.phaseList = phaseList;
    }

    public LCFSBuffer(int nClasses, int nServers) {
        this(nClasses, nServers, null);
    }

    public void addToBuffer(int classIdx) {
        if (this.totalInQueue < this.nServers) {
            this.inQueue[classIdx]++;
            this.totalInQueue++;

            this.serverQueue.addLast(classIdx);

            this.phaseList.addToService(classIdx);
            return;
        }

        this.inQueue[classIdx]++;
        this.totalInQueue++;

        this.deque.addFirst(classIdx);
    }

    public void addNToBuffer(int classIdx, int n) {
        for (int i = 0; i < n; i++) {
            this.addToBuffer(classIdx);
        }
    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {
        if (this.totalInQueue < this.nServers) {
            this.inQueue[classIdx]++;
            this.totalInQueue++;

            this.serverQueue.addLast(classIdx);

            this.phaseList.addToService(classIdx, position);
            return;
        }

        this.inQueue[classIdx]++;
        this.totalInQueue++;

        this.deque.addFirst(classIdx);
    }

    public int getInService(int classIdx) {
        Iterator<Integer> serverIterator = this.serverQueue.iterator();
        int acc = 0;
        while ((serverIterator.hasNext()) ) {
            if (serverIterator.next() == classIdx) {
                acc++;
            }
        }
        return acc;
    }

    public boolean isEmpty() {
        return this.totalInQueue == 0;
    }

    public void removeFirstOfClass(int classIdx) {
        Iterator<Integer> serverIterator = this.serverQueue.iterator();
        while (serverIterator.hasNext()) {
            if (serverIterator.next() == classIdx) {
                serverIterator.remove();
                this.totalInQueue--;
                this.inQueue[classIdx]--;

                this.moveToService();

                return;
            }
        }

        Iterator<Integer> dequeIterator = this.deque.iterator();
        while (dequeIterator.hasNext()) {
            if (dequeIterator.next() == classIdx) {
                dequeIterator.remove();
                this.totalInQueue--;
                this.inQueue[classIdx]--;
                return;
            }
        }

    }

    public void removeNClass(int n, int classIdx) {
        for (int i = 0; i < n; i++) {
            removeFirstOfClass(classIdx);
        }
    }

    public SSAStateCell createCopy() {
        LCFSBuffer copyBuffer = new LCFSBuffer(this.inQueue.length, this.nServers, this.phaseList.createCopy());
        copyBuffer.deque = new ArrayDeque<Integer>(this.deque);
        copyBuffer.serverQueue = new ArrayDeque<Integer>(this.serverQueue);
        copyBuffer.inQueue = Arrays.copyOf(this.inQueue, this.inQueue.length);
        copyBuffer.totalInQueue = this.totalInQueue;

        return copyBuffer;
    }

    public int getInQueue(int classIdx) {
        return this.inQueue[classIdx];
    }



    public boolean incrementPhase(int classIdx) {
        return this.phaseList.incrementPhase(classIdx, this.getInService(classIdx));
    }

    public boolean updatePhase(int classIdx, int startingPhase, int newPhase) {
        return this.phaseList.updatePhase(classIdx, startingPhase, newPhase);
    }

    public boolean updateGlobalPhase(int classIdx, int newPhase) {
        this.phaseList.updateGlobalPhase(classIdx, newPhase);
        return true;
    }

    public int incrementPhaseN(int n, int classIdx) {
        throw new RuntimeException("Not implemented");
    }

    public int getGlobalPhase(int classIdx) {
        return this.phaseList.getGlobalPhase(classIdx);
    }

    public PhaseList getPhaseList() {
        return this.phaseList;
    }

    public List<Integer> stateVector() {
        Stream<Integer> dequeArr = this.deque.stream();
        Stream<Integer> serverArr = this.serverQueue.stream();
        Stream<Integer> phaseArr = Arrays.stream(this.phaseList.getVector());
        return Stream.concat(Stream.concat(serverArr, dequeArr), phaseArr).collect(Collectors.toList());
    }
}
