package jline.solvers.ssa.state;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class LCFSPRBuffer extends SSAStateCell {
    public Deque<Integer> deque;
    protected int nServers;
    protected int[] inQueue;
    protected int totalInQueue;

    protected PhaseList phaseList;

    public LCFSPRBuffer(int nClasses, int nServers, PhaseList phaseList) {
        this.deque = new ArrayDeque<Integer>();
        this.nServers = nServers;
        this.inQueue = new int[nClasses];

        for (int i = 0; i < nClasses; i++) {
            this.inQueue[i] = 0;
        }
        this.totalInQueue = 0;
        this.phaseList = phaseList;
    }

    public LCFSPRBuffer(int nClasses, int nServers) {
        this(nClasses, nServers, null);
    }

    public void addToBuffer(int classIdx) {
        this.inQueue[classIdx]++;
        this.totalInQueue++;

        this.deque.addFirst(classIdx);
    }

    public void addNToBuffer(int classIdx, int n) {
        this.inQueue[classIdx] += n;
        this.totalInQueue += n;

        for (int i = 0; i < n; i++) {
            this.deque.addFirst(classIdx);
        }
    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {
        this.inQueue[classIdx]++;
        this.totalInQueue++;

        this.deque.addFirst(classIdx);
    }

    public int getInService(int classIdx) {
        Iterator<Integer> dequeIterator = deque.iterator();
        int nCt = 0;
        int acc = 0;
        while ((dequeIterator.hasNext()) && (nCt < this.nServers)) {
            if (dequeIterator.next() == classIdx) {
                acc++;
            }
            nCt++;
        }
        return acc;
    }

    public boolean isEmpty() {
        return this.totalInQueue == 0;
    }

    public void removeFirstOfClass(int classIdx) {
        Iterator<Integer> dequeIterator = this.deque.iterator();
        while (dequeIterator.hasNext()) {
            if (dequeIterator.next() == classIdx) {
                dequeIterator.remove();
                this.totalInQueue--;
                this.inQueue[classIdx]--;

                break;
            }
        }

        // replace lost service.. (we assume we removed classIdx from service)
        dequeIterator = deque.iterator();
        for (int i = 0; i < this.nServers-1; i++) {
            if (!dequeIterator.hasNext()) {
                return;
            }
            dequeIterator.next();
        }

        if (!dequeIterator.hasNext()) {
            return;
        }

        this.phaseList.addToService(dequeIterator.next());
    }

    public void removeNClass(int n, int classIdx) {
        Iterator<Integer> dequeIterator = deque.iterator();
        int nRemoved = 0;

        while ((dequeIterator.hasNext()) && (n > 0)) {
            if (dequeIterator.next() == classIdx) {
                dequeIterator.remove();
                n--;
                nRemoved++;
            }
        }

        this.totalInQueue -= nRemoved;
        this.inQueue[classIdx] -= nRemoved;

        int serviceReplacement = Math.min(this.getInService(classIdx),nRemoved);
        dequeIterator = deque.iterator();
        for (int i = 0; i < this.nServers-serviceReplacement; i++) {
            if (!dequeIterator.hasNext()) {
                return;
            }
            dequeIterator.next();
        }

        for (int i = 0; i < serviceReplacement; i++) {
            if (!dequeIterator.hasNext()) {
                return;
            }

            this.phaseList.addToService(dequeIterator.next());
        }
    }

    public SSAStateCell createCopy() {
        LCFSPRBuffer copyBuffer = new LCFSPRBuffer(this.inQueue.length, this.nServers, this.phaseList.createCopy());
        copyBuffer.deque = new ArrayDeque<Integer>(this.deque);
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
        Stream<Integer> phaseArr = Arrays.stream(this.phaseList.getVector());
        return Stream.concat(dequeArr, phaseArr).collect(Collectors.toList());
    }
}
