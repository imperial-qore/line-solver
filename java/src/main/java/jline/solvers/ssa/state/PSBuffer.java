package jline.solvers.ssa.state;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.lang.Math.min;

public class PSBuffer extends SSAStateCell {
    protected Integer[] classCounts;
    protected Random random;
    protected int nServers;
    protected int nClasses;

    protected PhaseList phaseList;

    protected int totalCt;

    public PSBuffer(Random random, int nClasses, int nServers, PhaseList phaseList) {
        this.random = random;
        this.nServers = nServers;
        this.nClasses = nClasses;

        this.phaseList = phaseList;

        this.classCounts = new Integer[this.nClasses];

        for (int i = 0; i < this.nClasses; i++) {
            this.classCounts[i] = 0;
        }

        this.totalCt = 0;
    }

    public void addToBuffer(int classIdx) {
        this.phaseList.addToService(classIdx);
        this.classCounts[classIdx] += 1;
        this.totalCt += 1;
    }

    public void addNToBuffer(int classIdx, int n) {
        this.phaseList.addToServiceN(classIdx, 0, n);
        this.classCounts[classIdx] += n;
        this.totalCt += n;
    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {
        this.phaseList.addToService(classIdx, position);
        this.classCounts[classIdx] += 1;
        this.totalCt += 1;
    }

    public int getInService(int classIdx) {
        double val = (this.nServers*((double)this.classCounts[classIdx]/(double)this.totalCt));
        int base = (int)Math.floor(val);
        if (this.random.nextDouble() <= val-base) {
            return base + 1;
        }
        return base;
    }

    public boolean isEmpty() {
        for (int i = 0; i < this.nClasses; i++) {
            if (this.classCounts[i] != 0) {
                return false;
            }
        }
        return true;
    }

    public void removeFirstOfClass(int classIdx) {
        if (this.classCounts[classIdx] == 0) {
            return;
        }
        this.classCounts[classIdx] -= 1;
        this.totalCt -= 1;
    }

    public void removeNClass(int n, int classIdx) {
        if (this.classCounts[classIdx] < n) {
            this.classCounts[classIdx] = 0;
        }

        this.classCounts[classIdx] -= n;
        this.totalCt -= n;
    }

    public SSAStateCell createCopy() {
        PSBuffer copyBuffer = new PSBuffer(this.random, this.nClasses, this.nServers, this.phaseList.createCopy());
        copyBuffer.classCounts = Arrays.copyOf(this.classCounts, this.nClasses);
        copyBuffer.nClasses = this.nClasses;
        copyBuffer.nServers = this.nServers;
        copyBuffer.random = this.random;
        copyBuffer.totalCt = this.totalCt;
        return copyBuffer;
    }

    public int getInQueue(int classIdx) {
        return this.classCounts[classIdx];
    }


    public boolean incrementPhase(int classIdx) {
        // everything in PS is in the queue.
        return this.phaseList.incrementPhase(classIdx, this.getInQueue(classIdx));
    }

    public boolean updatePhase(int classIdx, int startingPhase, int newPhase) {
        return this.phaseList.updatePhase(classIdx, startingPhase, newPhase);
    }

    public boolean updateGlobalPhase(int classIdx, int newPhase) {
        this.phaseList.updateGlobalPhase(classIdx, newPhase);
        return true;
    }

    public int incrementPhaseN(int n, int classIdx) {

        return this.phaseList.incrementPhaseN(n, classIdx, this.getInQueue(classIdx));
    }

    public int getGlobalPhase(int classIdx) {
        return this.phaseList.getGlobalPhase(classIdx);
    }

    public PhaseList getPhaseList() {
        return this.phaseList;
    }

    public List<Integer> stateVector() {
        return Stream.concat(Arrays.stream(this.classCounts), this.phaseList.getStream()).collect(Collectors.toList());
    }

    public int getTotalCapacity() {
        int acc = 0;
        for (int val : this.classCounts) {
            acc += val;
            if (acc >= this.nServers) {
                return this.nServers;
            }
        }
        return acc;
    }
}
