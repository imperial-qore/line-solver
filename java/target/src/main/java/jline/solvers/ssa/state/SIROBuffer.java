package jline.solvers.ssa.state;

import jline.lang.distributions.CumulativeDistribution;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SIROBuffer extends SSAStateCell {
    protected int[] inWaiting;
    protected int[] inService;
    protected Integer[] inQueue;

    protected int nServers;
    protected int nClasses;
    protected int totalInQueue;

    protected boolean isPreemptive;

    protected PhaseList phaseList;

    protected Random random;

    private void moveToService() {
        CumulativeDistribution<Integer> nextClassCumulativeDistribution = new CumulativeDistribution<Integer>(this.random);

        int totalInWaiting = 0;
        for (int i = 0; i < this.nClasses; i++) {
            nextClassCumulativeDistribution.addElement(i, this.inWaiting[i]);
            totalInWaiting += this.inWaiting[i];
        }

        if (totalInWaiting == 0) {
            return;
        }

        nextClassCumulativeDistribution.normalize(totalInWaiting);


        int classToAdd = nextClassCumulativeDistribution.sample(random);
        this.phaseList.addToService(classToAdd);
        this.inService[classToAdd] += 1;
        this.inWaiting[classToAdd] -= 1;
    }

    public SIROBuffer(Random random, int nClasses, int nServers, PhaseList phaseList, boolean isPreemptive) {
        this.inWaiting = new int[nClasses];
        this.inService = new int[nClasses];
        this.inQueue = new Integer[nClasses];

        this.nServers = nServers;
        this.nClasses = nClasses;

        for (int i = 0; i < nClasses; i++) {
            this.inWaiting[i] = 0;
            this.inService[i] = 0;
            this.inQueue[i] = 0;
        }

        this.totalInQueue = 0;

        this.phaseList = phaseList;
        this.isPreemptive = isPreemptive;

        this.random = random;
    }

    public SIROBuffer(Random random, int nClasses, int nServers, PhaseList phaseList) {
        this(random, nClasses, nServers, phaseList, false);
    }

    public void addToBuffer(int classIdx) {
        if (this.totalInQueue < this.nServers) {
            this.totalInQueue++;

            this.inService[classIdx] += 1;
            this.inQueue[classIdx] += 1;

            this.phaseList.addToService(classIdx);
            return;
        }

        if (this.isPreemptive) {
            double preemptionChance = ((double) this.nServers)/((double)(this.totalInQueue+1));
            if (this.random.nextDouble() <= preemptionChance) {
                this.totalInQueue++;

                this.inService[classIdx] += 1;
                this.inQueue[classIdx] += 1;

                this.phaseList.addToService(classIdx);
                return;
            }
        }

        this.totalInQueue++;

        this.inWaiting[classIdx] += 1;
        this.inQueue[classIdx] += 1;
    }

    public void addNToBuffer(int classIdx, int n) {
        int addToService = Math.min(Math.max(this.nServers-this.totalInQueue,0),n);
        int addToWaiting = n - addToService;
        this.inService[classIdx] += addToService;
        this.phaseList.addToServiceN(classIdx, 0, n);
        this.inWaiting[classIdx] += addToWaiting;
        this.totalInQueue += n;
        this.inQueue[classIdx] += n;
    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {
        if (this.totalInQueue < this.nServers) {
            this.totalInQueue++;

            this.inService[classIdx] += 1;
            this.inQueue[classIdx] += 1;

            this.phaseList.addToService(classIdx, position);
            return;
        }

        if (this.isPreemptive) {
            double preemptionChance = ((double) this.nServers)/((double)(this.totalInQueue+1));
            if (this.random.nextDouble() <= preemptionChance) {
                this.totalInQueue++;

                this.inService[classIdx] += 1;
                this.inQueue[classIdx] += 1;

                this.phaseList.addToService(classIdx, position);
                return;
            }
        }

        this.totalInQueue++;

        this.inWaiting[classIdx] += 1;
        this.inQueue[classIdx] += 1;
    }

    public int getInService(int classIdx) {
        return this.inService[classIdx];
    }

    public boolean isEmpty() {
        return this.totalInQueue == 0;
    }

    public void removeFirstOfClass(int classIdx) {
        if (this.inService[classIdx] > 0) {
            this.inService[classIdx] -= 1;
            this.inQueue[classIdx] -= 1;
            this.moveToService();
            this.totalInQueue -= 1;
        } else {
            this.inWaiting[classIdx] -= 1;
            this.inQueue[classIdx] -= 1;
            this.moveToService();
            this.totalInQueue -= 1;
            //System.out.println("failed departure");
        }
    }

    public void removeNClass(int n, int classIdx) {
        for (int i = 0; i < n; i++) {
            removeFirstOfClass(classIdx);
        }
    }

    public SSAStateCell createCopy() {
        SIROBuffer copyBuffer = new SIROBuffer(this.random, this.nClasses, this.nServers, this.phaseList.createCopy());
        copyBuffer.inWaiting = Arrays.copyOf(this.inWaiting, this.nClasses);
        copyBuffer.inService = Arrays.copyOf(this.inService, this.nClasses);
        copyBuffer.inQueue = Arrays.copyOf(this.inQueue, this.nClasses);
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
        return this.phaseList.incrementPhaseN(n, classIdx, this.getInService(classIdx));
    }

    public int getGlobalPhase(int classIdx) {
        return this.phaseList.getGlobalPhase(classIdx);
    }

    public PhaseList getPhaseList() {
        return this.phaseList;
    }

    public List<Integer> stateVector() {
        return Stream.concat(Arrays.stream(this.inQueue), this.phaseList.getStream()).collect(Collectors.toList());
    }
}
