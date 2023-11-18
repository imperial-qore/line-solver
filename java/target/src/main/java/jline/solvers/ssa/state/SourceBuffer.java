package jline.solvers.ssa.state;

import java.util.*;
import java.util.stream.Collectors;

public class SourceBuffer extends SSAStateCell {
    protected int nClasses;

    protected PhaseList phaseList;

    public SourceBuffer(int nClasses, PhaseList phaseList) {
        this.nClasses = nClasses;

        this.phaseList = phaseList;
    }

    public void addToBuffer(int classIdx) {

    }

    public void addNToBuffer(int classIdx, int n) {

    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {

    }

    public int getInService(int classIdx) {
        return 1;
    }

    public boolean isEmpty() {
        return false;
    }

    public void removeFirstOfClass(int classIdx) {
        this.phaseList.addToService(classIdx); // start afresh
    }

    public void removeNClass(int n, int classIdx) {
        // since this is a source, there should only be one job being created at each time.
        this.removeFirstOfClass(classIdx);
    }

    public SSAStateCell createCopy() {
        SourceBuffer copyBuffer = new SourceBuffer(this.nClasses, this.phaseList.createCopy());

        return copyBuffer;
    }

    public int getInQueue(int classIdx) {
        return 1;
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
        return Arrays.stream(this.phaseList.getVector()).collect(Collectors.toList());
    }
}
