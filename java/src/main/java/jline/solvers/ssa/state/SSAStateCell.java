package jline.solvers.ssa.state;

import java.util.List;

public abstract class SSAStateCell {

    public SSAStateCell() {
    }

//    public abstract int peakBuffer();
//    public abstract int peakBufferAt(int atPoint);
//    public abstract int popFromBuffer();
    public abstract void addToBuffer(int classIdx);
    public abstract void addNToBuffer(int classIdx, int n);
    public abstract void addToBufferAtPosition(int nodeIdx, int classIdx, int position);
    public abstract int getInService(int classIdx);
    public abstract boolean isEmpty();
    public abstract void removeFirstOfClass(int classIdx);
    public abstract void removeNClass(int n, int classIdx);
    public abstract SSAStateCell createCopy();
    public abstract int getInQueue(int classIdx);

    public abstract boolean incrementPhase(int classIdx);
    public abstract int     incrementPhaseN(int n, int classIdx);
    public abstract boolean updatePhase(int classIdx, int startingPhase, int newPhase);
    public abstract boolean updateGlobalPhase(int classIdx, int newPhase);
    public abstract int getGlobalPhase(int classIdx);
    public abstract PhaseList getPhaseList();
    public abstract List<Integer> stateVector();
}
