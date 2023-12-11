package jline.solvers.ssa.state;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;

import jline.util.Matrix;


import java.util.*;

public class SSAStateMatrix {
    /*
        In theory, this should be a one-stop point to handle all stateful information about the system.

        The system isn't quite there yet, e.g. some events track stateful info (e.g. JoinOutputEvent).
     */

    // configuration parameters
    public NetworkStruct sn;
    protected double cutoff = -1;

    // information on the state (jobs at each station, value at other StatefulNode objects, and phases)
    public int[][] state; // [node][class]
    protected SSAStateCell[] buffers;

    // caching, for TimeWarp
    protected int[][] stateCache;
    protected SSAStateCell[] bufferCache;

    protected Random random;

    public SSAStateMatrix(NetworkStruct sn, Random random) {
        this.sn = sn;
        this.random = random;

        this.state = new int[this.sn.nstateful][this.sn.nclasses];
        for (int i = 0; i < this.sn.nstateful; i++) {
            for (int j = 0; j < this.sn.nclasses; j++) {
                this.state[i][j] = 0;
            }
        }

        this.stateCache = new int[sn.nstateful][sn.nclasses];

        // build StateCell instances according to the scheduling strategy at each node.
        this.buffers = new SSAStateCell[sn.nstateful];
        this.bufferCache = new SSAStateCell[sn.nstateful];
        for (int i = 0; i < sn.nstations; i++) {
            int[] phases_r = new int[this.sn.nclasses];
            for (int r = 0; r < sn.nclasses; r++) {
                phases_r[r] = (int) this.sn.phases.get(i,r);
            }

            PhaseList phaseList = new PhaseList(phases_r, this.sn.nclasses, this.random);
            if (sn.pie.get(i) != null) {
                for (int j = 0; j < sn.nclasses; j++) {
                    if (sn.pie.get(i).containsKey(j)) {
                        phaseList.setPhaseStart(j, sn.pie.get(i).get(j));
                    }
                }
            }
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
                this.buffers[i] = new FCFSBuffer(sn.nclasses, (int) sn.nservers.get(i), phaseList);
            } else if(sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
                this.buffers[i] = new INFBuffer(sn.nclasses, phaseList);
            } else if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT) {
                this.buffers[i] = new SourceBuffer(sn.nclasses, phaseList);
            } else if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.LCFS) {
                this.buffers[i] = new LCFSBuffer(sn.nclasses, (int) sn.nservers.get(i), phaseList);
            } else if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.LCFSPR) {
                this.buffers[i] = new LCFSPRBuffer(sn.nclasses, (int) sn.nservers.get(i), phaseList);
            } else if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.PS) {
                this.buffers[i] = new PSBuffer(this.random, sn.nclasses, (int) sn.nservers.get(i), phaseList);
            } else if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.SIRO) {
                this.buffers[i] = new SIROBuffer(this.random, sn.nclasses, (int) sn.nservers.get(i), phaseList, false);
            /*} else if (sn.schedStrategies[i] == SchedStrategy.SIROPR) {
                this.buffers[i] = new SIROClassBuffer(this.random, sn.nclasses, sn.numberOfServers[i], phaseList, true);*/
            } else {
                System.out.println(sn.sched.get(sn.stations.get(i)));
                throw new RuntimeException("Unsupported Scheduling Strategy");
            }
        }
    }

    public SSAStateMatrix(SSAStateMatrix that) {
        this.sn = that.sn;
        this.sn.nstateful = that.sn.nstateful;
        this.sn.nclasses = that.sn.nclasses;
        this.cutoff = that.cutoff;
        this.state = new int[that.state.length][that.state[0].length];
        for(int i=0;i<state.length;i++){
            System.arraycopy(that.state[i], 0, this.state[i], 0, state[0].length);
        }

        this.random = that.random;
        this.buffers = new SSAStateCell[that.sn.nstateful];
        for (int i = 0 ; i < that.sn.nstateful; i++) {
            this.buffers[i] = that.buffers[i].createCopy();
        }
    }

    public void addToBuffer(int nodeIdx, int classIdx) {
        this.buffers[nodeIdx].addToBuffer(classIdx);
    }

    public void addToBuffer(int nodeIdx, int classIdx, int count) {
        this.buffers[nodeIdx].addNToBuffer(classIdx, count);
    }

    public void addToBufferAtPosition(int nodeIdx, int classIdx, int position) {
        this.buffers[nodeIdx].addToBufferAtPosition(nodeIdx, classIdx, position);
    }

    public boolean stateArrival(int nodeIdx, int classIdx) {
        // arrive 1 instance of [class] at [node]
        // returns: true if successful, false otherwise
        if (state[nodeIdx][classIdx] >= sn.classcap.get(nodeIdx,classIdx)) {
            return false;
        }

        this.addToBuffer(nodeIdx, classIdx);
        this.state[nodeIdx][classIdx]++;

        return true;
    }

    public int stateArrivalN(int n, int nodeIdx, int classIdx) {
        /*
            Try to arrive n instances of [class] at [node]
            returns: number of UNapplied arrivals (e.g. expect a 0 from this in normal cases).
         */
        int curState = this.state[nodeIdx][classIdx];
        int maxState = (int) this.sn.classcap.get(nodeIdx,classIdx);


        int nToApply = Math.min(n, maxState - curState);
        int rem = Math.min(n - nToApply, n);

        this.addToBuffer(nodeIdx, classIdx, nToApply);
        this.setState(nodeIdx, classIdx, this.getState(nodeIdx, classIdx) + nToApply);
        return rem;
    }

    public boolean stateArrivalAtPosition(int nodeIdx, int classIdx, int position) {
        // arrive 1 instance of [class] at [node]
        // returns: true if successful, false otherwise
        if (state[nodeIdx][classIdx] >= sn.classcap.get(nodeIdx,classIdx)) {
            return false;
        }

        this.addToBufferAtPosition(nodeIdx, classIdx, position);
        this.state[nodeIdx][classIdx]++;

        return true;
    }
    public boolean stateDeparture(int nodeIdx, int classIdx) {
        // depart 1 instance of [class] from [node]
        // returns: true if successful departure, false otherwise
        if ((state[nodeIdx][classIdx] == 0)/* && (!this.allowIllegalStates)*/) {
            return false;
        }

        this.buffers[nodeIdx].removeFirstOfClass(classIdx);
        this.state[nodeIdx][classIdx]--;
        return true;
    }

    public int stateDepartureN(int n, int nodeIdx, int classIdx) {
        // depart n instances of [class] from [node]
        // returns: number of UNapplied departures
        int curState = this.state[nodeIdx][classIdx];
        int nToApply = Math.min(curState, n);

        this.buffers[nodeIdx].removeNClass(nToApply, classIdx);
        this.state[nodeIdx][classIdx] -= nToApply;

        return n-nToApply;
    }

    public int totalStateAtNode(int nodeIdx) {
        // total jobs at a certain node
        int totalState = 0;
        for (int i = 0; i < sn.nclasses; i++) {
            totalState += this.state[nodeIdx][i];
        }

        return totalState;
    }

    public int getCapacity(int nodeIdx, int classIdx) {
        return (int) this.sn.classcap.get(nodeIdx,classIdx);
    }

    public boolean atEmpty(int nodeIdx, int classIdx) {
        return this.state[nodeIdx][classIdx] == 0;
    }

    public int getState(int nodeIdx, int classIdx) {
        return this.state[nodeIdx][classIdx];
    }

    public void setState(int nodeIdx, int classIdx, int state) {
        // mostly used for debugging
        this.state[nodeIdx][classIdx] = state;
    }

    public int inProcess(int nodeIdx, int classIdx) {
        if (atEmpty(nodeIdx, classIdx)) {
            return 0;
        }

        return this.buffers[nodeIdx].getInService(classIdx);
    }

    public int psTotalCapacity(int nodeIdx) {
        return ((PSBuffer)this.buffers[nodeIdx]).getTotalCapacity();
    }

    public boolean incrementPhase(int nodeIdx, int classIdx) {
        /*
            Signal a class-specific phase update
         */
        return this.buffers[nodeIdx].incrementPhase(classIdx);
    }

    public int incrementPhaseN(int n, int nodeIdx, int classIdx) {
        return this.buffers[nodeIdx].incrementPhaseN(n, classIdx);
    }

    public boolean updatePhase (int nodeIdx, int classIdx, int startingPhase, int endingPhase) {
        return this.buffers[nodeIdx].updatePhase(classIdx, startingPhase, endingPhase);
    }

    public boolean updateGlobalPhase(int nodeIdx, int classIdx, int newPhase) {
        /*
            Signal a global phase update
         */
        return this.buffers[nodeIdx].updateGlobalPhase(classIdx, newPhase);
    }

    public int getGlobalPhase(int nodeIdx, int classIdx) {
        return this.buffers[nodeIdx].getGlobalPhase(classIdx);
    }


    public int getInPhase(int nodeIdx, int classIdx, int phase) {
        return this.buffers[nodeIdx].getPhaseList().getNInPhase(classIdx, phase);
    }


    public void cacheState() {
        /*
            Create caches for each StateCell and state
         */
        this.stateCache = new int[this.sn.nstateful][this.sn.nclasses];
        this.bufferCache = new SSAStateCell[this.sn.nstateful];
        for (int i = 0; i < this.sn.nstateful; i++) {
            if (this.sn.nclasses >= 0) System.arraycopy(this.state[i], 0, this.stateCache[i], 0, this.sn.nclasses);
            this.bufferCache[i] = this.buffers[i].createCopy();
        }
    }

    public void revertToCache() {
        this.state = this.stateCache;
        this.buffers = this.bufferCache;
    }
    @SuppressWarnings("unchecked")
    public List<Integer>[] getStateVectors() {
        /*
            Return state vectors for transient analysis:

            Ext: [Inf, s11, ... S1K,...sR1, ...SRk]
            FCFS, HOL, LCFS: [cb,...c1, s11, ... S1K,...sR1, ...SRk]
         */
        List<Integer>[] outList = new List[this.sn.nstateful];

        for (int i = 0; i < this.sn.nstateful; i++) {
            outList[i] = this.buffers[i].stateVector();
        }

        return outList;
    }

    public static boolean sameState(List<Integer>[] state1, List<Integer>[] state2){
        assert state1.length==state2.length;
        for(int i=0;i<state1.length;i++){
            if(state1[i].size()!=state2[i].size()){
                return false;
            }
            for(int j = 0;j<state1[i].size();j++){
                if(state1[i].get(j)!=state2[i].get(j)){
                    return false;
                }
            }
        }
        return true;
    }

    public void printStateVector(){
        List<Integer>[] state = this.getStateVectors();
        for (List<Integer> list : state){
            for (int i : list){
                System.out.print(i + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public Random getRandom() {
        return this.random;
    }

    @Override
    public boolean equals(Object o){
        SSAStateMatrix that = (SSAStateMatrix) o;
        List<Integer>[] state1 = this.getStateVectors();
        List<Integer>[] state2 = that.getStateVectors();
        assert state1.length==state2.length;
        for(int i=0;i<state1.length;i++){
            if(state1[i].size()!=state2[i].size()){
                return false;
            }
            for(int j = 0;j<state1[i].size();j++){
                if(!Objects.equals(state1[i].get(j), state2[i].get(j))){
                    return false;
                }
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        return Arrays.deepHashCode(this.getStateVectors());
    }


    public int findPhaseChange(SSAStateMatrix newNetworkState, int statefulIndex, int classIndex) {
        PhaseList phaseList1 = this.buffers[statefulIndex].getPhaseList();
        PhaseList phaseList2 = newNetworkState.buffers[statefulIndex].getPhaseList();
        int n = phaseList1.getNPhases(classIndex);
        for(int i = 0; i < n; i++) {
            if(phaseList1.getNInPhase(classIndex, i) != phaseList2.getNInPhase(classIndex, i)) {
                return i;
            }
        }
        return -1;
    }

    public PhaseList getPhaseList(int statefulIndex) {
        return this.buffers[statefulIndex].getPhaseList();
    }

    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
    }

    public boolean exceedsCutoff() {
        if(cutoff == -1) {
            return false;
        }
        for(int k = 0; k < sn.nclasses; k++) {
            double sum = 0;
            for (int i = 0; i < sn.nstateful; i++) {
                sum += state[i][k];
            }
            if (sum > cutoff) {
                return true;
            }
        }
        return false;
    }
}
