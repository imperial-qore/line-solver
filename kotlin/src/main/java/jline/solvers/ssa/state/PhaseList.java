package jline.solvers.ssa.state;

import jline.lang.distributions.CumulativeDistribution;
import jline.util.Matrix;

import java.util.*;
import java.util.stream.Stream;

public class PhaseList {
    protected int[] nPhases;
    protected int nClasses;

    protected Random random;

    protected Integer[] nInPhase;
    protected int[] globalPhases;
    protected int[] phaseListStart;
    protected int[] totalInList;

    protected int phaseListOffset = 0;

    protected List<CumulativeDistribution<Integer>> startingPhaseProbabilities;
    protected Set<Integer> customStartClasses;

    public PhaseList(int[] nPhases, int nClasses, Random random) {
        this.nPhases = nPhases; // number of phases per each class
        this.nClasses = nClasses;

        this.random = random;

        this.globalPhases = new int[nClasses];

        this.phaseListStart = new int[nClasses];
        this.totalInList = new int[nClasses];

        int cumClassPhase = 0;

        for (int i = 0; i < nClasses; i++) {
            this.phaseListStart[i] = cumClassPhase;
            cumClassPhase += this.nPhases[i];

            this.totalInList[i] = 0;
        }

        this.nInPhase = new Integer[cumClassPhase];
        for (int i = 0; i < cumClassPhase; i++) {
            this.nInPhase[i] = 0;
        }
        this.phaseListOffset = 0;

        this.startingPhaseProbabilities = new ArrayList<CumulativeDistribution<Integer>>();
        for (int i = 0; i < nClasses; i++) {
            this.startingPhaseProbabilities.add(new CumulativeDistribution<Integer>(this.random));
            this.startingPhaseProbabilities.get(i).addElement(0, 1);
            for (int j = 0; j < nPhases[i]; j++) {
                this.startingPhaseProbabilities.get(i).addElement(j, 0);
            }
        }
        this.customStartClasses = new HashSet<Integer>();
    }

    public void setPhaseStart(int classIdx, Matrix classProbabilities) {
        this.startingPhaseProbabilities.set(classIdx, new CumulativeDistribution<Integer>(this.random));
        for (int i = 0; i < classProbabilities.length(); i++) {
            this.startingPhaseProbabilities.get(classIdx).addElement(i, classProbabilities.get(i));
        }
        this.customStartClasses.add(classIdx);
    }

    public int getPhaseStart(int classIdx) {
        if (!this.customStartClasses.contains(classIdx)) {
            return 0;
        }

        return this.startingPhaseProbabilities.get(classIdx).sample(random);
    }

    public void setPhaseVector(Integer[] nInPhase) {
        int prevOffset = this.phaseListOffset;
        this.phaseListOffset = nInPhase.length - this.nInPhase.length - this.phaseListOffset;
        this.nInPhase = nInPhase;
        for (int i = this.phaseListOffset; i < this.nInPhase.length; i++) {
            this.nInPhase[i] = 0;
        }

        for (int i = 0; i < nClasses; i++) {
            this.phaseListStart[i] += this.phaseListOffset - prevOffset;
        }
    }

    public int getNInPhase(int classIdx, int phaseIdx) {
        return this.nInPhase[this.phaseListStart[classIdx] + phaseIdx];
    }

    public int getNInClass(int classIdx) {
        int acc = 0;
        for (int i = 0; i < this.nPhases[classIdx]; i++) {
            acc += this.nInPhase[this.phaseListStart[classIdx] + i];
        }
        return acc;
    }

    public void addToService(int classIdx, int startingPhase) {
        if (this.nPhases[classIdx] == 0) {
            return;
        }
        int offset = this.phaseListStart[classIdx];

        this.totalInList[classIdx] += 1;
        this.nInPhase[offset+startingPhase] += 1;
    }

    public void addToService(int classIdx) {
        addToService(classIdx, this.getPhaseStart(classIdx));
    }

    public void addToServiceN(int classIdx, int n) {
        for (int i = 0; i < n; i++) {
            this.addToService(classIdx);
        }
    }

    public void addToServiceN(int classIdx, int startingPhase, int n) {
        if (this.nPhases[classIdx] == 0) {
            return;
        }
        int offset = this.phaseListStart[classIdx];

        this.totalInList[classIdx] += n;
        this.nInPhase[offset+startingPhase] += n;
    }

    public boolean incrementPhase(int classIdx, int nInService) {
        if (this.nPhases[classIdx] == 0) {
            return true;
        }
        int offset = this.phaseListStart[classIdx];

        if (nInService > this.totalInList[classIdx]) {
            // ensure buffer matches service
            int shortfall = nInService - this.totalInList[classIdx];
            this.nInPhase[offset] += shortfall;
            this.totalInList[classIdx] += shortfall;
        }

        int nSelected = this.random.nextInt(nInService);

        int nSeen = 0;
        int maxPhase = this.nPhases[classIdx];

        for (int i = 0; i < maxPhase; i++) {
            nSeen += this.nInPhase[offset + i];

            if (nSeen > nSelected) {
                if (i == (maxPhase-1)) {
                    this.nInPhase[offset + i] -= 1;
                    this.totalInList[classIdx] -= 1;
                    return true;
                }

                this.nInPhase[offset + i] -= 1;
                this.nInPhase[offset + i+1] += 1;
                return false;
            }
        }

        return false;
    }

    public int incrementPhaseN(int n, int classIdx, int nInService) {
        /*
            incrementPhaseN -
                Run n phase increments, assuming the number in service remains constant

                Returns: number of departures
         */
        int nDepartures = 0;
        for (int i = 0; i < n; i++) {
            if (this.incrementPhase(classIdx, nInService)) {
                nDepartures += 1;
            }

            if (this.nInPhase[classIdx] == 0) {
                return nDepartures;
            }
        }
        return nDepartures;
    }

    public boolean updatePhase(int classIdx, int startingPhase, int newPhase) {
        int offset = this.phaseListStart[classIdx];

        if (startingPhase != -1 && this.nInPhase[offset+startingPhase] == 0) {
            return false;
        }
        if(startingPhase != -1) {
            this.nInPhase[offset + startingPhase] -= 1;
        }
        if (newPhase == -1) {
            // absorbing departure phase
            return true;
        }

        this.nInPhase[offset+newPhase] += 1;
        return true;
    }

    public void updateGlobalPhase(int classIdx, int newPhase) {
        this.globalPhases[classIdx] = newPhase;
    }

    public int getGlobalPhase(int classIdx) {
        return this.globalPhases[classIdx];
    }

    public CumulativeDistribution<Integer> getStartingPhaseProbabilities(int classIdx) {
        return this.startingPhaseProbabilities.get(classIdx);
    }

    public PhaseList createCopy() {
        PhaseList outList = new PhaseList(this.nPhases, this.nClasses, this.random);
        outList.nInPhase = new Integer[this.nInPhase.length];
        for(int i=0;i< nInPhase.length;i++){
            Integer temp = this.nInPhase[i];
            outList.nInPhase[i]=temp;
        }
        outList.globalPhases = this.globalPhases.clone();

        List<CumulativeDistribution<Integer>> originalList = this.startingPhaseProbabilities;
        List<CumulativeDistribution<Integer>> clonedList = new ArrayList<>();

        for (CumulativeDistribution<Integer> item : originalList) {
            clonedList.add(item.clone());
        }

        outList.startingPhaseProbabilities = clonedList;
        outList.phaseListStart = this.phaseListStart.clone();
        outList.totalInList = this.totalInList.clone();

        return outList;
    }

    public Integer[] getVector() {
        return this.nInPhase;
    }

    public List<Integer> getArray() {
        return Arrays.asList(this.nInPhase);
    }

    public Stream<Integer> getStream() {
        return Arrays.stream(this.nInPhase);
    }

    public int getNPhases(int classIndex) {
        return this.nPhases[classIndex];
    }
}
