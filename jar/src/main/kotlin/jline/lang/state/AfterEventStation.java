package jline.lang.state;

import static jline.GlobalConstants.NegInf;
import static jline.io.InputOutputKt.line_warning;

import jline.io.InputOutputKt;
import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class AfterEventStation implements Serializable {
    static Ret.EventResult afterEventStation(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Double>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        // Declare variables at method level to avoid scope issues
        Matrix sir = null;
        List<Matrix> kir = null;
        // Declare variables outside switch to avoid scope issues
        Matrix ni = null;
        Matrix nir = null;
        
        switch (event) {
            case ARV:
                // return if no space to accept the arrival, otherwise check scheduling strategy
                State.StateMarginalStatistics stats = ToMarginal.toMarginalAggr(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
                ni = stats.ni;
                nir = stats.nir;
                sir = stats.sir;
                kir = stats.kir;
                Matrix pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
                // For Place nodes (INF scheduling with NaN service), use uniform entry probability
                if (pentry != null && pentry.hasNaN()) {
                    boolean allNaN = true;
                    for (int i = 0; i < pentry.length(); i++) {
                        if (!Double.isNaN(pentry.get(i))) {
                            allNaN = false;
                            break;
                        }
                    }
                    if (allNaN && pentry.length() > 0) {
                        pentry = new Matrix(pentry.getNumRows(), pentry.getNumCols());
                        pentry.fill(1.0 / pentry.length());
                    }
                }
                outprob = new Matrix(0, 0);
                Matrix outprobK = new Matrix(0, 0);
                k_loop:
                for (int kentry = 0; kentry < K.get(jobClass); kentry++) {
                    Matrix spaceVarK = spaceVar.copy();
                    Matrix spaceSrvK = spaceSrv.copy();
                    Matrix spaceBufK = spaceBuf.copy();
                    switch (sn.sched.get(sn.stations.get(ist))) {
                        case EXT: // source, can receive any "virtual" arrival from the sink as long as it is from an open class
                            if (Utils.isInf(sn.njobs.get(jobClass))) {
                                outspace = inspace.copy();
                                outrate = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                                outrate.zero();
                                outprob = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                                outprob.ones();
                                break k_loop; // must leave switch and for loop to go straight to sim code
                            }
                            break;
                        case INF:
                        case PS:
                        case DPS:
                        case GPS:
                        case LPS:
                        case PSPRIO:
                        case DPSPRIO:
                        case GPSPRIO:
                            // due to nature of these policies, a new job enters service immediately.
                            boolean exceedsClassCap = false;
                            double istCap = classcap.get(ist, jobClass);
                            // classcap = 0 means no per-class constraint (MATLAB convention)
                            // Only enforce limit if classcap > 0
                            if (istCap > 0) {
                                int col = (int) (Ks.get(jobClass) + kentry);
                                for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                    if (spaceSrvK.get(row, col) >= istCap) {
                                        exceedsClassCap = true;
                                    }
                                }
                            }
                            if (!exceedsClassCap) {
                                // increment spacesrvk by one
                                int col = (int) (Ks.get(jobClass) + kentry);
                                for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                    spaceSrvK.set(row, col, spaceSrvK.get(row, col) + 1);
                                }
                                outprobK = new Matrix(spaceSrvK.getNumRows(), spaceSrvK.getNumRows());
                                outprobK.fill(pentry.get(kentry));
                            } else {
                                outprobK = new Matrix(spaceSrvK.getNumRows(), spaceSrvK.getNumRows());
                                outprobK.zero();
                            }
                            break;
                        case SIRO:
                        case SEPT:
                        case LEPT:
                            outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                            for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                if (ni.get(row, 0) < S.get(ist)) {
                                    // Idle servers available - job enters service
                                    int colSrvK = (int) (Ks.get(jobClass) + kentry);
                                    spaceSrvK.set(row, colSrvK, spaceSrvK.get(row, colSrvK) + 1);
                                } else {
                                    // All servers busy - job goes to buffer
                                    spaceBufK.set(row, jobClass, spaceBufK.get(row, jobClass) + 1);
                                }
                            }
                            // Set output probability for the current entry phase
                            outprobK.fill(pentry.get(kentry));
                            break;
                        case FCFS:
                        case HOL:
                        case LCFS:
                            // find states with all servers busy
                            // if MAP service, when empty restart from the phase stored in spaceVar for this class
                            // MATLAB: sn.nvars(ind,1:class) extracts columns 1 to class (1-based)
                            // Java: columns 0 to jobClass (0-based), same logical data
                            // Extract from col 0 to col jobClass+1 (exclusive end)

                            Matrix nVarsAtInd = Matrix.extract(sn.nvars, ind, ind + 1, 0, jobClass + 1);
                            int sum = (int) nVarsAtInd.elementSum();
                            // MATLAB: kentry == space_var(sum(sn.nvars(ind,1:class)))
                            // MATLAB uses 1-based indexing, so sum gives position 1-based
                            // Java needs 0-based column index: sum - 1
                            // For multi-row spaceVar (RROBIN), check if condition holds for any row
                            int spaceVarCol = sum - 1;
                            boolean mapPhaseMatch = false;
                            if (ismkvmodclass.get(jobClass) == 0) {
                                mapPhaseMatch = true;
                            } else if (spaceVarCol >= 0 && spaceVarCol < spaceVar.getNumCols()) {
                                // Check if any row matches the kentry phase
                                // kentry is 0-based in Java, but MAP output var values
                                // in the state space are 1-based (initDefault sets to 1),
                                // so compare kentry+1 with spaceVar value
                                for (int row = 0; row < spaceVar.getNumRows(); row++) {
                                    if ((kentry + 1) == spaceVar.get(row, spaceVarCol)) {
                                        mapPhaseMatch = true;
                                        break;
                                    }
                                }
                            }
                            if (mapPhaseMatch) {
                                if (ismkvmodclass.get(jobClass) == 1) {
                                    pentry.zero();
                                    pentry.set(kentry, 1);
                                }
                                // construct all_busy_srv, a matrix with 0s where sum of that row in space_srv_k is >= S.get(ist) and 1s where its <
                                Matrix all_busy_srv = new Matrix(spaceSrvK.getNumRows(), 1);
                                for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                    Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                    int rowSum = (int) row.elementSum();
                                    if (rowSum >= S.get(ist)) {
                                        all_busy_srv.set(i, 0, 1);
                                    } else {
                                        all_busy_srv.set(i, 0, 0);
                                    }
                                }

                                // find and modify states with an idle server
                                Matrix idle_srv = new Matrix(spaceSrvK.getNumRows(), 1);
                                for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                    Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                    int rowSum = (int) row.elementSum();
                                    if (rowSum < S.get(ist)) {
                                        idle_srv.set(i, 0, 1);
                                    } else {
                                        idle_srv.set(i, 0, 0);
                                    }
                                }

                                // job enters service, increments idle server terms in space_srv_k
                                int colSrvK = (int) (spaceSrvK.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                                for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                    if (idle_srv.get(row, 0) == 1) {
                                        spaceSrvK.set(row, colSrvK, spaceSrvK.get(row, colSrvK) + 1);
                                    }
                                }
                                // this section dynamically grows the number of elements in the buffer
                                //  ni is an Mx1
                                // check if there is space in buffer
                                boolean spaceInBuffer = false;
                                for (int i = 0; i < ni.getNumRows(); i++) {
                                    if (ni.get(i, 0) < capacity.get(ist)) {
                                        spaceInBuffer = true;
                                    }
                                }
                                if (spaceInBuffer) {
                                    boolean spaceForClass = false;
                                    double classCapLimit = classcap.get(ist, jobClass);
                                    // classcap = 0 means no per-class constraint (MATLAB convention)
                                    // If classCapLimit == 0, there's always space for the class
                                    if (classCapLimit == 0) {
                                        spaceForClass = true;
                                    } else {
                                        for (int i = 0; i < nir.getNumRows(); i++) {
                                            if (nir.get(i, jobClass) < classCapLimit) {
                                                spaceForClass = true;
                                            }
                                        }
                                    }
                                    if (spaceForClass) {
                                        // there is room for the job, check if buffer has empty slots. If not, append job slot
                                        boolean emptySlots = false;
                                        for (int i = 0; i < spaceBufK.getNumRows(); i++) {
                                            if (spaceBufK.get(i, 0) == 0) {
                                                emptySlots = true;
                                            }
                                        }
                                        if (!emptySlots) {
                                            // append job slot
                                            Matrix left = new Matrix(spaceBufK.getNumRows(), 1);
                                            left.zero();
                                            spaceBufK = Matrix.concatColumns(left, spaceBufK, null);
                                        }
                                    }
                                }
                                // get position of first empty slot in buffer
                                Matrix empty_slots = new Matrix(all_busy_srv.getNumRows(), 1);
                                empty_slots.fill(-1);
                                int spaceBufKCols = spaceBufK.getNumCols();
                                if (spaceBufKCols == 0) {
                                    //set empty_slots(all_busy_srv) = false;
                                    for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                                        if (all_busy_srv.get(i, 0) == 1) {
                                            empty_slots.set(i, 0, 0);
                                        }
                                    }
                                } else if (spaceBufKCols == 1) {
                                    // set empty_slots(all_busy_srv) = space_buf_k(all_busy_srv,:)==0;
                                    for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                                        if (all_busy_srv.get(i, 0) == 1) {
                                            if (spaceBufK.get(i, 0) == 0) {
                                                empty_slots.set(i, 0, 1);
                                            } else {
                                                empty_slots.set(i, 0, 0);
                                            }
                                        }
                                    }

                                } else {
                                    Matrix spaceBufKAllBusySrv = new Matrix(0, 0);
                                    for (int i = 0; i < spaceBufK.getNumRows(); i++) {
                                        if (all_busy_srv.get(i, 0) == 1) {
                                            if (spaceBufKAllBusySrv.isEmpty()) {
                                                spaceBufKAllBusySrv = Matrix.extractRows(spaceBufK, i, i + 1, null);
                                            } else {
                                                spaceBufKAllBusySrv = Matrix.concatRows(spaceBufKAllBusySrv, Matrix.extractRows(spaceBufK, i, i + 1, null), null);
                                            }
                                        }
                                    }

                                    for (int row = 0; row < spaceBufKAllBusySrv.getNumRows(); row++) {
                                        for (int c = 0; c < spaceBufKAllBusySrv.getNumCols(); c++) {
                                            if (spaceBufKAllBusySrv.get(row, c) == 0) {
                                                spaceBufKAllBusySrv.set(row, c, 1);
                                            } else {
                                                spaceBufKAllBusySrv.set(row, c, 0);
                                            }
                                        }
                                    }
                                    Matrix sizeSpaceBufK = new Matrix(1, spaceBufK.getNumCols());
                                    for (int i = 0; i < sizeSpaceBufK.getNumCols(); i++) {
                                        sizeSpaceBufK.set(0, i, i + 1);
                                    }
                                    if (!spaceBufKAllBusySrv.isEmpty() && !sizeSpaceBufK.isEmpty()) {
                                        Matrix elementMult = spaceBufKAllBusySrv.elementMult(sizeSpaceBufK, null);
                                        Matrix max = new Matrix(elementMult.getNumRows(), 1);
                                        for (int i = 0; i < elementMult.getNumRows(); i++) {
                                            max.set(i, 0, Matrix.extractRows(elementMult, i, i + 1, null).elementMax());
                                        }
                                        int max_ind = 0;
                                        for (int i = 0; i < all_busy_srv.getNumRows(); i++) {
                                            if (all_busy_srv.get(i, 0) == 1) {
                                                empty_slots.set(i, 0, max.get(max_ind));
                                                max_ind++;
                                            }
                                        }
                                    }
                                }
                                // ignore states where buffer has no empty slots
                                // set wbuf_empty = empty_slots >0
                                Matrix wbuf_empty = new Matrix(empty_slots.getNumRows(), 1);
                                boolean space_available = false;
                                for (int i = 0; i < empty_slots.getNumRows(); i++) {
                                    if (empty_slots.get(i, 0) > 0) {
                                        wbuf_empty.set(i, 0, 1);
                                        space_available = true;
                                    } else {
                                        wbuf_empty.set(i, 0, 0);
                                    }
                                }
                                if (space_available) {
                                    // space_srv_k set to the rows of s[ace_srv_k where wbuf_empty = 1
                                    Matrix spaceSrvKTmp = new Matrix(0, 0);
                                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                                        if (wbuf_empty.get(i, 0) == 1) {
                                            if (spaceSrvKTmp.isEmpty()) {
                                                spaceSrvKTmp = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                            } else {
                                                spaceSrvKTmp = Matrix.concatRows(spaceSrvKTmp, Matrix.extractRows(spaceSrvK, i, i + 1, null), null);
                                            }
                                        }
                                    }
                                    spaceSrvK = spaceSrvKTmp;
                                    Matrix spaceBufKTmp = new Matrix(0, 0);
                                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                                        if (wbuf_empty.get(i, 0) == 1) {
                                            if (spaceBufKTmp.isEmpty()) {
                                                spaceBufKTmp = Matrix.extractRows(spaceBufK, i, i + 1, null);
                                            } else {
                                                spaceBufKTmp = Matrix.concatRows(spaceBufKTmp, Matrix.extractRows(spaceBufK, i, i + 1, null), null);
                                            }
                                        }
                                    }
                                    spaceBufK = spaceBufKTmp;
                                    Matrix spaceVarKTmp = new Matrix(0, 0);
                                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                                        if (wbuf_empty.get(i, 0) == 1) {
                                            if (spaceVarKTmp.isEmpty()) {
                                                spaceVarKTmp = Matrix.extractRows(spaceVarK, i, i + 1, null);
                                            } else {
                                                spaceVarKTmp = Matrix.concatRows(spaceVarKTmp, Matrix.extractRows(spaceVarK, i, i + 1, null), null);
                                            }
                                        }
                                    }
                                    spaceVarK = spaceVarKTmp;
                                    Matrix emptySlotsTmp = new Matrix(0, 0);
                                    for (int i = 0; i < wbuf_empty.getNumRows(); i++) {
                                        if (wbuf_empty.get(i, 0) == 1) {
                                            if (emptySlotsTmp.isEmpty()) {
                                                emptySlotsTmp = Matrix.extractRows(empty_slots, i, i + 1, null);
                                            } else {
                                                emptySlotsTmp = Matrix.concatRows(emptySlotsTmp, Matrix.extractRows(empty_slots, i, i + 1, null), null);
                                            }
                                        }
                                    }
                                    empty_slots = emptySlotsTmp;
                                    Matrix dims = new Matrix(1, 2);
                                    dims.set(0, 0, spaceBufK.getNumRows());
                                    dims.set(0, 1, spaceBufK.getNumCols());

                                    Matrix row_indices = new Matrix(1, spaceBufK.getNumRows());
                                    for (int i = 0; i < row_indices.getNumCols(); i++) {
                                        row_indices.set(0, i, i);
                                    }

                                    // need col indices to be = empty slots, but decrement as matrix 0 indexed
                                    Matrix col_indices = empty_slots.copy();
                                    for (int r = 0; r < col_indices.getNumRows(); r++) {
                                        for (int c = 0; c < col_indices.getNumCols(); c++) {
                                            col_indices.set(r, c, col_indices.get(r, c) - 1);
                                        }
                                    }
                                    col_indices = col_indices.transpose();
                                    List<Integer> indices = Maths.sub2ind(dims, row_indices, col_indices);
                                    for (Integer n : indices) {
                                        // jobClass + 1 since final form needs jobs to be 1 indexed
                                        spaceBufK.set(n, jobClass + 1);
                                    }

                                }
                                outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                                outprobK.fill(pentry.get(kentry));
                            } else {
                                outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                                outprobK.zero(); // zero probability event
                            }
                            break;
                        case LCFSPR: // LCFS with Preemption
                            // find states with all servers busy
                            Matrix allBusySrv = new Matrix(spaceSrvK.getNumRows(), 1);
                            Matrix idleSrv = new Matrix(spaceSrvK.getNumRows(), 1);

                            for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                int rowSum = (int) row.elementSum();
                                if (rowSum >= S.get(ist)) {
                                    allBusySrv.set(i, 0, 1);
                                    idleSrv.set(i, 0, 0);
                                } else {
                                    allBusySrv.set(i, 0, 0);
                                    idleSrv.set(i, 0, 1);
                                }
                            }

                            // Reorder states so that idle ones come first
                            Matrix spaceBufKReordLcfspr = new Matrix(0, 0);
                            Matrix spaceSrvKReordLcfspr = new Matrix(0, 0);
                            Matrix spaceVarKReordLcfspr = new Matrix(0, 0);

                            // Add idle states first
                            for (int row = 0; row < idleSrv.getNumRows(); row++) {
                                if (idleSrv.get(row, 0) == 1) {
                                    if (spaceBufKReordLcfspr.isEmpty()) {
                                        spaceBufKReordLcfspr = Matrix.extractRows(spaceBufK, row, row + 1, null);
                                        spaceSrvKReordLcfspr = Matrix.extractRows(spaceSrvK, row, row + 1, null);
                                        spaceVarKReordLcfspr = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                    } else {
                                        spaceBufKReordLcfspr = Matrix.concatRows(spaceBufKReordLcfspr, Matrix.extractRows(spaceBufK, row, row + 1, null), null);
                                        spaceSrvKReordLcfspr = Matrix.concatRows(spaceSrvKReordLcfspr, Matrix.extractRows(spaceSrvK, row, row + 1, null), null);
                                        spaceVarKReordLcfspr = Matrix.concatRows(spaceVarKReordLcfspr, Matrix.extractRows(spaceVarK, row, row + 1, null), null);
                                    }
                                }
                            }

                            // If idle, the job enters service in phase kentry
                            boolean anyIdle = false;
                            for (int row = 0; row < idleSrv.getNumRows(); row++) {
                                if (idleSrv.get(row, 0) == 1) {
                                    anyIdle = true;
                                    break;
                                }
                            }

                            if (anyIdle && !spaceSrvKReordLcfspr.isEmpty()) {
                                int colSrvK = (int) (spaceSrvKReordLcfspr.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                                for (int row = 0; row < spaceSrvKReordLcfspr.getNumRows(); row++) {
                                    spaceSrvKReordLcfspr.set(row, colSrvK, spaceSrvKReordLcfspr.get(row, colSrvK) + 1);
                                }
                            }

                            // If all busy, expand output states for all possible choices of job class to preempt
                            Matrix psentryLcfspr = new Matrix(spaceBufKReordLcfspr.getNumRows(), 1);
                            psentryLcfspr.ones(); // probability scaling due to preemption

                            for (int classpreempt = 0; classpreempt < R; classpreempt++) {
                                for (int phasepreempt = 0; phasepreempt < K.get(classpreempt); phasepreempt++) {
                                    // Check if there are jobs of this class/phase to preempt
                                    Matrix siPreempt = new Matrix(spaceSrvK.getNumRows(), 1);
                                    for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                        int colPreempt = (int) (spaceSrvK.getNumCols() - K.elementSum() + Ks.get(classpreempt) + phasepreempt);
                                        siPreempt.set(row, 0, spaceSrvK.get(row, colPreempt));
                                    }

                                    Matrix busyPreempt = new Matrix(spaceSrvK.getNumRows(), 1);
                                    boolean anyBusyPreempt = false;
                                    for (int row = 0; row < siPreempt.getNumRows(); row++) {
                                        if (siPreempt.get(row, 0) > 0) {
                                            busyPreempt.set(row, 0, 1);
                                            anyBusyPreempt = true;
                                        } else {
                                            busyPreempt.set(row, 0, 0);
                                        }
                                    }

                                    if (anyBusyPreempt) {
                                        // Update probability scaling - gather all busy states first
                                        Matrix busyPreemptProbs = new Matrix(0, 1);
                                        // Calculate row sums for ALL rows in spaceSrvK (MATLAB: sum(space_srv_k,2))
                                        Matrix allRowSums = new Matrix(spaceSrvK.getNumRows(), 1);
                                        for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                            Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                            allRowSums.set(i, 0, row.elementSum());
                                        }
                                        
                                        for (int row = 0; row < busyPreempt.getNumRows(); row++) {
                                            if (busyPreempt.get(row, 0) == 1) {
                                                double siPreemptVal = siPreempt.get(row, 0);
                                                double rowSum = allRowSums.get(row, 0);
                                                Matrix newPsentry = new Matrix(1, 1);
                                                newPsentry.set(0, 0, siPreemptVal / rowSum);
                                                busyPreemptProbs = Matrix.concatRows(busyPreemptProbs, newPsentry, null);
                                            }
                                        }
                                        psentryLcfspr = Matrix.concatRows(psentryLcfspr, busyPreemptProbs, null);

                                        // Create preempted states
                                        Matrix spaceSrvKPreempt = new Matrix(0, 0);
                                        Matrix spaceBufKPreempt = new Matrix(0, 0);
                                        Matrix spaceVarKPreempt = new Matrix(0, 0);

                                        for (int row = 0; row < busyPreempt.getNumRows(); row++) {
                                            if (busyPreempt.get(row, 0) == 1) {
                                                if (spaceSrvKPreempt.isEmpty()) {
                                                    spaceSrvKPreempt = Matrix.extractRows(spaceSrvK, row, row + 1, null);
                                                    spaceBufKPreempt = Matrix.extractRows(spaceBufK, row, row + 1, null);
                                                    spaceVarKPreempt = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                                } else {
                                                    spaceSrvKPreempt = Matrix.concatRows(spaceSrvKPreempt, Matrix.extractRows(spaceSrvK, row, row + 1, null), null);
                                                    spaceBufKPreempt = Matrix.concatRows(spaceBufKPreempt, Matrix.extractRows(spaceBufK, row, row + 1, null), null);
                                                    spaceVarKPreempt = Matrix.concatRows(spaceVarKPreempt, Matrix.extractRows(spaceVarK, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (!spaceSrvKPreempt.isEmpty()) {
                                            // Remove preempted job
                                            int colPreempt = (int) (spaceSrvKPreempt.getNumCols() - K.elementSum() + Ks.get(classpreempt) + phasepreempt);
                                            for (int row = 0; row < spaceSrvKPreempt.getNumRows(); row++) {
                                                spaceSrvKPreempt.set(row, colPreempt, spaceSrvKPreempt.get(row, colPreempt) - 1);
                                            }

                                            // Add new job to service
                                            int colNew = (int) (spaceSrvKPreempt.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                                            for (int row = 0; row < spaceSrvKPreempt.getNumRows(); row++) {
                                                spaceSrvKPreempt.set(row, colNew, spaceSrvKPreempt.get(row, colNew) + 1);
                                            }

                                            // Add preempted job to buffer with (class, phase) pairs
                                            // Check if buffer needs expansion - match MATLAB logic
                                            if (isSimulation) {
                                                // Check if there's room and no empty slots
                                                boolean needsExpansion = true;
                                                if (spaceBufKPreempt.getNumCols() > 0) {
                                                    for (int row = 0; row < spaceBufKPreempt.getNumRows(); row++) {
                                                        for (int colBuf = 0; colBuf < spaceBufKPreempt.getNumCols(); colBuf++) {
                                                            if (spaceBufKPreempt.get(row, colBuf) == 0) {
                                                                needsExpansion = false;
                                                                break;
                                                            }
                                                        }
                                                        if (!needsExpansion) break;
                                                    }
                                                }

                                                if (needsExpansion) {
                                                    // Append two columns for (class, phase) pair - prepend like MATLAB
                                                    Matrix expansion = new Matrix(spaceBufKPreempt.getNumRows(), 2);
                                                    expansion.zero();
                                                    spaceBufKPreempt = Matrix.concatColumns(expansion, spaceBufKPreempt, null);
                                                }
                                            }

                                            // Find position for first empty slot - match MATLAB logic exactly
                                            Matrix emptySlots = new Matrix(spaceBufKPreempt.getNumRows(), 1);
                                            emptySlots.fill(-1);

                                            if (spaceBufKPreempt.getNumCols() == 0) {
                                                // No buffer space  
                                                emptySlots.zero();
                                            } else if (spaceBufKPreempt.getNumCols() == 2) {
                                                // Only one pair slot - check if first column (class) is empty
                                                for (int row = 0; row < spaceBufKPreempt.getNumRows(); row++) {
                                                    if (spaceBufKPreempt.get(row, 0) == 0) {
                                                        emptySlots.set(row, 0, 1); // Position 1 in MATLAB indexing
                                                    } else {
                                                        emptySlots.set(row, 0, 0); // No empty slot
                                                    }
                                                }
                                            } else {
                                                // Multiple pair slots - find first empty pair using MATLAB logic
                                                for (int row = 0; row < spaceBufKPreempt.getNumRows(); row++) {
                                                    int maxPos = -1;
                                                    for (int colPair = 0; colPair < spaceBufKPreempt.getNumCols(); colPair++) {
                                                        if (spaceBufKPreempt.get(row, colPair) == 0) {
                                                            maxPos = Math.max(maxPos, colPair + 1); // 1-based indexing
                                                        }
                                                    }
                                                    // Subtract 1 for (class, preempt-phase) pairs like MATLAB
                                                    emptySlots.set(row, 0, maxPos - 1);
                                                }
                                            }

                                            // Filter states where buffer has empty slots
                                            Matrix wbuf_empty = new Matrix(emptySlots.getNumRows(), 1);
                                            for (int i = 0; i < emptySlots.getNumRows(); i++) {
                                                if (emptySlots.get(i, 0) > 0) {
                                                    wbuf_empty.set(i, 0, 1);
                                                } else {
                                                    wbuf_empty.set(i, 0, 0);
                                                }
                                            }
                                            
                                            // Only process states with empty buffer slots
                                            if (wbuf_empty.elementSum() > 0) {
                                                // Filter matrices to only include states with empty buffer slots
                                                Matrix spaceSrvKPreemptFiltered = new Matrix(0, 0);
                                                Matrix spaceBufKPreemptFiltered = new Matrix(0, 0);
                                                Matrix spaceVarKPreemptFiltered = new Matrix(0, 0);
                                                Matrix emptySlotsFiltered = new Matrix(0, 1);
                                                
                                                for (int row = 0; row < wbuf_empty.getNumRows(); row++) {
                                                    if (wbuf_empty.get(row, 0) == 1) {
                                                        if (spaceSrvKPreemptFiltered.isEmpty()) {
                                                            spaceSrvKPreemptFiltered = Matrix.extractRows(spaceSrvKPreempt, row, row + 1, null);
                                                            spaceBufKPreemptFiltered = Matrix.extractRows(spaceBufKPreempt, row, row + 1, null);
                                                            spaceVarKPreemptFiltered = Matrix.extractRows(spaceVarKPreempt, row, row + 1, null);
                                                            emptySlotsFiltered = Matrix.extractRows(emptySlots, row, row + 1, null);
                                                        } else {
                                                            spaceSrvKPreemptFiltered = Matrix.concatRows(spaceSrvKPreemptFiltered, Matrix.extractRows(spaceSrvKPreempt, row, row + 1, null), null);
                                                            spaceBufKPreemptFiltered = Matrix.concatRows(spaceBufKPreemptFiltered, Matrix.extractRows(spaceBufKPreempt, row, row + 1, null), null);
                                                            spaceVarKPreemptFiltered = Matrix.concatRows(spaceVarKPreemptFiltered, Matrix.extractRows(spaceVarKPreempt, row, row + 1, null), null);
                                                            emptySlotsFiltered = Matrix.concatRows(emptySlotsFiltered, Matrix.extractRows(emptySlots, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                
                                                // Store preempted job in buffer - match MATLAB sub2ind logic
                                                for (int row = 0; row < spaceBufKPreemptFiltered.getNumRows(); row++) {
                                                    int emptySlot = (int) emptySlotsFiltered.get(row, 0);
                                                    if (emptySlot > 0) { // MATLAB uses 1-based indexing
                                                        // Convert back to 0-based for Java
                                                        int zeroBasedSlot = emptySlot - 1;
                                                        spaceBufKPreemptFiltered.set(row, zeroBasedSlot, classpreempt + 1); // Store class (1-based)
                                                        spaceBufKPreemptFiltered.set(row, zeroBasedSlot + 1, phasepreempt + 1); // Store phase (1-based)
                                                    }
                                                }
                                                
                                                // Use filtered matrices
                                                spaceSrvKPreempt = spaceSrvKPreemptFiltered;
                                                spaceBufKPreempt = spaceBufKPreemptFiltered;
                                                spaceVarKPreempt = spaceVarKPreemptFiltered;
                                            }

                                            // Add to reordered states
                                            spaceBufKReordLcfspr = Matrix.concatRows(spaceBufKReordLcfspr, spaceBufKPreempt, null);
                                            spaceSrvKReordLcfspr = Matrix.concatRows(spaceSrvKReordLcfspr, spaceSrvKPreempt, null);
                                            spaceVarKReordLcfspr = Matrix.concatRows(spaceVarKReordLcfspr, spaceVarKPreempt, null);
                                        }
                                    }
                                }
                            }

                            // Set final output matrices
                            spaceBufK = spaceBufKReordLcfspr;
                            spaceSrvK = spaceSrvKReordLcfspr;
                            spaceVarK = spaceVarKReordLcfspr;

                            // Update probability
                            if (!psentryLcfspr.isEmpty()) {
                                for (int row = 0; row < psentryLcfspr.getNumRows(); row++) {
                                    if (outprobK.isEmpty()) {
                                        outprobK = new Matrix(1, 1);
                                        outprobK.set(0, 0, pentry.get(kentry) * psentryLcfspr.get(row, 0));
                                    } else {
                                        Matrix newProb = new Matrix(1, 1);
                                        newProb.set(0, 0, pentry.get(kentry) * psentryLcfspr.get(row, 0));
                                        outprobK = Matrix.concatRows(outprobK, newProb, null);
                                    }
                                }
                            } else {
                                outprobK = new Matrix(1, 1);
                                outprobK.set(0, 0, pentry.get(kentry));
                            }
                            break;
                        case LCFSPI: // LCFS with Preemption Independent (restart from phase 1)
                            // find states with all servers busy
                            Matrix allBusySrvPI = new Matrix(spaceSrvK.getNumRows(), 1);
                            Matrix idleSrvPI = new Matrix(spaceSrvK.getNumRows(), 1);

                            for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                int rowSum = (int) row.elementSum();
                                if (rowSum >= S.get(ist)) {
                                    allBusySrvPI.set(i, 0, 1);
                                    idleSrvPI.set(i, 0, 0);
                                } else {
                                    allBusySrvPI.set(i, 0, 0);
                                    idleSrvPI.set(i, 0, 1);
                                }
                            }

                            // Reorder states so that idle ones come first
                            Matrix spaceBufKReordLcfspi = new Matrix(0, 0);
                            Matrix spaceSrvKReordLcfspi = new Matrix(0, 0);
                            Matrix spaceVarKReordLcfspi = new Matrix(0, 0);

                            // Add idle states first
                            for (int row = 0; row < idleSrvPI.getNumRows(); row++) {
                                if (idleSrvPI.get(row, 0) == 1) {
                                    if (spaceBufKReordLcfspi.isEmpty()) {
                                        spaceBufKReordLcfspi = Matrix.extractRows(spaceBufK, row, row + 1, null);
                                        spaceSrvKReordLcfspi = Matrix.extractRows(spaceSrvK, row, row + 1, null);
                                        spaceVarKReordLcfspi = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                    } else {
                                        spaceBufKReordLcfspi = Matrix.concatRows(spaceBufKReordLcfspi, Matrix.extractRows(spaceBufK, row, row + 1, null), null);
                                        spaceSrvKReordLcfspi = Matrix.concatRows(spaceSrvKReordLcfspi, Matrix.extractRows(spaceSrvK, row, row + 1, null), null);
                                        spaceVarKReordLcfspi = Matrix.concatRows(spaceVarKReordLcfspi, Matrix.extractRows(spaceVarK, row, row + 1, null), null);
                                    }
                                }
                            }

                            // If idle, the job enters service in phase kentry
                            boolean anyIdlePI = false;
                            for (int row = 0; row < idleSrvPI.getNumRows(); row++) {
                                if (idleSrvPI.get(row, 0) == 1) {
                                    anyIdlePI = true;
                                    break;
                                }
                            }

                            if (anyIdlePI && !spaceSrvKReordLcfspi.isEmpty()) {
                                int colSrvK = (int) (spaceSrvKReordLcfspi.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                                for (int row = 0; row < spaceSrvKReordLcfspi.getNumRows(); row++) {
                                    spaceSrvKReordLcfspi.set(row, colSrvK, spaceSrvKReordLcfspi.get(row, colSrvK) + 1);
                                }
                            }

                            // If all busy, expand output states for all possible choices of job class to preempt
                            Matrix psentryLcfspi = new Matrix(spaceBufKReordLcfspi.getNumRows(), 1);
                            psentryLcfspi.ones(); // probability scaling due to preemption

                            for (int classpreempt = 0; classpreempt < R; classpreempt++) {
                                for (int phasepreempt = 0; phasepreempt < K.get(classpreempt); phasepreempt++) {
                                    // Check if there are jobs of this class/phase to preempt
                                    Matrix siPreemptPI = new Matrix(spaceSrvK.getNumRows(), 1);
                                    for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                        int colPreempt = (int) (spaceSrvK.getNumCols() - K.elementSum() + Ks.get(classpreempt) + phasepreempt);
                                        siPreemptPI.set(row, 0, spaceSrvK.get(row, colPreempt));
                                    }

                                    Matrix busyPreemptPI = new Matrix(spaceSrvK.getNumRows(), 1);
                                    boolean anyBusyPreemptPI = false;
                                    for (int row = 0; row < siPreemptPI.getNumRows(); row++) {
                                        if (siPreemptPI.get(row, 0) > 0) {
                                            busyPreemptPI.set(row, 0, 1);
                                            anyBusyPreemptPI = true;
                                        } else {
                                            busyPreemptPI.set(row, 0, 0);
                                        }
                                    }

                                    if (anyBusyPreemptPI) {
                                        // Update probability scaling - gather all busy states first
                                        Matrix busyPreemptProbsPI = new Matrix(0, 1);
                                        // Calculate row sums for ALL rows in spaceSrvK (MATLAB: sum(space_srv_k,2))
                                        Matrix allRowSumsPI = new Matrix(spaceSrvK.getNumRows(), 1);
                                        for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                            Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                            allRowSumsPI.set(i, 0, row.elementSum());
                                        }
                                        
                                        for (int row = 0; row < busyPreemptPI.getNumRows(); row++) {
                                            if (busyPreemptPI.get(row, 0) == 1) {
                                                double siPreemptVal = siPreemptPI.get(row, 0);
                                                double rowSum = allRowSumsPI.get(row, 0);
                                                Matrix newPsentry = new Matrix(1, 1);
                                                newPsentry.set(0, 0, siPreemptVal / rowSum);
                                                busyPreemptProbsPI = Matrix.concatRows(busyPreemptProbsPI, newPsentry, null);
                                            }
                                        }
                                        psentryLcfspi = Matrix.concatRows(psentryLcfspi, busyPreemptProbsPI, null);

                                        // Create preempted states
                                        Matrix spaceSrvKPreemptPI = new Matrix(0, 0);
                                        Matrix spaceBufKPreemptPI = new Matrix(0, 0);
                                        Matrix spaceVarKPreemptPI = new Matrix(0, 0);

                                        for (int row = 0; row < busyPreemptPI.getNumRows(); row++) {
                                            if (busyPreemptPI.get(row, 0) == 1) {
                                                if (spaceSrvKPreemptPI.isEmpty()) {
                                                    spaceSrvKPreemptPI = Matrix.extractRows(spaceSrvK, row, row + 1, null);
                                                    spaceBufKPreemptPI = Matrix.extractRows(spaceBufK, row, row + 1, null);
                                                    spaceVarKPreemptPI = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                                } else {
                                                    spaceSrvKPreemptPI = Matrix.concatRows(spaceSrvKPreemptPI, Matrix.extractRows(spaceSrvK, row, row + 1, null), null);
                                                    spaceBufKPreemptPI = Matrix.concatRows(spaceBufKPreemptPI, Matrix.extractRows(spaceBufK, row, row + 1, null), null);
                                                    spaceVarKPreemptPI = Matrix.concatRows(spaceVarKPreemptPI, Matrix.extractRows(spaceVarK, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (!spaceSrvKPreemptPI.isEmpty()) {
                                            // Remove preempted job
                                            int colPreempt = (int) (spaceSrvKPreemptPI.getNumCols() - K.elementSum() + Ks.get(classpreempt) + phasepreempt);
                                            for (int row = 0; row < spaceSrvKPreemptPI.getNumRows(); row++) {
                                                spaceSrvKPreemptPI.set(row, colPreempt, spaceSrvKPreemptPI.get(row, colPreempt) - 1);
                                            }

                                            // Add new job to service
                                            int colNew = (int) (spaceSrvKPreemptPI.getNumCols() - K.elementSum() + Ks.get(jobClass) + kentry);
                                            for (int row = 0; row < spaceSrvKPreemptPI.getNumRows(); row++) {
                                                spaceSrvKPreemptPI.set(row, colNew, spaceSrvKPreemptPI.get(row, colNew) + 1);
                                            }

                                            // Add preempted job to buffer with (class, phase) pairs
                                            // Check if buffer needs expansion - match MATLAB logic
                                            if (isSimulation) {
                                                // Check if there's room and no empty slots
                                                boolean needsExpansion = true;
                                                if (spaceBufKPreemptPI.getNumCols() > 0) {
                                                    for (int row = 0; row < spaceBufKPreemptPI.getNumRows(); row++) {
                                                        for (int colBuf = 0; colBuf < spaceBufKPreemptPI.getNumCols(); colBuf++) {
                                                            if (spaceBufKPreemptPI.get(row, colBuf) == 0) {
                                                                needsExpansion = false;
                                                                break;
                                                            }
                                                        }
                                                        if (!needsExpansion) break;
                                                    }
                                                }

                                                if (needsExpansion) {
                                                    // Append two columns for (class, phase) pair - prepend like MATLAB
                                                    Matrix expansion = new Matrix(spaceBufKPreemptPI.getNumRows(), 2);
                                                    expansion.zero();
                                                    spaceBufKPreemptPI = Matrix.concatColumns(expansion, spaceBufKPreemptPI, null);
                                                }
                                            }

                                            // Find position for first empty slot - match MATLAB logic exactly
                                            Matrix emptySlotsPI = new Matrix(spaceBufKPreemptPI.getNumRows(), 1);
                                            emptySlotsPI.fill(-1);

                                            if (spaceBufKPreemptPI.getNumCols() == 0) {
                                                // No buffer space  
                                                emptySlotsPI.zero();
                                            } else if (spaceBufKPreemptPI.getNumCols() == 2) {
                                                // Only one pair slot - check if first column (class) is empty
                                                for (int row = 0; row < spaceBufKPreemptPI.getNumRows(); row++) {
                                                    if (spaceBufKPreemptPI.get(row, 0) == 0) {
                                                        emptySlotsPI.set(row, 0, 1); // Position 1 in MATLAB indexing
                                                    } else {
                                                        emptySlotsPI.set(row, 0, 0); // No empty slot
                                                    }
                                                }
                                            } else {
                                                // Multiple pair slots - find first empty pair using MATLAB logic
                                                for (int row = 0; row < spaceBufKPreemptPI.getNumRows(); row++) {
                                                    int maxPos = -1;
                                                    for (int colPair = 0; colPair < spaceBufKPreemptPI.getNumCols(); colPair++) {
                                                        if (spaceBufKPreemptPI.get(row, colPair) == 0) {
                                                            maxPos = Math.max(maxPos, colPair + 1); // 1-based indexing
                                                        }
                                                    }
                                                    // Subtract 1 for (class, preempt-phase) pairs like MATLAB
                                                    emptySlotsPI.set(row, 0, maxPos - 1);
                                                }
                                            }

                                            // Filter states where buffer has empty slots
                                            Matrix wbuf_emptyPI = new Matrix(emptySlotsPI.getNumRows(), 1);
                                            for (int i = 0; i < emptySlotsPI.getNumRows(); i++) {
                                                if (emptySlotsPI.get(i, 0) > 0) {
                                                    wbuf_emptyPI.set(i, 0, 1);
                                                } else {
                                                    wbuf_emptyPI.set(i, 0, 0);
                                                }
                                            }
                                            
                                            // Only process states with empty buffer slots
                                            if (wbuf_emptyPI.elementSum() > 0) {
                                                // Filter matrices to only include states with empty buffer slots
                                                Matrix spaceSrvKPreemptFilteredPI = new Matrix(0, 0);
                                                Matrix spaceBufKPreemptFilteredPI = new Matrix(0, 0);
                                                Matrix spaceVarKPreemptFilteredPI = new Matrix(0, 0);
                                                Matrix emptySlotsFilteredPI = new Matrix(0, 1);
                                                
                                                for (int row = 0; row < wbuf_emptyPI.getNumRows(); row++) {
                                                    if (wbuf_emptyPI.get(row, 0) == 1) {
                                                        if (spaceSrvKPreemptFilteredPI.isEmpty()) {
                                                            spaceSrvKPreemptFilteredPI = Matrix.extractRows(spaceSrvKPreemptPI, row, row + 1, null);
                                                            spaceBufKPreemptFilteredPI = Matrix.extractRows(spaceBufKPreemptPI, row, row + 1, null);
                                                            spaceVarKPreemptFilteredPI = Matrix.extractRows(spaceVarKPreemptPI, row, row + 1, null);
                                                            emptySlotsFilteredPI = Matrix.extractRows(emptySlotsPI, row, row + 1, null);
                                                        } else {
                                                            spaceSrvKPreemptFilteredPI = Matrix.concatRows(spaceSrvKPreemptFilteredPI, Matrix.extractRows(spaceSrvKPreemptPI, row, row + 1, null), null);
                                                            spaceBufKPreemptFilteredPI = Matrix.concatRows(spaceBufKPreemptFilteredPI, Matrix.extractRows(spaceBufKPreemptPI, row, row + 1, null), null);
                                                            spaceVarKPreemptFilteredPI = Matrix.concatRows(spaceVarKPreemptFilteredPI, Matrix.extractRows(spaceVarKPreemptPI, row, row + 1, null), null);
                                                            emptySlotsFilteredPI = Matrix.concatRows(emptySlotsFilteredPI, Matrix.extractRows(emptySlotsPI, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                
                                                // Store preempted job in buffer - LCFSPI: use phase 1 placeholder (actual restart uses pie distribution)
                                                for (int row = 0; row < spaceBufKPreemptFilteredPI.getNumRows(); row++) {
                                                    int emptySlot = (int) emptySlotsFilteredPI.get(row, 0);
                                                    if (emptySlot > 0) { // MATLAB uses 1-based indexing
                                                        // Convert back to 0-based for Java
                                                        int zeroBasedSlot = emptySlot - 1;
                                                        spaceBufKPreemptFilteredPI.set(row, zeroBasedSlot, classpreempt + 1); // Store class (1-based)
                                                        // For LCFSPI: store phase 1 as placeholder, actual restart will use pie distribution
                                                        spaceBufKPreemptFilteredPI.set(row, zeroBasedSlot + 1, 1); // Phase placeholder for LCFSPI
                                                    }
                                                }
                                                
                                                // Use filtered matrices
                                                spaceSrvKPreemptPI = spaceSrvKPreemptFilteredPI;
                                                spaceBufKPreemptPI = spaceBufKPreemptFilteredPI;
                                                spaceVarKPreemptPI = spaceVarKPreemptFilteredPI;
                                            }

                                            // Add to reordered states
                                            spaceBufKReordLcfspi = Matrix.concatRows(spaceBufKReordLcfspi, spaceBufKPreemptPI, null);
                                            spaceSrvKReordLcfspi = Matrix.concatRows(spaceSrvKReordLcfspi, spaceSrvKPreemptPI, null);
                                            spaceVarKReordLcfspi = Matrix.concatRows(spaceVarKReordLcfspi, spaceVarKPreemptPI, null);
                                        }
                                    }
                                }
                            }

                            // Set final output matrices
                            spaceBufK = spaceBufKReordLcfspi;
                            spaceSrvK = spaceSrvKReordLcfspi;
                            spaceVarK = spaceVarKReordLcfspi;

                            // Update probability
                            if (!psentryLcfspi.isEmpty()) {
                                for (int row = 0; row < psentryLcfspi.getNumRows(); row++) {
                                    if (outprobK.isEmpty()) {
                                        outprobK = new Matrix(1, 1);
                                        outprobK.set(0, 0, pentry.get(kentry) * psentryLcfspi.get(row, 0));
                                    } else {
                                        Matrix newProb = new Matrix(1, 1);
                                        newProb.set(0, 0, pentry.get(kentry) * psentryLcfspi.get(row, 0));
                                        outprobK = Matrix.concatRows(outprobK, newProb, null);
                                    }
                                }
                            } else {
                                outprobK = new Matrix(1, 1);
                                outprobK.set(0, 0, pentry.get(kentry));
                            }
                            break;
                    }
                    // form the new state
                    Matrix outspaceKTmp = Matrix.concatColumns(spaceBufK, spaceSrvK, null);
                    Matrix outspaceK = Matrix.concatColumns(outspaceKTmp, spaceVarK, null);
                    // remove states where new arrival violates capacity or cutoff constraints
                    State.StateMarginalStatistics oi_oir = ToMarginal.toMarginalAggr(sn, ind, outspaceK, K, Ks, spaceBufK, spaceSrvK, spaceVarK);
                    Matrix oi = oi_oir.ni;
                    Matrix oir = oi_oir.nir;

                    Matrix en_o = new Matrix(oi.getNumRows(), 1);
                    for (int row = 0; row < oi.getNumRows(); row++) {
                        Matrix m = new Matrix(oi.getNumRows(), 1);
                        m.fill(capacity.get(ist));
                        boolean violates = false;
                        for (int col = 0; col < oi.getNumCols(); col++) {
                            if (m.get(row, 0) < oi.get(row, col)) {
                                violates = true;
                            }
                        }

                        double classCapLimit = classcap.get(ist, jobClass);
                        // classcap = 0 means no per-class constraint (MATLAB convention)
                        // If classCapLimit == 0, always pass the capacity check
                        boolean passesClassCapCheck = (classCapLimit == 0) || (classCapLimit >= oir.get(row, jobClass));
                        if (passesClassCapCheck && !violates) {
                            en_o.set(row, 0, 1);
                        }
                    }

                    // need to extract all rows of outspace_k where en_o is true
                    Matrix outspace_k_en_o = new Matrix(0, 0);
                    for (int row = 0; row < en_o.getNumRows(); row++) {
                        if (en_o.get(row, 0) == 1) {
                            if (outspace_k_en_o.isEmpty()) {
                                outspace_k_en_o = Matrix.extractRows(outspaceK, row, row + 1, null);
                            } else {
                                outspace_k_en_o = Matrix.concatRows(outspace_k_en_o, Matrix.extractRows(outspaceK, row, row + 1, null), null);
                            }
                        }
                    }

                    Matrix outprob_k_en_o = new Matrix(0, 0);
                    for (int row = 0; row < en_o.getNumRows(); row++) {
                        if (en_o.get(row, 0) == 1) {
                            if (outprob_k_en_o.isEmpty()) {
                                outprob_k_en_o = Matrix.extractRows(outprobK, row, row + 1, null);
                            } else {
                                outprob_k_en_o = Matrix.concatRows(outprob_k_en_o, Matrix.extractRows(outprobK, row, row + 1, null), null);
                            }
                        }
                    }

                    // Skip append when en_o filtered out all rows (matches MATLAB behavior where
                    // concatenating empty 0-row matrices is a no-op)
                    if (!outspace_k_en_o.isEmpty()) {
                        if (outspace.getNumCols() > outspace_k_en_o.getNumCols()) {
                            Matrix zeros = new Matrix(outspace_k_en_o.getNumRows(), outspace.getNumCols() - outspace_k_en_o.getNumCols());
                            zeros.zero();
                            Matrix bottom = Matrix.concatColumns(zeros, outspace_k_en_o, null);
                            outspace = Matrix.concatRows(outspace, bottom, null);
                        } else if (outspace.getNumCols() < outspace_k_en_o.getNumCols()) {
                            Matrix zeros = new Matrix(outspace.getNumRows(), outspace_k_en_o.getNumCols() - outspace.getNumCols());
                            zeros.zero();
                            Matrix top = Matrix.concatColumns(zeros, outspace, null);
                            outspace = Matrix.concatRows(top, outspace_k_en_o, null);
                        } else {
                            outspace = Matrix.concatRows(outspace, outspace_k_en_o, null);
                        }
                        Matrix newRates = new Matrix(outspace_k_en_o.getNumRows(), 1);
                        newRates.fill(-1);
                        outrate = Matrix.concatRows(outrate, newRates, null);
                        outprob = Matrix.concatRows(outprob, outprob_k_en_o, null);
                    }
                }

                // ARV results are NOT cached (matching MATLAB behavior which only caches DEP and PHASE)
                if (isSimulation) {
                    if (outprob.getNumRows() > 1) {
                        Matrix cum_sum = outprob.cumsumViaCol();
                        Matrix sum_by_col = outprob.sumCols();
                        Matrix cum_prob = Matrix.scaleMult(cum_sum, 1.0 / sum_by_col.value());

                        int firing_ctr = -1;
                        double rand = Maths.rand();
                        // we need the indicies where rand is bigger than cum_prob
                        for (int row = 0; row < cum_prob.getNumRows(); row++) {
                            if (rand > cum_prob.get(row, 0)) {
                                firing_ctr = row;
                            }
                        }
                        firing_ctr++;
                        outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                        outrate = new Matrix(1, 1);
                        outrate.set(0, 0, -1);
                        outprob = new Matrix(1, 1);
                        outprob.set(0, 0, 1);
                    }
                }
                break;
            case DEP:
                boolean busy = false;
                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                    for (int col = (int) Ks.get(jobClass); col < (Ks.get(jobClass) + K.get(jobClass)); col++) {
                        if (spaceSrv.get(row, col) > 0) {
                            busy = true;
                        }
                    }
                }
                if (busy) {
                    SchedStrategy strategy = sn.sched.get(sn.stations.get(ist));
                    sir = new Matrix(0, 0);
                    kir = new ArrayList<Matrix>();
                    if (hasOnlyExp && (strategy == SchedStrategy.PS || strategy == SchedStrategy.INF || strategy == SchedStrategy.DPS || strategy == SchedStrategy.GPS)) {
                        nir = spaceSrv.copy();
                        // set ni to sum of nir row-wise
                        ni = nir.sumRows();
                        sir = nir.copy();
                        kir.add(sir.copy());
                    } else {
                        State.StateMarginalStatistics dep_stats = ToMarginal.toMarginal(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
                        ni = dep_stats.ni;
                        nir = dep_stats.nir;
                        sir = dep_stats.sir;
                        kir = dep_stats.kir;
                    }

                    if (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass)) == RoutingStrategy.RROBIN) {
                        // Implement round-robin routing state update
                        // MATLAB: sn.nvars(ind,1:(R+class)) extracts columns 1 to R+class (1-based)
                        // Java: columns 0 to R+jobClass (0-based), same logical data
                        // Extract from col 0 to col R+jobClass+1 (exclusive end)
                        int nvarCols = R + jobClass + 1;
                        Matrix nvar_ind = new Matrix(1, nvarCols);
                        Matrix.extract(sn.nvars, ind, ind + 1, 0, nvarCols, nvar_ind, 0, 0);
                        int nvar_sum = (int) nvar_ind.elementSum();

                        // Get outlinks for this node and job class
                        Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(jobClass));

                        // MATLAB uses 1-based indexing: space_var(sum(...)) accesses element at position sum
                        // Java uses 0-based indexing: need to subtract 1 from nvar_sum to get column index
                        int spaceVarCol = nvar_sum - 1;

                        // Bounds check before accessing spaceVar
                        if (spaceVarCol >= 0 && spaceVarCol < spaceVar.getNumCols()) {
                            // Update RROBIN state for each row in spaceVar
                            // Each input state may have a different current RROBIN pointer
                            // outlinks is a row vector (1×N), use length() not getNumRows()
                            int numOutlinks = (int) outlinks.length();
                            for (int stateRow = 0; stateRow < spaceVar.getNumRows(); stateRow++) {
                                // Find current outlink index position for this state
                                int idx = -1;
                                double currentOutlink = spaceVar.get(stateRow, spaceVarCol);
                                for (int outlinkIdx = 0; outlinkIdx < numOutlinks; outlinkIdx++) {
                                    if (currentOutlink == outlinks.get(outlinkIdx)) {
                                        idx = outlinkIdx;
                                        break;
                                    }
                                }

                                // Update to next outlink (with wraparound)
                                if (idx >= 0) {
                                    if (idx < numOutlinks - 1) {
                                        spaceVar.set(stateRow, spaceVarCol, outlinks.get(idx + 1));
                                    } else {
                                        spaceVar.set(stateRow, spaceVarCol, outlinks.get(0));
                                    }
                                }
                            }
                        }
                    }

                    if (sir.get(jobClass) > 0) {
                        outprob = new Matrix(0, 0);
                        for (int k = 0; k < K.get(jobClass); k++) {
                            spaceSrv = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V)); // server state
                            int spaceBufCols = (int) (inspace.getNumCols() - K.elementSum() - V);
                            spaceBuf = Matrix.extract(inspace, 0, inspace.getNumRows(), 0, spaceBufCols); // buffer state
                            Matrix rate = new Matrix(spaceSrv.getNumRows(), 1);
                            rate.zero();
                            Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                            boolean en_set = false;
                            for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                if (spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) > 0) {
                                    en.set(row, 0, 1);
                                    en_set = true;
                                } else {
                                    en.set(row, 0, 0);
                                }
                            }
                            if (en_set) {
                                switch (sn.sched.get(sn.stations.get(ist))) {
                                    case EXT:
                                        // source, can produce an arrival from phase-k as long as it is from an open class
                                        if (Utils.isInf(sn.njobs.get(jobClass))) {
                                            pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
                                            for (int kentry = 0; kentry < K.get(jobClass); kentry++) {
                                                Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V), spaceSrv, 0, 0); // server state

                                                // record a departure
                                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                                    }
                                                }
                                                // record a new job arriving
                                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        spaceSrv.set(row, (int) (Ks.get(jobClass) + kentry), spaceSrv.get(row, (int) (Ks.get(jobClass) + kentry)) + 1);
                                                    }
                                                }
                                                // extract all rows of spaceBuf where en==1 rto space_buf_en
                                                Matrix inspaceEn = new Matrix(0, 0);
                                                Matrix spaceBufEn = new Matrix(0, 0);
                                                Matrix spaceSrvEn = new Matrix(0, 0);
                                                Matrix spaceVarEn = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (inspaceEn.isEmpty()) {
                                                            inspaceEn = Matrix.extractRows(inspace, row, row + 1, null);
                                                        } else {
                                                            inspaceEn = Matrix.concatRows(inspaceEn, Matrix.extractRows(inspace, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceBufEn.isEmpty()) {
                                                            spaceBufEn = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                        } else {
                                                            spaceBufEn = Matrix.concatRows(spaceBufEn, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                        }

                                                    }
                                                }
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceSrvEn.isEmpty()) {
                                                            spaceSrvEn = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        } else {
                                                            spaceSrvEn = Matrix.concatRows(spaceSrvEn, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceVarEn.isEmpty()) {
                                                            spaceVarEn = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        } else {
                                                            spaceVarEn = Matrix.concatRows(spaceVarEn, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        }
                                                    }
                                                }


                                                Matrix left_bottom = Matrix.concatColumns(spaceBufEn, spaceSrvEn, null);
                                                Matrix bottom = Matrix.concatColumns(left_bottom, spaceVarEn, null);
                                                outspace = Matrix.concatRows(outspace, bottom, null);
                                                if (ni.hasInfinite()) {
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    Matrix outrate_bottom = new Matrix(inspace.getNumRows(), 1);
                                                    outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, lldlimit - 1) * pentry.get(kentry) * mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k));
                                                    Matrix ones = new Matrix(inspaceEn.getNumRows(), 1);
                                                    ones.ones();
                                                    outrate_bottom = outrate_bottom.mult(ones);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                } else {
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    Matrix outrate_bottom = new Matrix(inspace.getNumRows(), 1);
                                                    // TODO: check usage of ni
                                                    outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit)) * pentry.get(kentry) * mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k));
                                                    Matrix ones = new Matrix(inspaceEn.getNumRows(), 1);
                                                    ones.ones();
                                                    outrate_bottom = outrate_bottom.mult(ones);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                }
                                                Matrix outprob_bottom = new Matrix(spaceBufEn.getNumRows(), 1);
                                                outprob_bottom.ones();
                                                outprob = Matrix.concatRows(outprob, outprob_bottom, null);
                                            }
                                        }
                                        break;
                                    case INF:
                                        // record a departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // get kir(en, class, k)
                                        Matrix kirEnClassK = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind, 0) == 1) {
                                                if (kirEnClassK.isEmpty()) {
                                                    kirEnClassK = new Matrix(1, 1);
                                                    kirEnClassK.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassK = Matrix.concatRows(kirEnClassK, new_elem, null);
                                                }
                                            }
                                        }
                                        for (int l_ind = 0; l_ind < kirEnClassK.getNumRows(); l_ind++) {
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassK.get(l_ind, 0));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassK.get(l_ind));
                                                if (l_ind < rate.getNumElements()) {
                                                    // replacing an existing element in rate
                                                    rate.set(l_ind, new_elem.value());
                                                } else {
                                                    // expand rate accordingly
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }
                                        // if state unchanged, add with rate 0
                                        Matrix spaceBufEn = new Matrix(0, 0);
                                        Matrix spaceSrvEn = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEn.isEmpty()) {
                                                    spaceBufEn = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                } else {
                                                    spaceBufEn = Matrix.concatRows(spaceBufEn, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                }
                                                if (spaceSrvEn.isEmpty()) {
                                                    spaceSrvEn = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                } else {
                                                    spaceSrvEn = Matrix.concatRows(spaceSrvEn, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                }
                                            }
                                        }
                                        Matrix space_var_last = Matrix.extractRows(spaceVar, spaceVar.getNumRows() - 1, spaceVar.getNumRows(), null);
                                        Matrix left_bottom = Matrix.concatColumns(spaceBufEn, spaceSrvEn, null);
                                        Matrix bottom = Matrix.concatColumns(left_bottom, space_var_last, null);
                                        outspace = Matrix.concatRows(outspace, bottom, null);
                                        Matrix rateEn = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEn.isEmpty()) {
                                                    rateEn = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEn = Matrix.concatRows(rateEn, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }
                                        if (ni.hasInfinite()) {
                                            // hit limited load-dependence
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEn, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEn, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom = new Matrix(rateEn.getNumRows(), 1);
                                        outprob_bottom.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom, null);
                                        break;

                                    case PS:
                                    case LPS:
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        Matrix kirEnClassKPs = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKPs.isEmpty()) {
                                                    kirEnClassKPs = new Matrix(1, 1);
                                                    kirEnClassKPs.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKPs = Matrix.concatRows(kirEnClassKPs, new_elem, null);
                                                }
                                            }
                                        }

                                        // assume active event
                                        for (int l_ind = 0; l_ind < kirEnClassKPs.getNumRows(); l_ind++) {
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * (kirEnClassKPs.get(l_ind) / ni.get(l_ind)) * Maths.min(ni.get(l_ind), S.get(ist)));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                new_elem.set(0, 0, mu_value * phi_value * (kirEnClassKPs.get(l_ind) / ni.get(l_ind)) * Maths.min(ni.get(l_ind), S.get(ist)));
                                                if (l_ind < rate.getNumElements()) {
                                                    // replacing an existing element in rate
                                                    rate.set(l_ind, new_elem.value());
                                                } else {
                                                    // expand rate accordingly
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        Matrix spaceBufEnPs = new Matrix(0, 0);
                                        Matrix spaceSrvEnPs = new Matrix(0, 0);
                                        Matrix spaceVarEn = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEnPs.isEmpty()) {
                                                    spaceBufEnPs = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                } else {
                                                    spaceBufEnPs = Matrix.concatRows(spaceBufEnPs, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                }
                                                if (spaceSrvEnPs.isEmpty()) {
                                                    spaceSrvEnPs = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                } else {
                                                    spaceSrvEnPs = Matrix.concatRows(spaceSrvEnPs, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                }
                                                if (spaceVarEn.isEmpty()) {
                                                    spaceVarEn = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceVarEn = Matrix.concatRows(spaceVarEn, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        Matrix left_bottom_ps = Matrix.concatColumns(spaceBufEnPs, spaceSrvEnPs, null);
                                        Matrix bottom_ps = Matrix.concatColumns(left_bottom_ps, spaceVarEn, null);
                                        outspace = Matrix.concatRows(outspace, bottom_ps, null);
                                        Matrix rateEnPs = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEnPs.isEmpty()) {
                                                    rateEnPs = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnPs = Matrix.concatRows(rateEnPs, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (ni.hasInfinite()) {
                                            // hit limited load-dependence
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnPs, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnPs, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_ps = new Matrix(rateEnPs.getNumRows(), 1);
                                        outprob_bottom_ps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_ps, null);
                                        break;
                                    case PSPRIO:
                                        int minPrio = Integer.MAX_VALUE; // min priority level = most urgent (lower value = higher priority in LINE)
                                        for (int r = 0; r < sn.nclasses; r++) {
                                            int rPrio = (int) sn.classprio.get(r);
                                            if (nir.get(0, r) > 0 && rPrio < minPrio) {
                                                minPrio = rPrio;
                                            }
                                        }
                                        // now check if the class is running or not
                                        if (sn.classprio.get(jobClass) == minPrio) {
                                            Matrix niprio = ni.copy();
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (sn.classprio.get(r) != minPrio) {
                                                    niprio.subEq(Matrix.singleton(nir.get(0, r)));
                                                }
                                            }
                                            // record departure after event
                                            for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                                }
                                            }

                                            Matrix kirEnClassKPsPrio = new Matrix(0, 0);
                                            for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                                if (en.get(l_ind) == 1) {
                                                    if (kirEnClassKPsPrio.isEmpty()) {
                                                        kirEnClassKPsPrio = new Matrix(1, 1);
                                                        kirEnClassKPsPrio.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    } else {
                                                        Matrix new_elem = new Matrix(1, 1);
                                                        new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                        kirEnClassKPsPrio = Matrix.concatRows(kirEnClassKPsPrio, new_elem, null);
                                                    }
                                                }
                                            }

                                            // assume active event
                                            for (int l_ind = 0; l_ind < kirEnClassKPsPrio.getNumRows(); l_ind++) {
                                                if (rate.isEmpty()) {
                                                    rate = new Matrix(1, 1);
                                                    rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * (kirEnClassKPsPrio.get(l_ind) / niprio.get(l_ind)) * Maths.min(niprio.get(l_ind), S.get(ist)));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                    double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                    new_elem.set(0, 0, mu_value * phi_value * (kirEnClassKPsPrio.get(l_ind) / niprio.get(l_ind)) * Maths.min(niprio.get(l_ind), S.get(ist)));
                                                    if (l_ind < rate.getNumElements()) {
                                                        // replacing an existing element in rate
                                                        rate.set(l_ind, new_elem.value());
                                                    } else {
                                                        // expand rate accordingly
                                                        rate = Matrix.concatRows(rate, new_elem, null);
                                                    }
                                                }
                                            }

                                            Matrix spaceBufEnPsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnPsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnPsPrio.isEmpty()) {
                                                        spaceBufEnPsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnPsPrio = Matrix.concatRows(spaceBufEnPsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnPsPrio.isEmpty()) {
                                                        spaceSrvEnPsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnPsPrio = Matrix.concatRows(spaceSrvEnPsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnPrio.isEmpty()) {
                                                        spaceVarEnPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnPrio = Matrix.concatRows(spaceVarEnPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            Matrix left_bottom_psprio = Matrix.concatColumns(spaceBufEnPsPrio, spaceSrvEnPsPrio, null);
                                            Matrix bottom_psprio = Matrix.concatColumns(left_bottom_psprio, spaceVarEnPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_psprio, null);
                                            Matrix rateEnPsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (rateEnPsPrio.isEmpty()) {
                                                        rateEnPsPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        rateEnPsPrio = Matrix.concatRows(rateEnPsPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            if (ni.hasInfinite()) {
                                                // hit limited load-dependence
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnPsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, (int) Maths.min(niprio.get(0), lldscaling.getNumCols() - 1));
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnPsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            }
                                            Matrix outprob_bottom_psprio = new Matrix(rateEnPsPrio.getNumRows(), 1);
                                            outprob_bottom_psprio.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_psprio, null);
                                        } else { // the class is blocked from executing
                                            Matrix spaceBufEnPsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnPsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnPsPrio.isEmpty()) {
                                                        spaceBufEnPsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnPsPrio = Matrix.concatRows(spaceBufEnPsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnPsPrio.isEmpty()) {
                                                        spaceSrvEnPsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnPsPrio = Matrix.concatRows(spaceSrvEnPsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnPrio.isEmpty()) {
                                                        spaceVarEnPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnPrio = Matrix.concatRows(spaceVarEnPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            Matrix left_bottom_psprio = Matrix.concatColumns(spaceBufEnPsPrio, spaceSrvEnPsPrio, null);
                                            Matrix bottom_psprio = Matrix.concatColumns(left_bottom_psprio, spaceVarEnPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_psprio, null);
                                            outrate = Matrix.concatRows(outrate, new Matrix(en.getNumRows(), 1), null);
                                            outprob = Matrix.concatRows(outprob, new Matrix(en.getNumRows(), 1), null);
                                        }

                                        break;

                                    case FCFS:
                                        // job departing
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }
                                        // set en_wbuf to states with jobs in buffer
                                        Matrix enWbuf = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbuf.set(row, 0, 1);
                                            } else {
                                                enWbuf.set(row, 0, 0);
                                            }
                                        }

                                        for (int kdest = 0; kdest < K.get(jobClass); kdest++) {
                                            Matrix space_buf_kd = spaceBuf.copy();
                                            Matrix space_var_kd = spaceVar.copy();
                                            if (ismkvmodclass.get(jobClass) == 1) {
                                                // set space_var_kd(en, sum(sn.nvars(ind,1:class)) = kdest
                                                // MATLAB: sn.nvars(ind,1:class) extracts columns 1 to class (1-based)
                                                // Java: columns 0 to jobClass (0-based), same logical data
                                                int nvarCols = jobClass + 1;
                                                Matrix nvar_ind = new Matrix(1, nvarCols);
                                                Matrix.extract(sn.nvars, ind, ind + 1, 0, nvarCols, nvar_ind, 0, 0);
                                                int nvar_sum = (int) nvar_ind.elementSum();
                                                // MATLAB uses 1-based indexing, Java needs 0-based
                                                int spaceVarCol = nvar_sum - 1;
                                                for (int row = 0; row < space_var_kd.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        // kdest is 0-based in Java, but MAP output var values
                                                        // in the state space are 1-based (initDefault sets to 1),
                                                        // so store kdest+1 to match MATLAB convention
                                                        space_var_kd.set(row, spaceVarCol, kdest + 1);
                                                    }
                                                }
                                            }
                                            Matrix rate_kd = rate.copy();

                                            // set rate_kd(en) = proc{ist}{class}{2}(k,kdest).*kir(en,class,k);


                                            Matrix kirEnClassKFcfs = new Matrix(0, 0);
                                            for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                                if (en.get(l_ind) == 1) {
                                                    if (kirEnClassKFcfs.isEmpty()) {
                                                        kirEnClassKFcfs = new Matrix(1, 1);
                                                        kirEnClassKFcfs.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    } else {
                                                        Matrix new_elem = new Matrix(1, 1);
                                                        new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                        kirEnClassKFcfs = Matrix.concatRows(kirEnClassKFcfs, new_elem, null);
                                                    }
                                                }
                                            }

                                            for (int row = 0; row < rate_kd.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    double D1val = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(1).get(k, kdest);
                                                    double kirVal = kirEnClassKFcfs.get(row);
                                                    double v = D1val * kirVal;
                                                    rate_kd.set(row, 0, v);
                                                }
                                            }
                                            Matrix en_wobuf = new Matrix(enWbuf.getNumRows(), enWbuf.getNumCols());
                                            // set all elems in en_wobuf to 1 where the same elem in enWBuf is 0 and vice versa
                                            boolean anyStateNoJobs = false;
                                            for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                if (enWbuf.get(row, 0) == 0) {
                                                    en_wobuf.set(row, 0, 1);
                                                    anyStateNoJobs = true;
                                                } else {
                                                    en_wobuf.set(row, 0, 0);
                                                }
                                            }
                                            Matrix rate_kd_no_jobs = new Matrix(0, 0);
                                            for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                                                if (en_wobuf.get(row, 0) == 1) {
                                                    if (rate_kd_no_jobs.isEmpty()) {
                                                        rate_kd_no_jobs = Matrix.extractRows(rate_kd, row, row + 1, null);
                                                    } else {
                                                        rate_kd_no_jobs = Matrix.concatRows(rate_kd_no_jobs, Matrix.extractRows(rate_kd, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            if (anyStateNoJobs) {
                                                // set outspace = [outspace; space_buf_kd(en_wobuf,:), space_srv(en_wobuf,:), space_var_kd(en_wobuf,:)];
                                                Matrix space_buf_kd_no_jobs = new Matrix(0, 0);
                                                for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                                                    if (en_wobuf.get(row, 0) == 1) {
                                                        if (space_buf_kd_no_jobs.isEmpty()) {
                                                            space_buf_kd_no_jobs = Matrix.extractRows(space_buf_kd, row, row + 1, null);
                                                        } else {
                                                            space_buf_kd_no_jobs = Matrix.concatRows(space_buf_kd_no_jobs, Matrix.extractRows(space_buf_kd, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                Matrix space_srv_no_jobs = new Matrix(0, 0);
                                                for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                                                    if (en_wobuf.get(row, 0) == 1) {
                                                        if (space_srv_no_jobs.isEmpty()) {
                                                            space_srv_no_jobs = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        } else {
                                                            space_srv_no_jobs = Matrix.concatRows(space_srv_no_jobs, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                Matrix space_var_kd_no_jobs = new Matrix(0, 0);
                                                for (int row = 0; row < en_wobuf.getNumRows(); row++) {
                                                    if (en_wobuf.get(row, 0) == 1) {
                                                        if (space_var_kd_no_jobs.isEmpty()) {
                                                            space_var_kd_no_jobs = Matrix.extractRows(space_var_kd, row, row + 1, null);
                                                        } else {
                                                            space_var_kd_no_jobs = Matrix.concatRows(space_var_kd_no_jobs, Matrix.extractRows(space_var_kd, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                Matrix left_bottom_fcfs = Matrix.concatColumns(space_buf_kd_no_jobs, space_srv_no_jobs, null);
                                                Matrix bottom_fcfs = Matrix.concatColumns(left_bottom_fcfs, space_var_kd_no_jobs, null);
                                                outspace = Matrix.concatRows(outspace, bottom_fcfs, null);
                                                // if all jobs (ni) are Infinite: hit limited load-dependence
                                                if (ni.hasInfinite()) {
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    // must multiply w rate_kd(en_wobuf, :)

                                                    Matrix outrate_bottom_fcfs = Matrix.scaleMult(rate_kd_no_jobs, cdscalingIst * lld);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_fcfs, null);
                                                } else {
                                                    // set outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_kd(en_wobuf,:)];
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lldscaling_ist = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);

                                                    //perform element wise multiplication between cdscalingIst, lldscaling_ist and rate_kd(en_wobuf,:)

                                                    Matrix outrate_bottom_fcfs = Matrix.scaleMult(rate_kd_no_jobs, cdscalingIst * lldscaling_ist);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_fcfs, null);

                                                }
                                            }
                                            // now process states with jobs in buffer
                                            Matrix outprob_bottom_fcfs = new Matrix(rate_kd_no_jobs.getNumRows(), 1);
                                            outprob_bottom_fcfs.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_fcfs, null);
                                            boolean any_jobs_in_buffer = false;
                                            for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                if (enWbuf.get(row, 0) == 1) {
                                                    any_jobs_in_buffer = true;
                                                }
                                            }
                                            if (any_jobs_in_buffer) { // if there is any state with jobs in the buffer
                                                // get class of job at head
                                                Matrix space_buf_kd_last = Matrix.extractColumn(space_buf_kd, space_buf_kd.getNumCols() - 1, null);
                                                // from space_buf_kd_last, extract all rows where en_wbuf = 1 into a matrix start_svc_class
                                                Matrix start_svc_class = new Matrix(0, 0);
                                                for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                    if (enWbuf.get(row, 0) == 1) {
                                                        if (start_svc_class.isEmpty()) {
                                                            start_svc_class = Matrix.extractRows(space_buf_kd_last, row, row + 1, null);
                                                        } else {
                                                            start_svc_class = Matrix.concatRows(start_svc_class, Matrix.extractRows(space_buf_kd_last, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                // if all elements of start_svc_lass are bigger than 0 set boolean x to true
                                                boolean all_elems_bigger_than_zero = true;
                                                for (int row = 0; row < start_svc_class.getNumRows(); row++) {
                                                    if (start_svc_class.get(row, 0) <= 0) {
                                                        all_elems_bigger_than_zero = false;
                                                    }
                                                }
                                                if (all_elems_bigger_than_zero) {
                                                    // update input buffer
                                                    // set space_buf_kd(en_wbuf,:) = [zeros(sum(en_wbuf),1),space_buf_kd(en_wbuf,1:end-1)];
                                                    Matrix left = new Matrix((int) enWbuf.elementSum(), 1);
                                                    left.zero();

                                                    Matrix space_buf_kd_end_removed = Matrix.extract(space_buf_kd, 0, space_buf_kd.getNumRows(), 0, space_buf_kd.getNumCols() - 1);

                                                    // extract into a matrix "right" all rows of space_buf_kd_end_removed where enWbuf = 1
                                                    Matrix right = new Matrix(0, 0);
                                                    for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                        if (enWbuf.get(row, 0) == 1) {
                                                            if (right.isEmpty()) {
                                                                right = Matrix.extractRows(space_buf_kd_end_removed, row, row + 1, null);
                                                            } else {
                                                                right = Matrix.concatRows(right, Matrix.extractRows(space_buf_kd_end_removed, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    Matrix new_space_buf_kd = Matrix.concatColumns(left, right, null);
                                                    // where en_wbuf = 1, take row of that index in "new_space_buf_kd" and assign it to that row in space_buf_kd
                                                    for (int row = 0; row < space_buf_kd.getNumRows(); row++) {
                                                        if (enWbuf.get(row) == 1) {
                                                            // write the row-th row in new_space_buf_kd into this row
                                                            for (int col = 0; col < space_buf_kd.getNumCols(); col++) {
                                                                space_buf_kd.set(row, col, new_space_buf_kd.get(row, col));
                                                            }
                                                        }
                                                    }


                                                    // Check if service class uses MAP/MMPP2 process
                                                    boolean start_svc_class_isMAP = false;
                                                    if (start_svc_class.value() > 0 && start_svc_class.value() <= ismkvmodclass.getNumRows()) {
                                                        start_svc_class_isMAP = ismkvmodclass.get((int) start_svc_class.value() - 1, 0) == 1;
                                                    }
                                                    
                                                    int kentry_start = 0;
                                                    int kentry_range = 0;
                                                    Matrix pentry_svc_class = new Matrix(0, 0);

                                                    if (start_svc_class_isMAP) {
                                                        // Markov-modulated case
                                                        Matrix klassPentryOriginal = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) start_svc_class.value() - 1));
                                                        pentry_svc_class = new Matrix(klassPentryOriginal.getNumRows(), klassPentryOriginal.getNumCols());
                                                        pentry_svc_class.zero();

                                                        if (start_svc_class.value() == jobClass + 1) {
                                                            // Successive service from the same class - new job enters in phase left by departing job
                                                            // Use kdest from the outer loop which tracks the phase of the departing job
                                                            // Use getNumElements() and single-index set() to handle both row and column vector pie orientations
                                                            if (kdest < pentry_svc_class.getNumElements()) {
                                                                pentry_svc_class.set(kdest, 1.0);
                                                                kentry_start = kdest;
                                                                kentry_range = 1;
                                                            } else {
                                                                kentry_start = 0;
                                                                kentry_range = klassPentryOriginal.getNumElements();
                                                                pentry_svc_class = klassPentryOriginal.copy();
                                                            }
                                                        } else {
                                                            // Resume phase from local variables
                                                            // MATLAB: sum(sn.nvars(ind,1:start_svc_class)) extracts columns 1 to start_svc_class (1-based)
                                                            // Java: columns 0 to start_svc_class-1 (0-based), same logical data
                                                            int nvarsSum = 0;
                                                            for (int r = 0; r < start_svc_class.value(); r++) {
                                                                nvarsSum += (int) sn.nvars.get(ind, r);
                                                            }
                                                            // MATLAB uses 1-based indexing, Java needs 0-based
                                                            int spaceVarCol = nvarsSum - 1;
                                                            if (spaceVarCol >= 0 && spaceVarCol < space_var_kd.getNumCols()) {
                                                                // MAP output var values are 1-based in state space,
                                                                // convert to 0-based for Java indexing
                                                                int kentry_resume = (int) space_var_kd.get(0, spaceVarCol) - 1;
                                                                // Use getNumElements() to handle both row and column vector pie orientations
                                                                if (kentry_resume >= 0 && kentry_resume < pentry_svc_class.getNumElements()) {
                                                                    pentry_svc_class.set(kentry_resume, 1.0);
                                                                    kentry_start = kentry_resume;
                                                                    kentry_range = 1;
                                                                } else {
                                                                    kentry_start = 0;
                                                                    kentry_range = klassPentryOriginal.getNumElements();
                                                                    pentry_svc_class = klassPentryOriginal.copy();
                                                                }
                                                            } else {
                                                                kentry_start = 0;
                                                                kentry_range = klassPentryOriginal.getNumElements();
                                                                pentry_svc_class = klassPentryOriginal.copy();
                                                            }
                                                        }
                                                    } else {
                                                        // I.i.d. case
                                                        kentry_start = 0;
                                                        pentry_svc_class = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get((int) start_svc_class.value() - 1));
                                                        kentry_range = (int) K.get((int) (start_svc_class.value() - 1));
                                                    }
                                                    for (int kentry = kentry_start; kentry < kentry_start + kentry_range; kentry++) {
                                                        // increment all values in space_srv at rows where enWbuf is 1 at the column = Ks(start_svc_class.value()+kentry)
                                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                spaceSrv.set(row, (int) (Ks.get((int) start_svc_class.value() - 1) + kentry), spaceSrv.get(row, (int) (Ks.get((int) start_svc_class.value() - 1) + kentry)) + 1);
                                                            }
                                                        }
                                                        // extract 3 matrices: space_buf_kd with only the rows where enWbuf is true, space_srv with only the rows where enWbuf is true, and space_var_kd with only the rows where enWbuf is true
                                                        Matrix space_buf_kd_en = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                if (space_buf_kd_en.isEmpty()) {
                                                                    space_buf_kd_en = Matrix.extractRows(space_buf_kd, row, row + 1, null);
                                                                } else {
                                                                    space_buf_kd_en = Matrix.concatRows(space_buf_kd_en, Matrix.extractRows(space_buf_kd, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }
                                                        Matrix space_srv_en = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                if (space_srv_en.isEmpty()) {
                                                                    space_srv_en = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                                } else {
                                                                    space_srv_en = Matrix.concatRows(space_srv_en, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }
                                                        Matrix space_var_kd_en = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                if (space_var_kd_en.isEmpty()) {
                                                                    space_var_kd_en = Matrix.extractRows(space_var_kd, row, row + 1, null);
                                                                } else {
                                                                    space_var_kd_en = Matrix.concatRows(space_var_kd_en, Matrix.extractRows(space_var_kd, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }
                                                        Matrix left_bottom_outspace = Matrix.concatColumns(space_buf_kd_en, space_srv_en, null);
                                                        Matrix bottom_outspace = Matrix.concatColumns(left_bottom_outspace, space_var_kd_en, null);
                                                        outspace = Matrix.concatRows(outspace, bottom_outspace, null);

                                                        Matrix rate_k = rate_kd.copy();
                                                        // multiply each element in rate_k in rows (across all columns) where enWbuf is one by pentry_svc_class(kentry)
                                                        for (int row = 0; row < rate_k.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                for (int col = 0; col < rate_k.getNumCols(); col++) {
                                                                    rate_k.set(row, col, rate_k.get(row, col) * pentry_svc_class.get(kentry));
                                                                }
                                                            }
                                                        }


                                                        Matrix rate_k_en = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                if (rate_k_en.isEmpty()) {
                                                                    rate_k_en = Matrix.extractRows(rate_k, row, row + 1, null);
                                                                } else {
                                                                    rate_k_en = Matrix.concatRows(rate_k_en, Matrix.extractRows(rate_k, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }

                                                        if (ni.hasInfinite()) {
                                                            // hit limited load-dependence
                                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                            Matrix outrate_bottom = Matrix.scaleMult(rate_k_en, cdscalingIst * lld);
                                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                        } else {
                                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                                            Matrix outrate_bottom = Matrix.scaleMult(rate_k_en, cdscalingIst * lld);
                                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                        }
                                                        // extract rows of rate_kd using enWbuf as indices
                                                        Matrix rate_kd_en = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbuf.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                if (rate_kd_en.isEmpty()) {
                                                                    rate_kd_en = Matrix.extractRows(rate_kd, row, row + 1, null);
                                                                } else {
                                                                    rate_kd_en = Matrix.concatRows(rate_kd_en, Matrix.extractRows(rate_kd, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }
                                                        Matrix outprob_cur = new Matrix(rate_kd_en.getNumRows(), 1);
                                                        outprob_cur.ones();

                                                        // at the indices where outrate == 0, set outprob_cur to 0. make outprob_cur bigger if needed
                                                        for (int row = 0; row < outrate.getNumRows(); row++) {
                                                            if (outrate.get(row) == 0) {
                                                                // may need to expand outprob_cur. add a new row to outprob_cur with one column containing value 0
                                                                if (row < outprob_cur.getNumElements()) {
                                                                    outprob_cur.set(row, 0);
                                                                } else {
                                                                    // extending outprob_cur
                                                                    Matrix new_elem = new Matrix(1, 1);
                                                                    new_elem.set(0, 0, 0);
                                                                    outprob_cur = Matrix.concatColumns(outprob_cur, new_elem, null);
                                                                }

                                                            }
                                                        }


                                                        outprob = Matrix.concatRows(outprob, outprob_cur.transpose(), null);
                                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                            if (enWbuf.get(row, 0) == 1) {
                                                                spaceSrv.set(row, (int) (Ks.get((int) start_svc_class.value() - 1) + kentry), spaceSrv.get(row, (int) (Ks.get((int) start_svc_class.value() - 1) + kentry)) - 1);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        // if state unchanged still add with rate 0
                                        break;

                                    case HOL: // FCFS priority - Head of Line
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        Matrix kirEnClassKHol = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKHol.isEmpty()) {
                                                    kirEnClassKHol = new Matrix(1, 1);
                                                    kirEnClassKHol.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKHol = Matrix.concatRows(kirEnClassKHol, new_elem, null);
                                                }
                                            }
                                        }

                                        for (int l_ind = 0; l_ind < kirEnClassKHol.getNumRows(); l_ind++) {
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKHol.get(l_ind));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKHol.get(l_ind));
                                                if (l_ind < rate.getNumElements()) {
                                                    rate.set(l_ind, new_elem.value());
                                                } else {
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        // set en_wbuf to states with jobs in buffer
                                        Matrix enWbufHol = new Matrix(en.getNumRows(), 1);
                                        Matrix enWobufHol = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufHol.set(row, 0, 1);
                                                enWobufHol.set(row, 0, 0);
                                            } else {
                                                enWbufHol.set(row, 0, 0);
                                                enWobufHol.set(row, 0, 1);
                                            }
                                        }

                                        // Process states without jobs in buffer first
                                        Matrix rateEnWobufHol = new Matrix(0, 0);
                                        for (int row = 0; row < enWobufHol.getNumRows(); row++) {
                                            if (enWobufHol.get(row, 0) == 1) {
                                                if (rateEnWobufHol.isEmpty()) {
                                                    rateEnWobufHol = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnWobufHol = Matrix.concatRows(rateEnWobufHol, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        Matrix spaceBufEnWobufHol = new Matrix(0, 0);
                                        Matrix spaceSrvEnWobufHol = new Matrix(0, 0);
                                        Matrix spaceVarEnWobufHol = new Matrix(0, 0);
                                        for (int row = 0; row < enWobufHol.getNumRows(); row++) {
                                            if (enWobufHol.get(row, 0) == 1) {
                                                if (spaceBufEnWobufHol.isEmpty()) {
                                                    spaceBufEnWobufHol = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    spaceSrvEnWobufHol = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    spaceVarEnWobufHol = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceBufEnWobufHol = Matrix.concatRows(spaceBufEnWobufHol, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    spaceSrvEnWobufHol = Matrix.concatRows(spaceSrvEnWobufHol, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    spaceVarEnWobufHol = Matrix.concatRows(spaceVarEnWobufHol, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (!spaceBufEnWobufHol.isEmpty()) {
                                            Matrix left_bottom_hol = Matrix.concatColumns(spaceBufEnWobufHol, spaceSrvEnWobufHol, null);
                                            Matrix bottom_hol = Matrix.concatColumns(left_bottom_hol, spaceVarEnWobufHol, null);
                                            outspace = Matrix.concatRows(outspace, bottom_hol, null);

                                            if (ni.hasInfinite()) {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom_hol = Matrix.scaleMult(rateEnWobufHol, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom_hol, null);
                                            } else {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                Matrix outrate_bottom_hol = Matrix.scaleMult(rateEnWobufHol, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom_hol, null);
                                            }

                                            Matrix outprob_bottom_hol = new Matrix(rateEnWobufHol.getNumRows(), 1);
                                            outprob_bottom_hol.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_hol, null);
                                        }

                                        // Process states with jobs in buffer using HOL priority logic
                                        boolean any_jobs_in_buffer_hol = false;
                                        for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                            if (enWbufHol.get(row, 0) == 1) {
                                                any_jobs_in_buffer_hol = true;
                                                break;
                                            }
                                        }

                                        if (any_jobs_in_buffer_hol) {
                                            // Create priority group matrix (0 for empty, classprio+1 for occupied)
                                            // Add 1 to classprio so that empty slots (0) are distinguishable from
                                            // highest priority jobs (classprio=0 -> priogroup=1)
                                            Matrix priogroup = new Matrix(1, R + 1);
                                            priogroup.set(0, 0, 0); // empty slot priority is 0
                                            for (int r = 0; r < R; r++) {
                                                priogroup.set(0, r + 1, sn.classprio.get(r) + 1);
                                            }

                                            // Transform buffer to priority groups
                                            Matrix spaceBufGroupg = new Matrix(spaceBuf.getNumRows(), spaceBuf.getNumCols());
                                            for (int row = 0; row < spaceBuf.getNumRows(); row++) {
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    int bufVal = (int) spaceBuf.get(row, col);
                                                    spaceBufGroupg.set(row, col, priogroup.get(0, bufVal));
                                                }
                                            }

                                            // Find minimum priority per row for states with jobs in buffer
                                            // LINE/JMT convention: lower priority value = higher priority (priority 0 is highest)
                                            Matrix startClassprio = new Matrix(0, 0);
                                            Matrix rightmostMinPos = new Matrix(0, 0);

                                            for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                                if (enWbufHol.get(row, 0) == 1) {
                                                    // Find min non-zero priority in this row (lower value = higher priority)
                                                    double minPrioHol = Double.POSITIVE_INFINITY;
                                                    for (int col = 0; col < spaceBufGroupg.getNumCols(); col++) {
                                                        double prio = spaceBufGroupg.get(row, col);
                                                        // Only consider non-zero priorities (0 means empty slot)
                                                        if (prio > 0 && prio < minPrioHol) {
                                                            minPrioHol = prio;
                                                        }
                                                    }

                                                    // Find rightmost position with min priority (for FCFS tie-breaking)
                                                    int rightmostPos = -1;
                                                    for (int col = spaceBufGroupg.getNumCols() - 1; col >= 0; col--) {
                                                        if (spaceBufGroupg.get(row, col) == minPrioHol) {
                                                            rightmostPos = col;
                                                            break;
                                                        }
                                                    }

                                                    if (startClassprio.isEmpty()) {
                                                        startClassprio = new Matrix(1, 1);
                                                        startClassprio.set(0, 0, minPrioHol);
                                                        rightmostMinPos = new Matrix(1, 1);
                                                        rightmostMinPos.set(0, 0, rightmostPos);
                                                    } else {
                                                        Matrix new_prio = new Matrix(1, 1);
                                                        new_prio.set(0, 0, minPrioHol);
                                                        startClassprio = Matrix.concatRows(startClassprio, new_prio, null);
                                                        Matrix new_pos = new Matrix(1, 1);
                                                        new_pos.set(0, 0, rightmostPos);
                                                        rightmostMinPos = Matrix.concatRows(rightmostMinPos, new_pos, null);
                                                    }
                                                }
                                            }

                                            // Get the actual class of the job to serve
                                            Matrix startSvcClass = new Matrix(0, 0);
                                            int startSvcClassRowIdx = 0;
                                            for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                                if (enWbufHol.get(row, 0) == 1) {
                                                    int pos = (int) rightmostMinPos.get(startSvcClassRowIdx, 0);
                                                    if (startSvcClass.isEmpty()) {
                                                        startSvcClass = new Matrix(1, 1);
                                                        startSvcClass.set(0, 0, spaceBuf.get(row, pos));
                                                    } else {
                                                        Matrix new_class = new Matrix(1, 1);
                                                        new_class.set(0, 0, spaceBuf.get(row, pos));
                                                        startSvcClass = Matrix.concatRows(startSvcClass, new_class, null);
                                                    }
                                                    startSvcClassRowIdx++;
                                                }
                                            }

                                            if (!startSvcClass.isEmpty() && startSvcClass.get(0, 0) > 0) {
                                                int svcClassIdx = (int) startSvcClass.get(0, 0) - 1; // Convert to 0-based index
                                                Matrix pentrySvcClass = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(svcClassIdx));
                                                if (pentrySvcClass == null) {
                                                    continue; // Skip if no entry probabilities are defined for this station-class combination
                                                }

                                                for (int kentry = 0; kentry < K.get(svcClassIdx); kentry++) {
                                                    Matrix spaceSrvKHol = spaceSrv.copy();
                                                    Matrix spaceBufKHol = spaceBuf.copy();

                                                    // Add job to service
                                                    int bufRowIdx = 0;
                                                    for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                                        if (enWbufHol.get(row, 0) == 1) {
                                                            spaceSrvKHol.set(row, (int) (Ks.get(svcClassIdx) + kentry), spaceSrvKHol.get(row, (int) (Ks.get(svcClassIdx) + kentry)) + 1);

                                                            // Remove job from buffer using MATLAB logic: [0, before_selected, after_selected]
                                                            int removePos = (int) rightmostMinPos.get(bufRowIdx, 0);
                                                            
                                                            // MATLAB: space_buf_k(j,:) = [0, space_buf_k(j,1:rightmostMinPos(j)-1), space_buf_k(j,(rightmostMinPos(j)+1):end)];
                                                            // Create temporary array to hold the new buffer configuration
                                                            double[] tempBuffer = new double[spaceBufKHol.getNumCols()];
                                                            tempBuffer[0] = 0; // Prepend 0
                                                            
                                                            // Copy jobs before the selected position (1:rightmostMinPos-1)
                                                            int destIdx = 1;
                                                            for (int col = 0; col < removePos; col++) {
                                                                tempBuffer[destIdx++] = spaceBufKHol.get(row, col);
                                                            }
                                                            
                                                            // Copy jobs after the selected position (rightmostMinPos+1:end)
                                                            for (int col = removePos + 1; col < spaceBufKHol.getNumCols(); col++) {
                                                                tempBuffer[destIdx++] = spaceBufKHol.get(row, col);
                                                            }
                                                            
                                                            // Fill remaining positions with 0
                                                            while (destIdx < spaceBufKHol.getNumCols()) {
                                                                tempBuffer[destIdx++] = 0;
                                                            }
                                                            
                                                            // Copy back to the matrix
                                                            for (int col = 0; col < spaceBufKHol.getNumCols(); col++) {
                                                                spaceBufKHol.set(row, col, tempBuffer[col]);
                                                            }
                                                            bufRowIdx++;
                                                        }
                                                    }

                                                    // Extract states with jobs in buffer
                                                    Matrix spaceBufKEnWbufHol = new Matrix(0, 0);
                                                    Matrix spaceSrvKEnWbufHol = new Matrix(0, 0);
                                                    Matrix spaceVarEnWbufHol = new Matrix(0, 0);
                                                    for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                                        if (enWbufHol.get(row, 0) == 1) {
                                                            if (spaceBufKEnWbufHol.isEmpty()) {
                                                                spaceBufKEnWbufHol = Matrix.extractRows(spaceBufKHol, row, row + 1, null);
                                                                spaceSrvKEnWbufHol = Matrix.extractRows(spaceSrvKHol, row, row + 1, null);
                                                                spaceVarEnWbufHol = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                            } else {
                                                                spaceBufKEnWbufHol = Matrix.concatRows(spaceBufKEnWbufHol, Matrix.extractRows(spaceBufKHol, row, row + 1, null), null);
                                                                spaceSrvKEnWbufHol = Matrix.concatRows(spaceSrvKEnWbufHol, Matrix.extractRows(spaceSrvKHol, row, row + 1, null), null);
                                                                spaceVarEnWbufHol = Matrix.concatRows(spaceVarEnWbufHol, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    Matrix left_bottom_hol_buf = Matrix.concatColumns(spaceBufKEnWbufHol, spaceSrvKEnWbufHol, null);
                                                    Matrix bottom_hol_buf = Matrix.concatColumns(left_bottom_hol_buf, spaceVarEnWbufHol, null);
                                                    outspace = Matrix.concatRows(outspace, bottom_hol_buf, null);

                                                    Matrix rateKHol = new Matrix(0, 0);
                                                    bufRowIdx = 0;
                                                    for (int row = 0; row < enWbufHol.getNumRows(); row++) {
                                                        if (enWbufHol.get(row, 0) == 1) {
                                                            if (rateKHol.isEmpty()) {
                                                                rateKHol = new Matrix(1, 1);
                                                                rateKHol.set(0, 0, rate.get(row, 0) * pentrySvcClass.get(kentry));
                                                            } else {
                                                                Matrix new_rate = new Matrix(1, 1);
                                                                new_rate.set(0, 0, rate.get(row, 0) * pentrySvcClass.get(kentry));
                                                                rateKHol = Matrix.concatRows(rateKHol, new_rate, null);
                                                            }
                                                            bufRowIdx++;
                                                        }
                                                    }

                                                    if (ni.hasInfinite()) {
                                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrate_bottom_hol_buf = Matrix.scaleMult(rateKHol, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_hol_buf, null);
                                                    } else {
                                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                        Matrix outrate_bottom_hol_buf = Matrix.scaleMult(rateKHol, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_hol_buf, null);
                                                    }

                                                    Matrix outprob_bottom_hol_buf = new Matrix(rateKHol.getNumRows(), 1);
                                                    outprob_bottom_hol_buf.ones();
                                                    outprob = Matrix.concatRows(outprob, outprob_bottom_hol_buf, null);
                                                }
                                            }
                                        }
                                        break;

                                    case DPS:
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        if (S.get(ist) > 1) {
                                            InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                            }), "Multi-server DPS stations are not supported yet.");
                                        }

                                        // in DPS, the scheduling parameter are the weights
                                        Matrix wDps = sn.schedparam.getRow(ist);
                                        wDps.scaleEq(1.0 / wDps.elementSum());

                                        Matrix kirEnClassKDps = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKDps.isEmpty()) {
                                                    kirEnClassKDps = new Matrix(1, 1);
                                                    kirEnClassKDps.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKDps = Matrix.concatRows(kirEnClassKDps, new_elem, null);
                                                }
                                            }
                                        }

                                        Matrix wDpsRepMatSum = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (wDpsRepMatSum.isEmpty()) {
                                                    wDpsRepMatSum = new Matrix(1, 1);
                                                    wDpsRepMatSum.set(0, 0, wDps.mult(nir.transpose()).sumRows().get(0, 0));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, wDps.mult(nir.transpose()).sumRows().get(0, 0));
                                                    wDpsRepMatSum = Matrix.concatRows(wDpsRepMatSum, new_elem, null);
                                                }
                                            }
                                        }

                                        for (int l_ind = 0; l_ind < kirEnClassKDps.getNumRows(); l_ind++) {
                                            double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                            double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu_value * phi_value * (
                                                        kirEnClassKDps.get(l_ind) / nir.get(jobClass)) * wDps.get(jobClass) *
                                                        nir.get(jobClass) / wDpsRepMatSum.get(l_ind));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu_value * phi_value * (
                                                        kirEnClassKDps.get(l_ind) / nir.get(jobClass)) * wDps.get(jobClass) *
                                                        nir.get(jobClass) / wDpsRepMatSum.get(l_ind));
                                                if (l_ind < rate.getNumElements()) {
                                                    // replacing an existing element in rate
                                                    rate.set(l_ind, new_elem.get(0, 0));
                                                } else {
                                                    // expand rate accordingly
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        Matrix spaceBufEnDps = new Matrix(0, 0);
                                        Matrix spaceSrvEnDps = new Matrix(0, 0);
                                        Matrix spaceVarEnDps = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEnDps.isEmpty()) {
                                                    spaceBufEnDps = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                } else {
                                                    spaceBufEnDps = Matrix.concatRows(spaceBufEnDps, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                }
                                                if (spaceSrvEnDps.isEmpty()) {
                                                    spaceSrvEnDps = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                } else {
                                                    spaceSrvEnDps = Matrix.concatRows(spaceSrvEnDps, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                }
                                                if (spaceVarEnDps.isEmpty()) {
                                                    spaceVarEnDps = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceVarEnDps = Matrix.concatRows(spaceVarEnDps, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        Matrix left_bottom_dps = Matrix.concatColumns(spaceBufEnDps, spaceSrvEnDps, null);
                                        Matrix bottom_dps = Matrix.concatColumns(left_bottom_dps, spaceVarEnDps, null);
                                        outspace = Matrix.concatRows(outspace, bottom_dps, null);
                                        Matrix rateEnDps = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEnDps.isEmpty()) {
                                                    rateEnDps = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnDps = Matrix.concatRows(rateEnDps, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (ni.hasInfinite()) {
                                            // hit limited load-dependence
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnDps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnDps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_dps = new Matrix(rateEnDps.getNumRows(), 1);
                                        outprob_bottom_dps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_dps, null);
                                        break;

                                    case DPSPRIO:
                                        int minPrioDps = Integer.MAX_VALUE; // min priority level = most urgent (lower value = higher priority in LINE)
                                        for (int r = 0; r < sn.nclasses; r++) {
                                            int rPrio = (int) sn.classprio.get(r);
                                            if (nir.get(0, r) > 0 && rPrio < minPrioDps) {
                                                minPrioDps = rPrio;
                                            }
                                        }
                                        // now check if the class is running or not
                                        if (sn.classprio.get(jobClass) == minPrioDps) {
                                            Matrix nirprio = nir.copy();
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (sn.classprio.get(r) != minPrioDps) {
                                                    nirprio.set(r, 0);
                                                }
                                            }
                                            Matrix niprio = nirprio.sumRows();

                                            // record departure
                                            for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                                }
                                            }

                                            if (S.get(ist) > 1) {
                                                InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                                }), "Multi-server DPS stations are not supported yet.");
                                            }

                                            // in DPS, the scheduling parameter are the weights
                                            Matrix wDpsPrio = sn.schedparam.getRow(ist);
                                            wDpsPrio.scaleEq(1.0 / wDpsPrio.elementSum());

                                            Matrix kirEnClassKDpsPrio = new Matrix(0, 0);
                                            for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                                if (en.get(l_ind) == 1) {
                                                    if (kirEnClassKDpsPrio.isEmpty()) {
                                                        kirEnClassKDpsPrio = new Matrix(1, 1);
                                                        kirEnClassKDpsPrio.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    } else {
                                                        Matrix new_elem = new Matrix(1, 1);
                                                        new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                        kirEnClassKDpsPrio = Matrix.concatRows(kirEnClassKDpsPrio, new_elem, null);
                                                    }
                                                }
                                            }

                                            Matrix wDpsPrioRepMatSum = new Matrix(0, 0);
                                            for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                                if (en.get(l_ind) == 1) {
                                                    if (wDpsPrioRepMatSum.isEmpty()) {
                                                        wDpsPrioRepMatSum = new Matrix(1, 1);
                                                        wDpsPrioRepMatSum.set(0, 0, wDpsPrio.mult(nirprio.transpose()).sumRows().get(0, 0));
                                                    } else {
                                                        Matrix new_elem = new Matrix(1, 1);
                                                        new_elem.set(0, 0, wDpsPrio.mult(nirprio.transpose()).sumRows().get(0, 0));
                                                        wDpsPrioRepMatSum = Matrix.concatRows(wDpsPrioRepMatSum, new_elem, null);
                                                    }
                                                }
                                            }

                                            for (int l_ind = 0; l_ind < kirEnClassKDpsPrio.getNumRows(); l_ind++) {
                                                double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                if (rate.isEmpty()) {
                                                    rate = new Matrix(1, 1);
                                                    rate.set(0, 0, mu_value * phi_value * (
                                                            kirEnClassKDpsPrio.get(l_ind) / nirprio.get(jobClass)) * wDpsPrio.get(jobClass) *
                                                            nirprio.get(jobClass) / wDpsPrioRepMatSum.get(l_ind));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, mu_value * phi_value * (
                                                            kirEnClassKDpsPrio.get(l_ind) / nirprio.get(jobClass)) * wDpsPrio.get(jobClass) *
                                                            nirprio.get(jobClass) / wDpsPrioRepMatSum.get(l_ind));
                                                    if (l_ind < rate.getNumElements()) {
                                                        // replacing an existing element in rate
                                                        rate.set(l_ind, new_elem.get(0, 0));
                                                    } else {
                                                        // expand rate accordingly
                                                        rate = Matrix.concatRows(rate, new_elem, null);
                                                    }
                                                }
                                            }

                                            Matrix spaceBufEnDpsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnDpsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnDpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnDpsPrio.isEmpty()) {
                                                        spaceBufEnDpsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnDpsPrio = Matrix.concatRows(spaceBufEnDpsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnDpsPrio.isEmpty()) {
                                                        spaceSrvEnDpsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnDpsPrio = Matrix.concatRows(spaceSrvEnDpsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnDpsPrio.isEmpty()) {
                                                        spaceVarEnDpsPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnDpsPrio = Matrix.concatRows(spaceVarEnDpsPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            Matrix left_bottom_dpsprio = Matrix.concatColumns(spaceBufEnDpsPrio, spaceSrvEnDpsPrio, null);
                                            Matrix bottom_dpsprio = Matrix.concatColumns(left_bottom_dpsprio, spaceVarEnDpsPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_dpsprio, null);
                                            Matrix rateEnDpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (rateEnDpsPrio.isEmpty()) {
                                                        rateEnDpsPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        rateEnDpsPrio = Matrix.concatRows(rateEnDpsPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            if (ni.hasInfinite()) {
                                                // hit limited load-dependence
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nirprio);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnDpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nirprio);
                                                double lld = lldscaling.get(ist, (int) Maths.min(niprio.get(0), lldscaling.getNumCols() - 1));
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnDpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            }
                                            Matrix outprob_bottom_dpsprio = new Matrix(rateEnDpsPrio.getNumRows(), 1);
                                            outprob_bottom_dpsprio.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_dpsprio, null);
                                        } else { // the class is blocked from executing
                                            Matrix spaceBufEnDpsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnDpsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnDpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnDpsPrio.isEmpty()) {
                                                        spaceBufEnDpsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnDpsPrio = Matrix.concatRows(spaceBufEnDpsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnDpsPrio.isEmpty()) {
                                                        spaceSrvEnDpsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnDpsPrio = Matrix.concatRows(spaceSrvEnDpsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnDpsPrio.isEmpty()) {
                                                        spaceVarEnDpsPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnDpsPrio = Matrix.concatRows(spaceVarEnDpsPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            Matrix left_bottom_dpsprio = Matrix.concatColumns(spaceBufEnDpsPrio, spaceSrvEnDpsPrio, null);
                                            Matrix bottom_dpsprio = Matrix.concatColumns(left_bottom_dpsprio, spaceVarEnDpsPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_dpsprio, null);
                                            outrate = Matrix.concatRows(outrate, new Matrix(en.getNumRows(), 1), null);
                                            outprob = Matrix.concatRows(outprob, new Matrix(en.getNumRows(), 1), null);
                                        }
                                        break;

                                    case GPS:
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        if (S.get(ist) > 1) {
                                            InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                            }), "Multi-server GPS stations are not supported yet.");
                                        }

                                        // in GPS, the scheduling parameter are the weights
                                        Matrix wGps = sn.schedparam.getRow(ist);
                                        wGps.scaleEq(1.0 / wGps.elementSum());

                                        Matrix kirEnClassKGps = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKGps.isEmpty()) {
                                                    kirEnClassKGps = new Matrix(1, 1);
                                                    kirEnClassKGps.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKGps = Matrix.concatRows(kirEnClassKGps, new_elem, null);
                                                }
                                            }
                                        }

                                        Matrix cirGps = new Matrix(nir.getNumRows(), nir.getNumCols());
                                        for (int row = 0; row < cirGps.getNumRows(); row++) {
                                            for (int col = 0; col < cirGps.getNumCols(); col++) {
                                                if (nir.get(row, col) < 1.0) {
                                                    cirGps.set(row, col, nir.get(row, col));
                                                } else {
                                                    cirGps.set(row, col, 1.0);
                                                }
                                            }
                                        }

                                        Matrix cirGps1D = new Matrix(0, 0);
                                        for (int col = 0; col < cirGps.getNumCols(); col++) {
                                            cirGps1D = Matrix.concatRows(cirGps1D, cirGps.getColumn(col), null);
                                        }
                                        double wcirGps = wGps.mult(cirGps1D).get(0);

                                        for (int l_ind = 0; l_ind < kirEnClassKGps.getNumRows(); l_ind++) {
                                            double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                            double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu_value * phi_value * (
                                                        kirEnClassKGps.get(l_ind) / nir.get(jobClass)) *
                                                        wGps.get(jobClass) / wcirGps);
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu_value * phi_value * (
                                                        kirEnClassKGps.get(l_ind) / nir.get(jobClass)) *
                                                        wGps.get(jobClass) / wcirGps);
                                                if (l_ind < rate.getNumElements()) {
                                                    // replacing an existing element in rate
                                                    rate.set(l_ind, new_elem.get(0, 0));
                                                } else {
                                                    // expand rate accordingly
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        Matrix spaceBufEnGps = new Matrix(0, 0);
                                        Matrix spaceSrvEnGps = new Matrix(0, 0);
                                        Matrix spaceVarEnGps = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEnGps.isEmpty()) {
                                                    spaceBufEnGps = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                } else {
                                                    spaceBufEnGps = Matrix.concatRows(spaceBufEnGps, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                }
                                                if (spaceSrvEnGps.isEmpty()) {
                                                    spaceSrvEnGps = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                } else {
                                                    spaceSrvEnGps = Matrix.concatRows(spaceSrvEnGps, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                }
                                                if (spaceVarEnGps.isEmpty()) {
                                                    spaceVarEnGps = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceVarEnGps = Matrix.concatRows(spaceVarEnGps, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        Matrix left_bottom_gps = Matrix.concatColumns(spaceBufEnGps, spaceSrvEnGps, null);
                                        Matrix bottom_gps = Matrix.concatColumns(left_bottom_gps, spaceVarEnGps, null);
                                        outspace = Matrix.concatRows(outspace, bottom_gps, null);
                                        Matrix rateEnGps = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEnGps.isEmpty()) {
                                                    rateEnGps = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnGps = Matrix.concatRows(rateEnGps, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (ni.hasInfinite()) {
                                            // hit limited load-dependence
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnGps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnGps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_gps = new Matrix(rateEnGps.getNumRows(), 1);
                                        outprob_bottom_gps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_gps, null);
                                        break;

                                    case GPSPRIO:
                                        int minPrioGps = Integer.MAX_VALUE; // min priority level = most urgent (lower value = higher priority in LINE)
                                        for (int r = 0; r < sn.nclasses; r++) {
                                            int rPrio = (int) sn.classprio.get(r);
                                            if (nir.get(0, r) > 0 && rPrio < minPrioGps) {
                                                minPrioGps = rPrio;
                                            }
                                        }
                                        // now check if the class is running or not
                                        if (sn.classprio.get(jobClass) == minPrioGps) {
                                            Matrix nirprio = nir.copy();
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (sn.classprio.get(r) != minPrioGps) {
                                                    nirprio.set(r, 0);
                                                }
                                            }
                                            Matrix niprio = nirprio.sumRows();

                                            // record departure
                                            for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                                }
                                            }

                                            if (S.get(ist) > 1) {
                                                InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                                }), "Multi-server GPS stations are not supported yet.");
                                            }

                                            // in GPS, the scheduling parameter are the weights
                                            Matrix wGpsPrio = sn.schedparam.getRow(ist);
                                            wGpsPrio.scaleEq(1.0 / wGpsPrio.elementSum());

                                            Matrix kirEnClassKGpsPrio = new Matrix(0, 0);
                                            for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                                if (en.get(l_ind) == 1) {
                                                    if (kirEnClassKGpsPrio.isEmpty()) {
                                                        kirEnClassKGpsPrio = new Matrix(1, 1);
                                                        kirEnClassKGpsPrio.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    } else {
                                                        Matrix new_elem = new Matrix(1, 1);
                                                        new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                        kirEnClassKGpsPrio = Matrix.concatRows(kirEnClassKGpsPrio, new_elem, null);
                                                    }
                                                }
                                            }

                                            Matrix cirGpsPrio = new Matrix(nirprio.getNumRows(), nirprio.getNumCols());
                                            for (int row = 0; row < cirGpsPrio.getNumRows(); row++) {
                                                for (int col = 0; col < cirGpsPrio.getNumCols(); col++) {
                                                    if (nirprio.get(row, col) < 1.0) {
                                                        cirGpsPrio.set(row, col, nirprio.get(row, col));
                                                    } else {
                                                        cirGpsPrio.set(row, col, 1.0);
                                                    }
                                                }
                                            }

                                            Matrix cirGpsPrio1D = new Matrix(0, 0);
                                            for (int col = 0; col < cirGpsPrio.getNumCols(); col++) {
                                                cirGpsPrio1D = Matrix.concatRows(cirGpsPrio1D, cirGpsPrio.getColumn(col), null);
                                            }
                                            double wcirGpsPrio = wGpsPrio.mult(cirGpsPrio1D).get(0);

                                            for (int l_ind = 0; l_ind < kirEnClassKGpsPrio.getNumRows(); l_ind++) {
                                                double mu_value = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double phi_value = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                if (rate.isEmpty()) {
                                                    rate = new Matrix(1, 1);
                                                    rate.set(0, 0, mu_value * phi_value * (
                                                            kirEnClassKGpsPrio.get(l_ind) / nirprio.get(jobClass)) *
                                                            wGpsPrio.get(jobClass) / wcirGpsPrio);
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, mu_value * phi_value * (
                                                            kirEnClassKGpsPrio.get(l_ind) / nirprio.get(jobClass)) *
                                                            wGpsPrio.get(jobClass) / wcirGpsPrio);
                                                    if (l_ind < rate.getNumElements()) {
                                                        // replacing an existing element in rate
                                                        rate.set(l_ind, new_elem.get(0, 0));
                                                    } else {
                                                        // expand rate accordingly
                                                        rate = Matrix.concatRows(rate, new_elem, null);
                                                    }
                                                }
                                            }

                                            Matrix spaceBufEnGpsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnGpsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnGpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnGpsPrio.isEmpty()) {
                                                        spaceBufEnGpsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnGpsPrio = Matrix.concatRows(spaceBufEnGpsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnGpsPrio.isEmpty()) {
                                                        spaceSrvEnGpsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnGpsPrio = Matrix.concatRows(spaceSrvEnGpsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnGpsPrio.isEmpty()) {
                                                        spaceVarEnGpsPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnGpsPrio = Matrix.concatRows(spaceVarEnGpsPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            Matrix left_bottom_gpsprio = Matrix.concatColumns(spaceBufEnGpsPrio, spaceSrvEnGpsPrio, null);
                                            Matrix bottom_gpsprio = Matrix.concatColumns(left_bottom_gpsprio, spaceVarEnGpsPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_gpsprio, null);
                                            Matrix rateEnGpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (rateEnGpsPrio.isEmpty()) {
                                                        rateEnGpsPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        rateEnGpsPrio = Matrix.concatRows(rateEnGpsPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            if (ni.hasInfinite()) {
                                                // hit limited load-dependence
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nirprio);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnGpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nirprio);
                                                double lld = lldscaling.get(ist, (int) Maths.min(niprio.get(0), lldscaling.getNumCols() - 1));
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnGpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            }
                                            Matrix outprob_bottom_gpsprio = new Matrix(rateEnGpsPrio.getNumRows(), 1);
                                            outprob_bottom_gpsprio.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_gpsprio, null);
                                        } else { // the class is blocked from executing
                                            Matrix spaceBufEnGpsPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnGpsPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnGpsPrio = new Matrix(0, 0);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    if (spaceBufEnGpsPrio.isEmpty()) {
                                                        spaceBufEnGpsPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnGpsPrio = Matrix.concatRows(spaceBufEnGpsPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    }
                                                    if (spaceSrvEnGpsPrio.isEmpty()) {
                                                        spaceSrvEnGpsPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    } else {
                                                        spaceSrvEnGpsPrio = Matrix.concatRows(spaceSrvEnGpsPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    }
                                                    if (spaceVarEnGpsPrio.isEmpty()) {
                                                        spaceVarEnGpsPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceVarEnGpsPrio = Matrix.concatRows(spaceVarEnGpsPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            Matrix left_bottom_gpsprio = Matrix.concatColumns(spaceBufEnGpsPrio, spaceSrvEnGpsPrio, null);
                                            Matrix bottom_gpsprio = Matrix.concatColumns(left_bottom_gpsprio, spaceVarEnGpsPrio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_gpsprio, null);
                                            outrate = Matrix.concatRows(outrate, new Matrix(en.getNumRows(), 1), null);
                                            outprob = Matrix.concatRows(outprob, new Matrix(en.getNumRows(), 1), null);
                                        }
                                        break;

                                    case LCFS: // Last Come First Served
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        Matrix kirEnClassKLcfs = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKLcfs.isEmpty()) {
                                                    kirEnClassKLcfs = new Matrix(1, 1);
                                                    kirEnClassKLcfs.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKLcfs = Matrix.concatRows(kirEnClassKLcfs, new_elem, null);
                                                }
                                            }
                                        }

                                        for (int l_ind = 0; l_ind < kirEnClassKLcfs.getNumRows(); l_ind++) {
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKLcfs.get(l_ind));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKLcfs.get(l_ind));
                                                if (l_ind < rate.getNumElements()) {
                                                    rate.set(l_ind, new_elem.value());
                                                } else {
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        // set en_wbuf to states with jobs in buffer
                                        Matrix enWbufLcfs = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufLcfs.set(row, 0, 1);
                                            } else {
                                                enWbufLcfs.set(row, 0, 0);
                                            }
                                        }

                                        // Find first non-zero column (leftmost job = most recent arrival)
                                        Matrix colFirstNnz = new Matrix(0, 0);
                                        Matrix startSvcClassLcfs = new Matrix(0, 0);

                                        for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                            if (enWbufLcfs.get(row, 0) == 1) {
                                                // Find first non-zero column in this row
                                                int firstCol = -1;
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    if (spaceBuf.get(row, col) != 0) {
                                                        firstCol = col;
                                                        break;
                                                    }
                                                }

                                                if (colFirstNnz.isEmpty()) {
                                                    colFirstNnz = new Matrix(1, 1);
                                                    colFirstNnz.set(0, 0, firstCol);
                                                    startSvcClassLcfs = new Matrix(1, 1);
                                                    if (firstCol >= 0) {
                                                        startSvcClassLcfs.set(0, 0, spaceBuf.get(row, firstCol));
                                                    } else {
                                                        startSvcClassLcfs.set(0, 0, 0);
                                                    }
                                                } else {
                                                    Matrix new_col = new Matrix(1, 1);
                                                    new_col.set(0, 0, firstCol);
                                                    colFirstNnz = Matrix.concatRows(colFirstNnz, new_col, null);
                                                    Matrix new_class = new Matrix(1, 1);
                                                    if (firstCol >= 0) {
                                                        new_class.set(0, 0, spaceBuf.get(row, firstCol));
                                                    } else {
                                                        new_class.set(0, 0, 0);
                                                    }
                                                    startSvcClassLcfs = Matrix.concatRows(startSvcClassLcfs, new_class, null);
                                                }
                                            }
                                        }

                                        // Remove job from buffer (set to 0)
                                        Matrix spaceBufLcfs = spaceBuf.copy();
                                        int bufRowIdx = 0;
                                        for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                            if (enWbufLcfs.get(row, 0) == 1) {
                                                int firstCol = (int) colFirstNnz.get(bufRowIdx, 0);
                                                if (firstCol >= 0) {
                                                    spaceBufLcfs.set(row, firstCol, 0);
                                                }
                                                bufRowIdx++;
                                            }
                                        }

                                        // Process states without buffer jobs first
                                        Matrix enWobufLcfs = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (enWbufLcfs.get(row, 0) == 0) {
                                                enWobufLcfs.set(row, 0, 1);
                                            } else {
                                                enWobufLcfs.set(row, 0, 0);
                                            }
                                        }

                                        // Handle states without buffer jobs (just departure, no buffer-to-service transition)
                                        boolean hasStatesWithoutBuffer = false;
                                        for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && enWbufLcfs.get(row, 0) == 0) {
                                                hasStatesWithoutBuffer = true;
                                                break;
                                            }
                                        }
                                        
                                        if (hasStatesWithoutBuffer) {
                                            // Extract states without buffer jobs
                                            Matrix spaceBufWobuf = new Matrix(0, 0);
                                            Matrix spaceSrvWobuf = new Matrix(0, 0);
                                            Matrix spaceVarWobuf = new Matrix(0, 0);
                                            Matrix rateWobuf = new Matrix(0, 0);
                                            
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1 && enWbufLcfs.get(row, 0) == 0) {
                                                    if (spaceBufWobuf.isEmpty()) {
                                                        spaceBufWobuf = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                        spaceSrvWobuf = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        spaceVarWobuf = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        rateWobuf = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        spaceBufWobuf = Matrix.concatRows(spaceBufWobuf, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                        spaceSrvWobuf = Matrix.concatRows(spaceSrvWobuf, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        spaceVarWobuf = Matrix.concatRows(spaceVarWobuf, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        rateWobuf = Matrix.concatRows(rateWobuf, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            
                                            // Add output for states without buffer jobs
                                            Matrix left_bottom_wobuf = Matrix.concatColumns(spaceBufWobuf, spaceSrvWobuf, null);
                                            Matrix bottom_wobuf = Matrix.concatColumns(left_bottom_wobuf, spaceVarWobuf, null);
                                            outspace = Matrix.concatRows(outspace, bottom_wobuf, null);
                                            
                                            // Calculate output rate
                                            if (ni.hasInfinite()) {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_wobuf = Matrix.scaleMult(rateWobuf, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_wobuf, null);
                                            } else {
                                                double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                Matrix outrate_wobuf = Matrix.scaleMult(rateWobuf, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_wobuf, null);
                                            }
                                            
                                            Matrix outprob_wobuf = new Matrix(rateWobuf.getNumRows(), 1);
                                            outprob_wobuf.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_wobuf, null);
                                        }
                                        
                                        if (!startSvcClassLcfs.isEmpty() && startSvcClassLcfs.getNumRows() > 0) {
                                            boolean hasValidStartClass = false;
                                            for (int i = 0; i < colFirstNnz.getNumRows(); i++) {
                                                if (colFirstNnz.get(i, 0) >= 0) {  // Check if we found a valid column, not if the class value is > 0
                                                    hasValidStartClass = true;
                                                    break;
                                                }
                                            }

                                            if (!hasValidStartClass) {
                                                // No valid job to start, just add current state
                                                Matrix spaceBufEnLcfs = new Matrix(0, 0);
                                                Matrix spaceSrvEnLcfs = new Matrix(0, 0);
                                                Matrix spaceVarEnLcfs = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceBufEnLcfs.isEmpty()) {
                                                            spaceBufEnLcfs = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                            spaceSrvEnLcfs = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                            spaceVarEnLcfs = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        } else {
                                                            spaceBufEnLcfs = Matrix.concatRows(spaceBufEnLcfs, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                            spaceSrvEnLcfs = Matrix.concatRows(spaceSrvEnLcfs, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                            spaceVarEnLcfs = Matrix.concatRows(spaceVarEnLcfs, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                Matrix left_bottom_lcfs = Matrix.concatColumns(spaceBufEnLcfs, spaceSrvEnLcfs, null);
                                                Matrix bottom_lcfs = Matrix.concatColumns(left_bottom_lcfs, spaceVarEnLcfs, null);
                                                outspace = Matrix.concatRows(outspace, bottom_lcfs, null);

                                                Matrix rateEnLcfs = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (rateEnLcfs.isEmpty()) {
                                                            rateEnLcfs = Matrix.extractRows(rate, row, row + 1, null);
                                                        } else {
                                                            rateEnLcfs = Matrix.concatRows(rateEnLcfs, Matrix.extractRows(rate, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                if (ni.hasInfinite()) {
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrate_bottom_lcfs = Matrix.scaleMult(rateEnLcfs, cdscalingIst * lld);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs, null);
                                                } else {
                                                    double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                    Matrix outrate_bottom_lcfs = Matrix.scaleMult(rateEnLcfs, cdscalingIst * lld);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs, null);
                                                }

                                                Matrix outprob_bottom_lcfs = new Matrix(rateEnLcfs.getNumRows(), 1);
                                                outprob_bottom_lcfs.ones();
                                                outprob = Matrix.concatRows(outprob, outprob_bottom_lcfs, null);
                                                break;
                                            }

                                            // Process each possible entry phase for the starting job
                                            int startClass = (int) startSvcClassLcfs.get(0, 0) - 1; // Convert to 0-based index
                                            if (startClass >= 0) {
                                                Matrix pentrySvcClass = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(startClass));
                                                if (pentrySvcClass == null) {
                                                    continue; // Skip if no entry probabilities are defined for this station-class combination
                                                }

                                                for (int kentry = 0; kentry < K.get(startClass); kentry++) {
                                                    Matrix spaceSrvKLcfs = spaceSrv.copy();

                                                    // Add job to service for states with buffer jobs
                                                    for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                                        if (enWbufLcfs.get(row, 0) == 1) {
                                                            spaceSrvKLcfs.set(row, (int) (Ks.get(startClass) + kentry), spaceSrvKLcfs.get(row, (int) (Ks.get(startClass) + kentry)) + 1);
                                                        }
                                                    }

                                                    // Extract states with enabled servers
                                                    Matrix spaceBufEnKLcfs = new Matrix(0, 0);
                                                    Matrix spaceSrvEnKLcfs = new Matrix(0, 0);
                                                    Matrix spaceVarEnKLcfs = new Matrix(0, 0);
                                                    for (int row = 0; row < en.getNumRows(); row++) {
                                                        if (en.get(row, 0) == 1) {
                                                            if (spaceBufEnKLcfs.isEmpty()) {
                                                                spaceBufEnKLcfs = Matrix.extractRows(spaceBufLcfs, row, row + 1, null);
                                                                spaceSrvEnKLcfs = Matrix.extractRows(spaceSrvKLcfs, row, row + 1, null);
                                                                spaceVarEnKLcfs = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                            } else {
                                                                spaceBufEnKLcfs = Matrix.concatRows(spaceBufEnKLcfs, Matrix.extractRows(spaceBufLcfs, row, row + 1, null), null);
                                                                spaceSrvEnKLcfs = Matrix.concatRows(spaceSrvEnKLcfs, Matrix.extractRows(spaceSrvKLcfs, row, row + 1, null), null);
                                                                spaceVarEnKLcfs = Matrix.concatRows(spaceVarEnKLcfs, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    Matrix left_bottom_lcfs_k = Matrix.concatColumns(spaceBufEnKLcfs, spaceSrvEnKLcfs, null);
                                                    Matrix bottom_lcfs_k = Matrix.concatColumns(left_bottom_lcfs_k, spaceVarEnKLcfs, null);
                                                    outspace = Matrix.concatRows(outspace, bottom_lcfs_k, null);

                                                    Matrix rateKLcfs = new Matrix(0, 0);
                                                    // Build rate for ALL enabled states, not just states with buffer
                                                    for (int row = 0; row < en.getNumRows(); row++) {
                                                        if (en.get(row, 0) == 1) {
                                                            double rateValue;
                                                            // Only multiply by entry probability for states with buffer jobs
                                                            if (enWbufLcfs.get(row, 0) == 1) {
                                                                rateValue = rate.get(row, 0) * pentrySvcClass.get(kentry);
                                                            } else {
                                                                rateValue = rate.get(row, 0);
                                                            }
                                                            
                                                            if (rateKLcfs.isEmpty()) {
                                                                rateKLcfs = new Matrix(1, 1);
                                                                rateKLcfs.set(0, 0, rateValue);
                                                            } else {
                                                                Matrix new_rate = new Matrix(1, 1);
                                                                new_rate.set(0, 0, rateValue);
                                                                rateKLcfs = Matrix.concatRows(rateKLcfs, new_rate, null);
                                                            }
                                                        }
                                                    }

                                                    if (ni.hasInfinite()) {
                                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrate_bottom_lcfs_k = Matrix.scaleMult(rateKLcfs, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs_k, null);
                                                    } else {
                                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                        Matrix outrate_bottom_lcfs_k = Matrix.scaleMult(rateKLcfs, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs_k, null);
                                                    }

                                                    Matrix outprob_bottom_lcfs_k = new Matrix(rateKLcfs.getNumRows(), 1);
                                                    outprob_bottom_lcfs_k.ones();
                                                    outprob = Matrix.concatRows(outprob, outprob_bottom_lcfs_k, null);

                                                    // Remove job from service to reset for next kentry
                                                    for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                                        if (enWbufLcfs.get(row, 0) == 1) {
                                                            spaceSrvKLcfs.set(row, (int) (Ks.get(startClass) + kentry), spaceSrvKLcfs.get(row, (int) (Ks.get(startClass) + kentry)) - 1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        break;

                                    case LCFSPR:
                                        // LCFSPR departure - record departure from service
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // Calculate rate for enabled states (simple calculation like MATLAB)
                                        rate = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                double muVal = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double phiVal = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double kirVal = kir.get(k).get(row, jobClass);
                                                rate.set(row, 0, muVal * phiVal * kirVal);
                                            } else {
                                                rate.set(row, 0, 0);
                                            }
                                        }

                                        Matrix enWbufLcfspr = new Matrix(en.getNumRows(), 1);
                                        // States with jobs in buffer (buffer is stored as pairs: class, phase)
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufLcfspr.set(row, 0, 1);
                                            } else {
                                                enWbufLcfspr.set(row, 0, 0);
                                            }
                                        }

                                        // Find first non-zero column (leftmost job = most recent)
                                        Matrix colFirstNnzLcfspr = new Matrix(0, 0);
                                        Matrix startSvcClassLcfspr = new Matrix(0, 0);
                                        Matrix kentryLcfspr = new Matrix(0, 0);

                                        for (int row = 0; row < enWbufLcfspr.getNumRows(); row++) {
                                            if (enWbufLcfspr.get(row, 0) == 1) {
                                                // Find first non-zero column in this row
                                                int firstCol = -1;
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    if (spaceBuf.get(row, col) != 0) {
                                                        firstCol = col;
                                                        break;
                                                    }
                                                }

                                                if (colFirstNnzLcfspr.isEmpty()) {
                                                    colFirstNnzLcfspr = new Matrix(1, 1);
                                                    colFirstNnzLcfspr.set(0, 0, firstCol);
                                                    startSvcClassLcfspr = new Matrix(1, 1);
                                                    kentryLcfspr = new Matrix(1, 1);
                                                    if (firstCol >= 0) {
                                                        startSvcClassLcfspr.set(0, 0, spaceBuf.get(row, firstCol)); // class
                                                        kentryLcfspr.set(0, 0, spaceBuf.get(row, firstCol + 1)); // phase
                                                    } else {
                                                        startSvcClassLcfspr.set(0, 0, 0);
                                                        kentryLcfspr.set(0, 0, 0);
                                                    }
                                                } else {
                                                    Matrix new_col = new Matrix(1, 1);
                                                    new_col.set(0, 0, firstCol);
                                                    colFirstNnzLcfspr = Matrix.concatRows(colFirstNnzLcfspr, new_col, null);
                                                    Matrix new_class = new Matrix(1, 1);
                                                    Matrix new_phase = new Matrix(1, 1);
                                                    if (firstCol >= 0) {
                                                        new_class.set(0, 0, spaceBuf.get(row, firstCol));
                                                        new_phase.set(0, 0, spaceBuf.get(row, firstCol + 1));
                                                    } else {
                                                        new_class.set(0, 0, 0);
                                                        new_phase.set(0, 0, 0);
                                                    }
                                                    startSvcClassLcfspr = Matrix.concatRows(startSvcClassLcfspr, new_class, null);
                                                    kentryLcfspr = Matrix.concatRows(kentryLcfspr, new_phase, null);
                                                }
                                            }
                                        }

                                        // Remove job from buffer (set both class and phase to 0)
                                        Matrix spaceBufLcfspr = spaceBuf.copy();
                                        int bufRowIdxLcfspr = 0;
                                        for (int row = 0; row < enWbufLcfspr.getNumRows(); row++) {
                                            if (enWbufLcfspr.get(row, 0) == 1) {
                                                int firstCol = (int) colFirstNnzLcfspr.get(bufRowIdxLcfspr, 0);
                                                if (firstCol >= 0) {
                                                    spaceBufLcfspr.set(row, firstCol, 0); // zero popped job class
                                                    spaceBufLcfspr.set(row, firstCol + 1, 0); // zero popped phase
                                                }
                                                bufRowIdxLcfspr++;
                                            }
                                        }

                                        // Check if we have valid jobs to start
                                        if (!startSvcClassLcfspr.isEmpty() && startSvcClassLcfspr.getNumRows() > 0) {
                                            boolean hasValidStartClassLcfspr = false;
                                            for (int i = 0; i < startSvcClassLcfspr.getNumRows(); i++) {
                                                if (startSvcClassLcfspr.get(i, 0) > 0) {
                                                    hasValidStartClassLcfspr = true;
                                                    break;
                                                }
                                            }

                                            if (!hasValidStartClassLcfspr) {
                                                // No valid job to start, just add current state
                                                Matrix spaceBufEnLcfspr = new Matrix(0, 0);
                                                Matrix spaceSrvEnLcfspr = new Matrix(0, 0);
                                                Matrix spaceVarEnLcfspr = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceBufEnLcfspr.isEmpty()) {
                                                            spaceBufEnLcfspr = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                            spaceSrvEnLcfspr = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                            spaceVarEnLcfspr = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        } else {
                                                            spaceBufEnLcfspr = Matrix.concatRows(spaceBufEnLcfspr, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                            spaceSrvEnLcfspr = Matrix.concatRows(spaceSrvEnLcfspr, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                            spaceVarEnLcfspr = Matrix.concatRows(spaceVarEnLcfspr, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                Matrix leftBottomLcfspr = Matrix.concatColumns(spaceBufEnLcfspr, spaceSrvEnLcfspr, null);
                                                Matrix bottomLcfspr = Matrix.concatColumns(leftBottomLcfspr, spaceVarEnLcfspr, null);
                                                outspace = Matrix.concatRows(outspace, bottomLcfspr, null);

                                                Matrix rateEnLcfspr = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (rateEnLcfspr.isEmpty()) {
                                                            rateEnLcfspr = Matrix.extractRows(rate, row, row + 1, null);
                                                        } else {
                                                            rateEnLcfspr = Matrix.concatRows(rateEnLcfspr, Matrix.extractRows(rate, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                if (ni.hasInfinite()) {
                                                    double cdscalingIstLcfspr = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lldLcfspr = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                                } else {
                                                    double cdscalingIstLcfspr = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lldLcfspr = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                    Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                                }

                                                Matrix outprobBottomLcfspr = new Matrix(rateEnLcfspr.getNumRows(), 1);
                                                outprobBottomLcfspr.ones();
                                                outprob = Matrix.concatRows(outprob, outprobBottomLcfspr, null);

                                                if (isSimulation && eventCache.isEnabled()) {
                                                    eventCache.put(key, new Ret.EventResult(outspace, outrate, outprob));
                                                }
                                                return new Ret.EventResult(outspace, outrate, outprob);
                                            }

                                            // Add job to service with preserved phase
                                            bufRowIdxLcfspr = 0;
                                            for (int row = 0; row < enWbufLcfspr.getNumRows(); row++) {
                                                if (enWbufLcfspr.get(row, 0) == 1) {
                                                    int startClass = (int) startSvcClassLcfspr.get(bufRowIdxLcfspr, 0) - 1; // Convert to 0-based index
                                                    int kentry = (int) kentryLcfspr.get(bufRowIdxLcfspr, 0);
                                                    // Use Ks.length() instead of Ks.getNumRows() since Ks is a row vector
                                                    if (startClass >= 0 && startClass < Ks.length()) {
                                                        // kentry from buffer is 1-based phase, convert to 0-based for column index
                                                        int colIndex = (int) (Ks.get(startClass) + kentry - 1);
                                                        if (colIndex >= 0 && colIndex < spaceSrv.getNumCols()) {
                                                            spaceSrv.set(row, colIndex, spaceSrv.get(row, colIndex) + 1);
                                                        }
                                                    }
                                                    bufRowIdxLcfspr++;
                                                }
                                            }
                                        }

                                        // Add state to output
                                        Matrix leftBottomLcfspr = Matrix.concatColumns(spaceBufLcfspr, spaceSrv, null);
                                        Matrix bottomLcfspr = Matrix.concatColumns(leftBottomLcfspr, spaceVar, null);

                                        Matrix spaceBufEnLcfspr = new Matrix(0, 0);
                                        Matrix spaceSrvEnLcfspr = new Matrix(0, 0);
                                        Matrix spaceVarEnLcfspr = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEnLcfspr.isEmpty()) {
                                                    spaceBufEnLcfspr = Matrix.extractRows(bottomLcfspr, row, row + 1, null);
                                                } else {
                                                    spaceBufEnLcfspr = Matrix.concatRows(spaceBufEnLcfspr, Matrix.extractRows(bottomLcfspr, row, row + 1, null), null);
                                                }
                                            }
                                        }
                                        outspace = Matrix.concatRows(outspace, spaceBufEnLcfspr, null);

                                        Matrix rateEnLcfspr = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEnLcfspr.isEmpty()) {
                                                    rateEnLcfspr = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnLcfspr = Matrix.concatRows(rateEnLcfspr, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (ni.hasInfinite()) {
                                            double cdscalingIstLcfspr = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lldLcfspr = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                            outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                        } else {
                                            double cdscalingIstLcfspr = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                            double lldLcfspr = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                            Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                            outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                        }

                                        Matrix outprobBottomLcfspr = new Matrix(rateEnLcfspr.getNumRows(), 1);
                                        outprobBottomLcfspr.ones();
                                        outprob = Matrix.concatRows(outprob, outprobBottomLcfspr, null);
                                        break;
                                    case LCFSPI:
                                        // LCFSPI departure - record departure from service and restart jobs from pie distribution
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }
                                        
                                        // Calculate basic rate for enabled states 
                                        Matrix rateLcfspi = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                double muVal = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double phiVal = phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k);
                                                double kirVal = kir.get(k).get(row, jobClass);
                                                rateLcfspi.set(row, 0, muVal * phiVal * kirVal);
                                            } else {
                                                rateLcfspi.set(row, 0, 0);
                                            }
                                        }
                                        
                                        // Handle states without buffer jobs (no job promotion)
                                        Matrix enWobufLcfspi = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) <= S.get(ist)) {
                                                enWobufLcfspi.set(row, 0, 1);
                                            } else {
                                                enWobufLcfspi.set(row, 0, 0);
                                            }
                                        }
                                        
                                        if (enWobufLcfspi.elementSum() > 0) {
                                            // Extract states without buffer jobs
                                            Matrix spaceBufEnWobufLcfspi = new Matrix(0, 0);
                                            Matrix spaceSrvEnWobufLcfspi = new Matrix(0, 0);
                                            Matrix spaceVarEnWobufLcfspi = new Matrix(0, 0);
                                            for (int row = 0; row < enWobufLcfspi.getNumRows(); row++) {
                                                if (enWobufLcfspi.get(row, 0) == 1) {
                                                    if (spaceBufEnWobufLcfspi.isEmpty()) {
                                                        spaceBufEnWobufLcfspi = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                        spaceSrvEnWobufLcfspi = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        spaceVarEnWobufLcfspi = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnWobufLcfspi = Matrix.concatRows(spaceBufEnWobufLcfspi, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                        spaceSrvEnWobufLcfspi = Matrix.concatRows(spaceSrvEnWobufLcfspi, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        spaceVarEnWobufLcfspi = Matrix.concatRows(spaceVarEnWobufLcfspi, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            
                                            if (!spaceBufEnWobufLcfspi.isEmpty()) {
                                                Matrix leftBottomWobufLcfspi = Matrix.concatColumns(spaceBufEnWobufLcfspi, spaceSrvEnWobufLcfspi, null);
                                                Matrix bottomWobufLcfspi = Matrix.concatColumns(leftBottomWobufLcfspi, spaceVarEnWobufLcfspi, null);
                                                outspace = Matrix.concatRows(outspace, bottomWobufLcfspi, null);
                                                
                                                Matrix rateEnWobufLcfspi = new Matrix(0, 0);
                                                for (int row = 0; row < enWobufLcfspi.getNumRows(); row++) {
                                                    if (enWobufLcfspi.get(row, 0) == 1) {
                                                        if (rateEnWobufLcfspi.isEmpty()) {
                                                            rateEnWobufLcfspi = Matrix.extractRows(rateLcfspi, row, row + 1, null);
                                                        } else {
                                                            rateEnWobufLcfspi = Matrix.concatRows(rateEnWobufLcfspi, Matrix.extractRows(rateLcfspi, row, row + 1, null), null);
                                                        }
                                                    }
                                                }
                                                
                                                if (ni.hasInfinite()) {
                                                    double cdscalingIstLcfspi = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lldLcfspi = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomWobufLcfspi = Matrix.scaleMult(rateEnWobufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomWobufLcfspi, null);
                                                } else {
                                                    double cdscalingIstLcfspi = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    double lldLcfspi = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                    Matrix outrateBottomWobufLcfspi = Matrix.scaleMult(rateEnWobufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomWobufLcfspi, null);
                                                }
                                                
                                                Matrix outprobBottomWobufLcfspi = new Matrix(rateEnWobufLcfspi.getNumRows(), 1);
                                                outprobBottomWobufLcfspi.ones();
                                                outprob = Matrix.concatRows(outprob, outprobBottomWobufLcfspi, null);
                                            }
                                        }
                                        
                                        // Handle states with buffer jobs - LCFSPI uses pie distribution for restart
                                        Matrix enWbufLcfspi = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufLcfspi.set(row, 0, 1);
                                            } else {
                                                enWbufLcfspi.set(row, 0, 0);
                                            }
                                        }
                                        
                                        if (enWbufLcfspi.elementSum() > 0) {
                                            // Find first non-zero column in buffer (last job to arrive - LCFS order)
                                            Matrix startSvcClassLcfspi = new Matrix(enWbufLcfspi.getNumRows(), 1);
                                            for (int row = 0; row < enWbufLcfspi.getNumRows(); row++) {
                                                if (enWbufLcfspi.get(row, 0) == 1) {
                                                    int firstCol = -1;
                                                    for (int col = 0; col < spaceBuf.getNumCols(); col += 2) { // Buffer stores (class, phase) pairs
                                                        if (spaceBuf.get(row, col) != 0) {
                                                            firstCol = col;
                                                            break;
                                                        }
                                                    }
                                                    if (firstCol >= 0) {
                                                        startSvcClassLcfspi.set(row, 0, spaceBuf.get(row, firstCol)); // Extract class
                                                    } else {
                                                        startSvcClassLcfspi.set(row, 0, 0);
                                                    }
                                                } else {
                                                    startSvcClassLcfspi.set(row, 0, 0);
                                                }
                                            }
                                            
                                            // Get unique start service classes that need processing
                                            java.util.Set<Integer> uniqueStartClasses = new java.util.HashSet<>();
                                            for (int row = 0; row < startSvcClassLcfspi.getNumRows(); row++) {
                                                int startClass = (int) startSvcClassLcfspi.get(row, 0);
                                                if (startClass > 0) {
                                                    uniqueStartClasses.add(startClass);
                                                }
                                            }
                                            
                                            // Process each unique start service class
                                            for (int startClass : uniqueStartClasses) {
                                                int startClassIdx = startClass - 1; // Convert to 0-based index
                                                Matrix pentrySvcClass = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(startClassIdx));
                                                if (pentrySvcClass == null) {
                                                    continue; // Skip if no entry probabilities are defined
                                                }
                                                
                                                // For each possible entry phase according to pie distribution
                                                for (int kentry = 0; kentry < K.get(startClassIdx); kentry++) {
                                                    // Create copies of state matrices for this phase
                                                    Matrix spaceBufLcfspi = spaceBuf.copy();
                                                    Matrix spaceSrvLcfspi = spaceSrv.copy();
                                                    
                                                    // Remove job from buffer and add to service
                                                    for (int row = 0; row < enWbufLcfspi.getNumRows(); row++) {
                                                        if (enWbufLcfspi.get(row, 0) == 1 && startSvcClassLcfspi.get(row, 0) == startClass) {
                                                            // Find first occurrence of this class in buffer
                                                            for (int col = 0; col < spaceBufLcfspi.getNumCols(); col += 2) {
                                                                if (spaceBufLcfspi.get(row, col) == startClass) {
                                                                    spaceBufLcfspi.set(row, col, 0); // Clear class
                                                                    spaceBufLcfspi.set(row, col + 1, 0); // Clear phase (ignored for LCFSPI)
                                                                    break;
                                                                }
                                                            }
                                                            // Add job to service in phase kentry (according to pie distribution)
                                                            spaceSrvLcfspi.set(row, (int) (Ks.get(startClassIdx) + kentry), 
                                                                             spaceSrvLcfspi.get(row, (int) (Ks.get(startClassIdx) + kentry)) + 1);
                                                        }
                                                    }
                                                    
                                                    // Build output states for this phase
                                                    Matrix enWbufClassLcfspi = new Matrix(enWbufLcfspi.getNumRows(), 1);
                                                    for (int row = 0; row < enWbufLcfspi.getNumRows(); row++) {
                                                        if (enWbufLcfspi.get(row, 0) == 1 && startSvcClassLcfspi.get(row, 0) == startClass) {
                                                            enWbufClassLcfspi.set(row, 0, 1);
                                                        } else {
                                                            enWbufClassLcfspi.set(row, 0, 0);
                                                        }
                                                    }
                                                    
                                                    if (enWbufClassLcfspi.elementSum() > 0) {
                                                        Matrix spaceBufEnWbufLcfspi = new Matrix(0, 0);
                                                        Matrix spaceSrvEnWbufLcfspi = new Matrix(0, 0);
                                                        Matrix spaceVarEnWbufLcfspi = new Matrix(0, 0);
                                                        
                                                        for (int row = 0; row < enWbufClassLcfspi.getNumRows(); row++) {
                                                            if (enWbufClassLcfspi.get(row, 0) == 1) {
                                                                if (spaceBufEnWbufLcfspi.isEmpty()) {
                                                                    spaceBufEnWbufLcfspi = Matrix.extractRows(spaceBufLcfspi, row, row + 1, null);
                                                                    spaceSrvEnWbufLcfspi = Matrix.extractRows(spaceSrvLcfspi, row, row + 1, null);
                                                                    spaceVarEnWbufLcfspi = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                                } else {
                                                                    spaceBufEnWbufLcfspi = Matrix.concatRows(spaceBufEnWbufLcfspi, Matrix.extractRows(spaceBufLcfspi, row, row + 1, null), null);
                                                                    spaceSrvEnWbufLcfspi = Matrix.concatRows(spaceSrvEnWbufLcfspi, Matrix.extractRows(spaceSrvLcfspi, row, row + 1, null), null);
                                                                    spaceVarEnWbufLcfspi = Matrix.concatRows(spaceVarEnWbufLcfspi, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }
                                                        
                                                        if (!spaceBufEnWbufLcfspi.isEmpty()) {
                                                            Matrix leftBottomWbufLcfspi = Matrix.concatColumns(spaceBufEnWbufLcfspi, spaceSrvEnWbufLcfspi, null);
                                                            Matrix bottomWbufLcfspi = Matrix.concatColumns(leftBottomWbufLcfspi, spaceVarEnWbufLcfspi, null);
                                                            outspace = Matrix.concatRows(outspace, bottomWbufLcfspi, null);
                                                            
                                                            // Apply pie distribution probability to rates
                                                            Matrix rateKLcfspi = rateLcfspi.copy();
                                                            for (int row = 0; row < enWbufClassLcfspi.getNumRows(); row++) {
                                                                if (enWbufClassLcfspi.get(row, 0) == 1) {
                                                                    rateKLcfspi.set(row, 0, rateKLcfspi.get(row, 0) * pentrySvcClass.get(kentry));
                                                                }
                                                            }
                                                            
                                                            Matrix rateEnWbufLcfspi = new Matrix(0, 0);
                                                            for (int row = 0; row < enWbufClassLcfspi.getNumRows(); row++) {
                                                                if (enWbufClassLcfspi.get(row, 0) == 1) {
                                                                    if (rateEnWbufLcfspi.isEmpty()) {
                                                                        rateEnWbufLcfspi = Matrix.extractRows(rateKLcfspi, row, row + 1, null);
                                                                    } else {
                                                                        rateEnWbufLcfspi = Matrix.concatRows(rateEnWbufLcfspi, Matrix.extractRows(rateKLcfspi, row, row + 1, null), null);
                                                                    }
                                                                }
                                                            }
                                                            
                                                            if (ni.hasInfinite()) {
                                                                double cdscalingIstLcfspi = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                                double lldLcfspi = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                                Matrix outrateBottomWbufLcfspi = Matrix.scaleMult(rateEnWbufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                                outrate = Matrix.concatRows(outrate, outrateBottomWbufLcfspi, null);
                                                            } else {
                                                                double cdscalingIstLcfspi = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                                double lldLcfspi = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                                Matrix outrateBottomWbufLcfspi = Matrix.scaleMult(rateEnWbufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                                outrate = Matrix.concatRows(outrate, outrateBottomWbufLcfspi, null);
                                                            }
                                                            
                                                            Matrix outprobBottomWbufLcfspi = new Matrix(rateEnWbufLcfspi.getNumRows(), 1);
                                                            outprobBottomWbufLcfspi.ones();
                                                            outprob = Matrix.concatRows(outprob, outprobBottomWbufLcfspi, null);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        break;

                                    case SIRO:
                                        // SIRO (Service In Random Order) - pick a job from buffer randomly by class
                                        rate = new Matrix(spaceSrv.getNumRows(), 1);
                                        rate.fill(0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                rate.set(row, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kir.get(k).get(row, jobClass));
                                            }
                                        }

                                        spaceSrv = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V));
                                        // Record departure for all states where en==1
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // First handle states where buffer is empty
                                        Matrix enWobufSiro = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                double bufSum = 0;
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    bufSum += spaceBuf.get(row, col);
                                                }
                                                if (bufSum == 0) {
                                                    enWobufSiro.set(row, 0, 1);
                                                } else {
                                                    enWobufSiro.set(row, 0, 0);
                                                }
                                            } else {
                                                enWobufSiro.set(row, 0, 0);
                                            }
                                        }

                                        // Add states without buffer to output
                                        Matrix spaceBufEnWobufSiro = new Matrix(0, 0);
                                        Matrix spaceSrvEnWobufSiro = new Matrix(0, 0);
                                        Matrix spaceVarEnWobufSiro = new Matrix(0, 0);
                                        for (int row = 0; row < enWobufSiro.getNumRows(); row++) {
                                            if (enWobufSiro.get(row, 0) == 1) {
                                                if (spaceBufEnWobufSiro.isEmpty()) {
                                                    spaceBufEnWobufSiro = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    spaceSrvEnWobufSiro = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    spaceVarEnWobufSiro = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceBufEnWobufSiro = Matrix.concatRows(spaceBufEnWobufSiro, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    spaceSrvEnWobufSiro = Matrix.concatRows(spaceSrvEnWobufSiro, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    spaceVarEnWobufSiro = Matrix.concatRows(spaceVarEnWobufSiro, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (!spaceBufEnWobufSiro.isEmpty()) {
                                            Matrix leftBottomWobufSiro = Matrix.concatColumns(spaceBufEnWobufSiro, spaceSrvEnWobufSiro, null);
                                            Matrix bottomWobufSiro = Matrix.concatColumns(leftBottomWobufSiro, spaceVarEnWobufSiro, null);
                                            outspace = Matrix.concatRows(outspace, bottomWobufSiro, null);

                                            Matrix rateEnWobufSiro = new Matrix(0, 0);
                                            for (int row = 0; row < enWobufSiro.getNumRows(); row++) {
                                                if (enWobufSiro.get(row, 0) == 1) {
                                                    if (rateEnWobufSiro.isEmpty()) {
                                                        rateEnWobufSiro = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        rateEnWobufSiro = Matrix.concatRows(rateEnWobufSiro, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            if (ni.hasInfinite()) {
                                                double cdscalingIstSiro = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lldSiro = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrateBottomWobufSiro = Matrix.scaleMult(rateEnWobufSiro, cdscalingIstSiro * lldSiro);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSiro, null);
                                            } else {
                                                double cdscalingIstSiro = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lldSiro = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                Matrix outrateBottomWobufSiro = Matrix.scaleMult(rateEnWobufSiro, cdscalingIstSiro * lldSiro);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSiro, null);
                                            }

                                            Matrix outprobBottomWobufSiro = new Matrix(rateEnWobufSiro.getNumRows(), 1);
                                            outprobBottomWobufSiro.ones();
                                            outprob = Matrix.concatRows(outprob, outprobBottomWobufSiro, null);
                                        }

                                        // Handle states with buffer - pick jobs randomly from each class
                                        for (int r = 0; r < R; r++) {
                                            Matrix rateR = rate.copy();
                                            Matrix spaceBufR = Matrix.extract(inspace, 0, inspace.getNumRows(), 0, (int) (inspace.getNumCols() - K.elementSum() - V));

                                            // Find states where class r has jobs in buffer
                                            Matrix enWbufSiro = new Matrix(en.getNumRows(), 1);
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1 && spaceBufR.get(row, r) > 0) {
                                                    enWbufSiro.set(row, 0, 1);
                                                } else {
                                                    enWbufSiro.set(row, 0, 0);
                                                }
                                            }

                                            // Remove job from buffer
                                            for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                if (enWbufSiro.get(row, 0) == 1) {
                                                    spaceBufR.set(row, r, spaceBufR.get(row, r) - 1);
                                                }
                                            }

                                            Matrix spaceSrvR = spaceSrv.copy();
                                            Matrix pentrySvcClass = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
                                            if (pentrySvcClass == null) {
                                                continue; // Skip if no entry probabilities are defined for this station-class combination
                                            }

                                            // Calculate pick probability for random selection
                                            for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                if (enWbufSiro.get(row, 0) == 1) {
                                                    double pickProb = (nir.get(row, r) - sir.get(row, r)) / (ni.get(row) - sir.getRow(row).elementSum());
                                                    if (pickProb >= 0) {
                                                        rateR.set(row, 0, rateR.get(row, 0) * pickProb);
                                                    }
                                                }
                                            }

                                            // For each entry phase
                                            for (int kentry = 0; kentry < K.get(r); kentry++) {
                                                // Add job to service in phase kentry
                                                for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                    if (enWbufSiro.get(row, 0) == 1) {
                                                        spaceSrvR.set(row, (int) (Ks.get(r) + kentry), spaceSrvR.get(row, (int) (Ks.get(r) + kentry)) + 1);
                                                    }
                                                }

                                                // Build output state
                                                Matrix spaceBufEnWbufSiro = new Matrix(0, 0);
                                                Matrix spaceSrvEnWbufSiro = new Matrix(0, 0);
                                                Matrix spaceVarEnWbufSiro = new Matrix(0, 0);
                                                for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                    if (enWbufSiro.get(row, 0) == 1) {
                                                        if (spaceBufEnWbufSiro.isEmpty()) {
                                                            spaceBufEnWbufSiro = Matrix.extractRows(spaceBufR, row, row + 1, null);
                                                            spaceSrvEnWbufSiro = Matrix.extractRows(spaceSrvR, row, row + 1, null);
                                                            spaceVarEnWbufSiro = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        } else {
                                                            spaceBufEnWbufSiro = Matrix.concatRows(spaceBufEnWbufSiro, Matrix.extractRows(spaceBufR, row, row + 1, null), null);
                                                            spaceSrvEnWbufSiro = Matrix.concatRows(spaceSrvEnWbufSiro, Matrix.extractRows(spaceSrvR, row, row + 1, null), null);
                                                            spaceVarEnWbufSiro = Matrix.concatRows(spaceVarEnWbufSiro, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                if (!spaceBufEnWbufSiro.isEmpty()) {
                                                    Matrix leftBottomWbufSiro = Matrix.concatColumns(spaceBufEnWbufSiro, spaceSrvEnWbufSiro, null);
                                                    Matrix bottomWbufSiro = Matrix.concatColumns(leftBottomWbufSiro, spaceVarEnWbufSiro, null);
                                                    outspace = Matrix.concatRows(outspace, bottomWbufSiro, null);

                                                    Matrix rateKSiro = rateR.copy();
                                                    for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                        if (enWbufSiro.get(row, 0) == 1) {
                                                            rateKSiro.set(row, 0, rateKSiro.get(row, 0) * pentrySvcClass.get(kentry));
                                                        }
                                                    }

                                                    Matrix rateEnWbufSiro = new Matrix(0, 0);
                                                    for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                        if (enWbufSiro.get(row, 0) == 1) {
                                                            if (rateEnWbufSiro.isEmpty()) {
                                                                rateEnWbufSiro = Matrix.extractRows(rateKSiro, row, row + 1, null);
                                                            } else {
                                                                rateEnWbufSiro = Matrix.concatRows(rateEnWbufSiro, Matrix.extractRows(rateKSiro, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    if (ni.hasInfinite()) {
                                                        double cdscalingIstSiro = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lldSiro = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrateBottomWbufSiro = Matrix.scaleMult(rateEnWbufSiro, cdscalingIstSiro * lldSiro);
                                                        outrate = Matrix.concatRows(outrate, outrateBottomWbufSiro, null);
                                                    } else {
                                                        double cdscalingIstSiro = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                        double lldSiro = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                        Matrix outrateBottomWbufSiro = Matrix.scaleMult(rateEnWbufSiro, cdscalingIstSiro * lldSiro);
                                                        outrate = Matrix.concatRows(outrate, outrateBottomWbufSiro, null);
                                                    }

                                                    Matrix outprobBottomWbufSiro = new Matrix(rateEnWbufSiro.getNumRows(), 1);
                                                    outprobBottomWbufSiro.ones();
                                                    outprob = Matrix.concatRows(outprob, outprobBottomWbufSiro, null);
                                                }

                                                // Reset server state for next kentry
                                                for (int row = 0; row < enWbufSiro.getNumRows(); row++) {
                                                    if (enWbufSiro.get(row, 0) == 1) {
                                                        spaceSrvR.set(row, (int) (Ks.get(r) + kentry), spaceSrvR.get(row, (int) (Ks.get(r) + kentry)) - 1);
                                                    }
                                                }
                                            }
                                        }
                                        break;

                                    case SEPT:
                                    case LEPT:
                                        // SEPT/LEPT (Shortest/Longest Expected Processing Time)
                                        rate = new Matrix(spaceSrv.getNumRows(), 1);
                                        rate.fill(0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                rate.set(row, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kir.get(k).get(row, jobClass));
                                            }
                                        }

                                        spaceSrv = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V));
                                        // Record departure for all states where en==1
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // First handle states where buffer is empty
                                        Matrix enWobufSeptLept = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                double bufSum = 0;
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    bufSum += spaceBuf.get(row, col);
                                                }
                                                if (bufSum == 0) {
                                                    enWobufSeptLept.set(row, 0, 1);
                                                } else {
                                                    enWobufSeptLept.set(row, 0, 0);
                                                }
                                            } else {
                                                enWobufSeptLept.set(row, 0, 0);
                                            }
                                        }

                                        // Handle states with empty buffer
                                        boolean enWobufSet = false;
                                        for (int row = 0; row < enWobufSeptLept.getNumRows(); row++) {
                                            if (enWobufSeptLept.get(row, 0) == 1) {
                                                enWobufSet = true;
                                                break;
                                            }
                                        }

                                        if (enWobufSet) {
                                            Matrix rateEnWobufSeptLept = new Matrix(0, 0);
                                            for (int row = 0; row < enWobufSeptLept.getNumRows(); row++) {
                                                if (enWobufSeptLept.get(row, 0) == 1) {
                                                    if (rateEnWobufSeptLept.isEmpty()) {
                                                        rateEnWobufSeptLept = Matrix.extractRows(rate, row, row + 1, null);
                                                    } else {
                                                        rateEnWobufSeptLept = Matrix.concatRows(rateEnWobufSeptLept, Matrix.extractRows(rate, row, row + 1, null), null);
                                                    }
                                                }
                                            }

                                            // Extract rows for enabled states
                                            Matrix spaceBufEnWobufSeptLept = new Matrix(0, 0);
                                            Matrix spaceSrvEnWobufSeptLept = new Matrix(0, 0);
                                            Matrix spaceVarEnWobufSeptLept = new Matrix(0, 0);
                                            
                                            for (int row = 0; row < enWobufSeptLept.getNumRows(); row++) {
                                                if (enWobufSeptLept.get(row, 0) == 1) {
                                                    if (spaceBufEnWobufSeptLept.isEmpty()) {
                                                        spaceBufEnWobufSeptLept = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                        spaceSrvEnWobufSeptLept = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        spaceVarEnWobufSeptLept = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnWobufSeptLept = Matrix.concatRows(spaceBufEnWobufSeptLept, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                        spaceSrvEnWobufSeptLept = Matrix.concatRows(spaceSrvEnWobufSeptLept, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        spaceVarEnWobufSeptLept = Matrix.concatRows(spaceVarEnWobufSeptLept, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            
                                            if (!spaceBufEnWobufSeptLept.isEmpty()) {
                                                Matrix leftBottomWobufSeptLept = Matrix.concatColumns(spaceBufEnWobufSeptLept, spaceSrvEnWobufSeptLept, null);
                                                Matrix bottomWobufSeptLept = Matrix.concatColumns(leftBottomWobufSeptLept, spaceVarEnWobufSeptLept, null);
                                                outspace = Matrix.concatRows(outspace, bottomWobufSeptLept, null);
                                            }

                                            if (ni.hasInfinite()) {
                                                double cdscalingIstSeptLept = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lldSeptLept = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrateBottomWobufSeptLept = Matrix.scaleMult(rateEnWobufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSeptLept, null);
                                            } else {
                                                double cdscalingIstSeptLept = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                double lldSeptLept = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldlimit));
                                                Matrix outrateBottomWobufSeptLept = Matrix.scaleMult(rateEnWobufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSeptLept, null);
                                            }

                                            Matrix outprobBottomWobufSeptLept = new Matrix(rateEnWobufSeptLept.getNumRows(), 1);
                                            outprobBottomWobufSeptLept.ones();
                                            outprob = Matrix.concatRows(outprob, outprobBottomWobufSeptLept, null);
                                        }

                                        // Handle states with non-empty buffer - need to select job based on expected processing time
                                        Matrix enWbufSeptLept = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && enWobufSeptLept.get(row, 0) == 0) {
                                                enWbufSeptLept.set(row, 0, 1);
                                            } else {
                                                enWbufSeptLept.set(row, 0, 0);
                                            }
                                        }

                                        boolean enWbufSet = false;
                                        for (int row = 0; row < enWbufSeptLept.getNumRows(); row++) {
                                            if (enWbufSeptLept.get(row, 0) == 1) {
                                                enWbufSet = true;
                                                break;
                                            }
                                        }

                                        if (enWbufSet) {
                                            // For SEPT/LEPT, we need to select which class enters service based on expected processing time
                                            // First, determine which classes have jobs in buffer
                                            for (int r = 0; r < R; r++) {
                                                boolean hasJobsInBuffer = false;
                                                for (int row = 0; row < spaceBuf.getNumRows(); row++) {
                                                    if (enWbufSeptLept.get(row, 0) == 1 && spaceBuf.get(row, r) > 0) {
                                                        hasJobsInBuffer = true;
                                                        break;
                                                    }
                                                }

                                                if (hasJobsInBuffer) {
                                                    // For states with jobs of class r in buffer
                                                    for (int kentry = 0; kentry < K.get(r); kentry++) {
                                                        Matrix spaceBufR = spaceBuf.copy();
                                                        Matrix spaceSrvR = spaceSrv.copy();

                                                        // Determine if this class should enter service based on SEPT/LEPT priority
                                                        boolean shouldEnterService = true;
                                                        double classRExpectedTime = 1.0 / sn.rates.get(ist, r);

                                                        // Check against other classes in buffer
                                                        for (int otherClass = 0; otherClass < R; otherClass++) {
                                                            if (otherClass != r) {
                                                                boolean hasOtherJobsInBuffer = false;
                                                                for (int row = 0; row < spaceBuf.getNumRows(); row++) {
                                                                    if (enWbufSeptLept.get(row, 0) == 1 && spaceBuf.get(row, otherClass) > 0) {
                                                                        hasOtherJobsInBuffer = true;
                                                                        break;
                                                                    }
                                                                }

                                                                if (hasOtherJobsInBuffer) {
                                                                    double otherClassExpectedTime = 1.0 / sn.rates.get(ist, otherClass);
                                                                    if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.SEPT) {
                                                                        // SEPT: Select shortest expected processing time
                                                                        if (otherClassExpectedTime < classRExpectedTime) {
                                                                            shouldEnterService = false;
                                                                            break;
                                                                        }
                                                                    } else { // LEPT
                                                                        // LEPT: Select longest expected processing time
                                                                        if (otherClassExpectedTime > classRExpectedTime) {
                                                                            shouldEnterService = false;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }

                                                        if (shouldEnterService) {
                                                            // Remove job from buffer and add to server
                                                            for (int row = 0; row < enWbufSeptLept.getNumRows(); row++) {
                                                                if (enWbufSeptLept.get(row, 0) == 1 && spaceBufR.get(row, r) > 0) {
                                                                    spaceBufR.set(row, r, spaceBufR.get(row, r) - 1);
                                                                    spaceSrvR.set(row, (int) (Ks.get(r) + kentry), spaceSrvR.get(row, (int) (Ks.get(r) + kentry)) + 1);
                                                                }
                                                            }

                                                            Matrix rateEnWbufSeptLept = new Matrix(0, 0);
                                                            for (int row = 0; row < enWbufSeptLept.getNumRows(); row++) {
                                                                if (enWbufSeptLept.get(row, 0) == 1 && spaceBuf.get(row, r) > 0) {
                                                                    if (rateEnWbufSeptLept.isEmpty()) {
                                                                        rateEnWbufSeptLept = Matrix.extractRows(rate, row, row + 1, null);
                                                                    } else {
                                                                        rateEnWbufSeptLept = Matrix.concatRows(rateEnWbufSeptLept, Matrix.extractRows(rate, row, row + 1, null), null);
                                                                    }
                                                                }
                                                            }

                                                            if (!rateEnWbufSeptLept.isEmpty()) {
                                                                // Extract rows for enabled states with jobs of class r in buffer
                                                                Matrix spaceBufREnWbuf = new Matrix(0, 0);
                                                                Matrix spaceSrvREnWbuf = new Matrix(0, 0);
                                                                Matrix spaceVarEnWbuf = new Matrix(0, 0);
                                                                
                                                                for (int row = 0; row < enWbufSeptLept.getNumRows(); row++) {
                                                                    if (enWbufSeptLept.get(row, 0) == 1 && spaceBuf.get(row, r) > 0) {
                                                                        if (spaceBufREnWbuf.isEmpty()) {
                                                                            spaceBufREnWbuf = Matrix.extractRows(spaceBufR, row, row + 1, null);
                                                                            spaceSrvREnWbuf = Matrix.extractRows(spaceSrvR, row, row + 1, null);
                                                                            spaceVarEnWbuf = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                                        } else {
                                                                            spaceBufREnWbuf = Matrix.concatRows(spaceBufREnWbuf, Matrix.extractRows(spaceBufR, row, row + 1, null), null);
                                                                            spaceSrvREnWbuf = Matrix.concatRows(spaceSrvREnWbuf, Matrix.extractRows(spaceSrvR, row, row + 1, null), null);
                                                                            spaceVarEnWbuf = Matrix.concatRows(spaceVarEnWbuf, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                                        }
                                                                    }
                                                                }
                                                                
                                                                if (!spaceBufREnWbuf.isEmpty()) {
                                                                    Matrix leftBottomWbufR = Matrix.concatColumns(spaceBufREnWbuf, spaceSrvREnWbuf, null);
                                                                    Matrix bottomWbufR = Matrix.concatColumns(leftBottomWbufR, spaceVarEnWbuf, null);
                                                                    outspace = Matrix.concatRows(outspace, bottomWbufR, null);
                                                                }

                                                                if (ni.hasInfinite()) {
                                                                    double cdscalingIstSeptLept = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                                    double lldSeptLept = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                                    Matrix outrateBottomWbufSeptLept = Matrix.scaleMult(rateEnWbufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                                    outrate = Matrix.concatRows(outrate, outrateBottomWbufSeptLept, null);
                                                                } else {
                                                                    double cdscalingIstSeptLept = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                                    double lldSeptLept = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                                    Matrix outrateBottomWbufSeptLept = Matrix.scaleMult(rateEnWbufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                                    outrate = Matrix.concatRows(outrate, outrateBottomWbufSeptLept, null);
                                                                }

                                                                Matrix outprobBottomWbufSeptLept = new Matrix(rateEnWbufSeptLept.getNumRows(), 1);
                                                                outprobBottomWbufSeptLept.ones();
                                                                outprob = Matrix.concatRows(outprob, outprobBottomWbufSeptLept, null);
                                                            }

                                                            // Reset server state for next kentry
                                                            for (int row = 0; row < enWbufSeptLept.getNumRows(); row++) {
                                                                if (enWbufSeptLept.get(row, 0) == 1) {
                                                                    spaceSrvR.set(row, (int) (Ks.get(r) + kentry), spaceSrvR.get(row, (int) (Ks.get(r) + kentry)) - 1);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        break;
                                    default:
                                        throw new RuntimeException(String.format("Scheduling strategy %s is not supported",
                                                sn.sched.get(sn.nodes.get(ind)).toString()));


                                }


                            }
                        }
                        Ret.EventResult result_d = new Ret.EventResult(outspace, outrate, outprob);
                        eventCache.put(key, result_d);
                        if (isSimulation) {
                            if (outspace.getNumRows() > 1) {
                                Matrix tot_rate = outrate.sumCols();
                                Matrix cum_sum = outrate.cumsumViaCol();
                                Matrix cum_rate = Matrix.scaleMult(cum_sum, 1.0 / tot_rate.value());
                                int firing_ctr = -1;
                                double rand = Maths.rand();
                                // we need the indicies where rand is bigger than cum_prob
                                for (int row = 0; row < cum_rate.getNumRows(); row++) {
                                    if (rand > cum_rate.get(row)) {
                                        firing_ctr = row;
                                    }
                                }
                                firing_ctr++;
                                outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                                double outrate_val = outrate.elementSum();
                                outrate = new Matrix(1, 1);
                                outrate.set(0, 0, outrate_val);
                                outprob = Matrix.extractRows(outprob, firing_ctr, firing_ctr + 1, null);
                            }

                        }
                    }
                } else {
                    Ret.EventResult result_d = new Ret.EventResult(outspace, outrate, outprob);
                    eventCache.put(key, result_d);
                }
                break;
            case PHASE:
                outspace = new Matrix(0, 0);
                outrate = new Matrix(0, 0);
                outprob = new Matrix(0, 0);
                State.StateMarginalStatistics stateMarginalStatistics = ToMarginal.toMarginal(sn, ind, inspace, K, Ks, spaceBuf, spaceSrv, spaceVar);
                ni = stateMarginalStatistics.ni;
                nir = stateMarginalStatistics.nir;
                kir = stateMarginalStatistics.kir;

                if (nir.get(jobClass) > 0) {
                    for (int k = 0; k < K.get(jobClass); k++) {
                        // set en = space_srv(:,Ks(class)+k) > 0;
                        // set en to a matrix which has a 1 if that row in space_srv in column Ks(class) + k is bigger than 0, and a 0 if not
                        Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                        en.zero();
                        boolean any_en = false;
                        for (int row = 0; row < en.getNumRows(); row++) {
                            if (spaceSrv.get(row, (int) Ks.get(jobClass) + k) > 0) {
                                en.set(row, 0, 1);
                                any_en = true;
                            }
                        }

                        if (any_en) {
                            for (int kdest = 0; kdest < K.get(jobClass); kdest++) {
                                if (kdest != k) {
                                    Matrix rate = new Matrix(1, 1);

                                    Matrix spaceSrvK = new Matrix(0, 0);
                                    for (int i = 0; i < en.getNumRows(); i++) {
                                        if (en.get(i, 0) == 1) {
                                            if (spaceSrvK.isEmpty()) {
                                                spaceSrvK = Matrix.extractRows(spaceSrv, i, i + 1, null);
                                            } else {
                                                spaceSrvK = Matrix.concatRows(spaceSrvK, Matrix.extractRows(spaceSrv, i, i + 1, null), null);
                                            }
                                        }
                                    }

                                    Matrix spaceBufK = new Matrix(0, 0);
                                    for (int i = 0; i < en.getNumRows(); i++) {
                                        if (en.get(i, 0) == 1) {
                                            if (spaceBufK.isEmpty()) {
                                                spaceBufK = Matrix.extractRows(spaceBuf, i, i + 1, null);
                                            } else {
                                                spaceBufK = Matrix.concatRows(spaceBufK, Matrix.extractRows(spaceBuf, i, i + 1, null), null);
                                            }
                                        }
                                    }

                                    Matrix spaceVarK = new Matrix(0, 0);
                                    for (int i = 0; i < en.getNumRows(); i++) {
                                        if (en.get(i, 0) == 1) {
                                            if (spaceVarK.isEmpty()) {
                                                spaceVarK = Matrix.extractRows(spaceVar, i, i + 1, null);
                                            } else {
                                                spaceVarK = Matrix.concatRows(spaceVarK, Matrix.extractRows(spaceVar, i, i + 1, null), null);
                                            }
                                        }
                                    }

                                    // markov-modulated case
                                    // MATLAB: space_var_k(sum(sn.nvars(ind,1:class))) = kdest
                                    // extracts columns 1 to class (1-based)
                                    // Java: columns 0 to jobClass (0-based), same logical data
                                    if (ismkvmodclass.get(jobClass) != 0) {
                                        int nvarsSum = 0;
                                        for (int i = 0; i <= jobClass; i++) {
                                            nvarsSum += (int) sn.nvars.get(ind, i);
                                        }
                                        // MATLAB uses 1-based indexing, Java needs 0-based for column
                                        int spaceVarCol = nvarsSum - 1;
                                        for (int row = 0; row < spaceVarK.getNumRows(); row++) {
                                            // kdest is 0-based in Java, but MAP output var values
                                            // in the state space are 1-based (initDefault sets to 1),
                                            // so store kdest+1 to match MATLAB convention
                                            spaceVarK.set(row, spaceVarCol, kdest + 1);
                                        }
                                    }

                                    for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                        spaceSrvK.set(row, (int) (Ks.get(jobClass) + k), spaceSrvK.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                    }
                                    for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                        spaceSrvK.set(row, (int) (Ks.get(jobClass) + kdest), spaceSrvK.get(row, (int) (Ks.get(jobClass) + kdest)) + 1);
                                    }

                                    switch (sn.sched.get(sn.stations.get(ist))) {
                                        case EXT:
                                            rate.set(0, 0, (int) proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest));
                                            break;
                                        case INF:
                                            double proc_value_inf = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value_inf = kir.get(k).get(jobClass);
                                            rate.set(0, 0, proc_value_inf * kir_value_inf);
                                            break;
                                        case PS:
                                        case LPS:
                                            double proc_value_ps = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value_ps = kir.get(k).get(jobClass);
                                            Matrix numerator = new Matrix(1, 1);
                                            numerator.set(0, 0, proc_value_ps * kir_value_ps);
                                            Matrix denom = new Matrix(1, 1);
                                            double ni_value = ni.get(0);
                                            denom.set(0, 0, ni_value * Maths.min(ni_value, S.get(ist)));
                                            rate = numerator.elementDivide(denom);
                                            break;
                                        case PSPRIO:
                                        case DPSPRIO:
                                        case GPSPRIO: {
                                            // Find minimum priority among present classes
                                            int minPrioPh = Integer.MAX_VALUE;
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (nir.get(0, r) > 0) {
                                                    int rPrio = (int) sn.classprio.get(r);
                                                    if (rPrio < minPrioPh) {
                                                        minPrioPh = rPrio;
                                                    }
                                                }
                                            }
                                            double ni_value_prio = ni.get(0);
                                            // If ni <= S (all jobs get service) or this class has highest priority
                                            if (ni_value_prio <= S.get(ist) || (int) sn.classprio.get(jobClass) == minPrioPh) {
                                                double proc_value_psprio = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                                double kir_value_psprio = kir.get(k).get(jobClass);
                                                Matrix numerator_prio = new Matrix(1, 1);
                                                numerator_prio.set(0, 0, proc_value_psprio * kir_value_psprio);
                                                Matrix denom_prio = new Matrix(1, 1);
                                                denom_prio.set(0, 0, ni_value_prio * Maths.min(ni_value_prio, S.get(ist)));
                                                rate = numerator_prio.elementDivide(denom_prio);
                                            } else {
                                                // Not in highest priority group - no service
                                                rate.set(0, 0, 0.0);
                                            }
                                            break;
                                        }
                                        case DPS:
                                            if (S.get(ist) > 1) {
                                                InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                                }), "Multi-server DPS not supported yet");
                                            }

                                            Matrix wDps = sn.schedparam.getRow(ist);
                                            wDps.scaleEq(1.0 / wDps.elementSum());

                                            Matrix wDpsRepMatSum = new Matrix(0, 0);
                                            for (int row = 0; row < nir.getNumRows(); row++) {
                                                if (wDpsRepMatSum.isEmpty()) {
                                                    wDpsRepMatSum = new Matrix(1, 1);
                                                    wDpsRepMatSum.set(0, 0, wDps.mult(nir.transpose()).sumRows().get(0, 0));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, wDps.mult(nir.transpose()).sumRows().get(0, 0));
                                                    wDpsRepMatSum = Matrix.concatRows(wDpsRepMatSum, new_elem, null);
                                                }
                                            }

                                            double proc_value_dps = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value_dps = kir.get(k).get(jobClass);
                                            rate.set(0, 0, proc_value_dps * kir_value_dps * wDps.get(jobClass) / wDpsRepMatSum.get(0));
                                            break;
                                        case GPS:
                                            if (S.get(ist) > 1) {
                                                InputOutputKt.line_error(InputOutputKt.mfilename(new Object() {
                                                }), "Multi-server GPS not supported yet");
                                            }

                                            Matrix cirGps = new Matrix(nir.getNumRows(), nir.getNumCols());
                                            for (int row = 0; row < cirGps.getNumRows(); row++) {
                                                for (int col = 0; col < cirGps.getNumCols(); col++) {
                                                    if (nir.get(row, col) < 1.0) {
                                                        cirGps.set(row, col, nir.get(row, col));
                                                    } else {
                                                        cirGps.set(row, col, 1.0);
                                                    }
                                                }
                                            }

                                            Matrix wGps = sn.schedparam.getRow(ist);
                                            wGps.scaleEq(1.0 / wGps.elementSum());

                                            double proc_value_gps = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value_gps = kir.get(k).get(jobClass);
                                            double wcirGps = wGps.mult(cirGps).get(0);
                                            rate.set(0, 0, proc_value_gps * kir_value_gps / nir.get(jobClass) * wGps.get(jobClass) / wcirGps);
                                            break;
                                        case FCFS:
                                        case HOL:
                                        case LCFS:
                                        case LCFSPR:
                                        case LCFSPI:
                                        case SIRO:
                                        case SEPT:
                                        case LEPT:
                                            double proc_value = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value = kir.get(k).get(jobClass);
                                            rate.set(0, 0, proc_value * kir_value);
                                            break;
                                    }

                                    // if class cannot be served locally, rate = NaN since mu{i, class} = NaN
                                    if (ni.hasInfinite()) {
                                        // hit limited load-dependence
                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                        Matrix outrate_bottom = Matrix.scaleMult(rate, cdscalingIst * lld);
                                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                    } else {
                                        double cdscalingIst = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                        double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                        Matrix outrate_bottom = Matrix.scaleMult(rate, cdscalingIst * lld);
                                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                    }
                                    Matrix outspace_bottom_left = Matrix.concatColumns(spaceBufK, spaceSrvK, null);
                                    Matrix bottom = Matrix.concatColumns(outspace_bottom_left, spaceVarK, null);
                                    outspace = Matrix.concatRows(outspace, bottom, null);
                                    Matrix outprob_bottom = new Matrix(rate.getNumRows(), 1);
                                    outprob_bottom.ones();
                                    outprob = Matrix.concatRows(outprob, outprob_bottom, null);
                                }
                            }
                        }
                        
                        // ADD MAP/MMPP2 phase transition logic
                        if (ismkvmodclass.get(jobClass, 0) == 1) {
                            // Process MAP phase transitions for departure events
                            Matrix origSpaceBuf = spaceBuf.copy();
                            Matrix origSpaceSrv = spaceSrv.copy();
                            Matrix origSpaceVar = spaceVar.copy();
                            
                            for (int mapK = 0; mapK < K.get(jobClass); mapK++) {
                                // Find states where this phase has jobs
                                Matrix mapEn = new Matrix(origSpaceSrv.getNumRows(), 1);
                                boolean hasJobs = false;
                                for (int row = 0; row < origSpaceSrv.getNumRows(); row++) {
                                    if (origSpaceSrv.get(row, (int)(Ks.get(jobClass) + mapK)) > 0) {
                                        mapEn.set(row, 0, 1);
                                        hasJobs = true;
                                    } else {
                                        mapEn.set(row, 0, 0);
                                    }
                                }
                                
                                if (hasJobs) {
                                    // Process phase transitions to all other phases
                                    for (int mapKdest = 0; mapKdest < K.get(jobClass); mapKdest++) {
                                        if (mapKdest != mapK) { // Don't transition to same phase
                                            Matrix spaceBufK = origSpaceBuf.copy();
                                            Matrix spaceSrvK = origSpaceSrv.copy();
                                            Matrix spaceVarK = origSpaceVar.copy();
                                            
                                            // Update local variables for MAP state
                                            // MATLAB: sn.nvars(ind,1:class) extracts columns 1 to class (1-based)
                                            // Java: columns 0 to jobClass (0-based), same logical data
                                            int nvarsSum = 0;
                                            for (int r = 0; r <= jobClass; r++) {
                                                nvarsSum += (int) sn.nvars.get(ind, r);
                                            }
                                            // MATLAB uses 1-based indexing, Java needs 0-based for column
                                            int spaceVarCol = nvarsSum - 1;
                                            if (spaceVarCol >= 0 && spaceVarCol < spaceVarK.getNumCols()) {
                                                for (int row = 0; row < spaceVarK.getNumRows(); row++) {
                                                    if (mapEn.get(row, 0) == 1) {
                                                        // mapKdest is 0-based in Java, but MAP output var values
                                                        // in the state space are 1-based (initDefault sets to 1),
                                                        // so store mapKdest+1 to match MATLAB convention
                                                        spaceVarK.set(row, spaceVarCol, mapKdest + 1);
                                                    }
                                                }
                                            }
                                            
                                            // Move job from phase mapK to phase mapKdest
                                            for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                                if (mapEn.get(row, 0) == 1) {
                                                    spaceSrvK.set(row, (int)(Ks.get(jobClass) + mapK), 
                                                                spaceSrvK.get(row, (int)(Ks.get(jobClass) + mapK)) - 1);
                                                    spaceSrvK.set(row, (int)(Ks.get(jobClass) + mapKdest),
                                                                spaceSrvK.get(row, (int)(Ks.get(jobClass) + mapKdest)) + 1);
                                                }
                                            }
                                            
                                            // Calculate transition rate based on scheduling strategy
                                            Matrix mapRate = new Matrix(mapEn.getNumRows(), 1);
                                            mapRate.zero();
                                            
                                            switch (sn.sched.get(sn.stations.get(ist))) {
                                                case EXT:
                                                    double procRate = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(mapK, mapKdest);
                                                    for (int row = 0; row < mapRate.getNumRows(); row++) {
                                                        if (mapEn.get(row, 0) == 1) {
                                                            mapRate.set(row, 0, procRate);
                                                        }
                                                    }
                                                    break;
                                                case INF:
                                                    double procRateInf = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(mapK, mapKdest);
                                                    for (int row = 0; row < mapRate.getNumRows(); row++) {
                                                        if (mapEn.get(row, 0) == 1 && mapK < kir.size()) {
                                                            double kirValue = kir.get(mapK).get(row, jobClass);
                                                            mapRate.set(row, 0, procRateInf * kirValue);
                                                        }
                                                    }
                                                    break;
                                                case PS:
                                                case LPS:
                                                    double procRatePs = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(mapK, mapKdest);
                                                    for (int row = 0; row < mapRate.getNumRows(); row++) {
                                                        if (mapEn.get(row, 0) == 1 && mapK < kir.size()) {
                                                            double kirValue = kir.get(mapK).get(row, jobClass);
                                                            // Bounds check for ni
                                                            double niValue = (row < ni.getNumRows()) ? ni.get(row, 0) : ni.get(0, 0);
                                                            double servers = Math.min(niValue, S.get(ist, 0));
                                                            if (niValue > 0) {
                                                                mapRate.set(row, 0, procRatePs * kirValue / niValue * servers);
                                                            }
                                                        }
                                                    }
                                                    break;
                                                case PSPRIO:
                                                case DPSPRIO:
                                                case GPSPRIO:
                                                    double procRatePsprio = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(mapK, mapKdest);
                                                    for (int row = 0; row < mapRate.getNumRows(); row++) {
                                                        if (mapEn.get(row, 0) == 1 && mapK < kir.size()) {
                                                            // Find minimum priority among present classes for this state
                                                            int minPrioMap = Integer.MAX_VALUE;
                                                            for (int r = 0; r < sn.nclasses; r++) {
                                                                if (nir.get(row, r) > 0) {
                                                                    int rPrio = (int) sn.classprio.get(r);
                                                                    if (rPrio < minPrioMap) {
                                                                        minPrioMap = rPrio;
                                                                    }
                                                                }
                                                            }
                                                            double kirValue = kir.get(mapK).get(row, jobClass);
                                                            double niValue = (row < ni.getNumRows()) ? ni.get(row, 0) : ni.get(0, 0);
                                                            double servers = Math.min(niValue, S.get(ist, 0));
                                                            // If ni <= S (all jobs get service) or this class has highest priority
                                                            if (niValue <= S.get(ist, 0) || (int) sn.classprio.get(jobClass) == minPrioMap) {
                                                                if (niValue > 0) {
                                                                    mapRate.set(row, 0, procRatePsprio * kirValue / niValue * servers);
                                                                }
                                                            } else {
                                                                // Not in highest priority group - no service
                                                                mapRate.set(row, 0, 0.0);
                                                            }
                                                        }
                                                    }
                                                    break;
                                                default:
                                                    // For other scheduling strategies (FCFS, LCFS, SIRO, etc.)
                                                    double procRateDefault = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(mapK, mapKdest);
                                                    for (int row = 0; row < mapRate.getNumRows(); row++) {
                                                        if (mapEn.get(row, 0) == 1 && mapK < kir.size()) {
                                                            double kirValue = kir.get(mapK).get(row, jobClass);
                                                            mapRate.set(row, 0, procRateDefault * kirValue);
                                                        }
                                                    }
                                                    break;
                                            }
                                            
                                            // Add enabled states to output
                                            for (int row = 0; row < mapEn.getNumRows(); row++) {
                                                if (mapEn.get(row, 0) == 1 && mapRate.get(row, 0) > 0) {
                                                    Matrix rowSpaceBuf = Matrix.extractRows(spaceBufK, row, row + 1, null);
                                                    Matrix rowSpaceSrv = Matrix.extractRows(spaceSrvK, row, row + 1, null);
                                                    Matrix rowSpaceVar = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                                    
                                                    Matrix newState = Matrix.concatColumns(rowSpaceBuf, rowSpaceSrv, null);
                                                    newState = Matrix.concatColumns(newState, rowSpaceVar, null);
                                                    outspace = Matrix.concatRows(outspace, newState, null);
                                                    
                                                    // Apply load and class dependent scaling
                                                    double cdscalingVal = cdscaling.get(sn.stations.get(ist)).apply(nir);
                                                    // Bounds check: ensure row is valid for ni matrix
                                                    double niValue = 0.0;
                                                    if (row < ni.getNumRows()) {
                                                        niValue = ni.get(row, 0);
                                                    } else {
                                                        // Use first row value as fallback and log warning
                                                        niValue = ni.get(0, 0);
                                                        line_warning("AfterEventStation", "MAP phase transition bounds mismatch - row=%d >= ni.getNumRows()=%d, mapEn.getNumRows()=%d, inspace.getNumRows()=%d", row, ni.getNumRows(), mapEn.getNumRows(), inspace.getNumRows());
                                                    }
                                                    double lldVal;
                                                    if (ni.hasInfinite()) {
                                                        lldVal = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    } else {
                                                        int lldCol = (int) Math.min(niValue, lldlimit);
                                                        // Ensure column is within bounds
                                                        lldCol = Math.max(0, Math.min(lldCol, lldscaling.getNumCols() - 1));
                                                        lldVal = lldscaling.get(ist, lldCol);
                                                    }
                                                    
                                                    Matrix rateRow = new Matrix(1, 1);
                                                    rateRow.set(0, 0, cdscalingVal * lldVal * mapRate.get(row, 0));
                                                    outrate = Matrix.concatRows(outrate, rateRow, null);
                                                    
                                                    Matrix probRow = new Matrix(1, 1);
                                                    probRow.set(0, 0, 1.0);
                                                    outprob = Matrix.concatRows(outprob, probRow, null);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    Ret.EventResult result_p = new Ret.EventResult(outspace, outrate, outprob);
                    eventCache.put(key, result_p);
                    if (isSimulation) {
                        if (outspace.getNumRows() > 1) {
                            Matrix tot_rate = outrate.sumCols();
                            Matrix cum_sum = outrate.cumsumViaCol();
                            Matrix cum_rate = Matrix.scaleMult(cum_sum, 1.0 / tot_rate.value());
                            int firing_ctr = -1;
                            double rand = Maths.rand();
                            // we need the indicies where rand is bigger than cum_prob
                            for (int row = 0; row < cum_rate.getNumRows(); row++) {
                                if (rand > cum_rate.get(row)) {
                                    firing_ctr = row;
                                }
                            }
                            firing_ctr++;
                            outspace = Matrix.extractRows(outspace, firing_ctr, firing_ctr + 1, null);
                            double outrate_val = outrate.elementSum();
                            outrate = new Matrix(1, 1);
                            outrate.set(0, 0, outrate_val);
                            outprob = Matrix.extractRows(outprob, firing_ctr, firing_ctr + 1, null);
                        }

                    }

                } else {
                    Ret.EventResult result_p = new Ret.EventResult(outspace, outrate, outprob);
                    eventCache.put(key, result_p);
                }
                break;
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}