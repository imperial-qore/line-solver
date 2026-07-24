package jline.lang.state;

import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.GlobalConstants;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Station;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class AfterEventCache implements Serializable {
    static Ret.EventResult afterEventCache(NetworkStruct sn, int ind, EventType event, int jobClass, boolean isSimulation,
                                           Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                           int M, int R,
                                           int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                           double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        // job arrives in class, then reads and moves into hit or miss class, then departs
        switch (event) {
            case ARV:
                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                    spaceSrv.set(row, jobClass, spaceSrv.get(row, jobClass) + 1);
                }
                // buf is empty
                outspace = Matrix.concatColumns(spaceSrv, spaceVar, null);
                // passive action, rate is unspecified
                outrate = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                outrate.ones();
                outrate.scaleEq(-1);
                break;
            case DEP:
                if (spaceSrv.get(jobClass) > 0) {
                    Set<Integer> retrievalClassIndices = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).retrievalClassIndices;
                    List<Double> p = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).pread.get(jobClass);
                    int totalCacheCapacity = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).totalCacheCapacity;

                    Matrix var = spaceVar.copy();
                    // A retrieval-class job departs the cache only to BEGIN a retrieval (cache -> queue).
                    if (retrievalClassIndices.contains(jobClass)) {
                        // Retrieval classes have one-hot read access representing the item the retrieval belongs to
                        int item = p.indexOf(1.0);
                        if (isInRetrievalSystem(var, item, totalCacheCapacity)) {
                            break;
                        }
                        // Beginning a retrieval: record the item as in-flight as the job leaves for the queue.
                        addToRetrievalSystem(var, item, totalCacheCapacity);
                    }

                    for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        spaceSrv.set(row, jobClass, spaceSrv.get(row, jobClass) - 1);
                    }

                    switch (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass))) {
                        case RROBIN:
                            // MATLAB: sn.nvars(ind,1:(R+class)) extracts columns 1 to R+class (1-based)
                            // Java: columns 0 to R+jobClass (0-based), same logical data
                            int nvarCols = R + jobClass + 1;
                            Matrix nvar_ind = new Matrix(1, nvarCols);
                            Matrix.extract(sn.nvars, ind, ind + 1, 0, nvarCols, nvar_ind, 0, 0);
                            int nvar_sum = (int) nvar_ind.elementSum();
                            // MATLAB uses 1-based indexing, Java needs 0-based
                            int spaceVarIdx = nvar_sum - 1;
                            int idx = -1;
                            Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(jobClass));
                            int numOutlinks = (int) outlinks.length();
                            for (int row = 0; row < numOutlinks; row++) {
                                if (spaceVar.get(spaceVarIdx) == outlinks.get(row)) {
                                    idx = row;
                                    break;
                                }
                            }
                            if (idx >= 0 && idx < numOutlinks - 1) {
                                spaceVar.set(spaceVarIdx, outlinks.get(idx + 1));
                            } else {
                                spaceVar.set(spaceVarIdx, outlinks.get(0));
                            }
                            break;
                    }
                    // buf is empty
                    outspace = Matrix.concatColumns(spaceSrv, var, null);
                    // immediate action
                    outrate = new Matrix(outspace.getNumRows(), outspace.getNumRows());
                    outrate.ones();
                    outrate.scaleEq(GlobalConstants.Immediate);
                    break;
                }
                break;
            case READ:
                int n = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).nitems;
                Matrix m = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).itemcap;
                int totalCacheCapacity = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).totalCacheCapacity;
                Matrix[][] ac = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).accost;
                Matrix hitclass = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).hitclass;
                Matrix missclass = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).missclass;
                Matrix retrievalClasses = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).retrievalClasses;
                int h = m.getNumCols();
                Set<Integer> retrievalClassIndices = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).retrievalClassIndices;
                ReplacementStrategy replacement = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).replacestrat;
                double qadm = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).qlru;

                if (spaceSrv.sumCols(jobClass) > 0 && (int) spaceSrv.elementSum() == 1) {
                    List<Double> p = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).pread.get(jobClass);
                    Matrix spaceSrvK = new Matrix(0, 0);
                    Matrix spaceVarK = new Matrix(0, 0);
                    outrate = new Matrix(0, 0);
                    Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                    en.zero();
                    boolean any_en = false;
                    for (int row = 0; row < en.getNumRows(); row++) {
                        if (spaceSrv.get(row, jobClass) > 0) {
                            en.set(row, 0, 1);
                            any_en = true;
                        }
                    }

                    if (any_en) {
                        // for e=find(en)'
                        for (int e = 0; e < en.getNumRows(); e++) {
                            if (en.get(e) == 1) {
                                int kset = -1;
                                int kend = -1;
                                int l = -1;
                                // If a retrieval is returning, it can only read one value
                                if (isSimulation || retrievalClassIndices.contains(jobClass)) {
                                    // pick one item
                                    List<Double> pcumsum = new ArrayList<Double>();
                                    pcumsum.add(p.get(0));
                                    for (int i = 1; i < p.size(); i++) {
                                        pcumsum.add(p.get(i) + pcumsum.get(i - 1));
                                    }
                                    double rand = Maths.rand();
                                    for (int row = 0; row < pcumsum.size(); row++) {
                                        if (rand > pcumsum.get(row)) {
                                            kset = row;
                                        }
                                    }
                                    kset++;
                                    kend = kset + 1;
                                    // pick one entry list
                                    Matrix accumsum = ac[jobClass][kset].getRow(0);
                                    accumsum = accumsum.cumsumViaRow();
                                    rand = Maths.rand();
                                    for (int col = 0; col < accumsum.getNumCols(); col++) {
                                        if (rand > accumsum.get(col)) {
                                            l = col;
                                        }
                                    }
                                    l++;
                                } else {
                                    // Check if varsparam is set for this node (used in testing)
                                    if (sn.varsparam != null && sn.varsparam.get(ind, 0) >= 0) {
                                        // Use the specific item specified in varsparam
                                        kset = (int) sn.varsparam.get(ind, 0);
                                        kend = kset + 1;
                                    } else {
                                        kset = 0;
                                        kend = n;
                                    }
                                }

                                // request to item k
                                for (int k = kset; k < kend; k++) {
                                    Matrix spaceSrvE = spaceSrv.getRow(e);
                                    spaceSrvE.set(jobClass, spaceSrvE.get(jobClass) - 1);
                                    Matrix var = spaceVar.getRow(e);

                                    int posk = -1;
                                    for (int col = 0; col < totalCacheCapacity; col++) {
                                        if (var.get(col) == k + 1) {
                                            posk = col;
                                            break;
                                        }
                                    }

                                    // CACHE MISS, either begin retrieval or retrieved item can enter any list based on
                                    // accessCost
                                    if (posk == -1) {
                                        Matrix retrievalClass = retrievalClasses.getRow(k);

                                        // This job should not continue, possible only as a consequence of the event
                                        // loop
                                        if (retrievalClassIndices.contains(jobClass) &&
                                                !isInRetrievalSystem(var, k, totalCacheCapacity)) {
                                            continue;
                                        }

                                        // If the item is not returning from retrieval, it must be retrieved or there is
                                        // a delayed hit, only begin a retrieval if the job class can switch to a
                                        // retrieval pending class for item k
                                        if (!retrievalClassIndices.contains(jobClass) &&
                                                jobClass < retrievalClasses.getNumCols() &&
                                                retrievalClass.get(jobClass) != -1) {
                                            // If retrieval has not started for the item, then begin it
                                            if (!isInRetrievalSystem(var, k, totalCacheCapacity)) {
                                                spaceSrvE.set((int) retrievalClass.get(jobClass), spaceSrvE.get((int) retrievalClass.get(jobClass)) + 1);
                                            }

                                            if (spaceSrvK.isEmpty()) {
                                                spaceSrvK = spaceSrvE.copy();
                                            } else {
                                                spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                            }
                                            if (spaceVarK.isEmpty()) {
                                                spaceVarK = var.copy();
                                            } else {
                                                spaceVarK = Matrix.concatRows(spaceVarK, var, null);
                                            }

                                            if (isSimulation) {
                                                if (outrate.isEmpty()) {
                                                    outrate = new Matrix(1, 1);
                                                    outrate.set(0, 0, GlobalConstants.Immediate);
                                                } else {
                                                    Matrix bottom_row = new Matrix(1, 1);
                                                    bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                    outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                }
                                                outprob = Matrix.concatRows(outprob, Matrix.singleton(p.get(k)), null);
                                                continue;
                                            }

                                            if (outrate.isEmpty()) {
                                                outrate = new Matrix(1, 1);
                                                outrate.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                            } else {
                                                Matrix bottom_row = new Matrix(1, 1);
                                                bottom_row.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                                outrate = Matrix.concatRows(outrate, bottom_row, null);
                                            }
                                            continue;
                                        }
                                        // Item has now been retrieved and can be marked as a miss
                                        spaceSrvE.set((int) missclass.get(jobClass), spaceSrvE.get((int) missclass.get(jobClass)) + 1);
                                        removeFromRetrievalSystem(var, k, totalCacheCapacity);
                                        Matrix varp;
                                        switch (replacement) {
                                            case FIFO:
                                            case LRU:
                                            case SFIFO:
                                            case HLRU:
                                                if (isSimulation) {
                                                    int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                    int headPos = cpos(m, listidx, 0);
                                                    int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                    varp = var.copy();
                                                    // Shift items in list listidx to make room at the head
                                                    if (m.get(listidx) > 1) {
                                                        Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
                                                    }
                                                    varp.set(headPos, k + 1); // head of list listidx
                                                    if (spaceSrvK.isEmpty()) {
                                                        spaceSrvK = spaceSrvE.copy();
                                                    } else {
                                                        spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    }
                                                    if (spaceVarK.isEmpty()) {
                                                        spaceVarK = varp.copy();
                                                    } else {
                                                        spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                    }
                                                    // no p(k) weighting since that goes in the outprob vec
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob.set(0, 0, ac[jobClass][kset].get(0, l) * p.get(kset));
                                                } else {
                                                    // Cache reject (column 0): pass through without caching
                                                    if (ac[jobClass][k].get(0, 0) > 0) {
                                                        varp = var.copy();
                                                        if (spaceSrvK.isEmpty()) {
                                                            spaceSrvK = spaceSrvE.copy();
                                                        } else {
                                                            spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        }
                                                        if (spaceVarK.isEmpty()) {
                                                            spaceVarK = varp.copy();
                                                        } else {
                                                            spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                        }
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, ac[jobClass][k].get(0, 0) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, ac[jobClass][k].get(0, 0) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                        outprob = Matrix.concatRows(outprob, Matrix.singleton(1.0), null);
                                                    }
                                                    // Iterate over all possible target lists (columns 1 to h)
                                                    for (l = 1; l <= h; l++) {
                                                        int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                        int headPos = cpos(m, listidx, 0);
                                                        int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                        varp = var.copy();
                                                        // Shift items in list listidx to make room at the head
                                                        if (m.get(listidx) > 1) {
                                                            Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
                                                        }
                                                        varp.set(headPos, k + 1); // head of list listidx
                                                        if (spaceSrvK.isEmpty()) {
                                                            spaceSrvK = spaceSrvE.copy();
                                                        } else {
                                                            spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        }
                                                        if (spaceVarK.isEmpty()) {
                                                            spaceVarK = varp.copy();
                                                        } else {
                                                            spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                        }
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, ac[jobClass][k].get(0, l) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, ac[jobClass][k].get(0, l) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                        outprob = Matrix.concatRows(outprob, Matrix.singleton(1.0), null);
                                                    }
                                                }
                                                break;
                                            case RR:
                                                if (isSimulation) {
                                                    int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                    int headPos = cpos(m, listidx, 0);
                                                    varp = var.copy();
                                                    // randi(m(listidx),1,1)
                                                    int r = (int) (Maths.rand() * m.get(listidx));
                                                    varp.set(headPos + r, k + 1);
                                                    if (spaceSrvK.isEmpty()) {
                                                        spaceSrvK = spaceSrvE.copy();
                                                    } else {
                                                        spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    }
                                                    if (spaceVarK.isEmpty()) {
                                                        spaceVarK = varp.copy();
                                                    } else {
                                                        spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                    }
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob.set(0, 0, ac[jobClass][kset].get(0, l) * p.get(kset));
                                                } else {
                                                    // Cache reject (column 0): pass through without caching
                                                    if (ac[jobClass][k].get(0, 0) > 0) {
                                                        varp = var.copy();
                                                        if (spaceSrvK.isEmpty()) {
                                                            spaceSrvK = spaceSrvE.copy();
                                                        } else {
                                                            spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        }
                                                        if (spaceVarK.isEmpty()) {
                                                            spaceVarK = varp.copy();
                                                        } else {
                                                            spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                        }
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, ac[jobClass][k].get(0, 0) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, ac[jobClass][k].get(0, 0) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                    }
                                                    // Iterate over all possible target lists
                                                    for (l = 1; l <= h; l++) {
                                                        int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                        int headPos = cpos(m, listidx, 0);
                                                        // random position in list listidx
                                                        for (int r = 0; r < (int) m.get(listidx); r++) {
                                                            varp = var.copy();
                                                            varp.set(headPos + r, k + 1);
                                                            if (spaceSrvK.isEmpty()) {
                                                                spaceSrvK = spaceSrvE.copy();
                                                            } else {
                                                                spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                            }
                                                            if (spaceVarK.isEmpty()) {
                                                                spaceVarK = varp.copy();
                                                            } else {
                                                                spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                            }
                                                            if (outrate.isEmpty()) {
                                                                outrate = new Matrix(1, 1);
                                                                outrate.set(0, 0, ac[jobClass][k].get(0, l) * p.get(k) / m.get(listidx) * GlobalConstants.Immediate);
                                                            } else {
                                                                Matrix bottom_row = new Matrix(1, 1);
                                                                bottom_row.set(0, 0, ac[jobClass][k].get(0, l) * p.get(k) / m.get(listidx) * GlobalConstants.Immediate);
                                                                outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                            }
                                                        }
                                                    }
                                                }
                                                break;
                                            case QLRU:
                                                // q-LRU: admit (LRU head insert) with probability q, else pass through.
                                                if (isSimulation) {
                                                    int listidx = l - 1;
                                                    varp = var.copy();
                                                    if (listidx >= 0 && Maths.rand() <= qadm) {
                                                        int headPos = cpos(m, listidx, 0);
                                                        int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                        if (m.get(listidx) > 1) {
                                                            Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
                                                        }
                                                        varp.set(headPos, k + 1);
                                                    }
                                                    spaceSrvK = spaceSrvK.isEmpty() ? spaceSrvE.copy() : Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    spaceVarK = spaceVarK.isEmpty() ? varp.copy() : Matrix.concatRows(spaceVarK, varp, null);
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob.set(0, 0, ac[jobClass][kset].get(0, l) * p.get(kset));
                                                } else {
                                                    double rejw = ac[jobClass][k].get(0, 0) + (1 - qadm) * (1 - ac[jobClass][k].get(0, 0));
                                                    if (rejw > 0) {
                                                        varp = var.copy();
                                                        spaceSrvK = spaceSrvK.isEmpty() ? spaceSrvE.copy() : Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        spaceVarK = spaceVarK.isEmpty() ? varp.copy() : Matrix.concatRows(spaceVarK, varp, null);
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, rejw * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, rejw * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                        outprob = Matrix.concatRows(outprob, Matrix.singleton(1.0), null);
                                                    }
                                                    for (l = 1; l <= h; l++) {
                                                        int listidx = l - 1;
                                                        int headPos = cpos(m, listidx, 0);
                                                        int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                        varp = var.copy();
                                                        if (m.get(listidx) > 1) {
                                                            Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
                                                        }
                                                        varp.set(headPos, k + 1);
                                                        spaceSrvK = spaceSrvK.isEmpty() ? spaceSrvE.copy() : Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        spaceVarK = spaceVarK.isEmpty() ? varp.copy() : Matrix.concatRows(spaceVarK, varp, null);
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, qadm * ac[jobClass][k].get(0, l) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, qadm * ac[jobClass][k].get(0, l) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                        outprob = Matrix.concatRows(outprob, Matrix.singleton(1.0), null);
                                                    }
                                                }
                                                break;
                                        }
                                    }
                                    // CACHE HIT in list i < h, move to list i+1
                                    else if (posk < (m.elementSum() - m.get(h - 1))) {
                                        spaceSrvE.set((int) hitclass.get(jobClass), spaceSrvE.get((int) hitclass.get(jobClass)) + 1);
                                        // i = min(find(posk <= cumsum(m)));
                                        int i = -1;
                                        Matrix mcumsum = m.cumsumViaRow();
                                        for (int col = 0; col < mcumsum.getNumCols(); col++) {
                                            if (posk < mcumsum.get(col)) {
                                                i = col;
                                                // quit for loop once get the first index that met the condition
                                                break;
                                            }
                                        }
                                        // j = posk - sum(m(1:i-1));
                                        int j = posk;
                                        for (int m_ind = 0; m_ind < i; m_ind++) {
                                            j -= m.get(m_ind);
                                        }

                                        switch (replacement) {
                                            case FIFO:
                                                if (isSimulation) {
                                                    Matrix varp = var.copy();

                                                    // probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))
                                                    Matrix aci = new Matrix(1, ac[jobClass][k].getNumCols() - (i + 1));
                                                    Matrix.extract(ac[jobClass][k], i + 1, i + 2,
                                                            i + 1, ac[jobClass][k].getNumCols(),
                                                            aci, 0, 0);
                                                    aci = aci.scale(1 / aci.elementSum());
                                                    int probchoose = Maths.probchoose(aci);
                                                    int inew = i + probchoose;

                                                    if (inew != i) {
                                                        varp.set(cpos(m, i, j), var.get(cpos(m, inew, (int) m.get(inew) - 1)));
                                                        Matrix.extract(var, 0, 1, cpos(m, inew, 0), cpos(m, inew, (int) m.get(inew) - 1), varp, 0, cpos(m, inew, 1));
                                                        varp.set(cpos(m, inew, 0), k + 1);
                                                    }

                                                    if (spaceSrvK.isEmpty()) {
                                                        spaceSrvK = spaceSrvE.copy();
                                                    } else {
                                                        spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    }
                                                    if (spaceVarK.isEmpty()) {
                                                        spaceVarK = varp.copy();
                                                    } else {
                                                        spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                    }
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                        outprob = Matrix.concatRows(outprob, Matrix.singleton(ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k).doubleValue()), null);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                } else {
                                                    // MATLAB: for inew = i:h (h = numel(m), 1-indexed)
                                                    // In Java (0-indexed): inew goes from i to h-1
                                                    for (int inew = i; inew < h; inew++) {
                                                        Matrix varp = var.copy();
                                                        varp.set(cpos(m, i, j), var.get(cpos(m, inew, (int) m.get(inew) - 1)));
                                                        if (m.get(inew) > 1) {
                                                            Matrix.extract(var, 0, 1, cpos(m, inew, 0), cpos(m, inew, (int) m.get(inew) - 1), varp, 0, cpos(m, inew, 1));
                                                        }
                                                        varp.set(cpos(m, inew, 0), k + 1);

                                                        if (spaceSrvK.isEmpty()) {
                                                            spaceSrvK = spaceSrvE.copy();
                                                        } else {
                                                            spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        }
                                                        if (spaceVarK.isEmpty()) {
                                                            spaceVarK = varp.copy();
                                                        } else {
                                                            spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                        }
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                    }
                                                }
                                                break;
                                            case RR:
                                                if (isSimulation) {
                                                    Matrix varp = var.copy();

                                                    // probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))
                                                    Matrix aci = new Matrix(1, ac[jobClass][k].getNumCols() - (i + 1));
                                                    Matrix.extract(ac[jobClass][k], i + 1, i + 2,
                                                            i + 1, ac[jobClass][k].getNumCols(),
                                                            aci, 0, 0);
                                                    aci = aci.scale(1 / aci.elementSum());
                                                    int probchoose = Maths.probchoose(aci);
                                                    int inew = i + probchoose;

                                                    int r = (int) (Maths.rand() * m.get(inew));
                                                    varp.set(cpos(m, i, j), var.get(cpos(m, inew, r)));
                                                    varp.set(cpos(m, inew, r), k + 1);
                                                    outprob = Matrix.concatRows(outprob, Matrix.singleton(ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k).doubleValue() / m.get(inew)), null);

                                                    if (spaceSrvK.isEmpty()) {
                                                        spaceSrvK = spaceSrvE.copy();
                                                    } else {
                                                        spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    }
                                                    if (spaceVarK.isEmpty()) {
                                                        spaceVarK = varp.copy();
                                                    } else {
                                                        spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                    }
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                } else {
                                                    for (int inew = i; inew < h; inew++) {
                                                        for (int r = 0; r < (int) m.get(inew); r++) {
                                                            Matrix varp = var.copy();
                                                            varp.set(cpos(m, i, j), var.get(cpos(m, inew, r)));
                                                            varp.set(cpos(m, inew, r), k + 1);

                                                            if (spaceSrvK.isEmpty()) {
                                                                spaceSrvK = spaceSrvE.copy();
                                                            } else {
                                                                spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                            }
                                                            if (spaceVarK.isEmpty()) {
                                                                spaceVarK = varp.copy();
                                                            } else {
                                                                spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                            }
                                                            if (outrate.isEmpty()) {
                                                                outrate = new Matrix(1, 1);
                                                                outrate.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) / m.get(inew) * GlobalConstants.Immediate);
                                                            } else {
                                                                Matrix bottom_row = new Matrix(1, 1);
                                                                bottom_row.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) / m.get(inew) * GlobalConstants.Immediate);
                                                                outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                            }
                                                        }
                                                    }
                                                }
                                                break;
                                            case LRU:
                                            case SFIFO:
                                            case HLRU:
                                            case QLRU:
                                                if (isSimulation) {
                                                    Matrix varp = var.copy();

                                                    // probchoose(ac{class,k}(1+i,(1+i):end)/sum(ac{class,k}(1+i,(1+i):end)))
                                                    Matrix aci = new Matrix(1, ac[jobClass][k].getNumCols() - (i + 1));
                                                    Matrix.extract(ac[jobClass][k], i + 1, i + 2,
                                                            i + 1, ac[jobClass][k].getNumCols(),
                                                            aci, 0, 0);
                                                    aci = aci.scale(1 / aci.elementSum());
                                                    int probchoose = Maths.probchoose(aci);
                                                    int inew = i + probchoose;

                                                    Matrix.extract(var, 0, 1, cpos(m, i, 0), cpos(m, i, j), varp, 0, cpos(m, i, 1));
                                                    varp.set(cpos(m, i, 0), var.get(cpos(m, inew, (int) m.get(inew) - 1)));
                                                    Matrix.extract(var, 0, 1, cpos(m, inew, 0), cpos(m, inew, (int) m.get(inew) - 1), varp, 0, cpos(m, inew, 1));
                                                    varp.set(cpos(m, inew, 0), k + 1);

                                                    if (spaceSrvK.isEmpty()) {
                                                        spaceSrvK = spaceSrvE.copy();
                                                    } else {
                                                        spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                    }
                                                    if (spaceVarK.isEmpty()) {
                                                        spaceVarK = varp.copy();
                                                    } else {
                                                        spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                    }
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob = Matrix.concatRows(outprob, Matrix.singleton(ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k).doubleValue()), null);
                                                } else {
                                                    for (int inew = i; inew < h; inew++) {
                                                        Matrix varp = var.copy();
                                                        if (j > 0) {
                                                            Matrix.extract(var, 0, 1, cpos(m, i, 0), cpos(m, i, j), varp, 0, cpos(m, i, 1));
                                                        }
                                                        varp.set(cpos(m, i, 0), var.get(cpos(m, inew, (int) m.get(inew) - 1)));
                                                        if (m.get(inew) > 1) {
                                                            Matrix.extract(var, 0, 1, cpos(m, inew, 0), cpos(m, inew, (int) m.get(inew) - 1), varp, 0, cpos(m, inew, 1));
                                                        }
                                                        varp.set(cpos(m, inew, 0), k + 1);

                                                        if (spaceSrvK.isEmpty()) {
                                                            spaceSrvK = spaceSrvE.copy();
                                                        } else {
                                                            spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                        }
                                                        if (spaceVarK.isEmpty()) {
                                                            spaceVarK = varp.copy();
                                                        } else {
                                                            spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                        }
                                                        if (outrate.isEmpty()) {
                                                            outrate = new Matrix(1, 1);
                                                            outrate.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) * GlobalConstants.Immediate);
                                                        } else {
                                                            Matrix bottom_row = new Matrix(1, 1);
                                                            bottom_row.set(0, 0, ac[jobClass][k].get(1 + i, 1 + inew) * p.get(k) * GlobalConstants.Immediate);
                                                            outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                        }
                                                    }
                                                }
                                                break;
                                        }
                                    }
                                    // CACHE HIT in list h
                                    else {
                                        spaceSrvE.set((int) hitclass.get(jobClass), spaceSrvE.get((int) hitclass.get(jobClass)) + 1);
                                        int i = h;
                                        // j = posk - sum(m(1:i-1));
                                        int j = posk;
                                        for (int m_ind = 0; m_ind < i - 1; m_ind++) {
                                            j -= m.get(m_ind);
                                        }

                                        switch (replacement) {
                                            case RR:
                                            case FIFO:
                                            case SFIFO:
                                                if (spaceSrvK.isEmpty()) {
                                                    spaceSrvK = spaceSrvE.copy();
                                                } else {
                                                    spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                }
                                                if (spaceVarK.isEmpty()) {
                                                    spaceVarK = var.copy();
                                                } else {
                                                    spaceVarK = Matrix.concatRows(spaceVarK, var, null);
                                                }

                                                if (isSimulation) {
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob = Matrix.concatRows(outprob, Matrix.singleton(p.get(k).doubleValue()), null);
                                                } else {
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                }
                                                break;
                                            case LRU:
                                            case HLRU:
                                            case QLRU:
                                                Matrix varp = var.copy();
                                                Matrix.extract(var, 0, 1, cpos(m, h - 1, 0), cpos(m, h - 1, j), varp, 0, cpos(m, h - 1, 1));
                                                varp.set(cpos(m, h - 1, 0), var.get(cpos(m, h - 1, j)));


                                                if (spaceSrvK.isEmpty()) {
                                                    spaceSrvK = spaceSrvE.copy();
                                                } else {
                                                    spaceSrvK = Matrix.concatRows(spaceSrvK, spaceSrvE, null);
                                                }
                                                if (spaceVarK.isEmpty()) {
                                                    spaceVarK = varp.copy();
                                                } else {
                                                    spaceVarK = Matrix.concatRows(spaceVarK, varp, null);
                                                }

                                                if (isSimulation) {
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                    outprob = Matrix.concatRows(outprob, Matrix.singleton(p.get(k).doubleValue()), null);
                                                } else {
                                                    if (outrate.isEmpty()) {
                                                        outrate = new Matrix(1, 1);
                                                        outrate.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                                    } else {
                                                        Matrix bottom_row = new Matrix(1, 1);
                                                        bottom_row.set(0, 0, p.get(k) * GlobalConstants.Immediate);
                                                        outrate = Matrix.concatRows(outrate, bottom_row, null);
                                                    }
                                                }
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                        // if state is unchanged, still add with rate 0
                        outspace = Matrix.concatColumns(spaceSrvK, spaceVarK, null);
                    }
                }
                break;
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    public static int cpos(Matrix matrix, int i, int j) {
        int pos = 0;

        for (int idx = 0; idx < i; idx++) {
            pos += (int) matrix.get(idx);
        }
        pos += j;
        return pos;
    }

    /**
     * Performs a check to see if an item is currently in the retrieval system.
     *
     * <p>
     *     The retrieval system is encoded as a per-item occupancy bitmap appended after the cache contents:
     *     column {@code totalCacheCapacity + item} is non-zero iff that item is currently being retrieved.
     * </p>
     *
     * @param var: Matrix representing the space [Cache contents | Retrieval-system bitmap]
     * @param item: Zero-indexed item identifier
     * @param totalCacheCapacity: Total capacity of the cache
     * @return true if item is in the retrieval system, false otherwise
     */
    private static boolean isInRetrievalSystem(Matrix var, int item, int totalCacheCapacity) {
        int col = totalCacheCapacity + item;
        // No retrieval-system bitmap present (no retrieval system configured)
        if (col >= var.getNumCols()) {
            return false;
        }
        return var.get(col) != 0;
    }

    /**
     * Adds an item to the retrieval system by setting its occupancy bit.
     *
     * @param var: Matrix representing the space [Cache contents | Retrieval-system bitmap]
     * @param item: Zero-indexed item identifier
     * @param totalCacheCapacity: Total capacity of the cache
     */
    private static void addToRetrievalSystem(Matrix var, int item, int totalCacheCapacity) {
        int col = totalCacheCapacity + item;
        // No retrieval-system bitmap present (no retrieval system configured)
        if (col >= var.getNumCols()) {
            return;
        }
        var.set(col, 1);
    }

    /**
     * Removes an item from the retrieval system by clearing its occupancy bit.
     *
     * <p>
     *     If the retrieval system is not set, then this function does nothing.
     * </p>
     *
     * @param var: Matrix representing the space [Cache contents | Retrieval-system bitmap]
     * @param item: Zero-indexed item identifier
     * @param totalCacheCapacity: Total capacity of the cache
     */
    private static void removeFromRetrievalSystem(Matrix var, int item, int totalCacheCapacity) {
        int col = totalCacheCapacity + item;
        // No retrieval-system bitmap present (no retrieval system configured)
        if (col >= var.getNumCols()) {
            return;
        }
        var.set(col, 0);
    }
}
