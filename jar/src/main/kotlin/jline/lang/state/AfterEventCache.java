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
                    for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                        spaceSrv.set(row, jobClass, spaceSrv.get(row, jobClass) - 1);
                    }
                    switch (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass))) {
                        case RROBIN:
                            Matrix nvar_ind = new Matrix(1, R + jobClass);
                            Matrix.extract(sn.nvars, ind, ind + 1, 0, R + jobClass, nvar_ind, 0, 0);
                            int nvar_sum = (int) nvar_ind.elementSum();
                            int idx = -1;
                            Matrix outlinks = sn.nodeparam.get(sn.nodes.get(ind)).outlinks.get(sn.jobclasses.get(jobClass));
                            for (int row = 0; row < outlinks.getNumRows(); row++) {
                                if (spaceVar.get(nvar_sum) == outlinks.get(row)) {
                                    idx = row;
                                }
                            }
                            if (idx < outlinks.getNumRows()) {
                                spaceVar.set(nvar_sum, outlinks.get(idx + 1));
                            } else {
                                spaceVar.set(nvar_sum, outlinks.get(0));
                            }
                            break;
                    }
                    // buf is empty
                    outspace = Matrix.concatColumns(spaceSrv, spaceVar, null);
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
                Matrix[][] ac = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).accost;
                Matrix hitclass = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).hitclass;
                Matrix missclass = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).missclass;
                int h = m.getNumCols();
                ReplacementStrategy replacement = ((CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind))).replacestrat;

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
                                if (isSimulation) {
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
                                    for (int col = 0; col < var.getNumCols(); col++) {
                                        if (var.get(col) == k + 1) {
                                            posk = col;
                                            break;
                                        }
                                    }

                                    // CACHE MISS, can enter any list based on accessCost
                                    if (posk == -1) {
                                        spaceSrvE.set((int) missclass.get(jobClass), spaceSrvE.get((int) missclass.get(jobClass)) + 1);
                                        Matrix varp;
                                        switch (replacement) {
                                            case FIFO:
                                            case LRU:
                                            case SFIFO:
                                                if (isSimulation) {
                                                    int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                    int headPos = cpos(m, listidx, 0);
                                                    int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                    varp = var.copy();
                                                    // Shift items in list listidx to make room at the head
                                                    Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
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
                                                    // Iterate over all possible target lists (columns 1 to h)
                                                    for (l = 1; l <= h; l++) {
                                                        int listidx = l - 1; // l is accessCost column index, listidx is actual list (0-indexed)
                                                        int headPos = cpos(m, listidx, 0);
                                                        int tailPos = cpos(m, listidx, (int) m.get(listidx) - 1);
                                                        varp = var.copy();
                                                        // Shift items in list listidx to make room at the head
                                                        Matrix.extract(var, 0, 1, headPos, tailPos, varp, 0, headPos + 1);
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
                                                    for (int inew = i; inew < h; inew++) {
                                                        Matrix varp = var.copy();
                                                        varp.set(cpos(m, i, j), var.get(cpos(m, inew, (int) m.get(inew) - 1)));
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
}