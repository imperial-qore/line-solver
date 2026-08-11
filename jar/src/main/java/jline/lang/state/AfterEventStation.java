package jline.lang.state;

import static jline.GlobalConstants.NegInf;
import static jline.io.InputOutput.line_warning;

import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.RetrialPolicy;
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
    /**
     * Evaluate the class-dependence function beta_{i,r}(n) of a station and reduce
     * it to the scalar rate scaling for the given class.
     *
     * <p>The function may return either a 1x1 matrix (a chain-independent scaling
     * shared by every class) or a length-R row vector of the per-class scalings,
     * of which the classIdx element is taken; the latter expresses Sauer's
     * chain-dependent service rates mu_{r,i}(n) (Sauer 1983, eq. (40)).</p>
     *
     * <p>A station with no class dependence is ABSENT from the map (it is not
     * mapped to a constant 1, since under Sauer a constant 1 would assert that
     * every class completes at rate 1, which is not load independence). The state
     * machinery indexes every station unconditionally, so the neutral scaling is
     * supplied here.</p>
     */
    private static double cdScalar(Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                   Station st, Matrix nir, int classIdx) {
        SerializableFunction<Matrix, Matrix> f = (cdscaling != null) ? cdscaling.get(st) : null;
        if (f == null) {
            return 1.0;
        }
        // nir may span several inspace rows (e.g. fromMarginal-generated cells);
        // those rows share the same per-class population and differ only in
        // phase/buffer detail, so evaluate the handle on the first row -- a
        // genuine 1xR population vector -- rather than on the column-major
        // flattening of the whole matrix. This matches the row-0 convention the
        // surrounding code uses for lldscaling (ni.get(0)) and mirrors
        // State.cdclassfactor.m, whose per-row factors coincide here.
        Matrix n = (nir.getNumRows() > 1) ? Matrix.extractRows(nir, 0, 1, null) : nir;
        Matrix v = f.apply(n);
        return v.get(Math.min(classIdx, v.length() - 1));
    }

    static Ret.EventResult afterEventStation(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        return afterEventStation(sn, ind, inspace, event, jobClass, isSimulation, outspace, outrate, outprob, eventCache, M, R, S, phasessz, phaseshift, pie, ismkvmodclass, lldscaling, lldlimit, cdscaling, hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, spaceBuf, spaceSrv, spaceVar, key, false);
    }

    /**
     * noPromote: when true, a DEP at an FCFS-family station does not promote a
     * waiting job into the vacated server (immediate feedback, see
     * State.afterEvent). Only the DEP branch consults it.
     */
    static Ret.EventResult afterEventStation(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key, boolean noPromote) {
        // Server breakdown. The status is the trailing local-variable column
        // (0 = down, 1 = up), exclusive with the BAS marker and the polling
        // controller. A down server does not serve, so DEP and PHASE are suppressed
        // unless a degraded down-server rate was configured, in which case the
        // completion rate is rescaled by downRate/upRate after the handler returns.
        // Gating here keeps every scheduling branch below unaware of the server
        // status.
        boolean isBreakdownStation = sn.hasbreakdown != null && ind < sn.hasbreakdown.length()
                && sn.hasbreakdown.get(ind) == 1;
        double downRateScale = 1.0;
        if (isBreakdownStation && inspace != null && !inspace.isEmpty()
                && (event == EventType.DEP || event == EventType.PHASE)) {
            if (inspace.get(0, inspace.getNumCols() - 1) == 0) { // server down
                double downRate = 0.0;
                if (sn.downServiceRates != null && ist < sn.downServiceRates.getNumRows()
                        && jobClass < sn.downServiceRates.getNumCols()) {
                    downRate = sn.downServiceRates.get(ist, jobClass);
                }
                if (downRate <= 0) {
                    return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0));
                }
                double upRate = sn.rates.get(ist, jobClass);
                if (!Double.isFinite(upRate) || upRate <= 0) {
                    throw new RuntimeException("Station '" + sn.nodenames.get(ind) + "' declares a down-server "
                            + "service rate for class '" + sn.jobclasses.get(jobClass).getName() + "' but has no "
                            + "finite up-server service rate to rescale.");
                }
                downRateScale = downRate / upRate;
            }
        }
        switch (event) {
            case FAILURE:
                // An up server fails at the memoryless rate breakdownMu, whether or not
                // it is serving. Only the status column changes: jobs in service are not
                // lost and, service being memoryless here, they resume on repair. The
                // passive half of the synchronization is LOCAL, so no job moves.
                return handleStatusFlip(sn, ind, ist, inspace, isBreakdownStation, 1, 0,
                        (sn.breakdownMu != null) ? sn.breakdownMu.get(ist, 0) : 0.0);
            case REPAIR:
                // A down server is restored at the memoryless rate repairMu.
                return handleStatusFlip(sn, ind, ist, inspace, isBreakdownStation, 0, 1,
                        (sn.repairMu != null) ? sn.repairMu.get(ist, 0) : 0.0);
            default:
                break;
        }
        Ret.EventResult res = afterEventStationDispatch(sn, ind, inspace, event, jobClass, isSimulation, outspace, outrate, outprob, eventCache, M, R, S, phasessz, phaseshift, pie, ismkvmodclass, lldscaling, lldlimit, cdscaling, hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, spaceBuf, spaceSrv, spaceVar, key, noPromote);
        // Degraded service while the server is down: the scheduling handlers computed
        // the completion rate from the up-server service process, so rescale it to the
        // configured down-server rate. downRateScale is 1 whenever the server is up or
        // no degraded rate was configured, so this is a no-op in every other model.
        if (downRateScale != 1.0 && res.outrate != null && !res.outrate.isEmpty()) {
            return new Ret.EventResult(res.outspace, Matrix.scaleMult(res.outrate, downRateScale), res.outprob);
        }
        return res;
    }

    /**
     * Flips the server-status column of a breakdown station from fromStatus to
     * toStatus at the given memoryless rate, leaving the job counts untouched.
     *
     * @param sn                 network structure
     * @param ind                node index
     * @param ist                station index
     * @param inspace            the incoming state rows
     * @param isBreakdownStation whether the station carries a status column
     * @param fromStatus         the status the transition fires from
     * @param toStatus           the status the transition leads to
     * @param rate               the memoryless transition rate
     * @return the successor state, or an empty result when the transition is disabled
     */
    private static Ret.EventResult handleStatusFlip(NetworkStruct sn, int ind, int ist, Matrix inspace,
                                                    boolean isBreakdownStation, int fromStatus, int toStatus,
                                                    double rate) {
        if (!isBreakdownStation || inspace == null || inspace.isEmpty()) {
            return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0));
        }
        int last = inspace.getNumCols() - 1;
        if (inspace.get(0, last) != fromStatus) {
            return new Ret.EventResult(new Matrix(0, 0), new Matrix(0, 0), new Matrix(0, 0));
        }
        Matrix os = inspace.copy();
        for (int row = 0; row < os.getNumRows(); row++) {
            os.set(row, last, toStatus);
        }
        Matrix orate = new Matrix(os.getNumRows(), 1);
        Matrix oprob = new Matrix(os.getNumRows(), 1);
        for (int row = 0; row < os.getNumRows(); row++) {
            orate.set(row, 0, rate);
            oprob.set(row, 0, 1.0);
        }
        return new Ret.EventResult(os, orate, oprob);
    }

    private static Ret.EventResult afterEventStationDispatch(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key, boolean noPromote) {
        switch (event) {
            case ARV:
                return handleArv(sn, ind, inspace, event, jobClass, isSimulation, outspace, outrate, outprob, eventCache, M, R, S, phasessz, phaseshift, pie, ismkvmodclass, lldscaling, lldlimit, cdscaling, hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, spaceBuf, spaceSrv, spaceVar, key);
            case DEP:
                return handleDep(sn, ind, inspace, event, jobClass, isSimulation, outspace, outrate, outprob, eventCache, M, R, S, phasessz, phaseshift, pie, ismkvmodclass, lldscaling, lldlimit, cdscaling, hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, spaceBuf, spaceSrv, spaceVar, key, noPromote);
            case PHASE:
                return handlePhase(sn, ind, inspace, event, jobClass, isSimulation, outspace, outrate, outprob, eventCache, M, R, S, phasessz, phaseshift, pie, ismkvmodclass, lldscaling, lldlimit, cdscaling, hasOnlyExp, ist, K, Ks, mu, phi, proc, capacity, classcap, V, spaceBuf, spaceSrv, spaceVar, key);
            case RENEGE:
                return handleRenege(sn, ind, inspace, jobClass, ist, R, K, Ks, phasessz, phaseshift, pie, spaceBuf, spaceSrv, spaceVar);
            case RETRY:
                return handleRetry(sn, ind, inspace, jobClass, ist, R, S, K, Ks, phasessz, phaseshift, pie, spaceBuf, spaceSrv, spaceVar, isSimulation);
            case SWITCH:
                return handleSwitch(sn, ind, inspace, jobClass, ist, R, K, Ks, pie, spaceBuf, spaceSrv, spaceVar);
            default:
                return new Ret.EventResult(outspace, outrate, outprob);
        }
    }

    // Remove the buffer job at column `slot` and re-pad a zero on the left so the
    // buffer stays in the canonical right-aligned form used by the arrival handler.
    private static Matrix removeBufJobPrependZero(Matrix spaceBuf, int slot) {
        int nc = spaceBuf.getNumCols();
        Matrix out = new Matrix(1, nc);
        out.zero();
        int w = 1; // leave column 0 as the padding zero
        for (int c = 0; c < nc; c++) {
            if (c == slot) continue;
            out.set(0, w, spaceBuf.get(0, c));
            w++;
        }
        return out;
    }

    private static Matrix catBufSrvVar(Matrix buf, Matrix srv, Matrix var) {
        Matrix out = buf;
        if (srv != null && srv.getNumCols() > 0) out = Matrix.concatColumns(out, srv, null);
        if (var != null && var.getNumCols() > 0) out = Matrix.concatColumns(out, var, null);
        return out;
    }

    // Switchover of a polling server walking towards buffer `jobClass`. Unlike PHASE,
    // which carries only the internal transitions of a phase-type and leaves the
    // absorption to DEP, this event carries both: a completed switchover moves no job
    // and so has no departure to attach the absorption to (see Network.refreshSync).
    // It is therefore also emitted for a single-phase switchover, where it consists of
    // the absorption alone.
    private static Ret.EventResult handleSwitch(NetworkStruct sn, int ind, Matrix inspace, int jobClass, int ist, int R,
                                                Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> pie,
                                                Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar) {
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(0, 0);
        Polling.Info pinfoS = Polling.info(sn, ind);
        if (pinfoS == null || !pinfoS.hasSw[jobClass]) {
            return new Ret.EventResult(outspace, outrate, outprob);
        }
        for (int row = 0; row < inspace.getNumRows(); row++) {
            Matrix varRowS = Matrix.extractRows(spaceVar, row, row + 1, null);
            int[] ctlS = Polling.get(pinfoS, varRowS, -1);
            if (ctlS[0] != jobClass || ctlS[1] == 0) {
                continue; // the server is not inside the switchover into `jobClass`
            }
            Matrix bufRowS = Matrix.extractRows(spaceBuf, row, row + 1, null);
            Matrix srvRowS = Matrix.extractRows(spaceSrv, row, row + 1, null);
            Matrix D0S = pinfoS.swD0[jobClass];
            Matrix D1S = pinfoS.swD1[jobClass];
            int swk = ctlS[1];
            // internal transitions of the switchover phase-type
            for (int kdest = 1; kdest <= pinfoS.Ksw[jobClass]; kdest++) {
                if (kdest == swk) {
                    continue;
                }
                double d0 = D0S.get(swk - 1, kdest - 1);
                if (d0 <= 0) {
                    continue;
                }
                Matrix varS = Polling.set(pinfoS, varRowS, jobClass, kdest, 0);
                outspace = Matrix.concatRows(outspace, catBufSrvVar(bufRowS, srvRowS, varS), null);
                Matrix orS = new Matrix(1, 1); orS.set(0, 0, d0);
                outrate = Matrix.concatRows(outrate, orS, null);
                Matrix opS = new Matrix(1, 1); opS.set(0, 0, 1.0);
                outprob = Matrix.concatRows(outprob, opS, null);
            }
            // absorption: the server arrives at buffer `jobClass` and either opens a
            // visit there or walks on
            double rateS = 0;
            for (int c = 0; c < D1S.getNumCols(); c++) {
                rateS += D1S.get(swk - 1, c);
            }
            if (rateS <= 0) {
                continue;
            }
            int[] nbufS = new int[R];
            for (int r = 0; r < R; r++) {
                nbufS[r] = (int) Math.round(bufRowS.get(0, r));
            }
            int[] resS = Polling.next(pinfoS, jobClass, nbufS, R, true);
            Polling.Landing landS = Polling.land(pinfoS, resS[0], resS[1], resS[2],
                    bufRowS, srvRowS, varRowS, K, Ks, pie.get(sn.stations.get(ist)), sn.jobclasses, R);
            for (int jS = 0; jS < landS.rows.size(); jS++) {
                Matrix cand = landS.rows.get(jS);
                // A switchover that completes over an empty buffer starts the next leg
                // at once, and when that leg re-enters the same phase of the same
                // switchover (the single-buffer case, or any memoryless switchover) the
                // landing state is the departure state. Such a self-loop is not a
                // transition: emitting it would inflate the exit rate of the row.
                boolean same = cand.getNumCols() == inspace.getNumCols();
                if (same) {
                    for (int c = 0; c < cand.getNumCols(); c++) {
                        if (Math.abs(cand.get(0, c) - inspace.get(row, c)) > 1e-12) { same = false; break; }
                    }
                }
                if (same) {
                    continue;
                }
                outspace = Matrix.concatRows(outspace, cand, null);
                Matrix orS = new Matrix(1, 1);
                orS.set(0, 0, rateS * landS.probs.get(jS).doubleValue());
                outrate = Matrix.concatRows(outrate, orS, null);
                Matrix opS = new Matrix(1, 1); opS.set(0, 0, 1.0);
                outprob = Matrix.concatRows(outprob, opS, null);
            }
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    // Exponential-patience reneging: a waiting (queued, not-in-service) class-r job
    // abandons the queue at aggregate rate (waiting count) * impatienceMu; one
    // waiting job is removed and leaves the system. Mirrors MATLAB afterEventStation.
    private static Ret.EventResult handleRenege(NetworkStruct sn, int ind, Matrix inspace, int jobClass, int ist, int R,
                                                Matrix K, Matrix Ks, Matrix phasessz, Matrix phaseshift,
                                                Map<Station, Map<JobClass, Matrix>> pie, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar) {
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(0, 0);
        State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, inspace, phasessz, phaseshift, spaceBuf, spaceSrv, spaceVar);
        double waiting = stats.nir.get(0, jobClass) - stats.sir.get(0, jobClass);
        if (waiting > 0 && spaceBuf.getNumCols() > 0) {
            int slot = -1;
            for (int c = 0; c < spaceBuf.getNumCols(); c++) {
                if ((int) Math.round(spaceBuf.get(0, c)) == jobClass + 1) { slot = c; break; }
            }
            if (slot >= 0) {
                Matrix bufK = removeBufJobPrependZero(spaceBuf, slot);
                outspace = catBufSrvVar(bufK, spaceSrv, spaceVar);
                double muRate = sn.impatienceMu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0);
                outrate = new Matrix(1, 1); outrate.set(0, 0, waiting * muRate);
                outprob = new Matrix(1, 1); outprob.set(0, 0, 1.0);
            }
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    // Exponential retrial: an orbiting (buffered) class-r job retries entry; it
    // succeeds only when a server is free, in which case one orbiting job enters
    // service (entry phase drawn from pie) at aggregate rate (orbit) * retrialMu.
    private static Ret.EventResult handleRetry(NetworkStruct sn, int ind, Matrix inspace, int jobClass, int ist, int R, Matrix S,
                                               Matrix K, Matrix Ks, Matrix phasessz, Matrix phaseshift,
                                               Map<Station, Map<JobClass, Matrix>> pie, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, boolean isSimulation) {
        Matrix outspace = new Matrix(0, 0);
        Matrix outrate = new Matrix(0, 0);
        Matrix outprob = new Matrix(0, 0);
        State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, inspace, phasessz, phaseshift, spaceBuf, spaceSrv, spaceVar);
        double orbit = stats.nir.get(0, jobClass) - stats.sir.get(0, jobClass);
        double inSrv = spaceSrv.elementSum();
        int Sist = (int) S.get(ist);
        if (orbit > 0 && inSrv < Sist && spaceBuf.getNumCols() > 0) {
            int slot = -1;
            for (int c = 0; c < spaceBuf.getNumCols(); c++) {
                if ((int) Math.round(spaceBuf.get(0, c)) == jobClass + 1) { slot = c; break; }
            }
            if (slot >= 0) {
                Matrix bufK = removeBufJobPrependZero(spaceBuf, slot);
                Matrix pentry = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass));
                double retrialMu = sn.retrialMu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0);
                int Kr = (int) K.get(jobClass);
                int ksBase = (int) Ks.get(jobClass);
                for (int kentry = 0; kentry < Kr; kentry++) {
                    double pe = pentry.get(kentry);
                    if (pe <= 0) continue;
                    Matrix srvK = spaceSrv.copy();
                    srvK.set(0, ksBase + kentry, srvK.get(0, ksBase + kentry) + 1);
                    Matrix os = catBufSrvVar(bufK, srvK, spaceVar);
                    // LINEAR policy: every orbiting job carries its own timer, so the
                    // aggregate retrial rate scales with the orbit size. CONSTANT
                    // policy: one controller retries on behalf of the whole orbit, so
                    // the rate does not depend on the orbit size.
                    double orbitWeight = orbit;
                    if (sn.retrialPolicy != null) {
                        Map<JobClass, Integer> pmap = sn.retrialPolicy.get(sn.stations.get(ist));
                        Integer pol = (pmap != null) ? pmap.get(sn.jobclasses.get(jobClass)) : null;
                        if (pol != null && pol.intValue() == RetrialPolicy.CONSTANT) {
                            orbitWeight = 1.0;
                        }
                    }
                    Matrix orow = new Matrix(1, 1); orow.set(0, 0, orbitWeight * retrialMu * pe);
                    Matrix prow = new Matrix(1, 1); prow.set(0, 0, 1.0);
                    if (outspace.isEmpty()) { outspace = os; outrate = orow; outprob = prow; }
                    else { outspace = Matrix.concatRows(outspace, os, null); outrate = Matrix.concatRows(outrate, orow, null); outprob = Matrix.concatRows(outprob, prow, null); }
                }
                if (isSimulation && outspace.getNumRows() > 1) {
                    Matrix cum = outrate.cumsumViaCol();
                    double tot = outrate.elementSum();
                    Matrix cumr = Matrix.scaleMult(cum, 1.0 / tot);
                    int fc = -1; double rnd = Maths.rand();
                    for (int row = 0; row < cumr.getNumRows(); row++) if (rnd > cumr.get(row)) fc = row;
                    fc++;
                    outspace = Matrix.extractRows(outspace, fc, fc + 1, null);
                    Matrix nr = new Matrix(1, 1); nr.set(0, 0, tot); outrate = nr;
                    Matrix np = new Matrix(1, 1); np.set(0, 0, 1.0); outprob = np;
                }
            }
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    private static Ret.EventResult handleArv(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        Matrix sir = null;
        List<Matrix> kir = null;
        Matrix ni = null;
        Matrix nir = null;
                // A Place holds a marking, not a service facility: it has no servers and
                // no service phases, so an arriving token only increments the class
                // marking. The scheduling branches below would instead write the entering
                // job into a phase slot that a Place state does not carry, widening the
                // row so that it can no longer be matched against the marking state
                // space; the successor lookup then fails and the arrival is dropped.
                // Closed SPNs never notice, because there tokens reach a Place through
                // FIRE (State.afterGlobalEvent), not through ARV.
                if (sn.nodetype.get(ind) == NodeType.Place) {
                    Matrix placeOut = inspace.copy();
                    for (int i = 0; i < placeOut.getNumRows(); i++) {
                        placeOut.set(i, jobClass, placeOut.get(i, jobClass) + 1);
                    }
                    // passive action: the rate is set by the active node
                    Matrix placeRate = new Matrix(placeOut.getNumRows(), 1);
                    Matrix placeProb = new Matrix(placeOut.getNumRows(), 1);
                    for (int i = 0; i < placeOut.getNumRows(); i++) {
                        placeRate.set(i, 0, -1.0);
                        placeProb.set(i, 0, 1.0);
                    }
                    return new Ret.EventResult(placeOut, placeRate, placeProb);
                }
                // Signal / catastrophe arrival (G-network negative customer): instead of
                // joining, the signal removes job(s) from this station and is annihilated.
                // The batch size, the eligible victim classes and the victim-selection
                // policy are all read from sn; see SignalRemoval, which ports MATLAB
                // State.afterEventStationSignal.
                if (sn.issignal != null && sn.issignal.get(jobClass) > 0) {
                    // A REPLY signal is not a negative customer: it completes a
                    // synchronous call, releasing the server this station holds for the
                    // caller, and then joins as an ordinary job. Only stations that
                    // actually hold a block for it take this path; elsewhere a REPLY
                    // class is a plain job class and falls through to the normal arrival
                    // handling below.
                    boolean isReplySignal = sn.signaltype != null && jobClass < sn.signaltype.size()
                            && sn.signaltype.get(jobClass) == jline.lang.constant.SignalType.REPLY;
                    if (isReplySignal) {
                        if (ReplyBlock.info(sn, ind).width > 0) {
                            Ret.EventResult replyRes = AfterEventStationReply.apply(sn, ind, ist, jobClass,
                                    K, Ks, S, pie, spaceBuf, spaceSrv, spaceVar);
                            if (isSimulation && replyRes.outprob.getNumRows() > 1) {
                                // Which entry phase the reply starts in is a random choice.
                                // The generator needs every destination; a simulation must
                                // pick exactly one, so sample it here.
                                double totProb = replyRes.outprob.elementSum();
                                double u = Math.random() * totProb;
                                double cum = 0;
                                int firing = replyRes.outprob.getNumRows() - 1;
                                for (int i = 0; i < replyRes.outprob.getNumRows(); i++) {
                                    cum += replyRes.outprob.get(i, 0);
                                    if (u <= cum) {
                                        firing = i;
                                        break;
                                    }
                                }
                                Matrix oneProb = new Matrix(1, 1);
                                oneProb.set(0, 0, 1.0);
                                replyRes = new Ret.EventResult(
                                        Matrix.extractRows(replyRes.outspace, firing, firing + 1, null),
                                        Matrix.extractRows(replyRes.outrate, firing, firing + 1, null),
                                        oneProb);
                            }
                            return replyRes;
                        }
                    } else {
                        return SignalRemoval.handleSignalArrival(sn, ind, ist, inspace, jobClass, isSimulation,
                                phasessz, phaseshift, K, Ks, S, spaceBuf, spaceSrv, spaceVar);
                    }
                }
                // Ordinary SPN Place arrival: see _kb/11-conventions-and-gotchas.md
                // ("An ordinary Place must be special-cased before the generic
                // scheduling switch") for the rationale. Mirrors MATLAB/Python
                // State.afterEventStation.
                if (sn.nodetype != null && ind < sn.nodetype.size()
                        && sn.nodetype.get(ind) == NodeType.Place) {
                    Matrix outsp = inspace.copy();
                    Matrix orate = new Matrix(outsp.getNumRows(), 1);
                    Matrix oprob = new Matrix(outsp.getNumRows(), 1);
                    double cap = (classcap != null) ? classcap.get(ist, jobClass) : Double.POSITIVE_INFINITY;
                    for (int i = 0; i < outsp.getNumRows(); i++) {
                        orate.set(i, 0, -1); // passive: rate set by the active source
                        if (outsp.get(i, jobClass) < cap) {
                            outsp.set(i, jobClass, outsp.get(i, jobClass) + 1);
                            oprob.set(i, 0, 1.0);
                        } else {
                            oprob.set(i, 0, 0.0); // place full: arrival blocked and lost
                        }
                    }
                    return new Ret.EventResult(outsp, orate, oprob);
                }
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
                            // classcap is a hard bound: classcap = 0 bars the class from this
                            // station outright, it does not mean "unconstrained". refreshCapacity
                            // sets 0 for a class the station does not serve, and Solver_ssa/
                            // Solver_ctmc set 0 for signal classes at non-EXT stations, precisely
                            // so those slots stay unoccupied. Unbounded is encoded as Inf.
                            // Admit only when every row is strictly below the cap, matching
                            // MATLAB afterEventStation.m "if space_srv_k(:,Ks(class)+kentry) <
                            // classcap(ist,class)" (an if over a vector requires all entries true).
                            int col = (int) (Ks.get(jobClass) + kentry);
                            for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                if (spaceSrvK.get(row, col) >= istCap) {
                                    exceedsClassCap = true;
                                }
                            }
                            if (!exceedsClassCap) {
                                // increment spacesrvk by one
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
                        case POLLING: {
                            // The controller, not the arrival, decides who is served: an
                            // arriving job joins its class buffer and waits for the server
                            // to walk to it, even when the facility is idle, because the
                            // server is then in a switchover. The exception is a parked
                            // server, which Polling.next only ever produces with an empty
                            // station and immediate switchovers: it therefore reaches the
                            // arriving job in zero time and starts a visit on it at once.
                            Polling.Info pinfoA = Polling.info(sn, ind);
                            outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                            for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                int srvclassA = -1;
                                for (int r = 0; r < R; r++) {
                                    double tot = 0;
                                    for (int k = 0; k < (int) K.get(r); k++) {
                                        tot += spaceSrvK.get(row, (int) Ks.get(r) + k);
                                    }
                                    if (tot > 0) { srvclassA = r; break; }
                                }
                                Matrix varRowA = Matrix.extractRows(spaceVarK, row, row + 1, null);
                                int[] ctlA = Polling.get(pinfoA, varRowA, srvclassA);
                                if (srvclassA < 0 && ctlA[1] == 0) {
                                    // Parked: the arriving job is the only work present, so
                                    // the walk necessarily resolves to a visit on its own
                                    // class; it enters service here and never occupies the buffer.
                                    int[] nbufA = new int[R];
                                    for (int r = 0; r < R; r++) {
                                        nbufA[r] = (int) Math.round(spaceBufK.get(row, r));
                                    }
                                    nbufA[jobClass] += 1;
                                    int[] resA = Polling.next(pinfoA, ctlA[0], nbufA, R, true);
                                    int colSrvA = (int) (Ks.get(jobClass) + kentry);
                                    spaceSrvK.set(row, colSrvA, spaceSrvK.get(row, colSrvA) + 1);
                                    Matrix setA = Polling.set(pinfoA, varRowA, resA[0], 0, resA[2]);
                                    for (int c = 0; c < setA.getNumCols(); c++) {
                                        spaceVarK.set(row, c, setA.get(0, c));
                                    }
                                } else {
                                    spaceBufK.set(row, jobClass, spaceBufK.get(row, jobClass) + 1);
                                }
                            }
                            outprobK.fill(pentry.get(kentry));
                            break;
                        }
                        case SIRO:
                        case SEPT:
                        case LEPT:
                            outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                            for (int row = 0; row < spaceSrvK.getNumRows(); row++) {
                                // Idle-server test on server occupancy, not total count ni:
                                // these agree for work-conserving states but an
                                // immediate-feedback self-loop transiently yields an idle
                                // server with a non-empty buffer, where the fed-back job must
                                // re-enter the vacated server (mirrors the FCFS case below).
                                double srvCountArv = 0;
                                for (int scol = 0; scol < spaceSrvK.getNumCols(); scol++) {
                                    srvCountArv += spaceSrvK.get(row, scol);
                                }
                                if (srvCountArv < S.get(ist)) {
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
                        case FCFSPRIO:
                        case LCFSPRIO:
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
                                // Servers held by synchronous calls awaiting a REPLY are NOT
                                // available to an arriving job: subtract them from the server
                                // count. Zero for every model without reply signals.
                                Matrix nbA = ReplyBlock.blockedTotal(sn, ind, spaceVarK);
                                // construct all_busy_srv, a matrix with 0s where sum of that row in space_srv_k is >= S.get(ist) and 1s where its <
                                Matrix all_busy_srv = new Matrix(spaceSrvK.getNumRows(), 1);
                                for (int i = 0; i < spaceSrvK.getNumRows(); i++) {
                                    Matrix row = Matrix.extractRows(spaceSrvK, i, i + 1, null);
                                    int rowSum = (int) row.elementSum();
                                    double SeffA = S.get(ist) - (i < nbA.getNumRows() ? nbA.get(i, 0) : 0);
                                    if (rowSum >= SeffA) {
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
                                    double SeffA = S.get(ist) - (i < nbA.getNumRows() ? nbA.get(i, 0) : 0);
                                    if (rowSum < SeffA) {
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
                                    // classcap is a hard bound: 0 bars the class outright rather
                                    // than meaning "unconstrained"; unbounded is encoded as Inf.
                                    // Mirrors MATLAB afterEventStation.m
                                    // "if any(nir(:,class) < classcap(ist,class))", which is false
                                    // for classcap = 0 because occupancies are non-negative.
                                    for (int i = 0; i < nir.getNumRows(); i++) {
                                        if (nir.get(i, jobClass) < classCapLimit) {
                                            spaceForClass = true;
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
                                // ignore states where buffer has no empty slots.
                                // A structurally free buffer column is NOT enough: the job
                                // may only be placed when the CAPACITY also permits it.
                                // Without this the job was written into the slot and the
                                // whole row was then vetoed by the capacity filter below
                                // (en_o), which returned an EMPTY outspace -- i.e. the
                                // arrival event never fired, so the upstream departure was
                                // never counted and the reported arrival rate collapsed
                                // from the OFFERED rate to the CARRIED one (BUG-85: M/M/1/1
                                // reported ArvR 0.444444 instead of 0.800000). Leaving the
                                // state unchanged instead makes the arrival a self-loop for
                                // a LOST (open) class -- see the refusal branch below and
                                // State.arrivalIsLost. classcap == 0 means no per-class
                                // constraint (MATLAB convention), matching the en_o filter.
                                //
                                // Gate the capacity check ONLY for a PHYSICAL finite
                                // capacity. When the bound is a state-space cutoff (open
                                // class, no physical cap; SSA folds the cutoff into the
                                // capacity/classcap arguments), place the job structurally
                                // and let the en_o filter below delete the beyond-cutoff row
                                // (pre-change truncation) -- firing the capacity gate at a
                                // cutoff would turn a truncation into a self-loop. See
                                // State.isPhysicalCapacity.
                                boolean physCap = State.isPhysicalCapacity(sn, ist, jobClass);
                                Matrix wbuf_empty = new Matrix(empty_slots.getNumRows(), 1);
                                boolean space_available = false;
                                double capLimit = capacity.get(ist);
                                double classCapLimit = classcap.get(ist, jobClass);
                                for (int i = 0; i < empty_slots.getNumRows(); i++) {
                                    int niRow = (ni != null && ni.getNumRows() > 1) ? i : 0;
                                    boolean hasRoom = (ni == null) || !physCap
                                            || (ni.get(niRow, 0) < capLimit
                                                && (classCapLimit == 0 || nir.get(niRow, jobClass) < classCapLimit));
                                    if (empty_slots.get(i, 0) > 0 && hasRoom) {
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

                                } else {
                                    // The arrival cannot be placed (all servers busy, no
                                    // room). Whether that is a LOSS or a BLOCK is decided by
                                    // the class type, not the drop rule: a CLOSED job has
                                    // nowhere to go so it cannot be lost, and an explicit
                                    // blocking rule (BAS/BBS/RSRD) asks for blocking too.
                                    // Removing the rows returns an EMPTY outspace, which
                                    // disables the upstream departure until space frees (and
                                    // is what the CTMC become-blocked edge tests for, so true
                                    // BAS fires). For a LOST arrival the rows are kept
                                    // unchanged instead, so the event still fires and ArvR
                                    // counts the offered job. See State.arrivalIsLost and the
                                    // hasRoom comment above.
                                    boolean anyBusy = all_busy_srv.elementMax() > 0;
                                    if (anyBusy && !State.arrivalIsLost(sn, ist, jobClass)) {
                                        Matrix srvTmp = new Matrix(0, spaceSrvK.getNumCols());
                                        Matrix bufTmp = new Matrix(0, spaceBufK.getNumCols());
                                        Matrix varTmp = new Matrix(0, spaceVarK.getNumCols());
                                        for (int i = 0; i < idle_srv.getNumRows(); i++) {
                                            if (idle_srv.get(i, 0) == 1) {
                                                srvTmp = srvTmp.isEmpty() ? Matrix.extractRows(spaceSrvK, i, i + 1, null)
                                                        : Matrix.concatRows(srvTmp, Matrix.extractRows(spaceSrvK, i, i + 1, null), null);
                                                bufTmp = bufTmp.isEmpty() ? Matrix.extractRows(spaceBufK, i, i + 1, null)
                                                        : Matrix.concatRows(bufTmp, Matrix.extractRows(spaceBufK, i, i + 1, null), null);
                                                varTmp = varTmp.isEmpty() ? Matrix.extractRows(spaceVarK, i, i + 1, null)
                                                        : Matrix.concatRows(varTmp, Matrix.extractRows(spaceVarK, i, i + 1, null), null);
                                            }
                                        }
                                        spaceSrvK = srvTmp;
                                        spaceBufK = bufTmp;
                                        spaceVarK = varTmp;
                                    }
                                }
                                outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                                outprobK.fill(pentry.get(kentry));
                            } else {
                                outprobK = new Matrix(spaceSrvK.getNumRows(), 1);
                                outprobK.zero(); // zero probability event
                            }
                            break;
                        case FCFSPRPRIO: // FCFS preempt-resume with priority groups
                        case LCFSPRPRIO: // LCFS preempt-resume with priority groups
                        case FCFSPIPRIO: // FCFS preempt-independent with priority groups
                        case LCFSPIPRIO: // LCFS preempt-independent with priority groups
                        case LCFSPR: // LCFS with Preemption
                            // Reset per-entry-phase probability
                            outprobK = new Matrix(0, 0);
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
                                // For priority variants, only higher-priority jobs can preempt
                                SchedStrategy schedIstArr = sn.sched.get(sn.stations.get(ist));
                                boolean isPrioSched = (schedIstArr == SchedStrategy.FCFSPRPRIO || schedIstArr == SchedStrategy.FCFSPIPRIO || schedIstArr == SchedStrategy.LCFSPRPRIO || schedIstArr == SchedStrategy.LCFSPIPRIO);
                                // Priority-awareness is a property of the DECLARED policy, never
                                // of the data. LCFSPRPRIO/FCFSPRPRIO are the priority-aware
                                // variants; inferring the discipline from "do the class
                                // priorities differ" turned plain LCFSPR/FCFSPR into something
                                // that is neither the base policy nor the PRIO variant.
                                boolean isPrioAware = isPrioSched;
                                if (isPrioAware) {
                                    if (sn.classprio.get(jobClass) >= sn.classprio.get(classpreempt)) continue; // arriving job has same or lower priority
                                }
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
                                                boolean isPreemptIndep = (schedIstArr == SchedStrategy.LCFSPI || schedIstArr == SchedStrategy.LCFSPIPRIO ||
                                                                          schedIstArr == SchedStrategy.FCFSPI || schedIstArr == SchedStrategy.FCFSPIPRIO);
                                                for (int row = 0; row < spaceBufKPreemptFiltered.getNumRows(); row++) {
                                                    int emptySlot = (int) emptySlotsFiltered.get(row, 0);
                                                    if (emptySlot > 0) { // MATLAB uses 1-based indexing
                                                        // Convert back to 0-based for Java
                                                        int zeroBasedSlot = emptySlot - 1;
                                                        spaceBufKPreemptFiltered.set(row, zeroBasedSlot, classpreempt + 1); // Store class (1-based)
                                                        if (isPreemptIndep) {
                                                            spaceBufKPreemptFiltered.set(row, zeroBasedSlot + 1, 1); // Preempt-independent: store phase 1 placeholder
                                                        } else {
                                                            spaceBufKPreemptFiltered.set(row, zeroBasedSlot + 1, phasepreempt + 1); // Preempt-resume: store actual phase (1-based)
                                                        }
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
                    Matrix oi;
                    Matrix oir;
                    if (outspaceK.isEmpty() || spaceSrvK.getNumRows() == 0 || spaceSrvK.getNumCols() == 0) {
                        // No arrival state was produced, e.g. a priority-preempt arrival
                        // that cannot preempt any same/lower-priority in-service job while
                        // all servers are busy: the arrival is disabled. Skip the marginal
                        // count (which would extract a column from an empty server matrix
                        // and throw) and leave the empty result as a no-op, matching MATLAB
                        // State.afterEventStation.
                        oi = new Matrix(0, 0);
                        oir = new Matrix(0, 0);
                    } else {
                        State.StateMarginalStatistics oi_oir = ToMarginal.toMarginalAggr(sn, ind, outspaceK, K, Ks, spaceBufK, spaceSrvK, spaceVarK);
                        oi = oi_oir.ni;
                        oir = oi_oir.nir;
                    }

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
                        // classcap is a hard bound: 0 admits only the zero-occupancy row rather
                        // than meaning "unconstrained"; unbounded is encoded as Inf. Mirrors
                        // MATLAB afterEventStation.m "en_o = classcap(ist,class) >= oir(:,class)".
                        boolean passesClassCapCheck = (classCapLimit >= oir.get(row, jobClass));
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

                // Balking (QUEUE_LENGTH): with probability balkProb the arriving
                // class-r job refuses to join based on the pre-arrival total
                // station population ni; the job is lost (state unchanged) and
                // the admitted branches are scaled by (1-balkProb). Only the
                // QUEUE_LENGTH strategy is a pure state function; EXPECTED_WAIT /
                // COMBINED are rejected in the analyzers. Mirrors MATLAB
                // afterEventStation.m.
                if (sn.balkingStrategy != null && !outspace.isEmpty()) {
                    jline.lang.nodes.Station stn = sn.stations.get(ist);
                    jline.lang.JobClass jc = sn.jobclasses.get(jobClass);
                    Map<jline.lang.JobClass, jline.lang.constant.BalkingStrategy> smap = sn.balkingStrategy.get(stn);
                    if (smap != null && smap.get(jc) == jline.lang.constant.BalkingStrategy.QUEUE_LENGTH) {
                        double balkProb = 0.0;
                        java.util.List<jline.lang.constant.BalkingThreshold> thr =
                                (sn.balkingThresholds != null && sn.balkingThresholds.get(stn) != null)
                                        ? sn.balkingThresholds.get(stn).get(jc) : null;
                        int qlen = (int) Math.round(ni.get(0));
                        if (thr != null) {
                            for (jline.lang.constant.BalkingThreshold t : thr) {
                                if (t.matches(qlen)) { balkProb = t.getProbability(); break; }
                            }
                        }
                        if (balkProb > 0.0) {
                            outprob = Matrix.scaleMult(outprob, 1.0 - balkProb);
                            Matrix inrow = inspace;
                            if (outspace.getNumCols() > inrow.getNumCols()) {
                                Matrix zeros = new Matrix(inrow.getNumRows(), outspace.getNumCols() - inrow.getNumCols());
                                zeros.zero();
                                inrow = Matrix.concatColumns(zeros, inrow, null);
                            }
                            outspace = Matrix.concatRows(outspace, inrow, null);
                            Matrix br = new Matrix(1, 1); br.set(0, 0, -1);
                            outrate = Matrix.concatRows(outrate, br, null);
                            Matrix bp = new Matrix(1, 1); bp.set(0, 0, balkProb);
                            outprob = Matrix.concatRows(outprob, bp, null);
                        }
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
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    private static Ret.EventResult handleDep(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key, boolean noPromote) {
        Matrix sir = null;
        List<Matrix> kir = null;
        Matrix ni = null;
        Matrix nir = null;
                // Marked (MMAP) source class: the shared modulating chain lives
                // in the carrier's phase block (mark index 1); this class's
                // departures fire from there using its per-mark matrix D1k
                // (M3A cell index 1+mark). Mirrors MATLAB afterEventStation.
                int markofclass = -1;
                int phclass = jobClass;
                if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT
                        && sn.markidx != null && ist < sn.markidx.getNumRows()
                        && sn.markidx.get(ist, jobClass) > 0) {
                    markofclass = (int) sn.markidx.get(ist, jobClass);
                    for (int r = 0; r < R; r++) {
                        if (sn.markidx.get(ist, r) == 1) {
                            phclass = r;
                            break;
                        }
                    }
                }
                boolean busy = false;
                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                    for (int col = (int) Ks.get(phclass); col < (Ks.get(phclass) + K.get(phclass)); col++) {
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

                    if (sn.routing.get(sn.nodes.get(ind)).get(sn.jobclasses.get(jobClass)) == RoutingStrategy.WRROBIN) {
                        // WRROBIN advances a POSITION pointer cyclically through the
                        // weighted-outlink cycle (destinations repeated by weight), not
                        // a node index. Without this the pointer never advances and only
                        // the initial destination is ever selected.
                        int nvarCols = R + jobClass + 1;
                        Matrix nvar_ind = new Matrix(1, nvarCols);
                        Matrix.extract(sn.nvars, ind, ind + 1, 0, nvarCols, nvar_ind, 0, 0);
                        int nvar_sum = (int) nvar_ind.elementSum();
                        int spaceVarCol = nvar_sum - 1;
                        Matrix wol = sn.nodeparam.get(sn.nodes.get(ind)).weightedOutlinks.get(sn.jobclasses.get(jobClass));
                        int cyc = (wol == null) ? 0 : (int) wol.length();
                        if (spaceVarCol >= 0 && spaceVarCol < spaceVar.getNumCols() && cyc > 0) {
                            for (int stateRow = 0; stateRow < spaceVar.getNumRows(); stateRow++) {
                                int pos = (int) spaceVar.get(stateRow, spaceVarCol);
                                int nxt = (pos >= 1 && pos < cyc) ? (pos + 1) : 1;
                                spaceVar.set(stateRow, spaceVarCol, nxt);
                            }
                        }
                    }

                    if (sir.get(phclass) > 0) {
                        outprob = new Matrix(0, 0);
                        for (int k = 0; k < K.get(phclass); k++) {
                            spaceSrv = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V)); // server state
                            int spaceBufCols = (int) (inspace.getNumCols() - K.elementSum() - V);
                            spaceBuf = Matrix.extract(inspace, 0, inspace.getNumRows(), 0, spaceBufCols); // buffer state
                            Matrix rate = new Matrix(spaceSrv.getNumRows(), 1);
                            rate.zero();
                            Matrix en = new Matrix(spaceSrv.getNumRows(), 1);
                            boolean en_set = false;
                            for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                if (spaceSrv.get(row, (int) (Ks.get(phclass) + k)) > 0) {
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
                                            Matrix D1_srv;
                                            if (markofclass > 0) {
                                                // per-mark arrival matrix over the shared chain
                                                D1_srv = (Matrix) proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(1 + markofclass);
                                            } else {
                                                D1_srv = (Matrix) proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(1);
                                            }
                                            for (int kentry = 0; kentry < K.get(phclass); kentry++) {
                                                double arvRate = D1_srv.get(k, kentry);
                                                if (arvRate <= 0) {
                                                    continue;
                                                }
                                                Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V), spaceSrv, 0, 0); // server state

                                                // record a departure
                                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        spaceSrv.set(row, (int) (Ks.get(phclass) + k), spaceSrv.get(row, (int) (Ks.get(phclass) + k)) - 1);
                                                    }
                                                }
                                                // record a new job arriving
                                                for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        spaceSrv.set(row, (int) (Ks.get(phclass) + kentry), spaceSrv.get(row, (int) (Ks.get(phclass) + kentry)) + 1);
                                                    }
                                                }
                                                // extract all rows of spaceBuf where en==1
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
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    Matrix outrate_bottom = new Matrix(inspaceEn.getNumRows(), 1);
                                                    outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, lldlimit - 1) * arvRate);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                } else {
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    Matrix outrate_bottom = new Matrix(inspaceEn.getNumRows(), 1);
                                                    outrate_bottom.fill(cdscalingIst * lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit)) * arvRate);
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
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEn, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnPs, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnPs, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_ps = new Matrix(rateEnPs.getNumRows(), 1);
                                        outprob_bottom_ps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_ps, null);
                                        break;
                                    case PSPRIO: {
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
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnPsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                    }

                                    case FCFS:
                                        // job departing
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }
                                        // set en_wbuf to states with jobs in buffer.
                                        // Servers held for a pending REPLY are unavailable, so a
                                        // job waits in the buffer already when ni exceeds the
                                        // REMAINING servers. Zero for models without replies.
                                        Matrix nbD = ReplyBlock.blockedTotal(sn, ind, spaceVar);
                                        Matrix enWbuf = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            double SeffD = S.get(ist) - (row < nbD.getNumRows() ? nbD.get(row, 0) : 0);
                                            if (en.get(row, 0) == 1 && ni.get(row) > SeffD) {
                                                enWbuf.set(row, 0, 1);
                                            } else {
                                                enWbuf.set(row, 0, 0);
                                            }
                                        }
                                        boolean isRetrialStationDep = false;
                                        if (sn.retrialProc != null) {
                                            Map<JobClass, jline.util.matrix.MatrixCell> rpD = sn.retrialProc.get(sn.stations.get(ist));
                                            if (rpD != null) {
                                                for (JobClass jcD : rpD.keySet()) {
                                                    jline.util.matrix.MatrixCell mcD = rpD.get(jcD);
                                                    if (mcD != null && mcD.size() > 0) { isRetrialStationDep = true; break; }
                                                }
                                            }
                                        }
                                        if (noPromote || isRetrialStationDep) { // immediate feedback / retrial orbit: hold server, do not promote a waiting job
                                            enWbuf.zero();
                                        }
                                        // Synchronous call: this departing job keeps its server
                                        // until its REPLY signal returns here, so the server is
                                        // NOT handed to a waiting job; it is recorded as held in
                                        // the reply block instead. Mirrors LDES, which omits the
                                        // markServerIdle call and records a pendingReply.
                                        if (sn.replyblock != null && !sn.replyblock.isEmpty()
                                                && ind < sn.replyblock.getNumRows()
                                                && sn.replyblock.get(ind, jobClass) > 0) {
                                            ReplyBlock.Info rinfoD = ReplyBlock.info(sn, ind);
                                            enWbuf.zero();
                                            int slotD = rinfoD.slot[jobClass];
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                if (en.get(row, 0) == 1) {
                                                    spaceVar.set(row, slotD, spaceVar.get(row, slotD) + 1);
                                                }
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
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    // must multiply w rate_kd(en_wobuf, :)

                                                    Matrix outrate_bottom_fcfs = Matrix.scaleMult(rate_kd_no_jobs, cdscalingIst * lld);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_fcfs, null);
                                                } else {
                                                    // set outrate = [outrate; cdscaling{ist}(nir).*lldscaling(ist,min(ni,lldlimit)).*rate_kd(en_wobuf,:)];
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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

                                                        Matrix outrate_bottom;
                                                        if (ni.hasInfinite()) {
                                                            // hit limited load-dependence
                                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                            outrate_bottom = Matrix.scaleMult(rate_k_en, cdscalingIst * lld);
                                                        } else {
                                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0) - 1, lldscaling.getNumCols() - 1));
                                                            outrate_bottom = Matrix.scaleMult(rate_k_en, cdscalingIst * lld);
                                                        }
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
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

                                                        // Zero the probability of the branches added *by this kentry* when their
                                                        // rate is zero. The mask must index outrate_bottom (this kentry's new
                                                        // rates), not the whole accumulated outrate: outrate grows by one entry
                                                        // per kentry while outprob_cur has one entry per enabled state, so
                                                        // masking on outrate expanded outprob_cur once any earlier branch had
                                                        // rate 0 -- which happens whenever pentry_svc_class has zeros, i.e. a PH
                                                        // whose entry vector pie does not reach every phase. That misaligned
                                                        // outprob against outspace/outrate, so the sampled branch read a bogus
                                                        // (often 0) probability and the departure rate under-counted.
                                                        for (int row = 0; row < outrate_bottom.getNumRows() && row < outprob_cur.getNumElements(); row++) {
                                                            if (outrate_bottom.get(row) == 0) {
                                                                outprob_cur.set(row, 0);
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
                                    case FCFSPRIO: // FCFS priority (alias)
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
                                        if (noPromote) { // immediate feedback: hold server, do not promote a waiting job
                                            enWbufHol.zero();
                                            enWobufHol.ones();
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
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom_hol = Matrix.scaleMult(rateEnWobufHol, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom_hol, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrate_bottom_hol_buf = Matrix.scaleMult(rateKHol, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_hol_buf, null);
                                                    } else {
                                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                            InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnDps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnDps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_dps = new Matrix(rateEnDps.getNumRows(), 1);
                                        outprob_bottom_dps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_dps, null);
                                        break;

                                    case DPSPRIO: {
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
                                                InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nirprio, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnDpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nirprio, jobClass);
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
                                    }

                                    case GPS:
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        if (S.get(ist) > 1) {
                                            InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnGps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        } else {
                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldscaling.getNumCols() - 1));
                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnGps, cdscalingIst * lld);
                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                        }
                                        Matrix outprob_bottom_gps = new Matrix(rateEnGps.getNumRows(), 1);
                                        outprob_bottom_gps.ones();
                                        outprob = Matrix.concatRows(outprob, outprob_bottom_gps, null);
                                        break;

                                    case GPSPRIO: {
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
                                                InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nirprio, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom = Matrix.scaleMult(rateEnGpsPrio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nirprio, jobClass);
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
                                    }

                                    case LCFSPRIO: {
                                        // LCFS priority - like HOL but LCFS order within priority groups
                                        // record departure
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        Matrix kirEnClassKLcfsprio = new Matrix(0, 0);
                                        for (int l_ind = 0; l_ind < en.getNumElements(); l_ind++) {
                                            if (en.get(l_ind) == 1) {
                                                if (kirEnClassKLcfsprio.isEmpty()) {
                                                    kirEnClassKLcfsprio = new Matrix(1, 1);
                                                    kirEnClassKLcfsprio.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                } else {
                                                    Matrix new_elem = new Matrix(1, 1);
                                                    new_elem.set(0, 0, kir.get(k).get(l_ind, jobClass));
                                                    kirEnClassKLcfsprio = Matrix.concatRows(kirEnClassKLcfsprio, new_elem, null);
                                                }
                                            }
                                        }

                                        for (int l_ind = 0; l_ind < kirEnClassKLcfsprio.getNumRows(); l_ind++) {
                                            if (rate.isEmpty()) {
                                                rate = new Matrix(1, 1);
                                                rate.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKLcfsprio.get(l_ind));
                                            } else {
                                                Matrix new_elem = new Matrix(1, 1);
                                                new_elem.set(0, 0, mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k) * kirEnClassKLcfsprio.get(l_ind));
                                                if (l_ind < rate.getNumElements()) {
                                                    rate.set(l_ind, new_elem.value());
                                                } else {
                                                    rate = Matrix.concatRows(rate, new_elem, null);
                                                }
                                            }
                                        }

                                        // set en_wbuf to states with jobs in buffer
                                        Matrix enWbufLcfsprio = new Matrix(en.getNumRows(), 1);
                                        Matrix enWobufLcfsprio = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufLcfsprio.set(row, 0, 1);
                                                enWobufLcfsprio.set(row, 0, 0);
                                            } else {
                                                enWbufLcfsprio.set(row, 0, 0);
                                                enWobufLcfsprio.set(row, 0, 1);
                                            }
                                        }
                                        if (noPromote) { // immediate feedback: hold server, do not promote a waiting job
                                            enWbufLcfsprio.zero();
                                            enWobufLcfsprio.ones();
                                        }

                                        // Process states without jobs in buffer first
                                        Matrix rateEnWobufLcfsprio = new Matrix(0, 0);
                                        for (int row = 0; row < enWobufLcfsprio.getNumRows(); row++) {
                                            if (enWobufLcfsprio.get(row, 0) == 1) {
                                                if (rateEnWobufLcfsprio.isEmpty()) {
                                                    rateEnWobufLcfsprio = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnWobufLcfsprio = Matrix.concatRows(rateEnWobufLcfsprio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        Matrix spaceBufEnWobufLcfsprio = new Matrix(0, 0);
                                        Matrix spaceSrvEnWobufLcfsprio = new Matrix(0, 0);
                                        Matrix spaceVarEnWobufLcfsprio = new Matrix(0, 0);
                                        for (int row = 0; row < enWobufLcfsprio.getNumRows(); row++) {
                                            if (enWobufLcfsprio.get(row, 0) == 1) {
                                                if (spaceBufEnWobufLcfsprio.isEmpty()) {
                                                    spaceBufEnWobufLcfsprio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                    spaceSrvEnWobufLcfsprio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                    spaceVarEnWobufLcfsprio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                } else {
                                                    spaceBufEnWobufLcfsprio = Matrix.concatRows(spaceBufEnWobufLcfsprio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                    spaceSrvEnWobufLcfsprio = Matrix.concatRows(spaceSrvEnWobufLcfsprio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                    spaceVarEnWobufLcfsprio = Matrix.concatRows(spaceVarEnWobufLcfsprio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (!spaceBufEnWobufLcfsprio.isEmpty()) {
                                            Matrix left_bottom_lcfsprio = Matrix.concatColumns(spaceBufEnWobufLcfsprio, spaceSrvEnWobufLcfsprio, null);
                                            Matrix bottom_lcfsprio = Matrix.concatColumns(left_bottom_lcfsprio, spaceVarEnWobufLcfsprio, null);
                                            outspace = Matrix.concatRows(outspace, bottom_lcfsprio, null);

                                            if (ni.hasInfinite()) {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_bottom_lcfsprio = Matrix.scaleMult(rateEnWobufLcfsprio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom_lcfsprio, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                Matrix outrate_bottom_lcfsprio = Matrix.scaleMult(rateEnWobufLcfsprio, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_bottom_lcfsprio, null);
                                            }

                                            Matrix outprob_bottom_lcfsprio = new Matrix(rateEnWobufLcfsprio.getNumRows(), 1);
                                            outprob_bottom_lcfsprio.ones();
                                            outprob = Matrix.concatRows(outprob, outprob_bottom_lcfsprio, null);
                                        }

                                        // Process states with jobs in buffer using LCFSPRIO logic
                                        boolean any_jobs_in_buffer_lcfsprio = false;
                                        for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                            if (enWbufLcfsprio.get(row, 0) == 1) {
                                                any_jobs_in_buffer_lcfsprio = true;
                                                break;
                                            }
                                        }

                                        if (any_jobs_in_buffer_lcfsprio) {
                                            // Create priority group matrix (Inf for empty, classprio for occupied)
                                            Matrix priogroupLcfsprio = new Matrix(1, R + 1);
                                            priogroupLcfsprio.set(0, 0, Double.POSITIVE_INFINITY); // empty slot
                                            for (int r = 0; r < R; r++) {
                                                priogroupLcfsprio.set(0, r + 1, sn.classprio.get(r));
                                            }

                                            // Transform buffer to priority groups
                                            Matrix spaceBufGroupLcfsprio = new Matrix(spaceBuf.getNumRows(), spaceBuf.getNumCols());
                                            for (int row = 0; row < spaceBuf.getNumRows(); row++) {
                                                for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                    int bufVal = (int) spaceBuf.get(row, col);
                                                    spaceBufGroupLcfsprio.set(row, col, priogroupLcfsprio.get(0, bufVal));
                                                }
                                            }

                                            // Find minimum priority per row (highest priority class in buffer)
                                            // Then find LEFTMOST position with that priority (LCFS order)
                                            Matrix startSvcClassLcfsprio = new Matrix(0, 0);
                                            Matrix leftmostMinPosLcfsprio = new Matrix(0, 0);

                                            for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                                if (enWbufLcfsprio.get(row, 0) == 1) {
                                                    // Find min priority in this row
                                                    double minPrioLcfsprio = Double.POSITIVE_INFINITY;
                                                    for (int col = 0; col < spaceBufGroupLcfsprio.getNumCols(); col++) {
                                                        double prio = spaceBufGroupLcfsprio.get(row, col);
                                                        if (prio < minPrioLcfsprio) {
                                                            minPrioLcfsprio = prio;
                                                        }
                                                    }

                                                    // Find leftmost position with min priority (LCFS = most recent)
                                                    int leftmostPos = -1;
                                                    for (int col = 0; col < spaceBufGroupLcfsprio.getNumCols(); col++) {
                                                        if (spaceBufGroupLcfsprio.get(row, col) == minPrioLcfsprio) {
                                                            leftmostPos = col;
                                                            break;
                                                        }
                                                    }

                                                    if (leftmostMinPosLcfsprio.isEmpty()) {
                                                        leftmostMinPosLcfsprio = new Matrix(1, 1);
                                                        leftmostMinPosLcfsprio.set(0, 0, leftmostPos);
                                                        startSvcClassLcfsprio = new Matrix(1, 1);
                                                        startSvcClassLcfsprio.set(0, 0, leftmostPos >= 0 ? spaceBuf.get(row, leftmostPos) : 0);
                                                    } else {
                                                        Matrix new_pos = new Matrix(1, 1);
                                                        new_pos.set(0, 0, leftmostPos);
                                                        leftmostMinPosLcfsprio = Matrix.concatRows(leftmostMinPosLcfsprio, new_pos, null);
                                                        Matrix new_class = new Matrix(1, 1);
                                                        new_class.set(0, 0, leftmostPos >= 0 ? spaceBuf.get(row, leftmostPos) : 0);
                                                        startSvcClassLcfsprio = Matrix.concatRows(startSvcClassLcfsprio, new_class, null);
                                                    }
                                                }
                                            }

                                            // Check if there is a valid job to start
                                            boolean hasValidStartLcfsprio = false;
                                            int startClassLcfsprio = 0;
                                            if (!startSvcClassLcfsprio.isEmpty()) {
                                                startClassLcfsprio = (int) startSvcClassLcfsprio.get(0, 0);
                                                if (startClassLcfsprio > 0) {
                                                    hasValidStartLcfsprio = true;
                                                }
                                            }

                                            if (hasValidStartLcfsprio) {
                                                int startClassIdxLcfsprio = startClassLcfsprio - 1; // Convert to 0-based
                                                Matrix pentrySvcClassLcfsprio = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(startClassIdxLcfsprio));

                                                for (int kentry = 0; kentry < K.get(startClassIdxLcfsprio); kentry++) {
                                                    Matrix spaceSrvKLcfsprio = spaceSrv.copy();
                                                    Matrix spaceBufKLcfsprio = spaceBuf.copy();

                                                    // Add job to service
                                                    int bufRowIdx = 0;
                                                    for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                                        if (enWbufLcfsprio.get(row, 0) == 1) {
                                                            spaceSrvKLcfsprio.set(row, (int) (Ks.get(startClassIdxLcfsprio) + kentry),
                                                                    spaceSrvKLcfsprio.get(row, (int) (Ks.get(startClassIdxLcfsprio) + kentry)) + 1);
                                                            // Remove from leftmost position, shift buffer
                                                            int pos = (int) leftmostMinPosLcfsprio.get(bufRowIdx, 0);
                                                            if (pos >= 0) {
                                                                // Shift: [0, buf[0..pos-1], buf[pos+1..end]]
                                                                Matrix newBufRow = new Matrix(1, spaceBufKLcfsprio.getNumCols());
                                                                newBufRow.set(0, 0, 0); // prepend empty slot
                                                                int destCol = 1;
                                                                for (int col = 0; col < spaceBufKLcfsprio.getNumCols(); col++) {
                                                                    if (col != pos) {
                                                                        if (destCol < newBufRow.getNumCols()) {
                                                                            newBufRow.set(0, destCol, spaceBufKLcfsprio.get(row, col));
                                                                            destCol++;
                                                                        }
                                                                    }
                                                                }
                                                                for (int col = 0; col < spaceBufKLcfsprio.getNumCols(); col++) {
                                                                    spaceBufKLcfsprio.set(row, col, newBufRow.get(0, col));
                                                                }
                                                            }
                                                            bufRowIdx++;
                                                        }
                                                    }

                                                    // Build output for buffer states
                                                    Matrix spaceBufEnWbufLcfsprio = new Matrix(0, 0);
                                                    Matrix spaceSrvEnWbufLcfsprio = new Matrix(0, 0);
                                                    Matrix spaceVarEnWbufLcfsprio = new Matrix(0, 0);
                                                    for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                                        if (enWbufLcfsprio.get(row, 0) == 1) {
                                                            if (spaceBufEnWbufLcfsprio.isEmpty()) {
                                                                spaceBufEnWbufLcfsprio = Matrix.extractRows(spaceBufKLcfsprio, row, row + 1, null);
                                                                spaceSrvEnWbufLcfsprio = Matrix.extractRows(spaceSrvKLcfsprio, row, row + 1, null);
                                                                spaceVarEnWbufLcfsprio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                            } else {
                                                                spaceBufEnWbufLcfsprio = Matrix.concatRows(spaceBufEnWbufLcfsprio, Matrix.extractRows(spaceBufKLcfsprio, row, row + 1, null), null);
                                                                spaceSrvEnWbufLcfsprio = Matrix.concatRows(spaceSrvEnWbufLcfsprio, Matrix.extractRows(spaceSrvKLcfsprio, row, row + 1, null), null);
                                                                spaceVarEnWbufLcfsprio = Matrix.concatRows(spaceVarEnWbufLcfsprio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    if (!spaceBufEnWbufLcfsprio.isEmpty()) {
                                                        Matrix left_bottom_wbuf_lcfsprio = Matrix.concatColumns(spaceBufEnWbufLcfsprio, spaceSrvEnWbufLcfsprio, null);
                                                        Matrix bottom_wbuf_lcfsprio = Matrix.concatColumns(left_bottom_wbuf_lcfsprio, spaceVarEnWbufLcfsprio, null);
                                                        outspace = Matrix.concatRows(outspace, bottom_wbuf_lcfsprio, null);

                                                        // Apply pie probability to rate
                                                        Matrix rateKLcfsprio = rate.copy();
                                                        for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                                            if (enWbufLcfsprio.get(row, 0) == 1) {
                                                                rateKLcfsprio.set(row, 0, rateKLcfsprio.get(row, 0) * pentrySvcClassLcfsprio.get(kentry));
                                                            }
                                                        }

                                                        Matrix rateEnWbufLcfsprio = new Matrix(0, 0);
                                                        for (int row = 0; row < enWbufLcfsprio.getNumRows(); row++) {
                                                            if (enWbufLcfsprio.get(row, 0) == 1) {
                                                                if (rateEnWbufLcfsprio.isEmpty()) {
                                                                    rateEnWbufLcfsprio = Matrix.extractRows(rateKLcfsprio, row, row + 1, null);
                                                                } else {
                                                                    rateEnWbufLcfsprio = Matrix.concatRows(rateEnWbufLcfsprio, Matrix.extractRows(rateKLcfsprio, row, row + 1, null), null);
                                                                }
                                                            }
                                                        }

                                                        if (ni.hasInfinite()) {
                                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                            double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnWbufLcfsprio, cdscalingIst * lld);
                                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                        } else {
                                                            double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                            double lld = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit));
                                                            Matrix outrate_bottom = Matrix.scaleMult(rateEnWbufLcfsprio, cdscalingIst * lld);
                                                            outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                                        }

                                                        Matrix outprob_bottom_wbuf_lcfsprio = new Matrix(rateEnWbufLcfsprio.getNumRows(), 1);
                                                        outprob_bottom_wbuf_lcfsprio.ones();
                                                        outprob = Matrix.concatRows(outprob, outprob_bottom_wbuf_lcfsprio, null);
                                                    }
                                                }
                                            }
                                        }
                                        break;
                                    }

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
                                        if (noPromote) { // immediate feedback: hold server, do not promote a waiting job
                                            enWbufLcfs.zero();
                                        }

                                        // Find column to promote from buffer.
                                        // Plain LCFS is not priority-aware: it always promotes the
                                        // most recent arrival. Priorities are honored by
                                        // SchedStrategy.LCFSPRIO, which has its own case, exactly
                                        // as FCFS relates to FCFSPRIO. Branching here on whether
                                        // the class priorities differ silently turned every LCFS
                                        // station with distinct priorities into an LCFSPRIO one.
                                        final boolean hasDiffPrioDepLcfs = false;
                                        Matrix colFirstNnz = new Matrix(0, 0);
                                        Matrix startSvcClassLcfs = new Matrix(0, 0);

                                        for (int row = 0; row < enWbufLcfs.getNumRows(); row++) {
                                            if (enWbufLcfs.get(row, 0) == 1) {
                                                int firstCol = -1;
                                                if (hasDiffPrioDepLcfs) {
                                                    // Priority-aware: find leftmost among highest-priority class
                                                    double bestPrio = Double.MAX_VALUE;
                                                    for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                        if (spaceBuf.get(row, col) != 0) {
                                                            int cls = (int) spaceBuf.get(row, col) - 1; // 0-based class
                                                            double prio = sn.classprio.get(cls);
                                                            if (prio < bestPrio) {
                                                                bestPrio = prio;
                                                                firstCol = col; // leftmost among best priority
                                                            }
                                                        }
                                                    }
                                                } else {
                                                    // Default: find first non-zero column (leftmost = most recent)
                                                    for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                        if (spaceBuf.get(row, col) != 0) {
                                                            firstCol = col;
                                                            break;
                                                        }
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
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrate_wobuf = Matrix.scaleMult(rateWobuf, cdscalingIst * lld);
                                                outrate = Matrix.concatRows(outrate, outrate_wobuf, null);
                                            } else {
                                                double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrate_bottom_lcfs = Matrix.scaleMult(rateEnLcfs, cdscalingIst * lld);
                                                    outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs, null);
                                                } else {
                                                    double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrate_bottom_lcfs_k = Matrix.scaleMult(rateKLcfs, cdscalingIst * lld);
                                                        outrate = Matrix.concatRows(outrate, outrate_bottom_lcfs_k, null);
                                                    } else {
                                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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

                                        // Find column to promote from buffer
                                        // Plain LCFSPR resumes the most recently preempted job;
                                        // LCFSPRPRIO is the priority-aware variant and carries its
                                        // own handling. See the arrival path for the same rule.
                                        final boolean hasDiffPrioDepLcfspr = false;
                                        Matrix colFirstNnzLcfspr = new Matrix(0, 0);
                                        Matrix startSvcClassLcfspr = new Matrix(0, 0);
                                        Matrix kentryLcfspr = new Matrix(0, 0);

                                        for (int row = 0; row < enWbufLcfspr.getNumRows(); row++) {
                                            if (enWbufLcfspr.get(row, 0) == 1) {
                                                int firstCol = -1;
                                                if (hasDiffPrioDepLcfspr) {
                                                    // Priority-aware: find leftmost pair among highest-priority class
                                                    double bestPrio = Double.MAX_VALUE;
                                                    for (int col = 0; col < spaceBuf.getNumCols() - 1; col += 2) {
                                                        if (spaceBuf.get(row, col) != 0) {
                                                            int cls = (int) spaceBuf.get(row, col) - 1; // 0-based class
                                                            double prio = sn.classprio.get(cls);
                                                            if (prio < bestPrio) {
                                                                bestPrio = prio;
                                                                firstCol = col; // leftmost among best priority
                                                            }
                                                        }
                                                    }
                                                } else {
                                                    // Default: find first non-zero column (leftmost = most recent)
                                                    for (int col = 0; col < spaceBuf.getNumCols(); col++) {
                                                        if (spaceBuf.get(row, col) != 0) {
                                                            firstCol = col;
                                                            break;
                                                        }
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
                                                    double cdscalingIstLcfspr = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldLcfspr = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                                } else {
                                                    double cdscalingIstLcfspr = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                            double cdscalingIstLcfspr = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lldLcfspr = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                            outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                        } else {
                                            double cdscalingIstLcfspr = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lldLcfspr = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                            Matrix outrateBottomLcfspr = Matrix.scaleMult(rateEnLcfspr, cdscalingIstLcfspr * lldLcfspr);
                                            outrate = Matrix.concatRows(outrate, outrateBottomLcfspr, null);
                                        }

                                        Matrix outprobBottomLcfspr = new Matrix(rateEnLcfspr.getNumRows(), 1);
                                        outprobBottomLcfspr.ones();
                                        outprob = Matrix.concatRows(outprob, outprobBottomLcfspr, null);
                                        break;

                                    case LCFSPRPRIO:
                                    case FCFSPRPRIO: {
                                        // LCFSPRPRIO/FCFSPRPRIO departure - like LCFSPR but with priority-based buffer selection
                                        SchedStrategy schedIstDep = sn.sched.get(sn.stations.get(ist));

                                        // Record departure from service
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // Calculate rate for enabled states
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

                                        Matrix enWbufPrio = new Matrix(en.getNumRows(), 1);
                                        // States with jobs in buffer
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufPrio.set(row, 0, 1);
                                            } else {
                                                enWbufPrio.set(row, 0, 0);
                                            }
                                        }

                                        // Build priority array
                                        double[] classprio = new double[R];
                                        for (int r = 0; r < R; r++) {
                                            classprio[r] = sn.classprio.get(r);
                                        }

                                        // Find highest-priority job in buffer for each enabled row
                                        Matrix colTargetPrio = new Matrix(0, 0);
                                        Matrix startSvcClassPrio = new Matrix(0, 0);
                                        Matrix kentryPrio = new Matrix(0, 0);

                                        for (int row = 0; row < enWbufPrio.getNumRows(); row++) {
                                            if (enWbufPrio.get(row, 0) == 1) {
                                                double minPrio = Double.POSITIVE_INFINITY;
                                                int targetBufCol = -1;

                                                if (schedIstDep == SchedStrategy.LCFSPRPRIO) {
                                                    // LCFSPRPRIO: scan left-to-right, use < so leftmost (most recently preempted) wins
                                                    for (int col = 0; col < spaceBuf.getNumCols(); col += 2) {
                                                        double bufVal = spaceBuf.get(row, col);
                                                        if (bufVal > 0) {
                                                            int cls = (int) bufVal - 1; // 0-based class index
                                                            double prio = classprio[cls];
                                                            if (prio < minPrio) {
                                                                minPrio = prio;
                                                                targetBufCol = col;
                                                            }
                                                        }
                                                    }
                                                } else {
                                                    // FCFSPRPRIO: scan right-to-left, first match = rightmost (longest waiting)
                                                    for (int col = spaceBuf.getNumCols() - 2; col >= 0; col -= 2) {
                                                        double bufVal = spaceBuf.get(row, col);
                                                        if (bufVal > 0) {
                                                            int cls = (int) bufVal - 1; // 0-based class index
                                                            double prio = classprio[cls];
                                                            if (prio < minPrio) {
                                                                minPrio = prio;
                                                                targetBufCol = col;
                                                            }
                                                        }
                                                    }
                                                }

                                                if (colTargetPrio.isEmpty()) {
                                                    colTargetPrio = new Matrix(1, 1);
                                                    colTargetPrio.set(0, 0, targetBufCol);
                                                    startSvcClassPrio = new Matrix(1, 1);
                                                    kentryPrio = new Matrix(1, 1);
                                                    if (targetBufCol >= 0) {
                                                        startSvcClassPrio.set(0, 0, spaceBuf.get(row, targetBufCol)); // class (1-based)
                                                        kentryPrio.set(0, 0, spaceBuf.get(row, targetBufCol + 1)); // phase (1-based)
                                                    } else {
                                                        startSvcClassPrio.set(0, 0, 0);
                                                        kentryPrio.set(0, 0, 0);
                                                    }
                                                } else {
                                                    Matrix new_col = new Matrix(1, 1);
                                                    new_col.set(0, 0, targetBufCol);
                                                    colTargetPrio = Matrix.concatRows(colTargetPrio, new_col, null);
                                                    Matrix new_class = new Matrix(1, 1);
                                                    Matrix new_phase = new Matrix(1, 1);
                                                    if (targetBufCol >= 0) {
                                                        new_class.set(0, 0, spaceBuf.get(row, targetBufCol));
                                                        new_phase.set(0, 0, spaceBuf.get(row, targetBufCol + 1));
                                                    } else {
                                                        new_class.set(0, 0, 0);
                                                        new_phase.set(0, 0, 0);
                                                    }
                                                    startSvcClassPrio = Matrix.concatRows(startSvcClassPrio, new_class, null);
                                                    kentryPrio = Matrix.concatRows(kentryPrio, new_phase, null);
                                                }
                                            }
                                        }

                                        // Remove job from buffer (set both class and phase to 0)
                                        Matrix spaceBufPrio = spaceBuf.copy();
                                        int bufRowIdxPrio = 0;
                                        for (int row = 0; row < enWbufPrio.getNumRows(); row++) {
                                            if (enWbufPrio.get(row, 0) == 1) {
                                                int targetCol = (int) colTargetPrio.get(bufRowIdxPrio, 0);
                                                if (targetCol >= 0) {
                                                    spaceBufPrio.set(row, targetCol, 0); // zero popped job class
                                                    spaceBufPrio.set(row, targetCol + 1, 0); // zero popped phase
                                                }
                                                bufRowIdxPrio++;
                                            }
                                        }

                                        // Check if we have valid jobs to start
                                        if (!startSvcClassPrio.isEmpty() && startSvcClassPrio.getNumRows() > 0) {
                                            boolean hasValidStartClassPrio = false;
                                            for (int i = 0; i < startSvcClassPrio.getNumRows(); i++) {
                                                if (startSvcClassPrio.get(i, 0) > 0) {
                                                    hasValidStartClassPrio = true;
                                                    break;
                                                }
                                            }

                                            if (!hasValidStartClassPrio) {
                                                // No valid job to start, just add current state
                                                Matrix spaceBufEnPrio = new Matrix(0, 0);
                                                Matrix spaceSrvEnPrio = new Matrix(0, 0);
                                                Matrix spaceVarEnPrio = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (spaceBufEnPrio.isEmpty()) {
                                                            spaceBufEnPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                            spaceSrvEnPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                            spaceVarEnPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                        } else {
                                                            spaceBufEnPrio = Matrix.concatRows(spaceBufEnPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                            spaceSrvEnPrio = Matrix.concatRows(spaceSrvEnPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                            spaceVarEnPrio = Matrix.concatRows(spaceVarEnPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                Matrix leftBottomPrio = Matrix.concatColumns(spaceBufEnPrio, spaceSrvEnPrio, null);
                                                Matrix bottomPrio = Matrix.concatColumns(leftBottomPrio, spaceVarEnPrio, null);
                                                outspace = Matrix.concatRows(outspace, bottomPrio, null);

                                                Matrix rateEnPrio = new Matrix(0, 0);
                                                for (int row = 0; row < en.getNumRows(); row++) {
                                                    if (en.get(row, 0) == 1) {
                                                        if (rateEnPrio.isEmpty()) {
                                                            rateEnPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                        } else {
                                                            rateEnPrio = Matrix.concatRows(rateEnPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                if (ni.hasInfinite()) {
                                                    double cdscalingIstPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldPrio = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomPrio = Matrix.scaleMult(rateEnPrio, cdscalingIstPrio * lldPrio);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomPrio, null);
                                                } else {
                                                    double cdscalingIstPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldPrio = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                    Matrix outrateBottomPrio = Matrix.scaleMult(rateEnPrio, cdscalingIstPrio * lldPrio);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomPrio, null);
                                                }

                                                Matrix outprobBottomPrio = new Matrix(rateEnPrio.getNumRows(), 1);
                                                outprobBottomPrio.ones();
                                                outprob = Matrix.concatRows(outprob, outprobBottomPrio, null);

                                                if (isSimulation && eventCache.isEnabled()) {
                                                    eventCache.put(key, new Ret.EventResult(outspace, outrate, outprob));
                                                }
                                                return new Ret.EventResult(outspace, outrate, outprob);
                                            }

                                            // Add job to service with preserved phase
                                            bufRowIdxPrio = 0;
                                            for (int row = 0; row < enWbufPrio.getNumRows(); row++) {
                                                if (enWbufPrio.get(row, 0) == 1) {
                                                    int startClass = (int) startSvcClassPrio.get(bufRowIdxPrio, 0) - 1; // Convert to 0-based index
                                                    int kentry = (int) kentryPrio.get(bufRowIdxPrio, 0);
                                                    if (startClass >= 0 && startClass < Ks.length()) {
                                                        int colIndex = (int) (Ks.get(startClass) + kentry - 1);
                                                        if (colIndex >= 0 && colIndex < spaceSrv.getNumCols()) {
                                                            spaceSrv.set(row, colIndex, spaceSrv.get(row, colIndex) + 1);
                                                        }
                                                    }
                                                    bufRowIdxPrio++;
                                                }
                                            }
                                        }

                                        // Add state to output
                                        Matrix leftBottomPrio = Matrix.concatColumns(spaceBufPrio, spaceSrv, null);
                                        Matrix bottomPrio = Matrix.concatColumns(leftBottomPrio, spaceVar, null);

                                        Matrix spaceBufEnPrio = new Matrix(0, 0);
                                        Matrix spaceSrvEnPrio = new Matrix(0, 0);
                                        Matrix spaceVarEnPrio = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (spaceBufEnPrio.isEmpty()) {
                                                    spaceBufEnPrio = Matrix.extractRows(bottomPrio, row, row + 1, null);
                                                } else {
                                                    spaceBufEnPrio = Matrix.concatRows(spaceBufEnPrio, Matrix.extractRows(bottomPrio, row, row + 1, null), null);
                                                }
                                            }
                                        }
                                        outspace = Matrix.concatRows(outspace, spaceBufEnPrio, null);

                                        Matrix rateEnPrio = new Matrix(0, 0);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                if (rateEnPrio.isEmpty()) {
                                                    rateEnPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                } else {
                                                    rateEnPrio = Matrix.concatRows(rateEnPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                }
                                            }
                                        }

                                        if (ni.hasInfinite()) {
                                            double cdscalingIstPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lldPrio = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                            Matrix outrateBottomPrio = Matrix.scaleMult(rateEnPrio, cdscalingIstPrio * lldPrio);
                                            outrate = Matrix.concatRows(outrate, outrateBottomPrio, null);
                                        } else {
                                            double cdscalingIstPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                            double lldPrio = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                            Matrix outrateBottomPrio = Matrix.scaleMult(rateEnPrio, cdscalingIstPrio * lldPrio);
                                            outrate = Matrix.concatRows(outrate, outrateBottomPrio, null);
                                        }

                                        Matrix outprobBottomPrio = new Matrix(rateEnPrio.getNumRows(), 1);
                                        outprobBottomPrio.ones();
                                        outprob = Matrix.concatRows(outprob, outprobBottomPrio, null);
                                        break;
                                    }

                                    case FCFSPIPRIO:
                                    case LCFSPIPRIO: {
                                        // FCFSPIPRIO/LCFSPIPRIO departure - like FCFSPRPRIO/LCFSPRPRIO but preempt-independent (restart from pie)
                                        SchedStrategy schedIstPiPrio = sn.sched.get(sn.stations.get(ist));

                                        // Record departure from service
                                        for (int row = 0; row < spaceSrv.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1) {
                                                spaceSrv.set(row, (int) (Ks.get(jobClass) + k), spaceSrv.get(row, (int) (Ks.get(jobClass) + k)) - 1);
                                            }
                                        }

                                        // Calculate rate for enabled states
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

                                        Matrix enWbufPiPrio = new Matrix(en.getNumRows(), 1);
                                        Matrix enWobufPiPrio = new Matrix(en.getNumRows(), 1);
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) == 1 && ni.get(row) > S.get(ist)) {
                                                enWbufPiPrio.set(row, 0, 1);
                                                enWobufPiPrio.set(row, 0, 0);
                                            } else {
                                                enWbufPiPrio.set(row, 0, 0);
                                                enWobufPiPrio.set(row, 0, 1);
                                            }
                                        }

                                        // Handle states without buffer jobs
                                        if (enWobufPiPrio.elementSum() > 0) {
                                            Matrix spaceBufEnWobufPiPrio = new Matrix(0, 0);
                                            Matrix spaceSrvEnWobufPiPrio = new Matrix(0, 0);
                                            Matrix spaceVarEnWobufPiPrio = new Matrix(0, 0);
                                            for (int row = 0; row < enWobufPiPrio.getNumRows(); row++) {
                                                if (enWobufPiPrio.get(row, 0) == 1) {
                                                    if (spaceBufEnWobufPiPrio.isEmpty()) {
                                                        spaceBufEnWobufPiPrio = Matrix.extractRows(spaceBuf, row, row + 1, null);
                                                        spaceSrvEnWobufPiPrio = Matrix.extractRows(spaceSrv, row, row + 1, null);
                                                        spaceVarEnWobufPiPrio = Matrix.extractRows(spaceVar, row, row + 1, null);
                                                    } else {
                                                        spaceBufEnWobufPiPrio = Matrix.concatRows(spaceBufEnWobufPiPrio, Matrix.extractRows(spaceBuf, row, row + 1, null), null);
                                                        spaceSrvEnWobufPiPrio = Matrix.concatRows(spaceSrvEnWobufPiPrio, Matrix.extractRows(spaceSrv, row, row + 1, null), null);
                                                        spaceVarEnWobufPiPrio = Matrix.concatRows(spaceVarEnWobufPiPrio, Matrix.extractRows(spaceVar, row, row + 1, null), null);
                                                    }
                                                }
                                            }
                                            if (!spaceBufEnWobufPiPrio.isEmpty()) {
                                                Matrix leftBottomWobufPiPrio = Matrix.concatColumns(spaceBufEnWobufPiPrio, spaceSrvEnWobufPiPrio, null);
                                                Matrix bottomWobufPiPrio = Matrix.concatColumns(leftBottomWobufPiPrio, spaceVarEnWobufPiPrio, null);
                                                outspace = Matrix.concatRows(outspace, bottomWobufPiPrio, null);

                                                Matrix rateEnWobufPiPrio = new Matrix(0, 0);
                                                for (int row = 0; row < enWobufPiPrio.getNumRows(); row++) {
                                                    if (enWobufPiPrio.get(row, 0) == 1) {
                                                        if (rateEnWobufPiPrio.isEmpty()) {
                                                            rateEnWobufPiPrio = Matrix.extractRows(rate, row, row + 1, null);
                                                        } else {
                                                            rateEnWobufPiPrio = Matrix.concatRows(rateEnWobufPiPrio, Matrix.extractRows(rate, row, row + 1, null), null);
                                                        }
                                                    }
                                                }

                                                if (ni.hasInfinite()) {
                                                    double cdscalingIstPiPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldPiPrio = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomPiPrio = Matrix.scaleMult(rateEnWobufPiPrio, cdscalingIstPiPrio * lldPiPrio);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomPiPrio, null);
                                                } else {
                                                    double cdscalingIstPiPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldPiPrio = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                    Matrix outrateBottomPiPrio = Matrix.scaleMult(rateEnWobufPiPrio, cdscalingIstPiPrio * lldPiPrio);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomPiPrio, null);
                                                }

                                                Matrix outprobBottomWobufPiPrio = new Matrix(rateEnWobufPiPrio.getNumRows(), 1);
                                                outprobBottomWobufPiPrio.ones();
                                                outprob = Matrix.concatRows(outprob, outprobBottomWobufPiPrio, null);
                                            }
                                        }

                                        // Handle states with buffer jobs
                                        if (enWbufPiPrio.elementSum() > 0) {
                                            // Build priority array
                                            double[] classprioPiPrio = new double[R];
                                            for (int r = 0; r < R; r++) {
                                                classprioPiPrio[r] = sn.classprio.get(r);
                                            }

                                            // Find highest-priority job in buffer for each enabled row
                                            Matrix startSvcClassPiPrio = new Matrix(0, 0);
                                            Matrix colTargetPiPrio = new Matrix(0, 0);

                                            for (int row = 0; row < enWbufPiPrio.getNumRows(); row++) {
                                                if (enWbufPiPrio.get(row, 0) == 1) {
                                                    double minPrio = Double.POSITIVE_INFINITY;
                                                    int targetBufCol = -1;

                                                    if (schedIstPiPrio == SchedStrategy.LCFSPIPRIO) {
                                                        // LCFSPIPRIO: scan left-to-right, leftmost highest-priority wins
                                                        for (int col = 0; col < spaceBuf.getNumCols(); col += 2) {
                                                            double bufVal = spaceBuf.get(row, col);
                                                            if (bufVal > 0) {
                                                                int cls = (int) bufVal - 1;
                                                                double prio = classprioPiPrio[cls];
                                                                if (prio < minPrio) {
                                                                    minPrio = prio;
                                                                    targetBufCol = col;
                                                                }
                                                            }
                                                        }
                                                    } else {
                                                        // FCFSPIPRIO: scan right-to-left, rightmost highest-priority wins
                                                        for (int col = spaceBuf.getNumCols() - 2; col >= 0; col -= 2) {
                                                            double bufVal = spaceBuf.get(row, col);
                                                            if (bufVal > 0) {
                                                                int cls = (int) bufVal - 1;
                                                                double prio = classprioPiPrio[cls];
                                                                if (prio < minPrio) {
                                                                    minPrio = prio;
                                                                    targetBufCol = col;
                                                                }
                                                            }
                                                        }
                                                    }

                                                    if (colTargetPiPrio.isEmpty()) {
                                                        colTargetPiPrio = new Matrix(1, 1);
                                                        colTargetPiPrio.set(0, 0, targetBufCol);
                                                        startSvcClassPiPrio = new Matrix(1, 1);
                                                        startSvcClassPiPrio.set(0, 0, targetBufCol >= 0 ? spaceBuf.get(row, targetBufCol) : 0);
                                                    } else {
                                                        Matrix new_col = new Matrix(1, 1);
                                                        new_col.set(0, 0, targetBufCol);
                                                        colTargetPiPrio = Matrix.concatRows(colTargetPiPrio, new_col, null);
                                                        Matrix new_class = new Matrix(1, 1);
                                                        new_class.set(0, 0, targetBufCol >= 0 ? spaceBuf.get(row, targetBufCol) : 0);
                                                        startSvcClassPiPrio = Matrix.concatRows(startSvcClassPiPrio, new_class, null);
                                                    }
                                                }
                                            }

                                            // Check if there is a valid job to start
                                            boolean hasValidStartPiPrio = false;
                                            int startClassPiPrio = 0;
                                            if (!startSvcClassPiPrio.isEmpty()) {
                                                startClassPiPrio = (int) startSvcClassPiPrio.get(0, 0);
                                                if (startClassPiPrio > 0) {
                                                    hasValidStartPiPrio = true;
                                                }
                                            }

                                            if (hasValidStartPiPrio) {
                                                int startClassIdxPiPrio = startClassPiPrio - 1; // Convert to 0-based
                                                Matrix pentrySvcClassPiPrio = pie.get(sn.stations.get(ist)).get(sn.jobclasses.get(startClassIdxPiPrio));

                                                // For preempt-independent, loop over kentry phases using pie distribution
                                                for (int kentry = 0; kentry < K.get(startClassIdxPiPrio); kentry++) {
                                                    Matrix spaceBufPiPrio = spaceBuf.copy();
                                                    Matrix spaceSrvPiPrio = spaceSrv.copy();

                                                    // Remove job from buffer (set both class and phase to 0) and add to service
                                                    int bufRowIdxPiPrio = 0;
                                                    for (int row = 0; row < enWbufPiPrio.getNumRows(); row++) {
                                                        if (enWbufPiPrio.get(row, 0) == 1) {
                                                            int targetCol = (int) colTargetPiPrio.get(bufRowIdxPiPrio, 0);
                                                            if (targetCol >= 0) {
                                                                spaceBufPiPrio.set(row, targetCol, 0); // zero class
                                                                spaceBufPiPrio.set(row, targetCol + 1, 0); // zero phase
                                                            }
                                                            // Add job to service in phase kentry (pie distribution)
                                                            spaceSrvPiPrio.set(row, (int) (Ks.get(startClassIdxPiPrio) + kentry),
                                                                    spaceSrvPiPrio.get(row, (int) (Ks.get(startClassIdxPiPrio) + kentry)) + 1);
                                                            bufRowIdxPiPrio++;
                                                        }
                                                    }

                                                    // Build output
                                                    Matrix leftBottomPiPrio = Matrix.concatColumns(spaceBufPiPrio, spaceSrvPiPrio, null);
                                                    Matrix bottomPiPrio = Matrix.concatColumns(leftBottomPiPrio, spaceVar, null);

                                                    Matrix spaceBufEnPiPrio = new Matrix(0, 0);
                                                    for (int row = 0; row < en.getNumRows(); row++) {
                                                        if (enWbufPiPrio.get(row, 0) == 1) {
                                                            if (spaceBufEnPiPrio.isEmpty()) {
                                                                spaceBufEnPiPrio = Matrix.extractRows(bottomPiPrio, row, row + 1, null);
                                                            } else {
                                                                spaceBufEnPiPrio = Matrix.concatRows(spaceBufEnPiPrio, Matrix.extractRows(bottomPiPrio, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }
                                                    outspace = Matrix.concatRows(outspace, spaceBufEnPiPrio, null);

                                                    // Apply pie probability to rates
                                                    Matrix rateKPiPrio = rate.copy();
                                                    for (int row = 0; row < enWbufPiPrio.getNumRows(); row++) {
                                                        if (enWbufPiPrio.get(row, 0) == 1) {
                                                            rateKPiPrio.set(row, 0, rateKPiPrio.get(row, 0) * pentrySvcClassPiPrio.get(kentry));
                                                        }
                                                    }

                                                    Matrix rateEnPiPrio = new Matrix(0, 0);
                                                    for (int row = 0; row < enWbufPiPrio.getNumRows(); row++) {
                                                        if (enWbufPiPrio.get(row, 0) == 1) {
                                                            if (rateEnPiPrio.isEmpty()) {
                                                                rateEnPiPrio = Matrix.extractRows(rateKPiPrio, row, row + 1, null);
                                                            } else {
                                                                rateEnPiPrio = Matrix.concatRows(rateEnPiPrio, Matrix.extractRows(rateKPiPrio, row, row + 1, null), null);
                                                            }
                                                        }
                                                    }

                                                    if (ni.hasInfinite()) {
                                                        double cdscalingIstPiPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                        double lldPiPrio = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrateBottomPiPrio = Matrix.scaleMult(rateEnPiPrio, cdscalingIstPiPrio * lldPiPrio);
                                                        outrate = Matrix.concatRows(outrate, outrateBottomPiPrio, null);
                                                    } else {
                                                        double cdscalingIstPiPrio = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                        double lldPiPrio = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                        Matrix outrateBottomPiPrio = Matrix.scaleMult(rateEnPiPrio, cdscalingIstPiPrio * lldPiPrio);
                                                        outrate = Matrix.concatRows(outrate, outrateBottomPiPrio, null);
                                                    }

                                                    Matrix outprobBottomPiPrio = new Matrix(rateEnPiPrio.getNumRows(), 1);
                                                    outprobBottomPiPrio.ones();
                                                    outprob = Matrix.concatRows(outprob, outprobBottomPiPrio, null);
                                                }
                                            }
                                        }
                                        break;
                                    }

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
                                                    double cdscalingIstLcfspi = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                    double lldLcfspi = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                    Matrix outrateBottomWobufLcfspi = Matrix.scaleMult(rateEnWobufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                    outrate = Matrix.concatRows(outrate, outrateBottomWobufLcfspi, null);
                                                } else {
                                                    double cdscalingIstLcfspi = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                                                double cdscalingIstLcfspi = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                                double lldLcfspi = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                                Matrix outrateBottomWbufLcfspi = Matrix.scaleMult(rateEnWbufLcfspi, cdscalingIstLcfspi * lldLcfspi);
                                                                outrate = Matrix.concatRows(outrate, outrateBottomWbufLcfspi, null);
                                                            } else {
                                                                double cdscalingIstLcfspi = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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

                                    case POLLING: {
                                        // A completion ends the visit unless the discipline
                                        // still allows another job of the same class to be
                                        // taken; when it ends, the server walks the cyclic
                                        // order to wherever the next tangible controller
                                        // state lies. Which of the two happens depends on the
                                        // buffer of each state row, so the rows are resolved
                                        // one by one rather than as a vectorized promotion.
                                        Polling.Info pinfoD = Polling.info(sn, ind);
                                        Matrix spaceSrvD = Matrix.extract(inspace, 0, inspace.getNumRows(), (int) (inspace.getNumCols() - K.elementSum() - V), (int) (inspace.getNumCols() - V));
                                        Matrix spaceBufD = Matrix.extract(inspace, 0, inspace.getNumRows(), 0, (int) (inspace.getNumCols() - K.elementSum() - V));
                                        // spaceVar is NOT re-sliced from inspace: it already
                                        // carries the RROBIN/WRROBIN pointer advance applied
                                        // upstream, and re-slicing would discard it.
                                        for (int row = 0; row < en.getNumRows(); row++) {
                                            if (en.get(row, 0) != 1) {
                                                continue;
                                            }
                                            double rateD = mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k)
                                                    * phi.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(k)
                                                    * kir.get(k).get(row, jobClass);
                                            if (rateD <= 0) {
                                                continue;
                                            }
                                            Matrix varRowD = Matrix.extractRows(spaceVar, row, row + 1, null);
                                            int[] ctlD = Polling.get(pinfoD, varRowD, jobClass);
                                            if (ctlD[1] != 0) {
                                                continue; // no job can complete while the server is walking
                                            }
                                            Matrix bufRowD = Matrix.extractRows(spaceBufD, row, row + 1, null);
                                            Matrix srvRowD = Matrix.extractRows(spaceSrvD, row, row + 1, null);
                                            srvRowD.set(0, (int) (Ks.get(jobClass) + k), srvRowD.get(0, (int) (Ks.get(jobClass) + k)) - 1);
                                            int[] nbufD = new int[R];
                                            for (int r = 0; r < R; r++) {
                                                nbufD[r] = (int) Math.round(bufRowD.get(0, r));
                                            }
                                            int ctrnextD;
                                            boolean goonD;
                                            switch (pinfoD.ptype) {
                                                case EXHAUSTIVE:
                                                    ctrnextD = 0;
                                                    goonD = nbufD[jobClass] > 0;
                                                    break;
                                                case GATED:
                                                    ctrnextD = ctlD[2] - 1; // one of the gated jobs completed
                                                    goonD = ctrnextD > 0;
                                                    break;
                                                case KLIMITED:
                                                    ctrnextD = ctlD[2] - 1; // one of the K permitted services used
                                                    goonD = ctrnextD > 0 && nbufD[jobClass] > 0;
                                                    break;
                                                case DECREMENTING:
                                                    ctrnextD = ctlD[2]; // the target level is fixed for the visit
                                                    goonD = nbufD[jobClass] > ctlD[2];
                                                    break;
                                                default:
                                                    throw new RuntimeException("Unsupported polling type: " + pinfoD.ptype);
                                            }
                                            int qD, modeD, budgetD;
                                            if (goonD) {
                                                qD = jobClass; modeD = Polling.MODE_VISIT; budgetD = ctrnextD;
                                            } else {
                                                int[] resD = Polling.next(pinfoD, jobClass, nbufD, R, false);
                                                qD = resD[0]; modeD = resD[1]; budgetD = resD[2];
                                            }
                                            Polling.Landing landD = Polling.land(pinfoD, qD, modeD, budgetD,
                                                    bufRowD, srvRowD, varRowD, K, Ks, pie.get(sn.stations.get(ist)), sn.jobclasses, R);
                                            for (int jD = 0; jD < landD.rows.size(); jD++) {
                                                outspace = Matrix.concatRows(outspace, landD.rows.get(jD), null);
                                                double cdD = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lldD;
                                                if (ni.hasInfinite()) {
                                                    lldD = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                } else {
                                                    lldD = lldscaling.get(ist, (int) Math.min(ni.get(row, 0), lldlimit) - 1);
                                                }
                                                Matrix orD = new Matrix(1, 1);
                                                orD.set(0, 0, cdD * lldD * rateD * landD.probs.get(jD).doubleValue());
                                                outrate = Matrix.concatRows(outrate, orD, null);
                                                Matrix opD = new Matrix(1, 1);
                                                opD.set(0, 0, 1.0);
                                                outprob = Matrix.concatRows(outprob, opD, null);
                                            }
                                        }
                                        break;
                                    }
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
                                        if (noPromote) { // immediate feedback: hold server, do not promote a waiting job
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                enWobufSiro.set(row, 0, en.get(row, 0));
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
                                                double cdscalingIstSiro = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lldSiro = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrateBottomWobufSiro = Matrix.scaleMult(rateEnWobufSiro, cdscalingIstSiro * lldSiro);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSiro, null);
                                            } else {
                                                double cdscalingIstSiro = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lldSiro = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
                                                Matrix outrateBottomWobufSiro = Matrix.scaleMult(rateEnWobufSiro, cdscalingIstSiro * lldSiro);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSiro, null);
                                            }

                                            Matrix outprobBottomWobufSiro = new Matrix(rateEnWobufSiro.getNumRows(), 1);
                                            outprobBottomWobufSiro.ones();
                                            outprob = Matrix.concatRows(outprob, outprobBottomWobufSiro, null);
                                        }

                                        // Handle states with buffer - pick jobs randomly from each class
                                        // (promotion suppressed under immediate feedback: bound 0)
                                        for (int r = 0; r < (noPromote ? 0 : R); r++) {
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
                                                        double cdscalingIstSiro = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                        double lldSiro = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                        Matrix outrateBottomWbufSiro = Matrix.scaleMult(rateEnWbufSiro, cdscalingIstSiro * lldSiro);
                                                        outrate = Matrix.concatRows(outrate, outrateBottomWbufSiro, null);
                                                    } else {
                                                        double cdscalingIstSiro = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                                        if (noPromote) { // immediate feedback: hold server, do not promote a waiting job
                                            for (int row = 0; row < en.getNumRows(); row++) {
                                                enWobufSeptLept.set(row, 0, en.get(row, 0));
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
                                                double cdscalingIstSeptLept = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                double lldSeptLept = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                Matrix outrateBottomWobufSeptLept = Matrix.scaleMult(rateEnWobufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                outrate = Matrix.concatRows(outrate, outrateBottomWobufSeptLept, null);
                                            } else {
                                                double cdscalingIstSeptLept = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                // Clamp as min(ni,lldlimit)-1 (0-based), matching MATLAB
                                                // lldscaling(ist,min(ni,lldlimit)) and the SIRO case; the prior
                                                // min(ni-1,lldlimit) could index lldlimit (out of bounds) once a
                                                // non-empty-buffer state reaches this freed-server path (immfeed).
                                                double lldSeptLept = lldscaling.get(ist, (int) Maths.min(ni.get(0), lldlimit) - 1);
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
                                                                    double cdscalingIstSeptLept = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                                                    double lldSeptLept = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                                                    Matrix outrateBottomWbufSeptLept = Matrix.scaleMult(rateEnWbufSeptLept, cdscalingIstSeptLept * lldSeptLept);
                                                                    outrate = Matrix.concatRows(outrate, outrateBottomWbufSeptLept, null);
                                                                } else {
                                                                    double cdscalingIstSeptLept = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
        // True BAS: when the front job is already blocked (completed, held at the server),
        // this DEP is the *instant transfer* of that job downstream — fire at rate 1e7
        // (effectively instant) and clear the blocked marker in the successor. The
        // complementary become-blocked transition (b:0->1 when the destination is full) is
        // added by the CTMC generator, since only it can see the destination's occupancy.
        // Gate on the dedicated sn.isbasblocking field (set for the blocking station
        // under BOTH declaration forms), not the shared marker column which a width-1
        // polling controller also sets. See BUG-83.
        if (sn.isbasblocking != null && ind < sn.isbasblocking.length() && sn.isbasblocking.get(ind) == 1
                && !inspace.isEmpty() && inspace.getNumCols() >= 1
                && inspace.get(0, inspace.getNumCols() - 1) == 1.0
                && !outspace.isEmpty()) {
            int obcol = outspace.getNumCols() - 1;
            for (int row = 0; row < outspace.getNumRows(); row++) {
                outspace.set(row, obcol, 0.0);
            }
            for (int row = 0; row < outrate.length(); row++) {
                outrate.set(row, 1.0e7);
            }
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    private static Ret.EventResult handlePhase(NetworkStruct sn, int ind, Matrix inspace, EventType event, int jobClass, boolean isSimulation,
                                             Matrix outspace, Matrix outrate, Matrix outprob, EventCache eventCache,
                                             int M, int R, Matrix S, Matrix phasessz, Matrix phaseshift, Map<Station, Map<JobClass, Matrix>> pie, Matrix ismkvmodclass,
                                             Matrix lldscaling, int lldlimit, Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                                             boolean hasOnlyExp, int ist, Matrix K, Matrix Ks, Map<Station, Map<JobClass, Matrix>> mu, Map<Station, Map<JobClass, Matrix>> phi,
                                             Map<Station, Map<JobClass, MatrixCell>> proc, Matrix capacity, Matrix classcap, double V, Matrix spaceBuf, Matrix spaceSrv, Matrix spaceVar, EventCacheKey key) {
        Matrix sir = null;
        List<Matrix> kir = null;
        Matrix ni = null;
        Matrix nir = null;
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
                                    if (ismkvmodclass.get(jobClass) != 0 && spaceVarK.getNumCols() > 0) {
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
                                            // D0 inter-phase (non-firing) rate of the arrival PH/MAP.
                                            // Must NOT be cast to int: fractional inter-phase rates
                                            // (e.g. Coxian/MMPP2 phase1->phase2 = 0.25) were truncated to 0,
                                            // freezing the source in phase 1 and yielding the phase-1 exit
                                            // rate as the (wrong) arrival rate. Keep the double value.
                                            rate.set(0, 0, proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest));
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
                                            double ni_value = ni.get(0);
                                            // proc*kir/ni*min(ni,S): the multiserver factor multiplies
                                            // the rate, matching afterEventStation.m case PHASE
                                            numerator.set(0, 0, proc_value_ps * kir_value_ps * Maths.min(ni_value, S.get(ist)));
                                            Matrix denom = new Matrix(1, 1);
                                            denom.set(0, 0, ni_value);
                                            rate = numerator.elementDivide(denom);
                                            break;
                                        case PSPRIO: {
                                            // Find minimum priority among present classes
                                            int minPrioPsPh = Integer.MAX_VALUE;
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (nir.get(0, r) > 0) {
                                                    int rPrio = (int) sn.classprio.get(r);
                                                    if (rPrio < minPrioPsPh) {
                                                        minPrioPsPh = rPrio;
                                                    }
                                                }
                                            }
                                            double ni_value_psprio = ni.get(0);
                                            // If ni <= S (all jobs get service) or this class has highest priority
                                            if (ni_value_psprio <= S.get(ist) || (int) sn.classprio.get(jobClass) == minPrioPsPh) {
                                                double proc_value_psprio = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                                double kir_value_psprio = kir.get(k).get(jobClass);
                                                Matrix numerator_prio = new Matrix(1, 1);
                                                // proc*kir/ni*min(ni,S), as in afterEventStation.m case PHASE
                                                numerator_prio.set(0, 0, proc_value_psprio * kir_value_psprio * Maths.min(ni_value_psprio, S.get(ist)));
                                                Matrix denom_prio = new Matrix(1, 1);
                                                denom_prio.set(0, 0, ni_value_psprio);
                                                rate = numerator_prio.elementDivide(denom_prio);
                                            } else {
                                                // Not in highest priority group - no service
                                                rate.set(0, 0, 0.0);
                                            }
                                            break;
                                        }
                                        case DPSPRIO: {
                                            int minPrioDpsPh = Integer.MAX_VALUE;
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (nir.get(0, r) > 0) {
                                                    int rPrio = (int) sn.classprio.get(r);
                                                    if (rPrio < minPrioDpsPh) minPrioDpsPh = rPrio;
                                                }
                                            }
                                            double ni_value_dpsprio = ni.get(0);
                                            if (ni_value_dpsprio <= S.get(ist) || (int) sn.classprio.get(jobClass) == minPrioDpsPh) {
                                                Matrix wDpsPh = sn.schedparam.getRow(ist);
                                                wDpsPh.scaleEq(1.0 / wDpsPh.elementSum());
                                                Matrix nirprioPh = nir.copy();
                                                if (ni_value_dpsprio > S.get(ist)) {
                                                    for (int r = 0; r < sn.nclasses; r++) {
                                                        if (sn.classprio.get(r) != sn.classprio.get(jobClass)) {
                                                            nirprioPh.set(r, 0);
                                                        }
                                                    }
                                                }
                                                double proc_val_dpsprio = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                                double kir_val_dpsprio = kir.get(k).get(jobClass);
                                                double wDenom = wDpsPh.mult(nirprioPh.transpose()).sumRows().get(0, 0);
                                                rate.set(0, 0, proc_val_dpsprio * kir_val_dpsprio * wDpsPh.get(jobClass) / wDenom);
                                            } else {
                                                rate.set(0, 0, 0.0);
                                            }
                                            break;
                                        }
                                        case GPSPRIO: {
                                            int minPrioGpsPh = Integer.MAX_VALUE;
                                            for (int r = 0; r < sn.nclasses; r++) {
                                                if (nir.get(0, r) > 0) {
                                                    int rPrio = (int) sn.classprio.get(r);
                                                    if (rPrio < minPrioGpsPh) minPrioGpsPh = rPrio;
                                                }
                                            }
                                            double ni_value_gpsprio = ni.get(0);
                                            if (ni_value_gpsprio <= S.get(ist) || (int) sn.classprio.get(jobClass) == minPrioGpsPh) {
                                                Matrix wGpsPh = sn.schedparam.getRow(ist);
                                                wGpsPh.scaleEq(1.0 / wGpsPh.elementSum());
                                                Matrix nirprioPh = nir.copy();
                                                if (ni_value_gpsprio > S.get(ist)) {
                                                    for (int r = 0; r < sn.nclasses; r++) {
                                                        if (sn.classprio.get(r) != sn.classprio.get(jobClass)) {
                                                            nirprioPh.set(r, 0);
                                                        }
                                                    }
                                                }
                                                Matrix cirPh = new Matrix(nirprioPh.getNumRows(), nirprioPh.getNumCols());
                                                for (int row = 0; row < cirPh.getNumRows(); row++) {
                                                    for (int col = 0; col < cirPh.getNumCols(); col++) {
                                                        cirPh.set(row, col, Math.min(nirprioPh.get(row, col), 1.0));
                                                    }
                                                }
                                                Matrix cirPh1D = new Matrix(0, 0);
                                                for (int col = 0; col < cirPh.getNumCols(); col++) {
                                                    cirPh1D = Matrix.concatRows(cirPh1D, cirPh.getColumn(col), null);
                                                }
                                                double wcirPh = wGpsPh.mult(cirPh1D).get(0);
                                                double proc_val_gpsprio = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                                double kir_val_gpsprio = kir.get(k).get(jobClass);
                                                rate.set(0, 0, proc_val_gpsprio * kir_val_gpsprio / nirprioPh.get(jobClass) * wGpsPh.get(jobClass) / wcirPh);
                                            } else {
                                                rate.set(0, 0, 0.0);
                                            }
                                            break;
                                        }
                                        case DPS:
                                            if (S.get(ist) > 1) {
                                                InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                                InputOutput.line_error(InputOutput.mfilename(new Object() {
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
                                        case FCFSPRIO:
                                        case LCFSPRIO:
                                        case LCFS:
                                        case LCFSPR:
                                        case LCFSPRPRIO:
                                        case FCFSPRPRIO:
                                        case LCFSPI:
                                        case FCFSPIPRIO:
                                        case LCFSPIPRIO:
                                        case SIRO:
                                        case SEPT:
                                        case LEPT:
                                        case POLLING:
                                            double proc_value = proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(jobClass)).get(0).get(k, kdest);
                                            double kir_value = kir.get(k).get(jobClass);
                                            rate.set(0, 0, proc_value * kir_value);
                                            break;
                                    }

                                    // if class cannot be served locally, rate = NaN since mu{i, class} = NaN
                                    if (ni.hasInfinite()) {
                                        // hit limited load-dependence
                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
                                        double lld = lldscaling.get(ist, lldscaling.getNumCols() - 1);
                                        Matrix outrate_bottom = Matrix.scaleMult(rate, cdscalingIst * lld);
                                        outrate = Matrix.concatRows(outrate, outrate_bottom, null);
                                    } else {
                                        double cdscalingIst = cdScalar(cdscaling, sn.stations.get(ist), nir, jobClass);
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
                        // NOTE: no separate MAP/MMPP2 phase block here. The generic
                        // (k,kdest) section above already performs the MAP local-variable
                        // update (spaceVarK) and emits each phase transition exactly once,
                        // matching MATLAB afterEventStation.m case EventType.PHASE. A
                        // duplicated MAP-specific block nested inside the k loop previously
                        // re-emitted every transition K times, inflating all MAP/MMPP2
                        // phase-switch rates by a factor K+1 in the CTMC generator.
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
        return new Ret.EventResult(outspace, outrate, outprob);
    }

    // -----------------------------------------------------------------------
    // Pass-and-swap (PAS) / order-independent station handler
    // -----------------------------------------------------------------------

    /**
     * Result of a pass-and-swap transition: the new ordered list and the departing class.
     * Public because the NRM engine (Solver_ssa_nrm) drives the same rewrite off its own
     * ordered list rather than off a state row, and the mechanism must not be duplicated:
     * a second copy could silently drift from the generator this solver is validated against.
     */
    public static class PasSwap {
        public final int[] cnew;
        public final int depClass; // 1-based class of the departing job
        PasSwap(int[] cnew, int depClass) { this.cnew = cnew; this.depClass = depClass; }
    }

    /**
     * Applies the pass-and-swap mechanism (Dorsman and Gardner 2024, Sect. 2.3)
     * triggered by the service completion of position p.
     *
     * @param c ordered list of 1-based class indices (c[0] oldest)
     * @param p 0-based position of the completing job
     * @param G (nclasses x nclasses) 0-based swapping-graph adjacency
     */
    public static PasSwap passAndSwap(int[] c, int p, Matrix G) {
        int n = c.length;
        List<Integer> chain = new ArrayList<Integer>();
        chain.add(p);
        int moving = c[p];   // 1-based
        int cur = p;
        while (true) {
            int q = -1;
            for (int j = cur + 1; j < n; j++) {
                if (G.get(moving - 1, c[j] - 1) != 0) { q = j; break; }
            }
            if (q < 0) break;
            chain.add(q);
            moving = c[q];
            cur = q;
        }
        int depClass = c[chain.get(chain.size() - 1)];
        int[] tmp = c.clone();
        for (int i = 0; i < chain.size() - 1; i++) {
            tmp[chain.get(i + 1)] = c[chain.get(i)];
        }
        int hole = chain.get(0);
        int[] cnew = new int[n - 1];
        int w = 0;
        for (int i = 0; i < n; i++) {
            if (i == hole) continue;
            cnew[w++] = tmp[i];
        }
        return new PasSwap(cnew, depClass);
    }

    /**
     * Event handler for pass-and-swap (PAS) / order-independent stations. The
     * local state is the ordered list of 1-based class indices (oldest first),
     * left-aligned in the first W = cap columns and right zero-padded; the
     * trailing V columns are routing variables (carried through). Service is
     * governed by the total rate function mu(c) (QueueNodeParam.svcRateFun,
     * called with a row Matrix of 0-based class indices) and the swapping graph
     * G (QueueNodeParam.swapGraph, 0-based).
     */
    public static Ret.EventResult afterEventStationPas(NetworkStruct sn, int ind, int ist, Matrix inspace,
                                                       EventType event, int jobClass, int R, int V, boolean isSimulation) {
        jline.lang.nodeparam.QueueNodeParam qp = (jline.lang.nodeparam.QueueNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
        SerializableFunction<Matrix, Double> muFun = (qp != null) ? qp.svcRateFun : null;
        Matrix G = (qp != null) ? qp.swapGraph : null;
        if (muFun == null) {
            throw new RuntimeException("PAS station has no service rate function mu(c).");
        }
        int nCols = inspace.getNumCols();
        int W = nCols - V;
        int cap = (int) sn.cap.get(ist);
        int nRows = inspace.getNumRows();
        int jobClass1 = jobClass + 1; // 1-based departing class to match list encoding

        List<double[]> outRows = new ArrayList<double[]>();
        List<Double> outRate = new ArrayList<Double>();
        List<Double> outProb = new ArrayList<Double>();

        for (int row = 0; row < nRows; row++) {
            // extract the ordered list (1-based, contiguous) and the var columns
            List<Integer> cl = new ArrayList<Integer>();
            for (int col = 0; col < W; col++) {
                int v = (int) inspace.get(row, col);
                if (v > 0) cl.add(v);
            }
            int n = cl.size();
            int[] c = new int[n];
            for (int i = 0; i < n; i++) c[i] = cl.get(i);

            if (event == EventType.ARV) {
                if (n >= cap) continue; // buffer full: arrival lost
                double[] outRow = new double[nCols];
                for (int i = 0; i < n; i++) outRow[i] = c[i];
                outRow[n] = jobClass1;
                for (int col = W; col < nCols; col++) outRow[col] = inspace.get(row, col);
                outRows.add(outRow);
                outRate.add(-1.0);
                outProb.add(1.0);
            } else if (event == EventType.DEP) {
                if (n == 0) continue;
                double muPrev = 0.0;
                for (int p = 0; p < n; p++) {
                    Matrix prefix = new Matrix(1, p + 1);
                    for (int i = 0; i <= p; i++) prefix.set(0, i, c[i] - 1); // 0-based for mu
                    double muCur = muFun.apply(prefix);
                    double ratep = muCur - muPrev;
                    muPrev = muCur;
                    if (ratep <= 0) continue;
                    PasSwap res = passAndSwap(c, p, G);
                    if (res.depClass != jobClass1) continue;
                    double[] outRow = new double[nCols];
                    for (int i = 0; i < res.cnew.length; i++) outRow[i] = res.cnew[i];
                    for (int col = W; col < nCols; col++) outRow[col] = inspace.get(row, col);
                    outRows.add(outRow);
                    outRate.add(ratep);
                    outProb.add(1.0);
                }
            }
            // PHASE: PAS service is exponential, no phase transitions
        }

        // In simulation mode the SSA engine expects a single sampled active
        // outcome per (station, event, class). The CTMC generator instead needs
        // every position-completion as a separate transition. For a DEP that can
        // complete at several positions of the requested class, collapse the rows
        // into one: keep the summed rate (the class-r departure rate, so the sum
        // over classes equals the total station rate mu(c)) and pick the next
        // state proportionally to the per-position marginal rates.
        if (isSimulation && event == EventType.DEP && outRows.size() > 1) {
            double totalRate = 0.0;
            for (Double rr : outRate) totalRate += rr.doubleValue();
            double u = Maths.rand() * totalRate;
            double acc = 0.0;
            int pick = outRows.size() - 1;
            for (int i = 0; i < outRate.size(); i++) {
                acc += outRate.get(i).doubleValue();
                if (u <= acc) { pick = i; break; }
            }
            double[] chosen = outRows.get(pick);
            outRows = new ArrayList<double[]>();
            outRows.add(chosen);
            outRate = new ArrayList<Double>();
            outRate.add(Double.valueOf(totalRate));
            outProb = new ArrayList<Double>();
            outProb.add(Double.valueOf(1.0));
        }

        int m = outRows.size();
        Matrix outspace = new Matrix(m, nCols);
        Matrix outrate = new Matrix(m, 1);
        Matrix outprob = new Matrix(m, 1);
        for (int i = 0; i < m; i++) {
            double[] r = outRows.get(i);
            for (int col = 0; col < nCols; col++) outspace.set(i, col, r[col]);
            outrate.set(i, 0, outRate.get(i));
            outprob.set(i, 0, outProb.get(i));
        }
        return new Ret.EventResult(outspace, outrate, outprob);
    }
}
