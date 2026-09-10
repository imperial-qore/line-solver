function updateMetricsMomentBased(self,it)
ensemble = self.ensemble;
lqn = self.lqn;
if ~self.hasconverged
    % first obtain servt of activities at hostlayers
    self.servt = zeros(lqn.nidx,1);
    self.residt = zeros(lqn.nidx,1);
    for r=1:size(self.servt_classes_updmap,1)
        idx = self.servt_classes_updmap(r,1);
        aidx = self.servt_classes_updmap(r,2);
        nodeidx = self.servt_classes_updmap(r,3);
        classidx = self.servt_classes_updmap(r,4);
        self.servt(aidx) = self.results{end,self.idxhash(idx)}.RN(nodeidx,classidx);
        self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);
        self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));

        % Compute residt from QN/TN_ref (matching updateMetricsDefault)
        layerIdx = self.idxhash(idx);
        layerSn = ensemble{layerIdx}.getStruct();
        % refclass is 0 on a chain with no reference class (open chain)
        [hasRef, refstat_k, refclass_c] = ln_layer_refcell(layerSn, classidx);
        if hasRef
            TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        else
            TN_ref = 0;
        end
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = self.results{end,layerIdx}.QN(nodeidx,classidx) / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
        end

        % see _kb/06-solver-catalog.md (LN section) for rationale
        zt_act = lqn_act_thinktime(lqn, aidx);
        if zt_act > 0
            self.servt(aidx) = self.servt(aidx) + zt_act;
            self.residt(aidx) = self.residt(aidx) + zt_act;
            self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));
        end

        % async-only targets carry no visit-ratio scaling (matching updateMetricsDefault)
        if aidx > lqn.ashift && aidx <= lqn.ashift + lqn.nacts
            for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
                if full(lqn.graph(eidx, aidx)) > 0
                    if full(any(lqn.isasynccaller(:, eidx))) && ~full(any(lqn.issynccaller(:, eidx)))
                        self.residt(aidx) = self.servt(aidx);
                    end
                    break;
                end
            end
        end
    end

    % estimate call response times at hostlayers
    self.callservt = zeros(lqn.ncalls,1);
    self.callresidt = zeros(lqn.ncalls,1);
    for r=1:size(self.call_classes_updmap,1)
        idx = self.call_classes_updmap(r,1);
        cidx = self.call_classes_updmap(r,2);
        nodeidx = self.call_classes_updmap(r,3);
        classidx = self.call_classes_updmap(r,4);
        if self.call_classes_updmap(r,3) > 1
            if nodeidx == 1
                self.callservt(cidx) = 0;
                self.callresidt(cidx) = 0;
            else
                % Include call multiplicity in callservt (matching updateMetricsDefault)
                self.callservt(cidx) = self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx) * self.lqn.callproc{cidx}.getMean;
                % callresidt uses WN which already includes visit multiplicity
                self.callresidt(cidx) = self.results{end, self.idxhash(idx)}.WN(nodeidx,classidx);
            end
        end
    end

    % then resolve the entry servt summing up these contributions; the terms are
    % residence times (Vtask=1), which the task/entry tput ratio below rescales to Ventry=1
    entry_servt = (eye(lqn.nidx+lqn.ncalls)-self.servtmatrix)\[self.residt;self.callresidt];
    entry_servt(1:lqn.eshift) = 0;

    % NO forwarding propagation here. lqn_fwd_rendezvous has already reconnected
    % every forwarding chain reachable from a synchronous call to the client that
    % issued the rendezvous (Franks 1999, Sec. 3.3.1), so the forwarded service is
    % in the caller's chain before this runs; adding it again inflated the caller
    % by exactly the forwarded entry's mean. An asynchronous call into a chain is
    % left untouched there by design -- a send-no-reply does not block -- so it
    % must not accumulate the forwarded service either. See BUGS.md BUG-91.

    self.servt(lqn.eshift+1:lqn.eshift+lqn.nentries) = entry_servt(lqn.eshift+1:lqn.eshift+lqn.nentries);
    entry_servt((lqn.ashift+1):end) = 0;

    % Compute entry-level residt using servtmatrix and activity residt
    % callresidt uses WN which already includes visit multiplicity
    entry_residt = self.servtmatrix*[self.residt;self.callresidt(:)];
    entry_residt(1:lqn.eshift) = 0;
    % Scale entry residt by task/entry throughput ratio (matching updateMetricsDefault)
    for eidx=(lqn.eshift+1):(lqn.eshift+lqn.nentries)
        tidx = lqn.parent(eidx);
        hidx = lqn.parent(tidx);
        if ~self.ignore(tidx) && ~self.ignore(hidx)
            hasSyncCallers = full(any(lqn.issynccaller(:, eidx)));
            if hasSyncCallers
                tidxclass = ensemble{self.idxhash(hidx)}.attribute.tasks(find(ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == tidx),1);
                eidxclass = ensemble{self.idxhash(hidx)}.attribute.entries(find(ensemble{self.idxhash(hidx)}.attribute.entries(:,2) == eidx),1);
                task_tput  = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,tidxclass));
                entry_tput = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,eidxclass));
                self.servt(eidx) = entry_servt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
                self.residt(eidx) = entry_residt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
            else
                self.residt(eidx) = entry_residt(eidx);
            end
        end
    end

    for r=1:size(self.call_classes_updmap,1)
        cidx = self.call_classes_updmap(r,2);
        eidx = lqn.callpair(cidx,2);
        if self.call_classes_updmap(r,3) > 1
            self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
        end
    end

    % determine call response times processes
    for r=1:size(self.call_classes_updmap,1)
        cidx = self.call_classes_updmap(r,2);
        eidx = lqn.callpair(cidx,2);
        if self.call_classes_updmap(r,3) > 1
            if it==1
                % note that respt is per visit, so number of calls is 1
                self.callservt(cidx) = self.servt(eidx);
                self.callservtproc{cidx} = self.servtproc{eidx};
            else
                % note that respt is per visit, so number of calls is 1
                self.callservtproc{cidx} = Exp.fitMean(self.callservt(cidx));
            end
        end
    end
else
    self.servtcdf = cell(lqn.nidx,1);
    repo = [];

    % first obtain servt of activities at hostlayers
    self.servt = zeros(lqn.nidx,1);
    self.residt = zeros(lqn.nidx,1);
    for r=1:size(self.servt_classes_updmap,1)
        idx = self.servt_classes_updmap(r,1);
        aidx = self.servt_classes_updmap(r,2);
        nodeidx = self.servt_classes_updmap(r,3);
        classidx = self.servt_classes_updmap(r,4);
        self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);

        % Compute residt from QN/TN_ref (matching updateMetricsDefault)
        layerIdx = self.idxhash(idx);
        layerSn = ensemble{layerIdx}.getStruct();
        % refclass is 0 on a chain with no reference class (open chain)
        [hasRef, refstat_k, refclass_c] = ln_layer_refcell(layerSn, classidx);
        if hasRef
            TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        else
            TN_ref = 0;
        end
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = self.results{end,layerIdx}.QN(nodeidx,classidx) / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
        end

        submodelidx = self.idxhash(idx);
        if submodelidx>length(repo)
            % see _kb/06-solver-catalog.md (LN section) for rationale
            try
                repo{submodelidx} = SolverFluid(ensemble{submodelidx}).getCdfRespT;
            catch
                repo{submodelidx} = self.solvers{submodelidx}.getCdfRespT;
            end
        end
        self.servtcdf{aidx} =  repo{submodelidx}{nodeidx,classidx};
    end

    self.callservtcdf = cell(lqn.ncalls,1);

    % estimate call response times at hostlayers
    self.callservt = zeros(lqn.ncalls,1);
    self.callresidt = zeros(lqn.ncalls,1);
    for r=1:size(self.call_classes_updmap,1)
        idx = self.call_classes_updmap(r,1);
        cidx = self.call_classes_updmap(r,2);
        nodeidx = self.call_classes_updmap(r,3);
        classidx = self.call_classes_updmap(r,4);
        if self.call_classes_updmap(r,3) > 1
            submodelidx = self.idxhash(idx);
            if submodelidx>length(repo)
                try
                    repo{submodelidx} = SolverFluid(ensemble{submodelidx}).getCdfRespT;
                catch
                    repo{submodelidx} = self.solvers{submodelidx}.getCdfRespT;
                end
            end
            try
                self.callservtcdf{cidx} =  repo{submodelidx}{nodeidx,classidx};
            catch
                self.callservtcdf{cidx} =  repo{submodelidx};
            end
            % Also set callresidt from WN (includes visit multiplicity)
            if nodeidx > 1
                self.callresidt(cidx) = self.results{end, self.idxhash(idx)}.WN(nodeidx,classidx);
            end
        end
    end
    cdf = [self.servtcdf;self.callservtcdf];

    % then resolve the entry servt summing up these contributions
    matrix = inv((eye(lqn.nidx+lqn.ncalls)-self.servtmatrix));
    for i = 1:1:lqn.nentries
        eidx = lqn.eshift+i;
        convolidx = find(matrix(eidx,:)>0);
        convolidx(find(convolidx<=lqn.eshift+lqn.nentries))=[];
        num = 0;
        ParamCell = {};
        while num<length(convolidx)
            fitidx = convolidx(num+1);
            [m1,m2,m3,~,~] = EmpiricalCDF(cdf{fitidx}).getMoments;

            % see _kb/06-solver-catalog.md (LN section) for rationale
            if fitidx <= lqn.nidx && lqn_act_thinktime(lqn, fitidx) > 0
                ztd = lqn.actthink{fitidx};
                t1 = ztd.getMean();
                sig2 = ztd.getSCV()*t1^2;
                t2 = sig2 + t1^2;
                t3 = ztd.getSkewness()*sig2^1.5 + 3*t1*t2 - 2*t1^3;
                m3 = m3 + 3*m2*t1 + 3*m1*t2 + t3;
                m2 = m2 + 2*m1*t1 + t2;
                m1 = m1 + t1;
            end
            % Use CoarseTol to skip near-zero mean CDFs (e.g., Immediate activities)
            % whose APH fitting produces extreme rate parameters that cause
            % matrix exponential evaluation to hang
            if m1>GlobalConstants.CoarseTol
                fitdist = APH.fitRawMoments(m1,m2,m3);
                aphparam{1} = fitdist.params{1}.paramValue; % alpha
                aphparam{2} = fitdist.params{2}.paramValue; % T

                % For call indices, multiply repetitions by mean number of calls
                % servtmatrix has 1.0 for calls, but we need callproc.getMean() repetitions
                reps = matrix(eidx,fitidx);
                if fitidx > lqn.nidx
                    cidx_local = fitidx - lqn.nidx;
                    reps = reps * lqn.callproc{cidx_local}.getMean();
                end
                integerRepetitions = floor(reps);
                fractionalPartRepetitions = reps - integerRepetitions;

                if fractionalPartRepetitions == 0
                    ParamCell = [ParamCell,repmat(aphparam,1,integerRepetitions)];
                elseif integerRepetitions>0 && fractionalPartRepetitions>0
                    ParamCell = [ParamCell,repmat(aphparam,1,integerRepetitions)];
                    % see _kb/06-solver-catalog.md (LN section) for rationale
                    zerodist = APH(1, -GlobalConstants.Immediate);
                    zeroparam{1} = zerodist.params{1}.paramValue;
                    zeroparam{2} = zerodist.params{2}.paramValue;
                    aph_pattern = 3; % branch structure
                    [aphparam{1},aphparam{2}] = aph_simplify(aphparam{1},aphparam{2},zeroparam{1},zeroparam{2},fractionalPartRepetitions,1-fractionalPartRepetitions, aph_pattern);
                    ParamCell = [ParamCell,aphparam];
                else
                    % see _kb/06-solver-catalog.md (LN section) for rationale
                    zerodist = APH(1, -GlobalConstants.Immediate);
                    zeroparam{1} = zerodist.params{1}.paramValue;
                    zeroparam{2} = zerodist.params{2}.paramValue;
                    aph_pattern = 3; % branch structure
                    [aphparam{1},aphparam{2}] = aph_simplify(aphparam{1},aphparam{2},zeroparam{1},zeroparam{2},fractionalPartRepetitions,1-fractionalPartRepetitions, aph_pattern);
                    ParamCell = [ParamCell,aphparam];
                end
                if fitidx <= lqn.nidx
                    self.servtproc{fitidx} = Exp.fitMean(m1);
                    self.servt(fitidx) = m1;
                else
                    self.callservtproc{fitidx-lqn.nidx} = Exp.fitMean(m1);
                    self.callservt(fitidx-lqn.nidx) = m1;
                end
            end
            num = num+1;
        end
        if isempty(ParamCell)
            self.servt(eidx) = 0;
        else
            [alpha,T] = aph_convseq(ParamCell); % convolution of sequential activities
            entry_dist = APH(alpha,T);
            self.entryproc{eidx-(lqn.nhosts+lqn.ntasks)} = entry_dist;
            self.servt(eidx) = entry_dist.getMean;
            self.servtproc{eidx} = Exp.fitMean(self.servt(eidx));
            self.entrycdfrespt{eidx-(lqn.nhosts+lqn.ntasks)} = entry_dist.evalCDF;
        end
    end

    % NO forwarding propagation here, for the reason given at the entry_servt
    % assembly above: lqn_fwd_rendezvous has already charged the forwarded
    % service to the caller. See BUGS.md BUG-91.

    % Compute entry-level residt using servtmatrix and activity residt
    % callresidt uses WN which already includes visit multiplicity
    entry_residt = self.servtmatrix*[self.residt;self.callresidt(:)];
    entry_residt(1:lqn.eshift) = 0;
    for eidx=(lqn.eshift+1):(lqn.eshift+lqn.nentries)
        tidx = lqn.parent(eidx);
        hidx = lqn.parent(tidx);
        if ~self.ignore(tidx) && ~self.ignore(hidx)
            hasSyncCallers = full(any(lqn.issynccaller(:, eidx)));
            if hasSyncCallers
                tidxclass = ensemble{self.idxhash(hidx)}.attribute.tasks(find(ensemble{self.idxhash(hidx)}.attribute.tasks(:,2) == tidx),1);
                eidxclass = ensemble{self.idxhash(hidx)}.attribute.entries(find(ensemble{self.idxhash(hidx)}.attribute.entries(:,2) == eidx),1);
                task_tput  = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,tidxclass));
                entry_tput = sum(self.results{end,self.idxhash(hidx)}.TN(ensemble{self.idxhash(hidx)}.attribute.clientIdx,eidxclass));
                self.residt(eidx) = entry_residt(eidx) * task_tput / max(GlobalConstants.Zero, entry_tput);
            else
                self.residt(eidx) = entry_residt(eidx);
            end
        end
    end

    % determine call response times processes
    for r=1:size(self.call_classes_updmap,1)
        cidx = self.call_classes_updmap(r,2);
        eidx = lqn.callpair(cidx,2);
        if self.call_classes_updmap(r,3) > 1
            if it==1
                self.callservt(cidx) = self.servt(eidx);
                self.callservtproc{cidx} = Exp.fitMean(self.servt(eidx));
            end
        end
    end

    % This pass IS the moment3 answer, and it is TERMINAL. Its entry laws are
    % convolutions of the activities' own response distributions; the branch
    % above instead reads QN/TN_ref, a residence per REFERENCE cycle, which the
    % entry assembly then treats as a per-entry-visit time. The two disagree by
    % the entry's visit ratio whenever it is not 1, so letting the iteration
    % fall back to that branch after this one has run DISCARDS the moment-based
    % laws and reports the other quantity. See BUGS.md BUG-97.
    self.momentPassDone = true;
end
self.ensemble = ensemble;
end
