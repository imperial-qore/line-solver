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
        c = find(layerSn.chains(:, classidx), 1);
        refclass_c = layerSn.refclass(c);
        refstat_k = layerSn.refstat(classidx);
        TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = self.results{end,layerIdx}.QN(nodeidx,classidx) / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
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

    % then resolve the entry servt summing up these contributions
    entry_servt = (eye(lqn.nidx+lqn.ncalls)-self.servtmatrix)\[self.servt;self.callservt];
    entry_servt(1:lqn.eshift) = 0;

    % Propagate forwarding calls: add target entry's service time to source entry
    for cidx = 1:lqn.ncalls
        if lqn.calltype(cidx) == CallType.FWD
            source_eidx = lqn.callpair(cidx, 1);
            target_eidx = lqn.callpair(cidx, 2);
            fwd_prob = lqn.callproc{cidx}.getMean();
            entry_servt(source_eidx) = entry_servt(source_eidx) + fwd_prob * entry_servt(target_eidx);
        end
    end

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
        c = find(layerSn.chains(:, classidx), 1);
        refclass_c = layerSn.refclass(c);
        refstat_k = layerSn.refstat(classidx);
        TN_ref = self.results{end,layerIdx}.TN(refstat_k, refclass_c);
        if TN_ref > GlobalConstants.FineTol
            self.residt(aidx) = self.results{end,layerIdx}.QN(nodeidx,classidx) / TN_ref;
        else
            self.residt(aidx) = self.results{end,layerIdx}.WN(nodeidx,classidx);
        end

        submodelidx = self.idxhash(idx);
        if submodelidx>length(repo)
            % Try SolverFluid first for actual response time CDFs with
            % higher-moment information; fall back to layer solver's
            % exponential approximation if Fluid fails on the layer model
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
                    zerodist = APH.fitMeanAndSCV(GlobalConstants.FineTol,0.99);
                    zeroparam{1} = zerodist.params{1}.paramValue;
                    zeroparam{2} = zerodist.params{2}.paramValue;
                    aph_pattern = 3; % branch structure
                    [aphparam{1},aphparam{2}] = aph_simplify(aphparam{1},aphparam{2},zeroparam{1},zeroparam{2},fractionalPartRepetitions,1-fractionalPartRepetitions, aph_pattern);
                    ParamCell = [ParamCell,aphparam];
                else
                    zerodist = APH.fitMeanAndSCV(GlobalConstants.FineTol,0.99);
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

    % Propagate forwarding calls: add target entry's service time to source entry
    for cidx = 1:lqn.ncalls
        if lqn.calltype(cidx) == CallType.FWD
            source_eidx = lqn.callpair(cidx, 1);
            target_eidx = lqn.callpair(cidx, 2);
            fwd_prob = lqn.callproc{cidx}.getMean();
            self.servt(source_eidx) = self.servt(source_eidx) + fwd_prob * self.servt(target_eidx);
            self.servtproc{source_eidx} = Exp.fitMean(self.servt(source_eidx));
        end
    end

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
end
self.ensemble = ensemble;
end
