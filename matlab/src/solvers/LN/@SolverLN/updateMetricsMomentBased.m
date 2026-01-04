function updateMetricsMomentBased(self,it)
ensemble = self.ensemble;
lqn = self.lqn;
if ~self.hasconverged
    % first obtain servt of activities at hostlayers
    self.servt = zeros(lqn.nidx,1);
    for r=1:size(self.servt_classes_updmap,1)
        idx = self.servt_classes_updmap(r,1);
        aidx = self.servt_classes_updmap(r,2);
        nodeidx = self.servt_classes_updmap(r,3);
        classidx = self.servt_classes_updmap(r,4);
        % this requires debugging, it should be WN but it has been
        % set temporarily to RN as bugs are left
        self.servt(aidx) = self.results{end,self.idxhash(idx)}.RN(nodeidx,classidx);
        %self.servt(aidx) = self.results{end,self.idxhash(idx)}.WN(nodeidx,classidx);
        self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);
        self.servtproc{aidx} = Exp.fitMean(self.servt(aidx));
    end

    % estimate call response times at hostlayers
    self.callservt = zeros(lqn.ncalls,1);
    for r=1:size(self.call_classes_updmap,1)
        idx = self.call_classes_updmap(r,1);
        cidx = self.call_classes_updmap(r,2);
        nodeidx = self.call_classes_updmap(r,3);
        classidx = self.call_classes_updmap(r,4);
        if self.call_classes_updmap(r,3) > 1

            %             self.callservt(cidx) = self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx);
            if nodeidx == 1
                self.callservt(cidx) = 0;
            else
                self.callservt(cidx) = self.results{end, self.idxhash(idx)}.RN(nodeidx,classidx);
                %self.callservt(cidx) = self.results{end, self.idxhash(idx)}.WN(nodeidx,classidx);
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
    for r=1:size(self.servt_classes_updmap,1)
        idx = self.servt_classes_updmap(r,1);
        aidx = self.servt_classes_updmap(r,2);
        nodeidx = self.servt_classes_updmap(r,3);
        classidx = self.servt_classes_updmap(r,4);
        self.tput(aidx) = self.results{end,self.idxhash(idx)}.TN(nodeidx,classidx);

        submodelidx = self.idxhash(idx);
        if submodelidx>length(repo)
            try
                repo{submodelidx} = self.solvers{submodelidx}.getCdfRespT;
            catch me
                switch me.identifier
                    case 'MATLAB:class:undefinedMethod'
                        line_warning(mfilename,'The solver for layer %d does not allow response time distribution calculation, switching to fluid solver.\n');
                        repo{submodelidx} = SolverFluid(ensemble{submodelidx}).getCdfRespT;
                end
            end
        end
        self.servtcdf{aidx} =  repo{submodelidx}{nodeidx,classidx};
    end

    self.callservtcdf = cell(lqn.ncalls,1);

    % estimate call response times at hostlayers
    self.callservt = zeros(lqn.ncalls,1);
    for r=1:size(self.call_classes_updmap,1)
        idx = self.call_classes_updmap(r,1);
        cidx = self.call_classes_updmap(r,2);
        nodeidx = self.call_classes_updmap(r,3);
        classidx = self.call_classes_updmap(r,4);
        if self.call_classes_updmap(r,3) > 1
            submodelidx = self.idxhash(idx);
            if submodelidx>length(repo)
                %             repo{submodelidx} = self.solvers{submodelidx}.getCdfRespT;
                try
                    repo{submodelidx} = self.solvers{submodelidx}.getCdfRespT;
                catch
                    repo{submodelidx} = [0,0;0.5,0;1,0]; % ??
                end
            end
            %         self.callservtcdf{cidx} =  repo{submodelidx}{nodeidx,classidx};
            try
                self.callservtcdf{cidx} =  repo{submodelidx}{nodeidx,classidx};
            catch
                self.callservtcdf{cidx} =  repo{submodelidx};
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
            if m1>GlobalConstants.FineTol
                fitdist = APH.fitRawMoments(m1,m2,m3);
                aphparam{1} = fitdist.params{1}.paramValue; % alpha
                aphparam{2} = fitdist.params{2}.paramValue; % T                
                integerRepetitions = floor(matrix(eidx,fitidx));
                fractionalPartRepetitions = matrix(eidx,fitidx)-integerRepetitions;
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