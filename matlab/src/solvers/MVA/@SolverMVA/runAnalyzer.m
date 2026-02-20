function runtime = runAnalyzer(self, options)
% RUNTIME = RUNANALYZER(OPTIONS)
% Run the solver

if nargin<2
    options = self.getOptions;
end

self.runAnalyzerChecks(options);

Solver.resetRandomGeneratorSeed(options.seed);

% Show library attribution if verbose and not yet shown
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    sn = self.getStruct();
    libs = SolverMVA.getLibrariesUsed(sn, options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

iter = 0;
%options.lang='java';

switch options.lang
    case 'java'
        sn = getStruct(self); % doesn't need initial state
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverMVA(jmodel, options);
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
        runtime = jsolver.result.runtime;
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        TN = reshape(TN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        lG = NaN;
        lastiter = NaN;
        for ind = 1:sn.nnodes
            if sn.nodetype(ind) == NodeType.Cache
                self.model.nodes{ind}.setResultHitProb(JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getHitRatio()));
                self.model.nodes{ind}.setResultMissProb(JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getMissRatio()));
            end
        end
        %self.model.refreshChains();
        self.model.refreshStruct(true);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        sn = getStruct(self); % doesn't need initial state
        forkLoop = true;
        forkIter = 0;
        % create artificial classes arrival rates
        forkLambda = GlobalConstants.FineTol * ones(1, 2*sn.nclasses*sum(sn.nodetype==NodeType.Fork));
        QN = GlobalConstants.Immediate * ones(1, sn.nclasses);
        QN_1 = 0*QN;
        UN = 0*QN;
        while (forkLoop && forkIter < options.iter_max)
            if self.model.hasFork
                forkIter = forkIter + 1;
                if forkIter == 1
                    switch options.config.fork_join
                        case {'heidelberger-trivedi', 'ht'}
                            [nonfjmodel, fjclassmap, fjforkmap, fj_auxiliary_delays] = ModelAdapter.ht(self.model);
                        case {'mmt', 'default', 'fjt'}
                            [nonfjmodel, fjclassmap, fjforkmap, fanout] = ModelAdapter.mmt(self.model, forkLambda);
                            [outer_forks, parent_forks] = ModelAdapter.sortForks(sn, fjforkmap, fjclassmap, nonfjmodel);
                    end
                elseif ~strcmp(options.config.fork_join, 'heidelberger-trivedi') & ~strcmp(options.config.fork_join, 'ht')
                    %line_printf('Fork-join iteration %d\n',forkIter);
                    nonfjSource = nonfjmodel.getSource;
                    for r=1:length(fjclassmap) % r is the auxiliary class
                        s = fjclassmap(r);
                        if s>0
                            if fanout(r)>0
                                if ~nonfjSource.arrivalProcess{r}.isDisabled
                                    nonfjSource.arrivalProcess{r}.setRate((fanout(r)-1)*forkLambda(r));
                                end
                            end
                        end
                    end
                    nonfjmodel.refreshRates();
                end
                sn = nonfjmodel.getStruct(false); % this ensures that we solve nonfjmodel instead of the original model
                % Convergence check: handle zero/NaN queue lengths from size mismatch and Immediate activities
                qn_ratio = QN_1 ./ QN;
                qn_ratio(isnan(qn_ratio)) = 1; % 0/0 = NaN → treat as converged (ratio=1)
                qn_ratio(isinf(qn_ratio)) = 0; % x/0 = Inf → treat as not converged (ratio=0, error=1)
                if isequal(size(QN_1), size(QN)) && max(abs(1 - qn_ratio(:))) < GlobalConstants.CoarseTol && (forkIter > 2)
                    forkLoop = false;
                else
                    if self.model.hasOpenClasses
                        sourceIndex = self.model.getSource.index;
                        UNnosource = UN; UNnosource(sourceIndex,:) = 0;
                        if any(find(sum(UNnosource(:,isinf(sn.njobs(1:size(QN,2)))),2)>0.99 * sn.nservers))
                            line_warning(mfilename,'The model may be unstable: the utilization of station %i for open classes exceeds 99 percent.\n',maxpos(sum(UNnosource,2)));
                        end
                    end
                    QN_1 = QN;
                end
            else
                forkLoop = false;
            end
            if strcmp(options.method,'exact')  && ~self.model.hasProductFormSolution
                line_error(mfilename,'The exact method requires the model to have a product-form solution. This model does not have one.\nYou can use Network.hasProductFormSolution() to check before running the solver.\n Run the ''mva'' method to obtain an approximation based on the exact MVA algorithm.\n');
            end
            if strcmp(options.method,'mva') && ~self.model.hasProductFormSolution
                line_warning(mfilename,'The exact method requires the model to have a product-form solution. This model does not have one.\nYou can use Network.hasProductFormSolution() to check before running the solver.\nSolverMVA will return an approximation generated by an exact MVA algorithm.');
            end

            method = options.method;

            % Check for size-based policies (SRPT, PSJF, FB, LRPT, SETF) in multiclass open systems
            queueIdx = find(sn.nodetype == NodeType.Queue);
            isSizeBasedPolicy = false;
            if ~isempty(queueIdx)
                schedType = sn.sched(sn.nodeToStation(queueIdx(1)));
                isSizeBasedPolicy = ismember(schedType, [SchedStrategy.SRPT, SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT, SchedStrategy.SETF]);
            end

            if sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink])) && isSizeBasedPolicy
                % Multiclass open system with size-based scheduling (SRPT, PSJF, FB, LRPT, SETF)
                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter] = solver_mva_qsys_sizebased_analyzer(sn, options, schedType);
            elseif sn.nclasses==1 && sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink])) % is an open queueing system
                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter] = solver_mva_qsys_analyzer(sn, options);
            elseif sn.nclasses>1 && sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Queue,NodeType.Sink])) && sn.sched(find(sn.nodetype==NodeType.Queue)) == SchedStrategy.POLLING % is an open polling system
                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter] = solver_mva_polling_analyzer(sn, options);
            elseif sn.nclosedjobs == 0 && length(sn.nodetype)==3 && all(sort(sn.nodetype)' == sort([NodeType.Source,NodeType.Cache,NodeType.Sink])) % is a non-rentrant cache
                % random initialization
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        prob = self.model.nodes{ind}.server.hitClass;
                        prob(prob>0) = 0.5;
                        self.model.nodes{ind}.setResultHitProb(prob);
                        self.model.nodes{ind}.setResultMissProb(1-prob);
                    end
                end
                self.model.refreshChains();
                % start iteration
                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_mva_cache_analyzer(sn, options);

                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        hitClass = self.model.nodes{ind}.getHitClass;
                        missClass = self.model.nodes{ind}.getMissClass;
                        hitprob = zeros(1,length(hitClass));
                        for k=1:length(self.model.nodes{ind}.getHitClass)
                            chain_k = sn.chains(:,k)>0;
                            inchain = sn.chains(chain_k,:)>0;
                            h = hitClass(k);
                            m = missClass(k);
                            if h>0 && m>0
                                hitprob(k) = XN(h) / sum(XN(inchain),"omitnan"); %#ok<NANSUM>
                            end
                        end
                        self.model.nodes{ind}.setResultHitProb(hitprob);
                        self.model.nodes{ind}.setResultMissProb(1-hitprob);
                    end
                end
                %self.model.refreshChains();
                self.model.refreshStruct(true);
            else % queueing network
                if any(sn.nodetype == NodeType.Cache) % if integrated caching-queueing
                    [QN,UN,RN,TN,CN,XN,lG,hitprob,missprob,runtime,lastiter] = solver_mva_cacheqn_analyzer(self, options);
                    for ind = 1:sn.nnodes
                        if sn.nodetype(ind) == NodeType.Cache
                            self.model.nodes{ind}.setResultHitProb(hitprob(ind,:));
                            self.model.nodes{ind}.setResultMissProb(missprob(ind,:));
                        end
                    end
                    %self.model.refreshChains();
                    self.model.refreshStruct(true);
                else % ordinary queueing network
                    switch method
                        case {'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'pb.upper', 'pb.lower', 'gb.upper', 'gb.lower', 'sb.upper', 'sb.lower'}
                            [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter] = solver_mva_bound_analyzer(sn, options);
                        otherwise
                            if ~isempty(sn.lldscaling) || ~isempty(sn.cdscaling)
                                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_mvald_analyzer(sn, options);
                            else
                                [QN,UN,RN,TN,CN,XN,lG,runtime,lastiter,actualmethod] = solver_mva_analyzer(sn, options);
                            end
                    end
                end
            end
            if self.model.hasFork
                nonfjstruct = sn;
                sn = self.getStruct;
                % Pre-compute linked routing matrix once (loop-invariant)
                switch options.config.fork_join
                    case {'mmt', 'default', 'fjt'}
                        Pcs = cell2mat(nonfjmodel.getLinkedRoutingMatrix);
                end
                for f=find(sn.nodetype == NodeType.Fork)'
                    switch options.config.fork_join
                        case {'mmt', 'default', 'fjt'}
                            TNfork = zeros(1,sn.nclasses);
                            for c=1:sn.nchains
                                inchain = find(sn.chains(c,:));
                                for r=inchain(:)'
                                    TNfork(r) =  (sn.nodevisits{c}(parent_forks(f),r) / sum(sn.visits{c}(sn.stationToStateful(sn.refstat(r)),inchain))) * sum(TN(sn.refstat(r),inchain));
                                end
                            end
                            % find join associated to fork f
                            joinIdx = find(sn.fj(f,:));
                            forkauxclasses = find(fjforkmap==f);
                            for s=forkauxclasses(:)'
                                r = fjclassmap(s); % original class associated to auxiliary class s
                                if isempty(joinIdx)
                                    forkLambda(s) = mean([forkLambda(s); TNfork(r)],1);
                                else
                                    joinStat = sn.nodeToStation(joinIdx);
                                    TN(sn.nodeToStation(joinIdx),r) = TN(sn.nodeToStation(joinIdx),r) + sum(TN(sn.nodeToStation(joinIdx), find(fjclassmap == r))) - TN(sn.nodeToStation(joinIdx), s);
                                    forkLambda(s) = mean([forkLambda(s); TN(sn.nodeToStation(joinIdx),r)],1);
                                end
                                if isempty(joinIdx) || ~outer_forks(f, r)
                                    % No join nodes for this fork, no synchronisation delay
                                    continue;
                                end
                                % Find the parallel paths coming out of the fork
                                ri = ModelAdapter.findPathsCS(sn, Pcs, f, joinIdx, r, [r,s], QN, TN, 0, fjclassmap, fjforkmap, nonfjmodel);
                                if isempty(ri)
                                    % No routing from fork for this class - set sync delay to 0 (Immediate)
                                    % Matches JAR behavior where empty Matrix gives syncDelay = 0
                                    syncDelay = 0;
                                else
                                    lambdai = 1./ri;
                                    d0 = 0;
                                    parallel_branches = length(ri);
                                    for pow=0:(parallel_branches - 1)
                                        current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                                        d0 = d0 + (-1)^pow * current_sum;
                                    end
                                    syncDelay = d0*sn.nodeparam{f}.fanOut - mean(ri);
                                end
                                % Set the synchronisation delays
                                nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{s}, Exp.fitMean(syncDelay));
                                if outer_forks(f, r)
                                    nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{r}, Exp.fitMean(syncDelay));
                                end
                            end
                        case {'heidelberger-trivedi', 'ht'}
                            joinIdx = find(sn.fj(f,:));
                            for c=1:sn.nchains
                                inchain = find(sn.chains(c,:));
                                for r=inchain(:)'
                                    if sn.nodevisits{c}(f,r) == 0
                                        continue;
                                    end
                                    % Obtain the response times on the parallel branches
                                    ri = RN(:, find(fjclassmap == r));
                                    ri(isnan(ri) | isinf(ri)) = 0;
                                    ri = sum(ri, 1, "omitnan") - RN(nonfjstruct.nodeToStation(fj_auxiliary_delays{joinIdx}), find(fjclassmap == r)) - RN(nonfjstruct.nodeToStation(joinIdx), find(fjclassmap == r));
                                    lambdai = 1./ri;
                                    d0 = 0;
                                    parallel_branches = length(self.model.nodes{f}.output.outputStrategy{r}{3});
                                    for pow=0:(parallel_branches - 1)
                                        current_sum = sum(1./sum(nchoosek(lambdai, pow + 1),2));
                                        d0 = d0 + (-1)^pow * current_sum;
                                    end
                                    di = d0*sn.nodeparam{f}.fanOut - ri;
                                    r0 = sum(RN(:, inchain), 2);
                                    r0(isnan(r0) | isinf(r0)) = 0;
                                    r0 = sum(r0, 1, "omitnan") - RN(nonfjstruct.nodeToStation(joinIdx), r);
                                    % Update the delays at the join node and at the auxiliary delay
                                    nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{r}, Exp.fitMean(d0*sn.nodeparam{f}.fanOut));
                                    idx = 1;
                                    for s=find(fjclassmap == r)
                                        nonfjmodel.nodes{joinIdx}.setService(nonfjmodel.classes{s}, Exp.fitMean(di(idx)));
                                        idx = idx + 1;
                                        nonfjmodel.nodes{fj_auxiliary_delays{joinIdx}}.setService(nonfjmodel.classes{s}, Exp.fitMean(r0));
                                    end

                                end
                            end
                    end
                end
                % Batch refreshRates after all sync delay updates (moved out of inner loops for performance)
                switch options.config.fork_join
                    case {'mmt', 'default', 'fjt'}
                        nonfjmodel.refreshRates();
                end
                switch options.config.fork_join
                    case {'heidelberger-trivedi', 'ht'}
                        nonfjmodel.refreshStruct();
                        % Delete the queue lengths, response times, throughputs and utilizations of the original classes at the join nodes
                        QN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap)) = 0;
                        RN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap)) = 0;
                        % Save the throughputs of the original classes at the join node
                        TN_orig = TN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap));
                        TN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap)) = 0;
                        UN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap)) = 0;

                        % Remove the times at the auxiliary delay
                        QN(nonfjstruct.nodeToStation(cell2mat(fj_auxiliary_delays)),:) = [];
                        UN(nonfjstruct.nodeToStation(cell2mat(fj_auxiliary_delays)),:) = [];
                        RN(nonfjstruct.nodeToStation(cell2mat(fj_auxiliary_delays)),:) = [];
                        TN(nonfjstruct.nodeToStation(cell2mat(fj_auxiliary_delays)),:) = [];
                        % merge back artificial classes into their original classes
                        for r=1:length(fjclassmap)
                            s = fjclassmap(r);
                            if s>0
                                QN(:,s) = QN(:,s) + QN(:,r);
                                UN(:,s) = UN(:,s) + UN(:,r);
                                % Add all throughputs of the auxiliary classes to facilitate the computation of the response times
                                TN(:,s) = TN(:,s) + TN(:,r);
                                RN(:,s) = QN(:,s) ./ TN(:,s);
                            end
                        end
                        % Re-set the throughputs for the original classes
                        TN(nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonzeros(fjclassmap)) = TN_orig;
                    case {'mmt', 'default', 'fjt'}
                        TN_orig = TN([nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Source))], nonzeros(fjclassmap));
                        % merge back artificial classes into their original classes
                        for r=1:length(fjclassmap)
                            s = fjclassmap(r);
                            if s>0
                                QN(:,s) = QN(:,s) + QN(:,r);
                                UN(:,s) = UN(:,s) + UN(:,r);
                                TN(:,s) = TN(:,s) + TN(:,r);
                                %RN(:,s) = RN(:,s) + RN(:,r);
                                % for i=find(snorig.nodetype == NodeType.Delay | snorig.nodetype == NodeType.Queue)'
                                %     TN(snorig.nodeToStation(i),s) = TN(snorig.nodeToStation(i),s) + TN(snorig.nodeToStation(i),r);
                                % end
                                RN(:,s) = QN(:,s) ./ TN(:,s);
                                %CN(:,s) = CN(:,s) + CN(:,r);
                                %XN(:,s) = XN(:,s) + XN(:,r);
                            end
                        end
                        TN([nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Join)), nonfjstruct.nodeToStation(find(sn.nodetype == NodeType.Source))], nonzeros(fjclassmap)) = TN_orig;
                end
                QN(:,fjclassmap>0) = [];
                UN(:,fjclassmap>0) = [];
                RN(:,fjclassmap>0) = [];
                TN(:,fjclassmap>0) = [];
                CN(:,fjclassmap>0) = [];
                XN(:,fjclassmap>0) = [];
            end
            iter = iter + lastiter;
            % Cap accumulated iterations for stability (Python parity)
            if iter > 10000
                iter = 10000;
                break; % Exit forkLoop early
            end
        end

        sn = self.model.getStruct();

        % Compute average residence time at steady-state
        AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
        WN = sn_get_residt_from_respt(sn, RN, self.getAvgResidTHandles());
        if strcmp(method,'default') && exist('actualmethod','var')
            self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,['default/' actualmethod],iter);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,method,iter);
        end
        self.result.Prob.logNormConstAggr = lG;
end
end


