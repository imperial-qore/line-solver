function out = fjFixedPoint(self, options, solveFcn)
% OUT = FJFIXEDPOINT(OPTIONS, SOLVEFCN)
%
% Solver-agnostic driver of the fork-join fixed point. The transformation
% (ModelAdapter.mmt, or ModelAdapter.ht under options.config.fork_join) turns
% the model into a plain network in which every fork is a router, every join a
% zero-service delay, and the parallelism is carried by auxiliary open classes
% of arrival rate (fanout-1)*forkLambda. Each pass solves that network,
% recomputes the synchronisation delays from the resulting metrics and updates
% forkLambda; the loop ends when the queue lengths stop moving.
%
% SOLVEFCN(sn, options) is the inner solve. It must return a struct with the
% fields QN, UN, RN, TN, CN, XN, lG, runtime, lastiter, method and
% actualmethod, i.e. the contract of @SolverMVA/mvaDispatch.m. Nothing in the
% loop reads a solver internal, so any NetworkSolver whose analyzer honours
% that contract can be driven here.
%
% On a model without forks the loop runs SOLVEFCN exactly once and returns its
% result unchanged.
%
% OUT carries QN, UN, RN, TN, CN, XN, lG, runtime, iter, method and
% actualmethod.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter = 0;
        sn = getStruct(self); % doesn't need initial state
        forkLoop = true;
        forkIter = 0;
        % create artificial classes arrival rates; warm-start from the retained
        % iterate. see _kb/05-solvers-overview.md for rationale
        forkLambda = GlobalConstants.FineTol * ones(1, 2*sn.nclasses*sum(sn.nodetype==NodeType.Fork));
        if options.config.fj_warmstart && ~isempty(self.fjForkLambda) ...
                && isequal(size(self.fjForkLambda), size(forkLambda))
            forkLambda = self.fjForkLambda;
        end
        QN = GlobalConstants.Immediate * ones(1, sn.nclasses);
        QN_1 = 0*QN;
        UN = 0*QN;
while (forkLoop && forkIter < options.iter_max)
            if self.model.hasFork
                forkIter = forkIter + 1;
                line_debug(options, 'Fork-join iteration %d', forkIter);
                if forkIter == 1
                    switch options.config.fork_join
                        case {'heidelberger-trivedi', 'ht'}
                            [nonfjmodel, fjclassmap, fjforkmap, fj_auxiliary_delays] = ModelAdapter.ht(self.model);
                            line_debug(options, 'Fork-join method: heidelberger-trivedi');
                        case {'mmt', 'default', 'fjt'}
                            % Reuse the MMT transformation across outer iterations
                            % while structVersion is unchanged; re-fed from the base
                            % model first. see _kb/05-solvers-overview.md for rationale
                            cacheUsable = ~isempty(self.mmtCache) && ...
                                self.mmtCache.structVersion == self.model.structVersion && ...
                                ModelAdapter.refreshServicesFromBase(self.mmtCache.nonfjmodel, self.mmtCache.prov);
                            if cacheUsable
                                nonfjmodel = self.mmtCache.nonfjmodel;
                                fjclassmap = self.mmtCache.fjclassmap;
                                fjforkmap = self.mmtCache.fjforkmap;
                                fanout = self.mmtCache.fanout;
                                outer_forks = self.mmtCache.outer_forks;
                                parent_forks = self.mmtCache.parent_forks;
                                nonfjmodel.refreshRates();
                                line_debug(options, 'Fork-join method: mmt (cached transformation)');
                            else
                                [nonfjmodel, fjclassmap, fjforkmap, fanout, prov] = ModelAdapter.mmt(self.model, forkLambda);
                                line_debug(options, 'Fork-join method: mmt');
                                [outer_forks, parent_forks] = ModelAdapter.sortForks(sn, fjforkmap, fjclassmap, nonfjmodel);
                                self.mmtCache = struct('nonfjmodel', nonfjmodel, 'prov', prov, ...
                                    'fjclassmap', fjclassmap, 'fjforkmap', fjforkmap, 'fanout', fanout, ...
                                    'outer_forks', outer_forks, 'parent_forks', parent_forks, ...
                                    'structVersion', self.model.structVersion);
                            end
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
                line_debug(options, 'Fork-join iter %d: rebuilt nonfjmodel struct (nstations=%d, nclasses=%d)', forkIter, sn.nstations, sn.nclasses);
                % Mixed absolute/relative convergence test on the MMT iterate,
                % using the SOLVER's iter_tol (not CoarseTol). see
                % _kb/05-solvers-overview.md for rationale
                qn_converged = abs(QN_1 - QN) <= GlobalConstants.Zero + options.iter_tol*abs(QN);
                if isequal(size(QN_1), size(QN)) && all(qn_converged(:)) && (forkIter > 2)
                    line_debug(options, 'Fork-join iter %d: converged (mixed abs/rel test)', forkIter);
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
            fjres = solveFcn(sn, options);
            QN = fjres.QN; UN = fjres.UN; RN = fjres.RN; TN = fjres.TN;
            CN = fjres.CN; XN = fjres.XN; lG = fjres.lG;
            runtime = fjres.runtime; lastiter = fjres.lastiter;
            method = fjres.method; actualmethod = fjres.actualmethod;
            if self.model.hasFork
                line_debug(options, 'Fork-join post-processing: method=%s, computing sync delays', options.config.fork_join);
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
                % Drop the auxiliary-class columns. CN is left untouched when not
                % class-indexed. see _kb/05-solvers-overview.md for rationale
                auxmask = fjclassmap>0;
                QN = dropAuxCols(QN, auxmask);
                UN = dropAuxCols(UN, auxmask);
                RN = dropAuxCols(RN, auxmask);
                TN = dropAuxCols(TN, auxmask);
                CN = dropAuxCols(CN, auxmask);
                XN = dropAuxCols(XN, auxmask);
            end
            iter = iter + lastiter;
            % Cap accumulated iterations for stability (Python parity)
            if iter > 10000
                iter = 10000;
                break; % Exit forkLoop early
            end
        end
        % see _kb/05-solvers-overview.md for rationale (non-convergence now warns)
        if self.model.hasFork
            if forkLoop && forkIter >= options.iter_max
                line_warning(mfilename,'The fork-join (%s) fixed point did not converge in options.iter_max=%d iterations; returning the interim solution.\n', options.config.fork_join, options.iter_max);
            end
            % Retain the MMT iterate so that a subsequent runAnalyzer call on
            % this solver (an outer LN iteration) resumes the fixed point here.
            self.fjForkLambda = forkLambda;
        end

out = struct('QN', QN, 'UN', UN, 'RN', RN, 'TN', TN, 'CN', CN, 'XN', XN, ...
    'lG', lG, 'runtime', runtime, 'iter', iter, 'method', method, ...
    'actualmethod', actualmethod);
end

function A = dropAuxCols(A, auxmask)
% A = DROPAUXCOLS(A, AUXMASK)
% Remove the auxiliary-class columns from a class-indexed metric. Arrays that
% are not indexed by the transformed model's classes (e.g. a per-chain system
% response time) are returned unchanged.
if ~isempty(A) && size(A,2) == numel(auxmask)
    A(:,auxmask) = [];
end
end
