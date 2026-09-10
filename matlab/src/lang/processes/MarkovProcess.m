classdef MarkovProcess < Process
    % A class for a continuous time Markov chain
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        infGen;
        stateSpace;
        isfinite;
    end

    methods
        function self = MarkovProcess(InfGen, isFinite, stateSpace)
            % SELF = MARKOVPROCESS(InfGen, isInfinite, stateSpace)
            self@Process('MarkovProcess', 1);
        
            self.infGen = ctmc_makeinfgen(InfGen);
            if nargin < 2
                self.isfinite = true;
            else
                self.isfinite = isFinite;
            end
            if nargin > 2
                self.stateSpace = stateSpace;
            else
                self.stateSpace = [];
            end
        end

        function A=toMarkovChain(self, q)
            if nargin==1
                q=(max(max(abs(self.infGen))))+rand;
            end
            P=self.infGen/q + eye(size(self.infGen));
            A=MarkovChain(P);
            A.setStateSpace(self.stateSpace);
        end
        
        function A=toDTMC(self, q)
            % TODTMC - Alias for toMarkovChain for backwards compatibility
            if nargin==1
                A = self.toMarkovChain();
            else
                A = self.toMarkovChain(q);
            end
        end

        function Qp = toTimeReversed(self)
            Qp = MarkovProcess(ctmc_timereverse(self.infGen));
        end

        function setStateSpace(self,stateSpace)
            self.stateSpace  = stateSpace;
        end

        function plot3(self)
            G = digraph(self.infGen-diag(diag(self.infGen)));
            nodeLbl = {};
            if ~isempty(self.stateSpace)
                for s=1:size(self.stateSpace,1)
                    if size(self.stateSpace,2)>1
                        nodeLbl{s} = sprintf('%s%d', sprintf('%d,', self.stateSpace(s,1:end-1)), self.stateSpace(s,end));
                    else
                        nodeLbl{s} = sprintf('%d', self.stateSpace(s,end));
                    end
                end
            end
            Q0 = self.infGen - diag(diag(self.infGen));
            [I,J,q]=find(Q0);
            edgeLbl = {};
            if ~isempty(self.stateSpace)
                for t=1:length(I)
                    edgeLbl{end+1,1} = nodeLbl{I(t)};
                    edgeLbl{end,2} = nodeLbl{J(t)};
                    edgeLbl{end,3} = strrep(strrep(sprintf('%.4f',q(t)),'000',''),'00','');
                end
            else
                for t=1:length(I)
                    edgeLbl{end+1,1} = num2str(I(t));
                    edgeLbl{end,2} = num2str(J(t));
                    edgeLbl{end,3} = strrep(strrep(sprintf('%.4f',q(t)),'000',''),'00','');
                end
            end
            h = plot(G,'Layout','force3','NodeLabel',nodeLbl,'EdgeLabel',edgeLbl(:,3));
        end

        function plot(self)
            if issym(self.infGen)
                line_error(mfilename,'CTMC.plot does not support symbolic CTMCs.');
            end
            G = digraph(self.infGen-diag(diag(self.infGen)));
            nodeLbl = {};
            if ~isempty(self.stateSpace)
                if iscell(self.stateSpace)
                    for s=1:size(self.stateSpace,1)
                        nodeLbl{s} = self.stateSpace{s};
                    end
                else
                    for s=1:size(self.stateSpace,1)
                        if size(self.stateSpace,2)>1
                            nodeLbl{s} = sprintf('%s%d', sprintf('%d,', self.stateSpace(s,1:end-1)), self.stateSpace(s,end));
                        else
                            nodeLbl{s} = sprintf('%d', self.stateSpace(s,end));
                        end
                    end
                end
            end
            Q0 = self.infGen - diag(diag(self.infGen));
            [I,J,q]=find(Q0);
            [~,sortIdx] = sortrows([I,J]);
            I=I(sortIdx);
            J=J(sortIdx);
            q=q(sortIdx);
            edgeLbl = {};
            if ~isempty(self.stateSpace)
                for t=1:length(I)
                    edgeLbl{end+1,1} = nodeLbl{I(t)};
                    edgeLbl{end,2} = nodeLbl{J(t)};
                    edgeLbl{end,3} = strrep(strrep(sprintf('%.4f',q(t)),'000',''),'00','');
                end
            else
                for t=1:length(I)
                    edgeLbl{end+1,1} = num2str(I(t));
                    edgeLbl{end,2} = num2str(J(t));
                    edgeLbl{end,3} = strrep(strrep(sprintf('%.4f',q(t)),'000',''),'00','');
                end
            end
            %             if length(nodeLbl) <= 6
            %                 colors = cell(1,length(nodeLbl)); for i=1:length(nodeLbl), colors{i}='w'; end
            %                 graphViz4Matlab('-adjMat',Q0,'-nodeColors',colors,'-nodeLabels',nodeLbl,'-edgeLabels',edgeLbl,'-layout',Circularlayout);
            %             else
            %                 graphViz4Matlab('-adjMat',Q0,'-nodeLabels',nodeLbl,'-edgeLabels',edgeLbl,'-layout',Springlayout);
            %             end
            h = plot(G,'Layout','force','NodeLabel',nodeLbl,'EdgeLabel',edgeLbl(:,3));
        end

        function Q = getGenerator(self)
            % Q = GETGENERATOR()

            % Get generator
            Q = self.infGen;
        end

        function [pi_i, num, den] = getProbState(self, state)
            % Use Cramer's rule to compute the probability of a single
            % state
            i = matchrow(self.stateSpace, state);
            Q = self.infGen; Q(:,1)=1;
            Q_i=Q; Q_i(i,:)=0; Q_i(i,1)=1;
            num=det(Q_i);
            den=det(Q);
            if issym(Q)
                pi_i=simplify(num/den);
            else
                pi_i=num/den;
            end
        end

        function pi = solve(self)
            if issym(self.infGen)
                pi = ctmc_solve(self.infGen);
            else
                pi = ctmc_solve_reducible(self.infGen);
            end
        end

        function [soujt, sts] = sample(self, n)
            if nargin<2
                n = 1;
            end
            [soujt, sts] = ctmc_simulate(self.infGen, [], n);
        end

        function [pi_t, kmax] = transient(self, pi0, t, method)
            % [PI_T, KMAX] = TRANSIENT(PI0, T, METHOD)
            % Distribution at time T from PI0 (uniform if empty). METHOD is
            % 'unif' (Jensen uniformization, the default) or 'foxglynn', whose
            % weights avoid evaluating the Poisson terms directly.
            if nargin<2 || isempty(pi0)
                pi0 = ones(1,length(self.infGen))/length(self.infGen);
            end
            if nargin<4 || isempty(method)
                method = 'unif';
            end
            switch lower(method)
                case 'foxglynn'
                    [pi_t, ~, right] = ctmc_foxglynn(reshape(pi0,1,[]), self.infGen, t);
                    kmax = right;
                otherwise
                    [pi_t, kmax] = ctmc_uniformization(reshape(pi0,1,[]), self.infGen, t);
            end
        end

        function p = solveRelative(self, refstate)
            % P = SOLVERELATIVE(REFSTATE)
            % Equilibrium distribution relative to REFSTATE, i.e. with
            % p(REFSTATE)=1. Unnormalized by construction, so it is defined
            % even where the normalizing constant is not.
            if nargin<2
                refstate = 1;
            end
            p = ctmc_relsolve(self.infGen, refstate);
        end

        function [p, eps, epsMAX] = aggregate(self, MS, method, param)
            % [P, EPS, EPSMAX] = AGGREGATE(MS, METHOD, PARAM)
            % Aggregation-disaggregation over the macrostate partition MS, a
            % cell array of state-index vectors. METHOD is 'courtois' (PARAM is
            % the randomization rate q), 'kms' or 'takahashi' (PARAM is the
            % iteration count, default 10), or 'multi' (PARAM is the
            % second-level partition MSS). EPS is the nearly-complete-
            % decomposability index of the partition and EPSMAX the largest
            % index for which the approximation is meant to hold.
            if nargin<3 || isempty(method)
                method = 'courtois';
            end
            if nargin<4
                param = [];
            end
            switch lower(method)
                case 'courtois'
                    if isempty(param)
                        [p,~,~,eps,epsMAX] = ctmc_courtois(self.infGen, MS);
                    else
                        [p,~,~,eps,epsMAX] = ctmc_courtois(self.infGen, MS, param);
                    end
                case 'kms'
                    if isempty(param), param = 10; end
                    [p,~,~,eps,epsMAX] = ctmc_kms(self.infGen, MS, param);
                case 'takahashi'
                    if isempty(param), param = 10; end
                    [p,~,~,~,eps,epsMAX] = ctmc_takahashi(self.infGen, MS, param);
                case 'multi'
                    if isempty(param)
                        line_error(mfilename,'The ''multi'' method requires the second-level partition MSS.');
                    end
                    [p,~,~,eps,epsMAX] = ctmc_multi(self.infGen, MS, param);
                otherwise
                    line_error(mfilename,'Unknown aggregation method ''%s''.', method);
            end
        end

        function [pi_t, kmax] = transientProb(self, pi0, t)
            % [PI_T, KMAX] = TRANSIENTPROB(PI0, T)
            % Alias of transient, under the name the JAR must use since
            % 'transient' is a Java keyword.
            if nargin<2
                pi0 = [];
            end
            [pi_t, kmax] = self.transient(pi0, t);
        end

        function [piTimeAvg, piExit] = timeAverage(self, pi0, t)
            % [PITIMEAVG, PIEXIT] = TIMEAVERAGE(PI0, T)
            % Time-averaged distribution over [0,T] and its endpoint.
            if nargin<2 || isempty(pi0)
                pi0 = ones(1,length(self.infGen))/length(self.infGen);
            end
            [piTimeAvg, piExit] = ctmc_timeaverage(reshape(pi0,1,[]), self.infGen, t);
        end

        function dpi = sens(self, dQ)
            % DPI = SENS(DQ)
            % Sensitivity of the stationary distribution to a scalar parameter,
            % given the derivative DQ of the generator.
            dpi = ctmc_sens(self.infGen, dQ, self.solve());
        end

        function S = stochComp(self, I)
            % S = STOCHCOMP(I)
            % Stochastic complement of the states I, a generator on that subset.
            % Use stochCompFull to also obtain the partitioned blocks.
            if nargin<2
                S = ctmc_stochcomp(self.infGen);
            else
                S = ctmc_stochcomp(self.infGen, I);
            end
        end

        function [S, Q11, Q12, Q21, Q22, T] = stochCompFull(self, I)
            % [S, Q11, Q12, Q21, Q22, T] = STOCHCOMPFULL(I)
            % Stochastic complement of the states I together with the blocks of
            % the generator partitioned by I and its complement, and the
            % return-path term T = Q12*inv(-Q22)*Q21, so that S = Q11 + T.
            if nargin<2
                [S, Q11, Q12, Q21, Q22, T] = ctmc_stochcomp(self.infGen);
            else
                [S, Q11, Q12, Q21, Q22, T] = ctmc_stochcomp(self.infGen, I);
            end
        end

        function h = hittingTime(self, targetStates)
            % H = HITTINGTIME(TARGETSTATES)
            % Mean TIME to reach any state in TARGETSTATES, zero on the target
            % set itself and Inf from a state that cannot reach it. The twin
            % MarkovChain.hittingTime counts STEPS instead, so the two answer
            % different questions about the same jump structure.
            h = ctmc_hitting_time(self.infGen, targetStates);
        end

        function bool = isFeasible(self)
            % BOOL = ISFEASIBLE()
            % True when the generator is a valid one.
            bool = ctmc_isfeasible(self.infGen);
        end

        function A = toEmbedded(self)
            % A = TOEMBEDDED()
            % Embedded jump chain, i.e. the DTMC of the states visited at
            % transition epochs. Unlike toDTMC (uniformization) it does not
            % preserve the stationary distribution, since it drops the holding
            % times; an absorbing state stays absorbing.
            Q = self.infGen;
            n = size(Q,1);
            exitRate = -diag(Q);
            P = Q - diag(diag(Q));
            for i = 1:n
                if exitRate(i) > 0
                    P(i,:) = P(i,:) / exitRate(i);
                else
                    P(i,i) = 1;
                end
            end
            A = MarkovChain(P);
            A.setStateSpace(self.stateSpace);
        end

    end

    methods (Static)
        function ctmcObj=rand(nStates) % creates a random CTMC
            ctmcObj = MarkovProcess(ctmc_rand(nStates));
        end

        function ctmcObj=fromSampleSysAggr(sa)
            isFinite = true;
            sampleState = sa.state{1};
            for r=2:length(sa.state)
                % per-node trajectories are time-aligned: join column-wise
                sampleState = [sampleState, sa.state{r}];
            end
            [stateSpace,~,stateHash] = unique(sampleState,'rows');
            dtmc = spalloc(length(stateSpace),length(stateSpace),length(stateSpace)); % assume O(n) elements with n states
            holdTime = zeros(length(stateSpace),1);
            for i=2:length(stateHash)
                if isempty(dtmc(stateHash(i-1),stateHash(i)))
                    dtmc(stateHash(i-1),stateHash(i)) = 0;
                end
                dtmc(stateHash(i-1),stateHash(i)) = dtmc(stateHash(i-1),stateHash(i)) + 1;
                holdTime(stateHash(i-1)) = holdTime(stateHash(i-1)) + sa.t(i) - sa.t(i-1);
            end
            % at this point, dtmc has absolute counts so not yet normalized
            holdTime = holdTime ./ sum(dtmc,2);
            infGen = ctmc_makeinfgen(dtmc_makestochastic(dtmc)./(holdTime*ones(1,length(stateSpace))));
            ctmcObj = MarkovProcess(infGen,isFinite,stateSpace);
        end
    end
end
