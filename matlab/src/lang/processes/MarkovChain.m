classdef MarkovChain < Process
    % An abstract class for a discrete time Markov chain
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        transMat;
        stateSpace;
        isfinite;
    end

    methods
        function self = MarkovChain(transMat, isFinite)
            % SELF = MARKOVCHAIN(transMat, isInfinite)
            self@Process('MarkovChain', 1);

            self.transMat = dtmc_makestochastic(transMat);
            self.stateSpace = [];
            if nargin < 2
                self.isfinite = true;
            else
                self.isfinite = isFinite;
            end
        end

        function A = toMarkovProcess(self)
            Q= self.transMat - eye(size(self.transMat));
            A=MarkovProcess(Q);
            A.setStateSpace(self.stateSpace);
        end
        
        function A = toCTMC(self)
            % TOCTMC - Alias for toMarkovProcess for backwards compatibility
            A = self.toMarkovProcess();
        end

        function Ap = toTimeReversed(self)
            Ap = MarkovChain(dtmc_timereverse(self.transMat));
        end

        function transMat = getTransMat(self)
            transMat = self.transMat;
        end

        function pi = solve(self)
            % PI = SOLVE()
            % Stationary distribution of the DTMC. Twin of MarkovProcess.solve.
            if issym(self.transMat)
                pi = dtmc_solve(self.transMat);
            else
                pi = dtmc_solve_reducible(self.transMat);
            end
        end

        function pi_t = transient(self, pi0, steps)
            % PI_T = TRANSIENT(PI0, STEPS)
            % Distribution at each step 0,...,STEPS from PI0 (uniform if empty).
            if nargin<2
                pi0 = [];
            end
            if nargin<3
                steps = 1;
            end
            pi_t = dtmc_transient(self.transMat, pi0, steps);
        end

        function pi_t = transientProb(self, pi0, steps)
            % PI_T = TRANSIENTPROB(PI0, STEPS)
            % Alias of transient, under the name the JAR must use since
            % 'transient' is a Java keyword.
            if nargin<2
                pi0 = [];
            end
            if nargin<3
                steps = 1;
            end
            pi_t = self.transient(pi0, steps);
        end

        function h = hittingTime(self, targetStates)
            % H = HITTINGTIME(TARGETSTATES)
            % Mean number of steps to reach any state in TARGETSTATES.
            h = dtmc_hitting_time(self.transMat, targetStates);
        end

        function S = stochComp(self, I)
            % S = STOCHCOMP(I)
            % Stochastic complement of the states I, a DTMC on that subset.
            % Use stochCompFull to also obtain the partitioned blocks.
            if nargin<2
                S = dtmc_stochcomp(self.transMat);
            else
                S = dtmc_stochcomp(self.transMat, I);
            end
        end

        function [S, P11, P12, P21, P22] = stochCompFull(self, I)
            % [S, P11, P12, P21, P22] = STOCHCOMPFULL(I)
            % Stochastic complement of the states I together with the blocks of
            % the transition matrix partitioned by I and its complement.
            if nargin<2
                [S, P11, P12, P21, P22] = dtmc_stochcomp(self.transMat);
            else
                [S, P11, P12, P21, P22] = dtmc_stochcomp(self.transMat, I);
            end
        end

        function [pi_t, kmax] = transientUnif(self, pi0, t)
            % [PI_T, KMAX] = TRANSIENTUNIF(PI0, T)
            % Distribution at time T of the DTMC seen through uniformization.
            % The chain is read as the randomized image of a CTMC, so T is
            % continuous here, unlike the step count taken by transient.
            n = size(self.transMat,1);
            if nargin<2 || isempty(pi0)
                pi0 = ones(1,n)/n;
            end
            if nargin<3
                t = 1;
            end
            [pi_t, kmax] = dtmc_uniformization(reshape(pi0,1,[]), self.transMat, t);
        end

        function bool = isFeasible(self)
            % BOOL = ISFEASIBLE()
            % True when the transition matrix is stochastic.
            bool = dtmc_isfeasible(self.transMat) > 0;
        end

        function sts = sample(self, n)
            % STS = SAMPLE(N) - Simulate n steps of the DTMC from a random initial state
            if nargin<2
                n = 1;
            end
            pi0 = rand(1, size(self.transMat,1)); pi0 = pi0/sum(pi0);
            sts = dtmc_simulate(self.transMat, pi0, n);
        end

        function setStateSpace(self,stateSpace)
            self.stateSpace  = stateSpace;
        end

        function plot(self)
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
            P0 = self.transMat;
            [I,J,q]=find(P0);
            edgeLbl = {};
            if ~isempty(self.stateSpace)
                for t=1:length(I)
                    edgeLbl{end+1,1} = nodeLbl{I(t)};
                    edgeLbl{end,2} = nodeLbl{J(t)};
                    edgeLbl{end,3} = sprintf('%.2f',(q(t)));
                end
            else
                for t=1:length(I)
                    edgeLbl{end+1,1} = num2str(I(t));
                    edgeLbl{end,2} = num2str(J(t));
                    edgeLbl{end,3} = sprintf('%.2f',(q(t)));
                end
            end
            if length(nodeLbl) <= 6
                colors = cell(1,length(nodeLbl)); for i=1:length(nodeLbl), colors{i}='w'; end
                graphViz4Matlab('-adjMat',P0,'-nodeColors',colors,'-nodeLabels',nodeLbl,'-edgeLabels',edgeLbl,'-layout',Circularlayout);
            else
                graphViz4Matlab('-adjMat',P0,'-nodeLabels',nodeLbl,'-edgeLabels',edgeLbl,'-layout',Springlayout);
            end
        end

    end

    methods (Static)
        function dtmcObj=rand(nStates) % creates a random DTMC
            dtmcObj = MarkovChain(dtmc_rand(nStates));
        end

        function dtmcObj=fromSampleSysAggr(sa)
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
            dtmc = dtmc_makestochastic(dtmc);
            dtmcObj = MarkovChain(dtmc, isFinite);
            dtmcObj.setStateSpace(stateSpace);
        end
    end
end
