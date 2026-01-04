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
                sampleState = State.cartesian(sampleState, sa.state{r});
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
        end
    end
end
