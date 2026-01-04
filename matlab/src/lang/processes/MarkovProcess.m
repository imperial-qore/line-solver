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
            end
        end

        function pi = solve(self)
            if issym(self.infGen)
                pi = ctmc_solve(self.infGen);
            else
                pi = ctmc_solve_reducible(self.infGen);
            end
        end

        function X = sample(self)
            X = ctmc_sample(self.infGen,1);
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
            holdTime = holdTime ./ sum(dtmc,2);
            infGen = ctmc_makeinfgen(dtmc_makestochastic(dtmc)./(holdTime*ones(1,length(stateSpace))));
            ctmcObj = MarkovProcess(infGen,isFinite,stateSpace);
        end
    end
end
