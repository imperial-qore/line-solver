classdef MarkedMarkovProcess < MarkovProcess
    % A class for continuous time Markov chain where transitions are labeled
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        eventFilt;
        eventList;
    end

    methods
        function self = MarkedMarkovProcess(InfGen, eventFilt, evs, isFinite, stateSpace)
            % SELF = MARKEDMARKOVPROCESS(InfGen, eventFilt, events, isFinite, stateSpace)
            if nargin<4
                isFinite = true;
            end
            self@MarkovProcess(InfGen, isFinite);
            self.name = 'MarkedMarkovProcess';
            self.eventFilt = eventFilt;
            self.eventList = evs;
            if nargin>=5
                self.setStateSpace(stateSpace);
            end
        end

        function eMAP = toMAP(self, varargin)
            % MAP = TOMAP(EVENT)

            if length(varargin)==1
                ev = varargin{1};
            else
                evtype = varargin{1};
                node = varargin{2};
                if ~isnumeric(node)
                    node = node.index;
                end
                class = varargin{3};
                if ~isnumeric(class)
                    class = class.index;
                end
                ev = Event(evtype, node, class);
            end

            for e=1:length(self.eventList) % for every event
                if self.eventList{e}.active{1}.node == ev.node && self.eventList{e}.active{1}.class == ev.class && self.eventList{e}.active{1}.event == ev.event
                    eMAP = MAP(self.infGen - self.eventFilt{e}, self.eventFilt{e});
                    break;
                end
            end
        end

        function pie = embeddedSolve(self, evset)
            % PIE = EMBEDDEDSOLVE(EVENTSET)

            if nargin<2
                evset = 1:length(self.eventFilt);
            end
            for ev=evset % for every event
                eMAP = MAP(self.infGen - self.eventFilt{ev}, self.eventFilt{ev});
                pie{ev} = eMAP.getEmbeddedProb;
                if issym(pie{ev})
                    pie{ev}=simplify(pie{ev})';
                end
            end
        end
    end

    methods (Static)
        function ctmcObj=fromSampleSysAggr(sa)
            isFinite = true;
            sampleState = sa.state{1};
            for r=2:length(sa.state)
                sampleState = State.cartesian(sampleState, sa.state{r});
            end
            [stateSpace,~,stateHash] = unique(sampleState,'rows');
            dtmc = spalloc(length(stateSpace),length(stateSpace),length(stateSpace));
            holdTime = zeros(length(stateSpace),1);
            activeEvent = cell2mat(cellfun(@(c) c.getProcess, sa.event, 'UniformOutput', false));
            [uniqueActiveEvents,~,evType] = unique(activeEvent,'rows');
            eventFilt = cell(1,size(uniqueActiveEvents,1));
            for e=1:size(uniqueActiveEvents,1)
                eventFilt{e} = spalloc(length(stateSpace),length(stateSpace),length(stateSpace));
            end
            for i=2:length(stateHash)
                e = evType(i-1);
                if isempty(dtmc(stateHash(i-1),stateHash(i)))
                    dtmc(stateHash(i-1),stateHash(i)) = 0;
                    eventFilt{e}(stateHash(i-1),stateHash(i)) = 0;
                end
                dtmc(stateHash(i-1),stateHash(i)) = dtmc(stateHash(i-1),stateHash(i)) + 1;
                eventFilt{e}(stateHash(i-1),stateHash(i)) = eventFilt{e}(stateHash(i-1),stateHash(i)) + 1;
                holdTime(stateHash(i-1)) = holdTime(stateHash(i-1)) + sa.t(i) - sa.t(i-1);
            end
            % at this point, dtmc has absolute counts so not yet normalized
            holdTime = holdTime ./ sum(dtmc,2);
            for e=1:size(uniqueActiveEvents,1)
                eventFilt{e} = eventFilt{e} ./ sum(dtmc,2);
            end
            infGen = ctmc_makeinfgen(dtmc_makestochastic(dtmc)./(holdTime*ones(1,length(stateSpace))));
            for e=1:size(uniqueActiveEvents,1)
                eventFilt{e} = eventFilt{e} .* (-diag(infGen)*ones(1,length(stateSpace)));
            end
            evList = cell(size(activeEvent,1),1);
            for e=1:length(evList)
                evList{e,1} = struct('active',sa.event{e},'passive',[]);                
            end
            ctmcObj = MarkedMarkovProcess(infGen,eventFilt,evList,isFinite,stateSpace);
        end
    end
end
