classdef Activity < LayeredNetworkElement
    % A stage of service in a Task of a LayeredNetwork.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
        
    properties
        hostDemand;
        hostDemandMean;             %double
        hostDemandSCV;              %double
        parent;
        parentName;                 %string
        boundToEntry;               %string
        callOrder;                  %string \in {'STOCHASTIC', 'DETERMINISTIC'}
        syncCallDests = cell(0);   %string array
        syncCallMeans = [];        %integer array
        asyncCallDests = cell(0);  %string array
        asyncCallMeans = [];       %integer array
        scheduling = [];
        thinkTime;                 %Distribution object
        thinkTimeMean;             %double
        thinkTimeSCV;              %double
        phase = 1;                 %integer: phase number (1 or 2)
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function obj = Activity(model, name, hostDemand, boundToEntry, callOrder)
            % OBJ = ACTIVITY(MODEL, NAME, HOSTDEMAND, BOUNDTOENTRY, CALLORDER)
            name = char(name);
            if nargin<2 %~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            obj@LayeredNetworkElement(name);
            
            if nargin<3 %~exist('hostDemand','var')
                hostDemand = GlobalConstants.FineTol;
            elseif isnumeric(hostDemand) && hostDemand == 0
                hostDemand = GlobalConstants.FineTol;
            end
            if nargin<4 %~exist('boundToEntry','var')
                boundToEntry = '';
            end
            if nargin<5 %~exist('callOrder','var')
                callOrder = 'STOCHASTIC';
            end
            
            obj.setHostDemand(hostDemand);
            obj.boundToEntry = boundToEntry;
            obj.setCallOrder(callOrder);

            % Initialize default think time
            obj.thinkTime = Immediate.getInstance();
            obj.thinkTimeMean = GlobalConstants.FineTol;
            obj.thinkTimeSCV = GlobalConstants.FineTol;

            if isa(model,'LayeredNetwork')
                model.activities{end+1} = obj;
                obj.model = model;
            elseif isa(model,'JLayeredNetwork')
                % JLayeredNetwork support would go here if it exists
                obj.obj = jline.lang.layered.Activity(model.obj, name, hostDemand, boundToEntry, callOrder);
                model.addActivity(obj);
            end
        end
        
        function obj = setParent(obj, parent)
            % OBJ = SETPARENT(OBJ, PARENT)
            
            if isa(parent,'Entry') || isa(parent,'Task')
                obj.parentName = parent.name;
                obj.parent = parent;
            else
                obj.parentName = parent;
                obj.parent = [];
            end
        end
        
        function self = on(self, parent)
            % OBJ = ON(OBJ, PARENT)
            if ~isa(parent,'Task')
                line_error(mfilename,'Invalid .on() argument: expected a task.')
            end
            if isempty(self.parent)
                parent.addActivity(self);
                self.parent = parent;
            else
                line_error(mfilename,'Parent task already defined.')
            end
        end
        
        function obj = setHostDemand(obj, hostDemand)
            % OBJ = SETHOSTDEMAND(OBJ, HOSTDEMAND)
            
            if isnumeric(hostDemand)
                if hostDemand <= GlobalConstants.FineTol
                    obj.hostDemand = Immediate.getInstance();
                    obj.hostDemandMean = GlobalConstants.FineTol;
                    obj.hostDemandSCV = GlobalConstants.FineTol;
                else
                    obj.hostDemand = Exp(1/hostDemand);
                    obj.hostDemandMean = hostDemand;
                    obj.hostDemandSCV = 1.0;
                end
            elseif isa(hostDemand,'Distribution')
                obj.hostDemand = hostDemand;
                obj.hostDemandMean = hostDemand.getMean();
                obj.hostDemandSCV = hostDemand.getSCV();
            end
        end
        
        function obj = repliesTo(obj, entry)
            % OBJ = REPLIESTO(OBJ, ENTRY)
            
            if ~isempty(obj.parent)
                switch SchedStrategy.fromText(obj.parent.scheduling)
                    case SchedStrategy.REF
                        line_error(mfilename,'Activities in reference tasks cannot reply.');
                    otherwise
                        entry.replyActivity{end+1} = obj.name;
                end
            else
                entry.replyActivity{end+1} = obj.name;
            end
        end
        
        function obj = boundTo(obj, entry)
            % OBJ = BOUNDTO(OBJ, ENTRY)
            
            if isa(entry,'Entry')
                obj.boundToEntry = entry.name;
            elseif ischar(entry)
                obj.boundToEntry = entry;
            else
                line_error(mfilename,'Wrong entry parameter for boundTo method.');
            end
        end
        
        function obj = setCallOrder(obj, callOrder)
            % OBJ = SETCALLORDER(OBJ, CALLORDER)

            if strcmpi(callOrder,'STOCHASTIC') || strcmpi(callOrder,'DETERMINISTIC')
                obj.callOrder = upper(callOrder);
            else
                obj.callOrder = 'STOCHASTIC';
            end
        end

        function obj = setThinkTime(obj, thinkTime)
            % OBJ = SETTHINKTIME(OBJ, THINKTIME)

            if isnumeric(thinkTime)
                if thinkTime <= GlobalConstants.FineTol
                    obj.thinkTime = Immediate.getInstance();
                    obj.thinkTimeMean = GlobalConstants.FineTol;
                    obj.thinkTimeSCV = GlobalConstants.FineTol;
                else
                    obj.thinkTime = Exp(1/thinkTime);
                    obj.thinkTimeMean = thinkTime;
                    obj.thinkTimeSCV = 1.0;
                end
            elseif isa(thinkTime,'Distribution')
                obj.thinkTime = thinkTime;
                obj.thinkTimeMean = thinkTime.getMean();
                obj.thinkTimeSCV = thinkTime.getSCV();
            end
        end

        function obj = setPhase(obj, phaseNum)
            % OBJ = SETPHASE(OBJ, PHASENUM)
            % Set the phase number for this activity.
            % Phase 1: activities before the reply is sent
            % Phase 2: activities after the reply is sent (post-reply processing)

            if ~isnumeric(phaseNum) || phaseNum < 1 || phaseNum > 2
                line_error(mfilename, 'Phase must be 1 or 2.');
            end
            obj.phase = phaseNum;
        end

        %synchCall
       
        function obj = synchCall(obj, synchCallDest, synchCallMean)
            % OBJ = SYNCHCALL(OBJ, SYNCHCALLDEST, SYNCHCALLMEAN)
           
            if nargin<3 %~exist('synchCallMean','var')
                synchCallMean = 1.0;
            end
            
            % Get the destination entry name
            if ischar(synchCallDest)
                destName = synchCallDest;
            else % object
                destName = synchCallDest.name;
            end
            
            % Check if this entry is already in the synchronous call destinations
            if ~isempty(obj.syncCallDests)
                for i = 1:length(obj.syncCallDests)
                    if strcmp(obj.syncCallDests{i}, destName)
                        error_msg = sprintf('Activity "%s" already has a synchronous call to entry "%s". ', ...
                            obj.name, destName);
                        error_msg = [error_msg, 'If you intend multiple calls, change the mean number of calls using synchCall(entry, mean_calls) or define an activity graph with explicit precedences.'];
                        line_error(mfilename, error_msg);
                    end
                end
            end
            
            % Add the new call destination
            obj.syncCallDests{length(obj.syncCallDests)+1} = destName;
            obj.syncCallMeans = [obj.syncCallMeans; synchCallMean];
        end
        
        %asynchCall
        function obj = asynchCall(obj, asynchCallDest, asynchCallMean)
            % OBJ = ASYNCHCALL(OBJ, ASYNCHCALLDEST, ASYNCHCALLMEAN)
            
            if nargin<3 %~exist('asynchCallMean','var')
                asynchCallMean = 1.0;
            end
            
            % Get the destination entry name
            if ischar(asynchCallDest)
                destName = asynchCallDest;
            else % object
                destName = asynchCallDest.name;
            end
            
            % Check if this entry is already in the asynchronous call destinations
            if ~isempty(obj.asyncCallDests)
                for i = 1:length(obj.asyncCallDests)
                    if strcmp(obj.asyncCallDests{i}, destName)
                        error_msg = sprintf('Activity "%s" already has an asynchronous call to entry "%s". ', ...
                            obj.name, destName);
                        error_msg = [error_msg, 'If you intend multiple calls, change the mean number of calls using asynchCall(entry, mean_calls) or define an activity graph with explicit precedences.'];
                        line_error(mfilename, error_msg);
                    end
                end
            end
            
            % Add the new call destination
            obj.asyncCallDests{length(obj.asyncCallDests)+1} = destName;
            obj.asyncCallMeans = [obj.asyncCallMeans; asynchCallMean];
        end

        % Getter methods for API consistency with Java/Python
        function val = getHostDemand(obj)
            % GETHOSTDEMAND Get the host demand distribution
            val = obj.hostDemand;
        end

        function val = getHostDemandMean(obj)
            % GETHOSTDEMANDMEAN Get the mean host demand
            val = obj.hostDemandMean;
        end

        function val = getHostDemandSCV(obj)
            % GETHOSTDEMANDSCV Get the SCV of host demand
            val = obj.hostDemandSCV;
        end

        function val = getCallOrder(obj)
            % GETCALLORDER Get the call order (STOCHASTIC or DETERMINISTIC)
            val = obj.callOrder;
        end

        function val = getBoundToEntry(obj)
            % GETBOUNDTOENTRY Get the entry this activity is bound to
            val = obj.boundToEntry;
        end

        function val = getParent(obj)
            % GETPARENT Get the parent task
            val = obj.parent;
        end

        function val = getSyncCallDests(obj)
            % GETSYNCCALLDESTS Get the synchronous call destinations
            val = obj.syncCallDests;
        end

        function val = getSyncCallMeans(obj)
            % GETSYNCCALLMEANS Get the synchronous call means
            val = obj.syncCallMeans;
        end

        function val = getAsyncCallDests(obj)
            % GETASYNCCALLDESTS Get the asynchronous call destinations
            val = obj.asyncCallDests;
        end

        function val = getAsyncCallMeans(obj)
            % GETASYNCCALLMEANS Get the asynchronous call means
            val = obj.asyncCallMeans;
        end

        function val = getThinkTimeMean(obj)
            % GETTHINKTIMEMEAN Get the mean think time
            val = obj.thinkTimeMean;
        end

        function val = getPhase(obj)
            % GETPHASE Get the phase number (1 or 2)
            val = obj.phase;
        end

    end

end