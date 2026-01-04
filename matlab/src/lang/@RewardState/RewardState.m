classdef RewardState < handle
    % REWARDSTATE Smart state accessor for reward functions
    %
    % Provides intuitive access to aggregated state populations at specific
    % stations and job classes, with support for aggregation operations.
    %
    % Usage:
    %   % Access specific station/class population (returns scalar)
    %   pop = state.at(queue1, class1);
    %
    %   % Access all classes at a station (returns RewardStateView)
    %   view = state.at(queue1);
    %   total = view.total();  % Sum across classes
    %
    %   % Access class across all stations (returns RewardStateView)
    %   view = state.forClass(class1);
    %   total = view.total();  % Sum across stations
    %
    % See also: RewardStateView, Reward, setReward

    properties (Access = private)
        stateVector         % Row vector: aggregated state [1 x (M*K)]
        sn                  % NetworkStruct reference
        nclasses            % Number of job classes
        nstations           % Number of stations
        nodeToStationMap    % containers.Map: node.index -> station index
        classToIndexMap     % containers.Map: jobclass.index -> class index
    end

    methods
        function self = RewardState(stateVec, sn, nodesToStation, classesToIdx)
            % REWARDSTATE Constructor
            %
            % SELF = REWARDSTATE(STATEVEC, SN, NODESTOSTATION, CLASSESTOIDX)
            %
            % STATEVEC        - Row vector of aggregated state [1 x (M*K)]
            % SN              - NetworkStruct with model information
            % NODESTOSTATION  - containers.Map: node.index -> station index
            % CLASSESTOIDX    - containers.Map: jobclass.index -> class index

            self.stateVector = stateVec;
            self.sn = sn;
            self.nclasses = sn.nclasses;
            self.nstations = sn.nstations;
            self.nodeToStationMap = nodesToStation;
            self.classToIndexMap = classesToIdx;
        end

        function result = at(self, node, jobclass)
            % AT Access population at station or station+class
            %
            % RESULT = AT(SELF, NODE) returns RewardStateView for all classes
            %          at the station containing NODE
            %
            % RESULT = AT(SELF, NODE, JOBCLASS) returns scalar population count
            %          of JOBCLASS at the station containing NODE
            %
            % Examples:
            %   view = state.at(queue1);              % All classes at queue1
            %   pop = state.at(queue1, class1);       % Class1 jobs at queue1
            %   total = state.at(queue1).total();     % Total jobs at queue1

            if nargin == 2
                % Return view for all classes at this station
                stationIdx = self.nodeToStationMap(node.index);
                startIdx = (stationIdx - 1) * self.nclasses + 1;
                endIdx = stationIdx * self.nclasses;

                if startIdx < 1 || endIdx > length(self.stateVector)
                    error('RewardState:InvalidIndex', ...
                        'Invalid state indices for station %d: [%d:%d]', ...
                        stationIdx, startIdx, endIdx);
                end

                subvector = self.stateVector(startIdx:endIdx);
                result = RewardStateView(subvector, self, stationIdx, []);
            elseif nargin == 3
                % Return scalar for specific station/class
                if ~isa(jobclass, 'JobClass')
                    error('RewardState:InvalidArgument', ...
                        'Second argument must be a JobClass object');
                end

                stationIdx = self.nodeToStationMap(node.index);
                classIdx = self.classToIndexMap(jobclass.index);

                if classIdx < 1 || classIdx > self.nclasses
                    error('RewardState:InvalidClass', ...
                        'Invalid class index %d for %s', classIdx, jobclass.name);
                end

                idx = (stationIdx - 1) * self.nclasses + classIdx;

                if idx < 1 || idx > length(self.stateVector)
                    error('RewardState:InvalidIndex', ...
                        'Invalid state index %d for station %d, class %s', ...
                        idx, stationIdx, jobclass.name);
                end

                result = self.stateVector(idx);
            else
                error('RewardState:InvalidNumberOfArguments', ...
                    'at() takes 2 or 3 arguments (self, node, [jobclass])');
            end
        end

        function result = forClass(self, jobclass)
            % FORCLASS Get populations for class across all stations
            %
            % RESULT = FORCLASS(SELF, JOBCLASS) returns RewardStateView
            %          showing population of JOBCLASS at each station
            %
            % Example:
            %   view = state.forClass(class1);
            %   total = view.total();  % Total class1 jobs in system
            %   view.atStation(queue1); % Class1 population at queue1

            if ~isa(jobclass, 'JobClass')
                error('RewardState:InvalidArgument', ...
                    'Argument must be a JobClass object');
            end

            classIdx = self.classToIndexMap(jobclass.index);

            if classIdx < 1 || classIdx > self.nclasses
                error('RewardState:InvalidClass', ...
                    'Invalid class index %d for %s', classIdx, jobclass.name);
            end

            values = zeros(1, self.nstations);
            for ist = 1:self.nstations
                idx = (ist - 1) * self.nclasses + classIdx;
                values(ist) = self.stateVector(idx);
            end

            result = RewardStateView(values, self, [], classIdx);
        end

        function result = forStation(self, station)
            % FORSTATION Alias for at(station)
            %
            % RESULT = FORSTATION(SELF, STATION) is equivalent to AT(SELF, STATION)
            %
            % Provided for semantic clarity when accessing all classes at a station

            result = self.at(station);
        end

        function sn = get_sn(self)
            % GET_SN Return NetworkStruct (for advanced usage)
            %
            % SN = GET_SN(SELF) returns the NetworkStruct, allowing advanced
            % users to access rates, capacities, and other parameters.
            %
            % Example (advanced):
            %   model.setReward('Advanced', @(state, sn) ...
            %       state.at(queue1) * sn.rates(sn.nodeToStation(queue1.index), 1));

            sn = self.sn;
        end
    end

    methods
        % Property getter for backward compatibility with @(state, sn) syntax
        function varargout = subsref(self, s)
            % Allow accessing sn property with dot notation
            if strcmp(s(1).type, '.')
                switch s(1).subs
                    case 'sn'
                        varargout{1} = self.get_sn();
                    otherwise
                        % Default behavior
                        varargout{1} = builtin('subsref', self, s);
                end
            else
                varargout{1} = builtin('subsref', self, s);
            end
        end
    end
end
