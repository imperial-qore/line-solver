classdef RewardStateView < handle
    % REWARDSTATEVIEW View of reward state subset with aggregation operations
    %
    % Provides aggregation operations (sum, max, min, count) on a subset
    % of the state vector. Returned by RewardState.at() and forClass().
    %
    % Aggregation Methods:
    %   total()   - Sum all values in this view
    %   max()     - Maximum value in this view
    %   min()     - Minimum value in this view
    %   count()   - Number of non-zero entries in this view
    %
    % Filter Methods (to transition between views):
    %   atClass(jobclass)   - Filter station view to specific class
    %   atStation(station)  - Filter class view to specific station
    %
    % Usage:
    %   % Station view (all classes)
    %   view = state.at(queue1);
    %   total = view.total();                    % Total jobs at queue1
    %   jobs_class1 = view.atClass(class1);     % Class1 jobs at queue1
    %
    %   % Class view (all stations)
    %   view = state.forClass(class1);
    %   total = view.total();                    % Total class1 jobs
    %   jobs_queue1 = view.atStation(queue1);   % Class1 jobs at queue1
    %
    % See also: RewardState, Reward, setReward

    properties (Access = private)
        values              % Vector of state subset values [1 x N]
        parentState         % RewardState reference (for filtering)
        stationIdx          % Station index if this is a station view, [] otherwise
        classIdx            % Class index if this is a class view, [] otherwise
    end

    methods
        function self = RewardStateView(vals, parent, stIdx, clIdx)
            % REWARDSTATEVIEW Constructor
            %
            % SELF = REWARDSTATEVIEW(VALS, PARENT, STIDX, CLIDX)
            %
            % VALS    - Vector of values [1 x N]
            % PARENT  - RewardState parent object
            % STIDX   - Station index if station view, [] if class view
            % CLIDX   - Class index if class view, [] if station view

            self.values = vals;
            self.parentState = parent;
            self.stationIdx = stIdx;
            self.classIdx = clIdx;
        end

        function result = total(self)
            % TOTAL Sum all values in this view
            %
            % RESULT = TOTAL(SELF) returns the sum of all state values
            %
            % Example:
            %   view = state.at(queue1);
            %   total_jobs = view.total();

            result = sum(self.values);
        end

        function result = max(self)
            % MAX Maximum value in this view
            %
            % RESULT = MAX(SELF) returns the maximum state value
            %
            % Example:
            %   view = state.at(queue1);
            %   max_class = max(view);

            if isempty(self.values)
                result = 0;
            else
                result = max(self.values);
            end
        end

        function result = min(self)
            % MIN Minimum value in this view
            %
            % RESULT = MIN(SELF) returns the minimum state value
            %
            % Example:
            %   view = state.at(queue1);
            %   min_class = min(view);

            if isempty(self.values)
                result = 0;
            else
                result = min(self.values);
            end
        end

        function result = count(self)
            % COUNT Number of non-zero entries
            %
            % RESULT = COUNT(SELF) returns the number of non-zero values
            %
            % Example:
            %   view = state.at(queue1);
            %   nonempty_classes = count(view);

            result = sum(self.values > 0);
        end

        function result = atClass(self, jobclass)
            % ATCLASS Filter station view to specific class
            %
            % RESULT = ATCLASS(SELF, JOBCLASS) returns the scalar population
            % of JOBCLASS at the station in this view.
            %
            % Only valid if this is a station view (created by state.at(station))
            %
            % Example:
            %   view = state.at(queue1);
            %   class1_at_queue1 = view.atClass(class1);

            if isempty(self.stationIdx)
                error('RewardStateView:InvalidOperation', ...
                    'atClass() is only valid for station views (created by state.at(station))');
            end

            % Return scalar for this station and class
            result = self.parentState.at(self.stationIdx, jobclass);
        end

        function result = atStation(self, station)
            % ATSTATION Filter class view to specific station
            %
            % RESULT = ATSTATION(SELF, STATION) returns the scalar population
            % of the class in this view at the specified STATION.
            %
            % Only valid if this is a class view (created by state.forClass(class))
            %
            % Example:
            %   view = state.forClass(class1);
            %   class1_at_queue1 = view.atStation(queue1);

            if isempty(self.classIdx)
                error('RewardStateView:InvalidOperation', ...
                    'atStation() is only valid for class views (created by state.forClass(class))');
            end

            % Return scalar for this class and station
            result = self.parentState.at(station, self.classIdx);
        end
    end

    % Allow arithmetic operations on views
    methods
        function result = plus(self, other)
            % PLUS Add scalar or compatible view
            if isa(other, 'RewardStateView')
                result = sum(self.values) + sum(other.values);
            else
                result = sum(self.values) + other;
            end
        end

        function result = minus(self, other)
            % MINUS Subtract scalar or compatible view
            if isa(other, 'RewardStateView')
                result = sum(self.values) - sum(other.values);
            else
                result = sum(self.values) - other;
            end
        end

        function result = times(self, other)
            % TIMES Multiply view values by scalar
            if isa(other, 'RewardStateView')
                result = sum(self.values .* other.values);
            else
                result = self.values * other;
            end
        end

        function result = mtimes(self, other)
            % MTIMES Matrix multiply (scalar multiplication)
            if isscalar(other)
                result = sum(self.values) * other;
            else
                error('RewardStateView:InvalidOperation', ...
                    'Matrix multiplication not supported');
            end
        end

        function result = ge(self, threshold)
            % GE Greater than or equal (for conditional rewards)
            result = double(sum(self.values) >= threshold);
        end

        function result = le(self, threshold)
            % LE Less than or equal (for conditional rewards)
            result = double(sum(self.values) <= threshold);
        end

        function result = gt(self, threshold)
            % GT Greater than (for conditional rewards)
            result = double(sum(self.values) > threshold);
        end

        function result = lt(self, threshold)
            % LT Less than (for conditional rewards)
            result = double(sum(self.values) < threshold);
        end
    end

    % Display
    methods
        function disp(self)
            % DISP Display RewardStateView
            fprintf('RewardStateView with values: [');
            fprintf('%g ', self.values);
            fprintf(']\n');
            fprintf('  Total: %g\n', self.total());
        end
    end
end
