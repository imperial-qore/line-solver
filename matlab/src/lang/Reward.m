classdef Reward
    % REWARD Factory class for common reward function templates
    %
    % This class provides static methods for creating common reward functions.
    % Used with model.setReward() to define metrics on queueing networks.
    %
    % Static Methods:
    %   queueLength(node, [jobclass]) - Queue length reward
    %   utilization(node, [jobclass]) - Server utilization reward
    %   blocking(node)                 - Blocking probability reward
    %   custom(fn)                     - Wrap a custom function
    %
    % Usage:
    %   model.setReward('QLen_Q1', Reward.queueLength(queue1));
    %   model.setReward('Util_Q1', Reward.utilization(queue1));
    %   model.setReward('Block_Q1', Reward.blocking(queue1));
    %   model.setReward('Cost', Reward.custom(@(state) state.at(q1).total()^2));
    %
    % See also: RewardState, RewardStateView, setReward

    methods (Static)
        function fn = queueLength(node, varargin)
            % QUEUELENGTH Queue length reward function
            %
            % FN = QUEUELENGTH(NODE) returns function for total jobs at NODE
            % FN = QUEUELENGTH(NODE, JOBCLASS) returns function for JOBCLASS jobs
            %
            % The returned function is suitable for use with model.setReward()
            %
            % Examples:
            %   % Total jobs at queue1
            %   model.setReward('QLen', Reward.queueLength(queue1));
            %
            %   % Class1 jobs at queue1
            %   model.setReward('QLen_C1', Reward.queueLength(queue1, class1));

            if nargin == 1
                % Return total jobs at node (all classes)
                fn = RewardDescriptor('QLen', node, [], @(state) state.at(node).total());
            else
                % Return jobs of specific class at node
                jobclass = varargin{1};
                fn = RewardDescriptor('QLen', node, jobclass, @(state) state.at(node, jobclass));
            end
        end

        function fn = utilization(node, varargin)
            % UTILIZATION Server utilization reward function
            %
            % FN = UTILIZATION(NODE) returns function for utilization at NODE
            % FN = UTILIZATION(NODE, JOBCLASS) returns function for class utilization
            %
            % Utilization is computed as min(jobs, nservers), representing the
            % fraction of servers in use.
            %
            % Note: For M/M/1 queues, this simplifies to min(jobs, 1)
            %
            % Examples:
            %   model.setReward('Util', Reward.utilization(queue1));
            %   model.setReward('Util_C1', Reward.utilization(queue1, class1));

            if nargin == 1
                % Total utilization at node (all classes)
                fn = RewardDescriptor('Util', node, [], @(state) compute_utilization(state, node, []));
            else
                % Utilization of specific class at node
                jobclass = varargin{1};
                fn = RewardDescriptor('Util', node, jobclass, @(state) compute_utilization(state, node, jobclass));
            end

            function util = compute_utilization(state, nd, jc)
                % Helper function for utilization computation
                sn = state.sn;
                stationIdx = sn.nodeToStation(nd.index);
                if isempty(jc)
                    % Total jobs at station
                    jobs = state.at(nd).total();
                else
                    % Jobs of specific class
                    jobs = state.at(nd, jc);
                end
                nservers = sn.nservers(stationIdx);
                util = min(jobs, nservers);
            end
        end

        function fn = blocking(node)
            % BLOCKING Blocking probability reward function
            %
            % FN = BLOCKING(NODE) returns 1 if NODE is at capacity, 0 otherwise
            %
            % This is useful for measuring congestion or capacity violations.
            %
            % Example:
            %   model.setReward('Block', Reward.blocking(queue1));

            fn = RewardDescriptor('Blocking', node, [], @(state) compute_blocking(state, node));

            function block = compute_blocking(state, nd)
                % Helper function for blocking computation
                sn = state.sn;
                stationIdx = sn.nodeToStation(nd.index);
                jobs = state.at(nd).total();
                capacity = sn.cap(stationIdx);
                block = double(jobs >= capacity);
            end
        end

        function fn = custom(userFn)
            % CUSTOM Wrap a custom reward function
            %
            % FN = CUSTOM(USERFN) wraps USERFN as a custom reward descriptor
            %
            % The returned value is callable exactly like USERFN. It is marked
            % as kind 'Custom', which is deliberately NOT serializable: an
            % arbitrary user function cannot be reproduced from JSON, so
            % linemodel_save warns and omits it rather than emitting a wrong
            % reward.
            %
            % Example:
            %   myReward = @(state) state.at(q1).total()^2 + state.at(q2).total();
            %   model.setReward('Custom', Reward.custom(myReward));

            fn = RewardDescriptor('Custom', [], [], userFn);
        end
    end
end
