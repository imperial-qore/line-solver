classdef (Sealed) HeteroSchedPolicy
    % HeteroSchedPolicy Enumeration of scheduling policies for heterogeneous multiserver queues
    %
    % HeteroSchedPolicy defines constants for scheduling policies that determine
    % how jobs are assigned to server types in heterogeneous multiserver queues.
    % These policies are used when a job's class is compatible with multiple
    % server types.
    %
    % @brief Scheduling policies for heterogeneous server assignment
    %
    % Available policies:
    % - ORDER: Assign to first available compatible server type (in definition order)
    % - ALIS: Assign Longest Idle Server (round-robin with busy servers at back)
    % - ALFS: Assign Longest Free Server (fairness sorting by coverage)
    % - FAIRNESS: Fair distribution across compatible server types
    % - FSF: Fastest Server First (based on expected service time)
    % - RAIS: Random Available Idle Server
    %
    % Example:
    % @code
    % queue = Queue(model, 'HeteroQueue', SchedStrategy.FCFS);
    % st1 = ServerType('Fast', 2);
    % queue.addServerType(st1);
    % queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        % ORDER - First available compatible server type in definition order (default)
        ORDER = 0;

        % ALIS - Assign Longest Idle Server (round-robin)
        ALIS = 1;

        % ALFS - Assign Longest Free Server with fairness sorting
        ALFS = 2;

        % FAIRNESS - Fair distribution across compatible server types
        FAIRNESS = 3;

        % FSF - Fastest Server First (minimum expected service time)
        FSF = 4;

        % RAIS - Random Available Idle Server
        RAIS = 5;
    end

    methods (Static)
        function policy = fromText(text)
            % FROMTEXT Convert text to HeteroSchedPolicy constant
            %
            % policy = FROMTEXT(text) converts a string to the corresponding
            % HeteroSchedPolicy constant value.
            %
            % @param text String representation of the policy
            % @return policy The corresponding HeteroSchedPolicy constant

            switch upper(text)
                case 'ORDER'
                    policy = HeteroSchedPolicy.ORDER;
                case 'ALIS'
                    policy = HeteroSchedPolicy.ALIS;
                case 'ALFS'
                    policy = HeteroSchedPolicy.ALFS;
                case 'FAIRNESS'
                    policy = HeteroSchedPolicy.FAIRNESS;
                case 'FSF'
                    policy = HeteroSchedPolicy.FSF;
                case 'RAIS'
                    policy = HeteroSchedPolicy.RAIS;
                otherwise
                    line_error(mfilename, 'Unknown HeteroSchedPolicy: %s', text);
            end
        end

        function text = toText(policy)
            % TOTEXT Convert HeteroSchedPolicy constant to text
            %
            % text = TOTEXT(policy) converts a HeteroSchedPolicy constant
            % to its string representation.
            %
            % @param policy The HeteroSchedPolicy constant
            % @return text String representation of the policy

            switch policy
                case HeteroSchedPolicy.ORDER
                    text = 'ORDER';
                case HeteroSchedPolicy.ALIS
                    text = 'ALIS';
                case HeteroSchedPolicy.ALFS
                    text = 'ALFS';
                case HeteroSchedPolicy.FAIRNESS
                    text = 'FAIRNESS';
                case HeteroSchedPolicy.FSF
                    text = 'FSF';
                case HeteroSchedPolicy.RAIS
                    text = 'RAIS';
                otherwise
                    text = 'UNKNOWN';
            end
        end
    end
end
