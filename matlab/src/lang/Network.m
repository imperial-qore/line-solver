classdef Network < MNetwork
    % Main queueing network model class for LINE analysis
    %
    % Provides methods for adding nodes, job classes, and links to create queueing networks.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    % PUBLIC METHODS
    methods (Access=public)
        %Constructor
        function self = Network(name, varargin)
            % NETWORK Create a new queueing network model
            %
            % @brief Creates a Network instance for queueing model construction
            % @param name String identifier for the network model
            % @param varargin Optional implementation parameter (ignored for performance)
            % @return self Network instance ready for model construction
            %
            % For compatibility, accepts but ignores the implementation argument.
            % Always uses MNetwork (MATLAB) implementation for optimal performance.
            
            self@MNetwork(name); % Always use MNetwork for performance
            
            % Parse optional implementation argument for compatibility
            % but ignore it - always use MNetwork
            if nargin >= 2 && ischar(varargin{1})
                implementation = lower(varargin{1});
                if strcmp(implementation, 'java') || strcmp(implementation, 'j') || strcmp(implementation, 'jnetwork')
                    warning('Network:JavaNotSupported', ...
                        'Java implementation requested but not supported in performance-optimized mode. Using MATLAB implementation.');
                end
            end
        end
    end
    
    % STATIC METHODS
    methods (Static)
        function model = cyclic(N, D, strategy, S)
            % MODEL = CYCLIC(N, D, STRATEGY, S)
            %
            % Generates a cyclic queueing network
            model = MNetwork.cyclic(N, D, strategy, S);
        end
        
        function model = tandem(lambda, D, strategy)
            % MODEL = TANDEM(LAMBDA, D, STRATEGY)
            %
            % Generates a tandem queueing network
            model = MNetwork.tandem(lambda, D, strategy);
        end

        function model = cluster(lambda, D, strategy, S, dispatching)
            % MODEL = SERVERFARM(LAMBDA, D, STRATEGY, S, DISPATCHING)
            %
            % Generates an open server-farm queueing network:
            % Source -> Dispatcher (Router) -> Server[1..M] -> Sink
            if nargin < 4 || isempty(S), S = ones(size(D,1),1); end
            if nargin < 5, dispatching = RoutingStrategy.RAND; end
            model = MNetwork.cluster(lambda, D, strategy, S, dispatching);
        end

        function model = clusterPs(lambda, D, dispatching)
            % MODEL = SERVERFARMPS(LAMBDA, D, DISPATCHING)
            %
            % Open PS cluster with one server per queue
            if nargin < 3, dispatching = RoutingStrategy.RAND; end
            M = size(D,1);
            strategy = cell(M,1);
            for i = 1:M, strategy{i} = SchedStrategy.PS; end
            model = MNetwork.cluster(lambda, D, strategy, ones(M,1), dispatching);
        end

        function model = clusterFcfs(lambda, D, S, dispatching)
            % MODEL = SERVERFARMFCFS(LAMBDA, D, S, DISPATCHING)
            %
            % Open FCFS cluster
            if nargin < 3 || isempty(S), S = ones(size(D,1),1); end
            if nargin < 4, dispatching = RoutingStrategy.RAND; end
            M = size(D,1);
            strategy = cell(M,1);
            for i = 1:M, strategy{i} = SchedStrategy.FCFS; end
            model = MNetwork.cluster(lambda, D, strategy, S, dispatching);
        end

        function model = clusterClosed(N, Z, D, strategy, S, dispatching)
            % MODEL = SERVERFARMCLOSED(N, Z, D, STRATEGY, S, DISPATCHING)
            %
            % Generates a closed server-farm queueing network:
            % Think (Delay) -> Dispatcher (Router) -> Server[1..M] -> Think
            if nargin < 5 || isempty(S), S = ones(size(D,1),1); end
            if nargin < 6, dispatching = RoutingStrategy.RAND; end
            model = MNetwork.clusterClosed(N, Z, D, strategy, S, dispatching);
        end
    end
end