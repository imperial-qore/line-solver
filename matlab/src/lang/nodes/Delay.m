classdef Delay < Queue
    % Delay Infinite server station with no queueing delays
    %
    % Delay is a specialized Queue with infinite servers, meaning jobs
    % experience service time but never wait in queue. Each arriving job
    % is immediately served, making this ideal for modeling think times,
    % processing delays, and other non-competing service processes.
    %
    % @brief Infinite server station for modeling non-queueing delays
    %
    % Key characteristics:
    % - Infinite number of servers (no queueing)
    % - Jobs experience service time but no waiting
    % - Suitable for modeling think times and processing delays
    % - Inherits all Queue functionality with modified scheduling
    % - Equivalent to M/M/∞ queueing system
    %
    % Delay stations are commonly used for:
    % - User think time modeling
    % - Network propagation delays  
    % - Processing time without resource contention
    % - Timeout and waiting periods
    % - Background processing tasks
    %
    % Example:
    % @code
    % delay = Delay(model, 'ThinkTime');
    % delay.setService(jobClass, Exp(1.0));  % Exponential service time
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Delay(model, name)
            % DELAY Create a Delay station instance
            %
            % @brief Creates a Delay station with infinite servers (no queueing)
            % @param model Network model to add the delay station to
            % @param name String identifier for the delay station
            % @return self Delay instance configured with infinite servers
            %
            % The constructor creates a delay station by calling the parent Queue
            % constructor with infinite scheduling and sets numberOfServers to Inf.
            
            self@Queue(model, name, SchedStrategy.INF);
            self.numberOfServers = Inf;
        end
    end
    
end

