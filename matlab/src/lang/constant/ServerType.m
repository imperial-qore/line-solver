classdef ServerType < Element
    % ServerType Represents a type of server within a heterogeneous multiserver queue
    %
    % ServerType defines a group of identical servers with:
    % - A unique name identifying this server type
    % - A count of servers of this type
    % - A list of job classes that are compatible with (can be served by) this type
    %
    % Server types enable modeling of heterogeneous multiserver queues where different
    % servers may have different service rates and serve different subsets of job classes.
    %
    % Example:
    % @code
    % fastServer = ServerType('Fast', 2);
    % fastServer.setCompatible([classA, classB]);
    % queue.addServerType(fastServer);
    % queue.setService(classA, fastServer, Exp(2.0));
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        id;                    % Unique identifier within the queue
        numOfServers;          % Number of servers of this type
        compatibleClasses;     % Cell array of compatible JobClass objects
        parentQueue;           % The Queue this server type belongs to
    end

    methods
        function self = ServerType(name, numOfServers, compatibleClasses)
            % SERVERTYPE Create a new server type
            %
            % self = SERVERTYPE(name, numOfServers) creates a server type with
            % the specified name and number of servers. Compatible classes can
            % be added later using setCompatible() or addCompatible().
            %
            % self = SERVERTYPE(name, numOfServers, compatibleClasses) creates
            % a server type with initial compatible classes.
            %
            % @param name String name identifying this server type
            % @param numOfServers Number of servers of this type (must be >= 1)
            % @param compatibleClasses (optional) Cell array or array of JobClass objects

            self@Element(name);

            if numOfServers < 1
                line_error(mfilename, 'Number of servers must be at least 1');
            end

            self.id = -1;  % Will be set when added to a queue
            self.numOfServers = numOfServers;
            self.parentQueue = [];

            if nargin < 3
                self.compatibleClasses = {};
            else
                if iscell(compatibleClasses)
                    self.compatibleClasses = compatibleClasses;
                else
                    % Convert array to cell array
                    self.compatibleClasses = cell(1, length(compatibleClasses));
                    for i = 1:length(compatibleClasses)
                        self.compatibleClasses{i} = compatibleClasses(i);
                    end
                end
            end
        end

        function id = getId(self)
            % GETID Get the server type ID
            %
            % id = GETID() returns the unique identifier of this server type
            % within its queue, or -1 if not yet added to a queue.

            id = self.id;
        end

        function self = setId(self, id)
            % SETID Set the server type ID
            %
            % self = SETID(id) sets the unique identifier. This is typically
            % called by the Queue when the server type is added.

            self.id = id;
        end

        function n = getNumOfServers(self)
            % GETNUMOFSERVERS Get the number of servers of this type
            %
            % n = GETNUMOFSERVERS() returns the number of servers.

            n = self.numOfServers;
        end

        function self = setNumOfServers(self, n)
            % SETNUMOFSERVERS Set the number of servers of this type
            %
            % self = SETNUMOFSERVERS(n) sets the number of servers (must be >= 1).

            if n < 1
                line_error(mfilename, 'Number of servers must be at least 1');
            end
            self.numOfServers = n;
        end

        function classes = getCompatibleClasses(self)
            % GETCOMPATIBLECLASSES Get the list of compatible job classes
            %
            % classes = GETCOMPATIBLECLASSES() returns a cell array of
            % compatible JobClass objects.

            classes = self.compatibleClasses;
        end

        function self = setCompatible(self, classes)
            % SETCOMPATIBLE Set the list of compatible job classes
            %
            % self = SETCOMPATIBLE(classes) sets the list of job classes that
            % can be served by this server type.
            %
            % @param classes Array or cell array of JobClass objects

            if iscell(classes)
                self.compatibleClasses = classes;
            else
                % Convert array to cell array
                self.compatibleClasses = cell(1, length(classes));
                for i = 1:length(classes)
                    self.compatibleClasses{i} = classes(i);
                end
            end
        end

        function self = addCompatible(self, jobClass)
            % ADDCOMPATIBLE Add a job class to the compatible list
            %
            % self = ADDCOMPATIBLE(jobClass) adds a job class to the list of
            % classes that can be served by this server type.
            %
            % @param jobClass The JobClass to add

            if ~self.isCompatible(jobClass)
                self.compatibleClasses{end+1} = jobClass;
            end
        end

        function self = removeCompatible(self, jobClass)
            % REMOVECOMPATIBLE Remove a job class from the compatible list
            %
            % self = REMOVECOMPATIBLE(jobClass) removes a job class from the
            % list of compatible classes.
            %
            % @param jobClass The JobClass to remove

            idx = [];
            for i = 1:length(self.compatibleClasses)
                if self.compatibleClasses{i} == jobClass
                    idx = i;
                    break;
                end
            end
            if ~isempty(idx)
                self.compatibleClasses(idx) = [];
            end
        end

        function result = isCompatible(self, jobClass)
            % ISCOMPATIBLE Check if a job class is compatible with this server type
            %
            % result = ISCOMPATIBLE(jobClass) returns true if the job class
            % can be served by this server type.
            %
            % @param jobClass The JobClass to check
            % @return result Boolean indicating compatibility

            result = false;
            for i = 1:length(self.compatibleClasses)
                if self.compatibleClasses{i} == jobClass
                    result = true;
                    return;
                end
            end
        end

        function n = getNumCompatibleClasses(self)
            % GETNUMCOMPATIBLECLASSES Get the number of compatible classes
            %
            % n = GETNUMCOMPATIBLECLASSES() returns the count of compatible classes.

            n = length(self.compatibleClasses);
        end

        function result = hasCompatibleClasses(self)
            % HASCOMPATIBLECLASSES Check if any compatible classes are defined
            %
            % result = HASCOMPATIBLECLASSES() returns true if at least one
            % compatible class is defined.

            result = ~isempty(self.compatibleClasses);
        end

        function queue = getParentQueue(self)
            % GETPARENTQUEUE Get the parent queue
            %
            % queue = GETPARENTQUEUE() returns the Queue this server type
            % belongs to, or [] if not yet added to a queue.

            queue = self.parentQueue;
        end

        function self = setParentQueue(self, queue)
            % SETPARENTQUEUE Set the parent queue
            %
            % self = SETPARENTQUEUE(queue) sets the parent queue. This is
            % typically called by the Queue when the server type is added.

            self.parentQueue = queue;
        end

        function summary(self)
            % SUMMARY Print a summary of this server type
            %
            % SUMMARY() displays information about this server type.

            line_printf('ServerType: <strong>%s</strong>', self.name);
            line_printf('  ID: %d', self.id);
            line_printf('  Number of servers: %d', self.numOfServers);
            line_printf('  Compatible classes: %d', length(self.compatibleClasses));
            for i = 1:length(self.compatibleClasses)
                line_printf('    - %s', self.compatibleClasses{i}.getName());
            end
        end
    end
end
