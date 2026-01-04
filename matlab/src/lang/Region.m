classdef Region < handle
    % A finite capacity region
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        UNBOUNDED = -1;
    end

    properties
        name;
        classes;
        nodes;
        globalMaxJobs;
        globalMaxMemory;
        classMaxJobs;
        classMaxMemory;
        dropRule;
        classWeight;
        classSize;
    end

    methods
        function self = Region(nodes,classes)
            self.name = '';
            self.nodes = {nodes{:}}';
            self.classes = classes;
            self.globalMaxJobs = Region.UNBOUNDED;
            self.globalMaxMemory = Region.UNBOUNDED;
            self.classWeight = ones(1,length(self.classes));
            self.dropRule = DropStrategy.WAITQ * ones(1,length(self.classes));
            self.classSize = ones(1,length(self.classes));
            self.classMaxJobs = Region.UNBOUNDED * ones(1,length(self.classes));
            self.classMaxMemory = Region.UNBOUNDED * ones(1,length(self.classes));
        end

        function self = setGlobalMaxJobs(self, njobs)
            self.globalMaxJobs = njobs;
        end

        function self = setGlobalMaxMemory(self, memlim)
            self.globalMaxMemory = memlim;
        end

        function self = setClassMaxJobs(self, class, njobs)
            self.classMaxJobs(class.index) = njobs;
        end

        function self = setClassWeight(self, class, weight)
            self.classWeight(class.index) = weight;
        end

        function self = setDropRule(self, class, dropStrategy)
            % SELF = SETDROPRULE(CLASS, DROPSTRATEGY)
            % Set the drop rule for a class.
            % dropStrategy can be:
            %   - A boolean: true = DROP, false = WAITQ (for backwards compatibility)
            %   - A DropStrategy enum value: DROP, WAITQ, BAS, BBS, RSRD
            if islogical(dropStrategy)
                if dropStrategy
                    self.dropRule(class.index) = DropStrategy.DROP;
                else
                    self.dropRule(class.index) = DropStrategy.WAITQ;
                end
            else
                self.dropRule(class.index) = dropStrategy;
            end
        end

        function strategy = getDropRule(self, class)
            % STRATEGY = GETDROPRULE(CLASS)
            % Get the drop strategy for a class.
            strategy = self.dropRule(class.index);
        end

        function self = setClassSize(self, class, size)
            self.classSize(class.index) = size;
        end

        function self = setClassMaxMemory(self, class, memlim)
            self.classMaxMemory(class.index) = memlim;
        end

        function self = setName(self, name)
            self.name = name;
        end

        function name = getName(self)
            name = self.name;
        end
    end

end
