classdef FiniteCapacityRegion < handle
    % A finite capacity region
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties (Constant)
        UNBOUNDED = -1;
    end

    properties
        classes;
        nodes;
        globalMaxJobs;
        globalMaxMemory;
        classMaxJobs;
        classMaxMemory;
        dropRule;
        %classWeight;
        classSize;
    end

    methods
        function self = FiniteCapacityRegion(nodes,classes)
            self.nodes = {nodes{:}}';
            self.classes = classes;
            self.globalMaxJobs = FiniteCapacityRegion.UNBOUNDED;
            self.globalMaxMemory = FiniteCapacityRegion.UNBOUNDED;
            %self.classWeight = ones(1,length(self.classes));
            self.dropRule = false(1,length(self.classes));
            self.classSize = ones(1,length(self.classes));
            self.classMaxJobs = FiniteCapacityRegion.UNBOUNDED * ones(1,length(self.classes));
            self.classMaxMemory = FiniteCapacityRegion.UNBOUNDED * ones(1,length(self.classes));
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

        function self = setDropRule(self, class, isDropEnabled)
            self.dropRule(class.index) = isDropEnabled;
        end

        function self = setClassSize(self, class, size)
            self.classSize(class.index) = size;
        end

        function self = setClassMaxMemory(self, class, memlim)
            self.classMaxMemory(class.index) = memlim;
        end
    end

end
