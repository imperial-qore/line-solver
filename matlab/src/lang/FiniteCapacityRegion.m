classdef FiniteCapacityRegion < handle
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
        constraintA;        % Linear constraint matrix (C_f x K) or empty
        constraintB;        % Linear constraint capacity (C_f x 1) or empty
    end

    methods
        function self = FiniteCapacityRegion(nodes,classes)
            self.name = '';
            self.nodes = {nodes{:}}';
            self.classes = classes;
            self.globalMaxJobs = FiniteCapacityRegion.UNBOUNDED;
            self.globalMaxMemory = FiniteCapacityRegion.UNBOUNDED;
            self.classWeight = ones(1,length(self.classes));
            self.dropRule = DropStrategy.WAITQ * ones(1,length(self.classes));
            self.classSize = ones(1,length(self.classes));
            self.classMaxJobs = FiniteCapacityRegion.UNBOUNDED * ones(1,length(self.classes));
            self.classMaxMemory = FiniteCapacityRegion.UNBOUNDED * ones(1,length(self.classes));
            self.constraintA = [];
            self.constraintB = [];
        end

        function self = setConstraint(self, A, b)
            % SETCONSTRAINT(A, B) — Set the linear constraint A * x <= B
            % where x is the per-class job count vector. A must be (C x K),
            % B must be (C x 1) for some number of constraints C.
            A = double(A);
            b = double(b(:));
            if size(A, 1) ~= length(b)
                line_error(mfilename, 'A and b must have matching number of rows.');
            end
            if size(A, 2) ~= length(self.classes)
                line_error(mfilename, 'A column count must equal number of classes (%d).', length(self.classes));
            end
            self.constraintA = A;
            self.constraintB = b;
        end

        function [A, b] = getLinearConstraints(self)
            % [A, B] = GETLINEARCONSTRAINTS() — Linear admission constraint
            % pair (A,b), or empties if not set.
            A = self.constraintA;
            b = self.constraintB;
        end

        function tf = hasLinearConstraints(self)
            tf = ~isempty(self.constraintA) && ~isempty(self.constraintB);
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
