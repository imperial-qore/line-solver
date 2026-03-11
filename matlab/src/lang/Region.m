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
        constraintA;  % Linear constraint matrix A (C x K), An <= b
        constraintB;  % Linear constraint vector b (C x 1)
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
            self.constraintA = [];
            self.constraintB = [];
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

        function self = setLinearConstraints(self, A, b)
            % SELF = SETLINEARCONSTRAINTS(A, B)
            % Set general linear admission constraints An <= b.
            %   A: Constraint matrix (C x K) where C is the number of
            %      constraints and K is the number of classes.
            %   b: Capacity vector (C x 1) or (1 x C).
            self.constraintA = A;
            self.constraintB = b(:);  % Ensure column vector
        end

        function tf = hasLinearConstraints(self)
            % TF = HASLINEARCONSTRAINTS()
            % Returns true if linear constraints have been set.
            tf = ~isempty(self.constraintA) && ~isempty(self.constraintB);
        end
    end

end
