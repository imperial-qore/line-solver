classdef ClassPriority < opt.DecisionVariable
    % ClassPriority  Optimize the priority of job classes. In 'levels' mode each
    % class gets an integer priority in [minPriority, maxPriority] (one dim per
    % class); in 'permutation' mode encoded keys induce a priority ordering
    % (n-1 dims). Mirrors native-Python ClassPriority.

    properties
        jobclasses    % cell array of JobClass
        mode
        minPriority
        maxPriority
    end

    methods
        function obj = ClassPriority(jobclasses, mode, priorityRange, name)
            if nargin < 2 || isempty(mode), mode = 'levels'; end
            if nargin < 3 || isempty(priorityRange), priorityRange = [1 10]; end
            if nargin < 4 || isempty(name), name = 'class_priorities'; end
            obj@opt.DecisionVariable(name);
            obj.jobclasses = jobclasses;
            obj.mode = mode;
            obj.minPriority = priorityRange(1);
            obj.maxPriority = priorityRange(2);
            if strcmp(mode, 'levels')
                obj.dimension = numel(jobclasses);
            else
                obj.dimension = max(1, numel(jobclasses) - 1);
            end
        end

        function c = getJobClasses(obj), c = obj.jobclasses; end
        function m = getMode(obj), m = obj.mode; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(obj.dimension); end

        function v = decode(obj, x)
            if strcmp(obj.mode, 'levels')
                v = zeros(1, numel(x));
                for i = 1:numel(x)
                    p = obj.minPriority + x(i) * (obj.maxPriority - obj.minPriority);
                    v(i) = round(p);
                end
            else
                n = numel(obj.jobclasses);
                if n == 1
                    v = 0; return;
                end
                keys = [x(:).' 0.0];
                % descending order (argsort of -keys), stable ties by index
                [~, order] = sortrows([-keys(:), (1:n).']);
                v = (order.' - 1);   % 0-based class indices
            end
        end

        function apply(obj, model, value)
            if strcmp(obj.mode, 'levels')
                for i = 1:numel(obj.jobclasses)
                    cls = obj.resolveOrSelf(model, obj.jobclasses{i});
                    cls.setPriority(value(i));
                end
            else
                m = numel(value);
                for rank = 1:m
                    classIdx = value(rank) + 1;   % 0-based -> 1-based
                    priority = m - (rank - 1);
                    cls = obj.resolveOrSelf(model, obj.jobclasses{classIdx});
                    cls.setPriority(priority);
                end
            end
        end

        function t = getVariableType(obj), t = 'class_priority'; end
    end

    methods (Access = private)
        function c = resolveOrSelf(~, model, jc)
            c = opt.DecisionVariable.resolveClass(model, jc);
            if isempty(c), c = jc; end
        end
    end
end
