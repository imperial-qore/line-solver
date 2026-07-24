classdef RoutingProbabilities < opt.DecisionVariable
    % RoutingProbabilities  Optimize routing of a job class from a source node
    % to target nodes. Stick-breaking encoding (dim = targets-1). On apply,
    % default routing is rebuilt from the connection matrix for every class,
    % then the overridden (class, source) row is set to the decoded
    % probabilities. Mirrors native-Python RoutingProbabilities.

    properties
        jobclass
        source
        targets      % cell array of Node
    end

    methods
        function obj = RoutingProbabilities(jobclass, source, targets, name)
            if nargin < 4 || isempty(name)
                name = [jobclass.getName() '_routing_from_' source.getName()];
            end
            obj@opt.DecisionVariable(name);
            obj.jobclass = jobclass;
            obj.source = source;
            obj.targets = targets;
            obj.dimension = max(1, numel(targets) - 1);
        end

        function c = getJobClass(obj), c = obj.jobclass; end
        function s = getSource(obj), s = obj.source; end
        function t = getTargets(obj), t = obj.targets; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(obj.dimension); end

        function probs = decode(obj, x)
            n = numel(obj.targets);
            if n == 1
                probs = 1.0; return;
            end
            probs = zeros(1, n);
            remaining = 1.0;
            for i = 1:(n-1)
                probs(i) = remaining * x(i);
                remaining = remaining - probs(i);
            end
            probs(n) = remaining;
        end

        function apply(obj, model, value)
            % Override only the source node's outgoing routing for this class,
            % via setProbRouting, leaving all other routes intact. This is the
            % MATLAB idiom that works whether the model was built with link()
            % or addLink() (model.link() is rejected after addLink()).
            cls = opt.DecisionVariable.resolveClass(model, obj.jobclass);
            if isempty(cls), return; end
            src = opt.DecisionVariable.resolveNode(model, obj.source);
            if isempty(src), return; end
            for i = 1:numel(obj.targets)
                target = opt.DecisionVariable.resolveNode(model, obj.targets{i});
                if ~isempty(target)
                    src.setProbRouting(cls, target, value(i));
                end
            end
        end

        function t = getVariableType(obj), t = 'routing'; end
    end
end
