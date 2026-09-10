classdef ProcessorMultiplicity < opt.DecisionVariable
    % ProcessorMultiplicity  Optimize the multiplicity (core count) of an LQN
    % Processor (integer). Mirrors native-Python ProcessorMultiplicity.

    properties
        processor
        minValue
        maxValue
    end

    methods
        function obj = ProcessorMultiplicity(processor, bounds, name)
            procName = opt.Layered.elemName(processor);
            if nargin < 3 || isempty(name)
                name = [procName '_multiplicity'];
            end
            obj@opt.DecisionVariable(name);
            obj.processor = procName;
            obj.minValue = round(bounds(1));
            obj.maxValue = round(bounds(2));
        end

        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            continuous = obj.minValue + x(1) * (obj.maxValue - obj.minValue);
            v = min(max(round(continuous), obj.minValue), obj.maxValue);
        end

        function apply(obj, model, value)
            proc = opt.Layered.resolveProcessor(model, obj.processor);
            if ~isempty(proc)
                proc.multiplicity = round(value);
            end
        end

        function t = getVariableType(obj), t = 'processor_multiplicity'; end

        function layers = getLayer(obj, ~)
            % A processor owns its own host layer, named after the processor.
            layers = {obj.processor};
        end

        function v = currentValue(obj, model)
            v = [];
            proc = opt.Layered.resolveProcessor(model, obj.processor);
            if isempty(proc), return; end
            m = proc.multiplicity;
            if isnumeric(m) && isfinite(m), v = round(m); end
        end
    end
end
