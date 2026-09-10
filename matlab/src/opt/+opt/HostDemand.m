classdef HostDemand < opt.DecisionVariable
    % HostDemand  Optimize the mean host demand D of an LQN Activity
    % (continuous). The processor-layer service rate is mu = 1/D; this is the
    % primary LQN tuning knob (analogous to ServiceRate for a flat station).
    % Exposes the partial-sensitivity gradient hooks
    % (sensKey/sensMetricTargets/rateJacobian/decodeJacobian). Mirrors
    % native-Python HostDemand.

    properties
        activity    % activity name (char)
        minDemand
        maxDemand
    end

    methods
        function obj = HostDemand(activity, bounds, name)
            actName = opt.Layered.elemName(activity);
            if nargin < 3 || isempty(name)
                name = [actName '_hostdemand'];
            end
            obj@opt.DecisionVariable(name);
            obj.activity = actName;
            obj.minDemand = bounds(1);
            obj.maxDemand = bounds(2);
        end

        function a = getActivity(obj), a = obj.activity; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            v = obj.minDemand + x(1) * (obj.maxDemand - obj.minDemand);
        end

        function apply(obj, model, value)
            act = opt.Layered.resolveActivity(model, obj.activity);
            if ~isempty(act)
                act.setHostDemand(double(value));
            end
        end

        function t = getVariableType(obj), t = 'host_demand'; end

        function layers = getLayer(obj, model)
            proc = opt.Layered.activityProcessorName(model, obj.activity);
            if isempty(proc), layers = {}; else, layers = {proc}; end
        end

        function v = currentValue(obj, model)
            v = [];
            act = opt.Layered.resolveActivity(model, obj.activity);
            if ~isempty(act)
                v = act.getHostDemandMean();
            end
        end

        % ---- partial-sensitivity gradient hooks ---------------------------

        function k = sensKey(obj, model)
            % Row key 'Station||JobClass' in the per-layer sensitivity table.
            % Host-layer rows are (Layer=processor, Station=processor,
            % JobClass=activity); the value is d(metric)/d(service rate).
            k = '';
            proc = opt.Layered.activityProcessorName(model, obj.activity);
            if isempty(proc), return; end
            k = [proc '||' obj.activity];
        end

        function targets = sensMetricTargets(obj, model)
            % Map each layer-row metric to the EvaluationResult key it
            % approximates: the host-layer row utilization tracks the processor
            % node's utilization (keyed by node name); its throughput/queue-
            % length/response-time track the activity node's (keyed 'act||act').
            proc = opt.Layered.activityProcessorName(model, obj.activity);
            act = obj.activity;
            targets = struct('Util', proc, 'Tput', [act '||' act], ...
                'QLen', [act '||' act], 'RespT', [act '||' act]);
        end

        function j = rateJacobian(~, value)
            % d(service rate)/d(demand) = d(1/D)/dD = -1/D^2 at D=value.
            v = double(value);
            if v <= 0, j = 0.0; else, j = -1.0 / (v * v); end
        end

        function j = decodeJacobian(obj, ~)
            j = obj.maxDemand - obj.minDemand;
        end
    end
end
