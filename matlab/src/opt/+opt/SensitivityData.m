classdef SensitivityData < handle
    % SensitivityData  Analytic performance sensitivities for a product-form
    % model, mirroring native-Python compute_model_sensitivities: metric kind
    % ('RespT'|'QLen'|'Tput'|'Util') -> metric key -> parameter key ->
    % d(metric)/d(parameter). Keys are canonical strings: a metric key is
    % 'station' (Util) or 'station||class'; a parameter key is
    % 'rate||station||class'. Backed by nested containers.Map.

    properties
        data     % Map kind -> (Map metricKey -> (Map paramKey -> value))
    end

    methods
        function obj = SensitivityData()
            obj.data = containers.Map('KeyType', 'char', 'ValueType', 'any');
        end

        function add(obj, kind, metricKey, paramKey, value)
            if ~isKey(obj.data, kind)
                obj.data(kind) = containers.Map('KeyType', 'char', 'ValueType', 'any');
            end
            byMetric = obj.data(kind);
            if ~isKey(byMetric, metricKey)
                byMetric(metricKey) = containers.Map('KeyType', 'char', 'ValueType', 'double'); %#ok<NASGU>
            end
            byParam = byMetric(metricKey);
            if isKey(byParam, paramKey)
                byParam(paramKey) = byParam(paramKey) + value;
            else
                byParam(paramKey) = value;
            end
        end

        function m = forKind(obj, kind)
            if isKey(obj.data, kind), m = obj.data(kind); else, m = []; end
        end

        function tf = isempty(obj)
            tf = (obj.data.Count == 0);
        end
    end

    methods (Static)
        function k = metricKey(station, jobclass)
            k = [station '||' jobclass];
        end
        function k = paramKey(station, jobclass)
            k = ['rate||' station '||' jobclass];
        end
    end
end
