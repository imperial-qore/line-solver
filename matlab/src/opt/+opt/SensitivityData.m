classdef SensitivityData < handle
    % SensitivityData  Analytic performance sensitivities for a product-form
    % model, mirroring native-Python compute_model_sensitivities: metric kind
    % ('RespT'|'QLen'|'Tput'|'Util') -> metric key -> parameter key ->
    % d(metric)/d(parameter). Keys are canonical strings: a metric key is
    % 'station' (Util) or 'station||class'; a parameter key is
    % 'rate||station||class'. Backed by nested dictionaries.

    properties
        data     % dict kind -> (dict metricKey -> (dict paramKey -> value))
    end

    methods
        function obj = SensitivityData()
            obj.data = configureDictionary('string', 'cell');
        end

        function add(obj, kind, metricKey, paramKey, value)
            % dictionary is a value type, so each nested level is written back
            if isKey(obj.data, kind)
                byMetric = obj.data{kind};
            else
                byMetric = configureDictionary('string', 'cell');
            end
            if isKey(byMetric, metricKey)
                byParam = byMetric{metricKey};
            else
                byParam = configureDictionary('string', 'double');
            end
            if isKey(byParam, paramKey)
                byParam(paramKey) = byParam(paramKey) + value;
            else
                byParam(paramKey) = value;
            end
            byMetric{metricKey} = byParam;
            obj.data{kind} = byMetric;
        end

        function m = forKind(obj, kind)
            if isKey(obj.data, kind), m = obj.data{kind}; else, m = []; end
        end

        function tf = isempty(obj)
            tf = (numEntries(obj.data) == 0);
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
