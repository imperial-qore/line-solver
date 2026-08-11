function sens = computeModelSensitivities(model, useCTMC)
% computeModelSensitivities  Analytic d(metric)/d(rate) for a model as an
% opt.SensitivityData, mirroring native-Python compute_model_sensitivities.
% Returns [] when unavailable (the optimizer then falls back to finite
% differences).
%
% The open and closed product-form branches are implemented in
% opt.sens.openSensitivities and opt.sens.closedSensitivities. Both are fast
% but narrow: they return [] for open, multiserver, or non-unit-visit
% networks. When they do, and USECTMC is true, the generator-derivative
% fallback opt.sens.ctmcSensitivities is tried instead, which is exact but
% pays for state-space generation.
%
% @param model Network model
% @param useCTMC Enable the CTMC fallback, default false
% @return sens opt.SensitivityData or []

if nargin < 2 || isempty(useCTMC)
    useCTMC = false;
end

sens = [];
try
    s = opt.sens.openSensitivities(model);
    if ~isempty(s)
        sens = s;
        return;
    end
    sens = opt.sens.closedSensitivities(model);
catch
    sens = [];
end

if isempty(sens) && useCTMC
    try
        sens = opt.sens.ctmcSensitivities(model);
    catch
        sens = [];
    end
end
end
