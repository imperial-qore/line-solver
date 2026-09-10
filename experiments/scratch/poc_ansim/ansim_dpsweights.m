function w = ansim_dpsweights(cfg)
% ANSIM_DPSWEIGHTS  DPS weights that emulate the priority order,
% w_r = W^(maxprio - prio_r); with W large the top group takes almost the whole
% server. This is the only preemption surrogate SolverFLD and SolverNC can take,
% and its residual error is the dps-gap floor.
w = cfg.Wdps .^ (max(cfg.prio) - cfg.prio);
end
