function assertPhaseTypeStates(self, what)
% ASSERTPHASETYPESTATES(WHAT)
% Refuse a query whose answer is a per-state probability under an ME.
%
% A matrix-exponential service embeds in the generator with negative
% off-diagonal entries, so the stationary vector is a SIGNED measure: only its
% aggregates over each phase block are probabilities. Mean measures stay exact
% (they are linear in that vector), but a per-state or transient answer is not a
% probability at all, and uniformization -- a Poisson mixture of powers of
% I + Q/lambda -- diverges on a signed generator. Such queries are refused
% rather than answered with a number that looks like a probability.
%
% Mirrors SolverCTMC._assert_phasetype_states in native Python and
% jline.solvers.ctmc.SolverCTMC.assertPhaseTypeStates.

sn = self.model.getStruct();
if isfield(sn,'isph') && ~isempty(sn.isph) && ~all(sn.isph(:))
    line_error(mfilename, sprintf(['%s is unavailable: the model has a matrix-exponential (ME) ' ...
        'service or arrival process, so the stationary vector of the generator is a signed ' ...
        'measure and per-state probabilities and uniformization-based transients do not exist. ' ...
        'Mean measures (getAvg, getAvgTable) remain exact.'], what));
end
end
