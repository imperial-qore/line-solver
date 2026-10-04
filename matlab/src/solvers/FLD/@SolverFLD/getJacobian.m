function [J, rhs, vars, equilibria] = getJacobian(self, options)
% [J, RHS, VARS, EQUILIBRIA] = GETJACOBIAN(OPTIONS)
%
% Jacobian of the mean-field ODE right-hand side, d f_i / d x_j, as a cell
% matrix of expression strings, computed exactly by the computer algebra
% backend (see SAGE.m).
%
% The Jacobian is what tells a fixed point apart from a limit cycle and gives
% the local convergence rate of the fluid approximation, neither of which a
% numerical integration reports. EQUILIBRIA, returned only when asked for as a
% fourth output, are the solutions of f(x) = 0; they can be empty when the
% system is beyond what the backend solves in closed form, which is a
% limitation of the solve and not an assertion that none exist.
%
% Only smooth drifts have a Jacobian: see getSymbolicDrift, which refuses the
% min-scaled methods by name rather than returning a one-sided derivative.
%
% @param self The SolverFLD instance
% @param options Solver options (optional, defaults to the solver's own)
% @return J Cell matrix of expressions, J{i,j} = d f_i / d x_j
% @return rhs The drift itself, one expression per state variable
% @return vars The state variable names
% @return equilibria Struct array of solutions of f(x) = 0 (optional)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


% lang='cpp' takes the Jacobian from line-cli (-s fluid -a jacobian). Both
% sides are SYMBOLIC and reach the same computer-algebra backend: the C++ arm
% builds the drift with fluid_symodes and differentiates it there, exactly as
% the SAGE call below does, so the entries are expressions on both paths.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [J, rhs, vars, equilibria] = CPPLINE.jacobian(self.name, self.model, self.options, nargout >= 4);
    return
end

if nargin < 2 || isempty(options)
    options = self.getOptions();
end
[rhs, vars] = self.getSymbolicDrift(options);

backend = 'auto';
if isfield(options, 'config') && isfield(options.config, 'symbolic')
    backend = options.config.symbolic;
end
url = SAGE.resolve(backend);
if isempty(url)
    url = SAGE.require();
end

want = {'jacobian'};
if nargout >= 4
    want{end+1} = 'equilibria';
end
timeout = 300;
if isfield(options, 'config') && isfield(options.config, 'symbolic_timeout')
    timeout = options.config.symbolic_timeout;
end
[J, ~, equilibria] = SAGE.fluidODEs(rhs, vars, want, url, timeout);
end
