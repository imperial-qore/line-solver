function [pi, num, den, stateSpace] = getSymbolicSolution(self)
% [PI, NUM, DEN, STATESPACE] = GETSYMBOLICSOLUTION()
%
% Symbolic stationary distribution of the CTMC as a function of the event rate
% symbols x1, ..., xE, i.e. the solution of pi*Q = 0 with sum(pi) = 1 over the
% field of rational functions in those symbols.
%
% The generator is assembled by getSymbolicGenerator, which needs no computer
% algebra because it is linear in the symbols. Solving with it does, and is
% delegated to the backend named by options.config.symbolic: the Symbolic Math
% Toolbox when it is licensed and the backend is 'auto' or 'matlab', otherwise
% the line-sage-rest service (see SAGE.m). NUM and DEN give the same vector
% over one common denominator.
%
% The expressions are not comparable with another codebase's by text: symbol
% numbering follows event enumeration order and the printed normal form
% depends on the engine. Substitute rates and compare numbers instead, as
% SAGE.eval does.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


% lang='cpp' cannot serve this getter; the reason is named, not blanket.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.cppUnsupported(self.name, 'getSymbolicSolution', ...
        ['the C++ port has no symbolic arithmetic backend; --arith exact is rational and not symbolic']);
end

backend = SolverCTMC.symbolicBackend(self);
timeout = 300;
if isprop(self, 'options') && isfield(self.options, 'config') && ...
        isfield(self.options.config, 'symbolic_timeout')
    timeout = self.options.config.symbolic_timeout;
end

[infGen, ~, ~, stateSpace] = self.getSymbolicGenerator();

useToolbox = SAGE.hasSymbolicToolbox() && ~strcmpi(backend, 'sage') && ...
    ~strncmpi(backend, 'http', 4);
if useToolbox
    pi = ctmc_solve(infGen);
    pi = simplify(pi(:).');
    [num, den] = numden(pi);
    num = simplify(num);
    den = simplify(den);
    if numel(symvar(den)) == 0 && isscalar(unique(den))
        den = den(1);
    end
    return
end

% The service takes the generator as expression strings; a sym generator is
% rendered on the way out, a cell one is already in that form.
symbols = {};
if isa(infGen, 'sym')
    symbols = arrayfun(@char, symvar(infGen), 'UniformOutput', false);
else
    % Symbols are x1..xE by construction; collect the ones that occur.
    joined = strjoin(reshape(infGen, 1, []), ' ');
    tokens = unique(regexp(joined, 'x\d+', 'match'));
    symbols = tokens;
end
url = SAGE.resolve(backend);
if isempty(url)
    url = SAGE.require();
end
[pi, num, den] = SAGE.solveCTMC(infGen, symbols, url, timeout);
end
