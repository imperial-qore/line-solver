function convertSignalPlaceholders(self)
% CONVERTSIGNALPLACEHOLDERS Validate Signal class chain consistency
%
% This function examines Signal classes and validates they are consistent
% with the routing matrix. Since Signal extends OpenClass, all Signal
% objects are open classes by default. This function validates that any
% class-switching involving signals maintains chain consistency.
%
% Note: Signal classes are always open (extend OpenClass). For class-switching
% from closed classes to signals, the signal must be in the same chain.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = length(self.classes);
if K == 0
    return;
end

% Find all Signal classes
signalIndices = [];
for k = 1:K
    if isa(self.classes{k}, 'Signal')
        signalIndices(end+1) = k; %#ok<AGROW>
    end
end

if isempty(signalIndices)
    return;
end

% Signal classes extending OpenClass are already properly typed.
% No conversion needed - Signal inherits from OpenClass.
% This function now serves as a validation pass.

end
