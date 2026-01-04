function resolveSignals(self)
% RESOLVESIGNALS Resolve Signal placeholders to OpenSignal or ClosedSignal
%
% This method is called during model finalization (refreshStruct) to
% convert Signal placeholder objects to their concrete types (OpenSignal
% or ClosedSignal) based on the network structure.
%
% For open networks (with Source node): Signal -> OpenSignal
% For closed networks (no Source node): Signal -> ClosedSignal
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Check if there are any Signal placeholders to resolve
hasSignals = false;
for r = 1:length(self.classes)
    % Check for Signal that is NOT already OpenSignal or ClosedSignal
    if isa(self.classes{r}, 'Signal') && ...
       ~isa(self.classes{r}, 'OpenSignal') && ...
       ~isa(self.classes{r}, 'ClosedSignal')
        hasSignals = true;
        break;
    end
end

if ~hasSignals
    return;  % No Signal placeholders to resolve
end

% Determine if the network is open (has Source node) or closed
isOpen = false;
for i = 1:length(self.nodes)
    if isa(self.nodes{i}, 'Source')
        isOpen = true;
        break;
    end
end

% For closed networks, find a default reference station
defaultRefstat = [];
if ~isOpen
    % Look for first Delay node, then first Queue node
    for i = 1:length(self.stations)
        if isa(self.stations{i}, 'Delay')
            defaultRefstat = self.stations{i};
            break;
        end
    end
    if isempty(defaultRefstat)
        for i = 1:length(self.stations)
            if isa(self.stations{i}, 'Queue')
                defaultRefstat = self.stations{i};
                break;
            end
        end
    end
    if isempty(defaultRefstat) && ~isempty(self.stations)
        % Fallback to first station
        defaultRefstat = self.stations{1};
    end
end

% Collect Signal placeholders to resolve (indices will shift during resolution)
signalsToResolve = {};
for r = 1:length(self.classes)
    if isa(self.classes{r}, 'Signal') && ...
       ~isa(self.classes{r}, 'OpenSignal') && ...
       ~isa(self.classes{r}, 'ClosedSignal')
        signalsToResolve{end+1} = self.classes{r}; %#ok<AGROW>
    end
end

% Resolve each Signal placeholder
% Note: resolve() removes the placeholder and creates concrete signal
for s = 1:length(signalsToResolve)
    sig = signalsToResolve{s};

    % Determine reference station for closed signals
    if ~isOpen
        % Prefer the reference station of the targetJobClass if available
        if ~isempty(sig.targetJobClass) && isprop(sig.targetJobClass, 'refstat')
            refstat = sig.targetJobClass.refstat;
        else
            refstat = defaultRefstat;
        end
    else
        refstat = [];  % Not used for open signals
    end

    % Resolve the Signal to its concrete type
    % This removes the placeholder from model.classes and adds the concrete
    sig.resolve(isOpen, refstat);
end

end
