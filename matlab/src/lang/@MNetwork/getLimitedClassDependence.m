function gamma = getLimitedClassDependence(self)
% gamma = GETLIMITEDCLASSDEPENDENCE()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
gamma = cell(M,0);
for i=1:M
    if ~isempty(self.stations{i}.lcdScaling)
        gamma{i} = self.stations{i}.lcdScaling; % function handle
    end
end
if ~isempty(gamma)
    % if the model has at least on class-dep station then fill-in the
    % others
    for i=1:M
        if isempty(self.stations{i}.lcdScaling)
            gamma{i} = @(nvec) 1;
        end
    end
end    
end
