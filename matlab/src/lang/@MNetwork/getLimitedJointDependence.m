function [gamma, cutoffs] = getLimitedJointDependence(self)
% [GAMMA, CUTOFFS] = GETLIMITEDJOINTDEPENDENCE()
%
% Returns the joint-dependence scaling tables and cutoffs for all stations
%
% GAMMA: cell array where gamma{ist} is the linearized scaling table for station ist
% CUTOFFS: M x K matrix where cutoffs(ist,:) gives per-class cutoffs

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);

gamma = cell(M, 1);
cutoffs = zeros(M, K);

for ist = 1:M
    if ~isempty(self.stations{ist}.ljdScaling)
        gamma{ist} = self.stations{ist}.ljdScaling;
        cutoffs(ist, :) = self.stations{ist}.ljdCutoffs;
    end
end
end
