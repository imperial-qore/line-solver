function [gamma, cutoffs] = getLimitedJointClassDependence(self)
% [GAMMA, CUTOFFS] = GETLIMITEDJOINTCLASSDEPENDENCE()
%
% Returns the per-class joint-dependence scaling tables and cutoffs for all stations
%
% GAMMA: cell array where gamma{ist}{r} is the linearized scaling table
%        for station ist and class r
% CUTOFFS: cell array where cutoffs{ist} gives per-class cutoffs vector
%          for station ist

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);

gamma = cell(M, 1);
cutoffs = cell(M, 1);

for ist = 1:M
    if ~isempty(self.stations{ist}.ljcdScaling)
        gamma{ist} = self.stations{ist}.ljcdScaling;
        cutoffs{ist} = self.stations{ist}.ljcdCutoffs;
    end
end
end
