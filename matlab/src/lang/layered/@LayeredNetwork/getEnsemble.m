function ensemble = getEnsemble(self, force)
% ENSEMBLE = GETENSEMBLE()
if isempty(self.ensemble)
    SolverLN(self).buildLayers;
end
ensemble = self.ensemble;
end