classdef Ensemble < Model
    % A model defined by a collection of sub-models.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        ensemble;
    end
    
    methods
        function self = Ensemble(models)
            % SELF = ENSEMBLE(MODELS)
            
            self@Model('Ensemble');
            self.ensemble = reshape(models,1,numel(models)); % flatten
        end
        
        function self = setEnsemble(self,ensemble)
            % SELF = SETENSEMBLE(SELF,ENSEMBLE)
            
            self.ensemble = ensemble;
        end
        
        function ensemble = getEnsemble(self)
            % ENSEMBLE = GETENSEMBLE()
            
            ensemble = self.ensemble;
        end
        
        function model = getModel(self, modelIdx)
            model = self.ensemble{modelIdx};
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each ensemble object
            for e=1:length(self.ensemble)
                clone.ensemble{e} = copy(self.ensemble{e});
            end
        end
    end
end
