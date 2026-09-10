classdef Forker < OutputSection
    % An output section forking jobs into sibling tasks
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        % Nominal tasks emitted on each outgoing link, shared by every link and
        % every class. For a fork whose degree is random this is the MEAN, so a
        % consumer that only knows this field gets E[tasks per link] rather than
        % a number the fork never emits.
        tasksPerLink;
        % Per-destination, per-class overrides, each a containers-free struct
        % array of (destNodeName, class index, value). Empty means "uniform,
        % use tasksPerLink". They are kept by DESTINATION NAME rather than by
        % link ordinal because the link order is an artefact of connmatrix
        % traversal and silently renumbers when the model is relinked.
        tasksPerLinkByDest = struct('dest',{},'class',{},'value',{});
        % Per-destination, per-class jobs-per-link distributions (DiscreteSampler
        % over a non-negative integer support). Empty means degenerate at the
        % corresponding tasksPerLink entry.
        tasksPerLinkDist = struct('dest',{},'class',{},'dist',{});
        % Per-destination, per-class branch activation probability. Empty means
        % 1.0: the branch always fires, which is the classic fork.
        branchProb = struct('dest',{},'class',{},'value',{});
    end

    methods
        %Constructor
        function self = Forker(customerClasses)
            % SELF = FORKER(CUSTOMERCLASSES)

            self@OutputSection('Forker');
            self.tasksPerLink = 1.0;
            initDispatcherJobClasses(self, customerClasses);
        end
    end
    
    methods (Access = 'private')
        function initDispatcherJobClasses(self, customerClasses)
            % INITDISPATCHERJOBCLASSES(CUSTOMERCLASSES)
            
            for i = 1 : length(customerClasses)
                self.outputStrategy{i} = {customerClasses{i}.name, RoutingStrategy.RAND};
            end
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i = 1 : length(self.outputStrategy)
                if ishandle(clone.outputStrategy{i}{1})
                    % this is a problem if one modifies the classes in the
                    % model because the one below is not an handle so it
                    % will not be modified
                    clone.outputStrategy{i}{1} = self.outputStrategy{i}{1}.copy;
                end
            end
        end
    end
end
