classdef WorkflowActivity < Element
    % A computational activity in a Workflow.
    %
    % WorkflowActivity represents a stage of computation in a standalone
    % workflow model. Unlike Activity in LayeredNetwork, this class is
    % designed for pure computational workflows without external calls.
    %
    % Each activity has a service time distribution that can be any
    % phase-type compatible distribution (Exp, Erlang, APH, HyperExp, etc.).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        hostDemand;         % Distribution object (Exp, APH, Erlang, etc.)
        hostDemandMean;     % double - mean service time
        hostDemandSCV;      % double - squared coefficient of variation
        workflow;           % Reference to parent Workflow
        index;              % Index in the workflow's activity list
        metadata;           % Optional struct for external metadata (e.g., WfCommons)
    end

    methods
        function obj = WorkflowActivity(workflow, name, hostDemand)
            % WORKFLOWACTIVITY Create a workflow activity
            %
            % OBJ = WORKFLOWACTIVITY(WORKFLOW, NAME, HOSTDEMAND)
            %
            % Parameters:
            %   workflow   - Parent Workflow object
            %   name       - Activity name (string)
            %   hostDemand - Service time distribution or numeric mean
            %
            % If hostDemand is numeric, it is treated as the mean of an
            % exponential distribution.

            if nargin < 2
                line_error(mfilename, 'Constructor requires at least workflow and name.');
            end

            obj@Element(name);
            obj.workflow = workflow;

            if nargin < 3
                hostDemand = GlobalConstants.FineTol;
            elseif isnumeric(hostDemand) && hostDemand == 0
                hostDemand = GlobalConstants.FineTol;
            end

            obj.setHostDemand(hostDemand);
        end

        function obj = setHostDemand(obj, hostDemand)
            % SETHOSTDEMAND Set the service time distribution
            %
            % OBJ = SETHOSTDEMAND(OBJ, HOSTDEMAND)
            %
            % Parameters:
            %   hostDemand - Distribution object or numeric mean

            if isnumeric(hostDemand)
                if hostDemand <= GlobalConstants.FineTol
                    obj.hostDemand = Immediate.getInstance();
                    obj.hostDemandMean = GlobalConstants.FineTol;
                    obj.hostDemandSCV = GlobalConstants.FineTol;
                else
                    obj.hostDemand = Exp(1/hostDemand);
                    obj.hostDemandMean = hostDemand;
                    obj.hostDemandSCV = 1.0;
                end
            elseif isa(hostDemand, 'Distribution')
                obj.hostDemand = hostDemand;
                obj.hostDemandMean = hostDemand.getMean();
                obj.hostDemandSCV = hostDemand.getSCV();
            else
                line_error(mfilename, 'hostDemand must be a Distribution or numeric value.');
            end
        end

        function [alpha, T] = getPHRepresentation(obj)
            % GETPHREPRESENTATION Get the phase-type representation
            %
            % [ALPHA, T] = GETPHREPRESENTATION(OBJ)
            %
            % Returns:
            %   alpha - Initial probability vector (1 x n)
            %   T     - Subgenerator matrix (n x n)
            %
            % For Markovian distributions, extracts the PH representation.
            % For other distributions, fits an APH from the first two moments.

            if isa(obj.hostDemand, 'Immediate')
                % Immediate activity: single absorbing state
                alpha = 1;
                T = -1e10;  % Very high rate (essentially immediate)
                return;
            end

            if isa(obj.hostDemand, 'Markovian')
                % Extract from Markovian distribution
                if ismethod(obj.hostDemand, 'getInitProb') && ismethod(obj.hostDemand, 'getSubgenerator')
                    alpha = obj.hostDemand.getInitProb();
                    T = obj.hostDemand.getSubgenerator();
                else
                    % Use process representation
                    proc = obj.hostDemand.getProcess();
                    T = proc{1};
                    % Extract alpha from D1 = -T*e*alpha
                    n = size(T, 1);
                    e = ones(n, 1);
                    absRate = -T * e;
                    D1 = proc{2};
                    % D1(i,j) = absRate(i) * alpha(j), so alpha = D1(k,:) / absRate(k) for any k with absRate(k) > 0
                    idx = find(absRate > GlobalConstants.FineTol, 1);
                    if ~isempty(idx)
                        alpha = D1(idx, :) / absRate(idx);
                    else
                        alpha = zeros(1, n);
                        alpha(1) = 1;
                    end
                end
            else
                % Fit APH from moments
                mean_val = obj.hostDemandMean;
                scv_val = obj.hostDemandSCV;
                if scv_val < GlobalConstants.FineTol
                    scv_val = 1.0;  % Default to exponential
                end
                aph = APH.fitMeanAndSCV(mean_val, scv_val);
                alpha = aph.getInitProb();
                T = aph.getSubgenerator();
            end

            % Ensure alpha is a row vector
            alpha = reshape(alpha, 1, []);
        end

        function n = getNumberOfPhases(obj)
            % GETNUMBEROFPHASES Get the number of phases
            %
            % N = GETNUMBEROFPHASES(OBJ)
            %
            % Returns the number of phases in the PH representation.

            if isa(obj.hostDemand, 'Immediate')
                n = 1;
            elseif isa(obj.hostDemand, 'Markovian')
                if isprop(obj.hostDemand, 'nPhases')
                    n = obj.hostDemand.nPhases;
                else
                    [~, T] = obj.getPHRepresentation();
                    n = size(T, 1);
                end
            else
                % Fit to get phase count
                [~, T] = obj.getPHRepresentation();
                n = size(T, 1);
            end
        end
    end
end
