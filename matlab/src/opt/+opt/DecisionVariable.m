classdef DecisionVariable < handle
    % DecisionVariable  Abstract base for line-opt decision variables. Mirrors
    % native-Python line_solver.opt.variables.DecisionVariable: each variable
    % encodes a tunable model parameter as continuous values in [0,1]
    % (getBounds), decodes them to the native domain (decode), and applies the
    % decoded value to a per-evaluation model copy (apply). Objects are
    % re-resolved by name in the target model because models are copied per
    % evaluation.

    properties
        name
        dimension = 1;
    end

    methods
        function obj = DecisionVariable(name)
            obj.name = name;
        end
        function n = getName(obj)
            n = obj.name;
        end
        function d = getDimension(obj)
            d = obj.dimension;
        end

        function layers = getLayer(~, ~)
            % LQN layer name(s) this variable perturbs, or {} for flat models.
            % Consumed by layer freezing (explicit frozenLayers and adaptive
            % auto-freeze): a variable whose layer set intersects the frozen
            % set is held fixed. LQN variable subclasses override this.
            layers = {};
        end

        function v = currentValue(~, ~)
            % The variable's current (decoded) value in the given model, or []
            % when not introspectable. Used by layer freezing to hold a
            % variable at the model's existing parameter value. LQN variable
            % subclasses override this.
            v = [];
        end
    end

    methods (Abstract)
        b = getBounds(obj)            % dimension x 2
        v = decode(obj, x)            % x is 1 x dimension in [0,1]
        apply(obj, model, value)
        t = getVariableType(obj)
    end

    methods (Static)
        function c = resolveClass(model, jobclass)
            c = [];
            target = jobclass.getName();
            classes = model.getClasses();
            for i = 1:numel(classes)
                if strcmp(classes{i}.getName(), target)
                    c = classes{i};
                    return;
                end
            end
        end

        function nd = resolveNode(model, node)
            nd = [];
            target = node.getName();
            nodes = model.getNodes();
            for i = 1:numel(nodes)
                if strcmp(nodes{i}.getName(), target)
                    nd = nodes{i};
                    return;
                end
            end
        end

        function conn = connectionMatrix(model)
            sn = model.getStruct();
            conn = full(sn.connmatrix);
        end

        function b = unitBounds(dim)
            b = [zeros(dim, 1), ones(dim, 1)];
        end

        function idx = indexOfNode(nodes, nameToFind)
            idx = -1;
            for i = 1:numel(nodes)
                if strcmp(nodes{i}.getName(), nameToFind)
                    idx = i;
                    return;
                end
            end
        end
    end
end
