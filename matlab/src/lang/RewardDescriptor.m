classdef RewardDescriptor
    % REWARDDESCRIPTOR Callable reward function carrying its own metadata
    %
    % A RewardDescriptor wraps a reward function handle together with a
    % structural description of what the reward measures. It is CALLABLE
    % exactly like the bare function handle it replaces, so that
    %
    %   fn = Reward.queueLength(queue1);
    %   value = fn(state);
    %
    % keeps working unchanged, while additionally exposing
    %
    %   fn.kind      - 'QLen' | 'Util' | 'Blocking' | 'Custom'
    %   fn.node      - Node object the reward refers to ([] if none)
    %   fn.jobclass  - JobClass object the reward refers to ([] if none)
    %   fn.fn        - the underlying function handle
    %
    % The metadata is what allows setReward/linemodel_save to serialize the
    % reward declaratively. A descriptor of kind 'Custom' wraps an arbitrary
    % user function and is deliberately NOT serializable: the writer warns and
    % omits it rather than emitting a reward it cannot reproduce.
    %
    % See also: Reward, setReward, RewardState
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        kind;     % char: 'QLen' | 'Util' | 'Blocking' | 'Custom'
        node;     % Node object, or [] when the reward is not node-scoped
        jobclass; % JobClass object, or [] when the reward is not class-scoped
        fn;       % function_handle actually evaluated
    end

    methods
        function self = RewardDescriptor(kind, node, jobclass, fn)
            % SELF = REWARDDESCRIPTOR(KIND, NODE, JOBCLASS, FN)
            if ~ischar(kind) && ~isstring(kind)
                line_error(mfilename, 'Reward kind must be a string.');
            end
            if ~isa(fn, 'function_handle')
                line_error(mfilename, 'Reward descriptor must wrap a function handle.');
            end
            self.kind = char(kind);
            self.node = node;
            self.jobclass = jobclass;
            self.fn = fn;
        end

        function varargout = subsref(self, s)
            % SUBSREF Make the descriptor callable as fn(state[, sn])
            switch s(1).type
                case '()'
                    out = self.fn(s(1).subs{:});
                    if numel(s) > 1
                        [varargout{1:max(nargout,1)}] = subsref(out, s(2:end));
                    else
                        varargout{1} = out;
                    end
                otherwise
                    [varargout{1:max(nargout,1)}] = builtin('subsref', self, s);
            end
        end

        function n = numArgumentsFromSubscript(self, s, indexingContext) %#ok<INUSD>
            % NUMARGUMENTSFROMSUBSCRIPT A reward always produces exactly one value
            n = 1;
        end
    end
end
