function varargout = delegate(self, method, nretout, varargin)
this_model = self.model;
lastError = [];

switch class(this_model)
    case 'Network'
        % first try with chosen solver, if the method is not available
        %     or fails keep going with the other candidates        
        if length(self.candidates)>1
            % Solver console: AUTO opens no run of its own, so that the
            % solver it picks narrates its own analysis rather than being
            % nested and silenced. It announces the choice instead.
            LineConsole.step('AUTO selecting a solver for %s among %d candidates', ...
                method, length(self.candidates));
            chosenSolver = chooseSolver(self, method);
            LineConsole.substep('chose %s (supports this model: %d)', ...
                class(chosenSolver), chosenSolver.supports(self.model));
            line_debug('AUTO solver analysis: chosen=%s, supports_model=%d', class(chosenSolver), chosenSolver.supports(self.model));
            if chosenSolver.supports(self.model)
                proposedSolvers = {chosenSolver, self.candidates{:}}; %#ok<CCAT>
            else
                proposedSolvers = self.candidates;
            end
        else % if the user wishes to use a precise solver
            proposedSolvers = {self.solvers{:}};
        end
        for s=1:length(proposedSolvers)
            try
                switch nretout
                    case 0
                        proposedSolvers{s}.(method)(varargin{:});
                    case 1
                        varargout{1} = proposedSolvers{s}.(method)(varargin{:});
                    case 2
                        [r1,r2] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                    case 3
                        [r1,r2,r3] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                    case 4
                        [r1,r2,r3,r4] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                    case 5
                        [r1,r2,r3,r4,r5] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                    case 6
                        [r1,r2,r3,r4,r5,r6] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                        varargout{6} = r6;
                    case 7
                        [r1,r2,r3,r4,r5,r6,r7] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                        varargout{6} = r6;
                        varargout{7} = r7;
                end
                line_debug('Successful method execution completed by %s', proposedSolvers{s}.getName);
                if self.options.verbose
                    %line_printf('Successful method execution completed by %s.\n',proposedSolvers{s}.getName);
                end
                return
            catch ME
                lastError = ME;
                if ~isempty(strfind(ME.message,'Unrecognized method'))
                    if self.options.verbose
                        line_printf('Method unsupported by %s.\n',proposedSolvers{s}.getName);
                    end
                else
                    %keyboard
                    line_warning(mfilename,[ME.message,'\n']);
                end
            end
        end
    case {'LayeredNetwork','Environment'}
        % first try with chosen solver, if the method is not available
        %     or fails keep going with the other candidates. Environment models
        %     take the same path: their candidates are the blending solvers,
        %     and without this case no call on them would reach a delegate.
        if numel(self.candidates) > 1
            chosenSolver = chooseSolver(self,method);
            % LayeredNetwork solvers use static-style supports(model) signature
            try
                [supportsModel, ~] = chosenSolver.supports(self.model);
            catch
                % Fallback: assume solver supports the model if supports() signature differs
                supportsModel = true;
            end
            if supportsModel
                proposedSolvers = {chosenSolver, self.candidates{:}}; %#ok<CCAT>
            else
                proposedSolvers = self.candidates;
            end
        else % a forced method name leaves a single solver, which owns the diagnostic
            proposedSolvers = {self.solvers{:}}; %#ok<CCAT>
        end
        for s=1:length(proposedSolvers)
            try
                switch nretout
                    case 0
                        proposedSolvers{s}.(method)(varargin{:});
                    case 1
                        varargout{1} = proposedSolvers{s}.(method)(varargin{:});
                    case 2
                        [r1,r2] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                    case 3
                        [r1,r2,r3] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                    case 4
                        [r1,r2,r3,r4] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                    case 5
                        [r1,r2,r3,r4,r5] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                    case 6
                        [r1,r2,r3,r4,r5,r6] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                        varargout{6} = r6;
                    case 7
                        [r1,r2,r3,r4,r5,r6,r7] = proposedSolvers{s}.(method)(varargin{:});
                        varargout{1} = r1;
                        varargout{2} = r2;
                        varargout{3} = r3;
                        varargout{4} = r4;
                        varargout{5} = r5;
                        varargout{6} = r6;
                        varargout{7} = r7;
                end
                if self.options.verbose
                    line_printf('Successful method execution completed by %s.\n',proposedSolvers{s}.getName);
                end
                return
            catch ME
                lastError = ME;
                if self.options.verbose
                    line_printf('Switching %s.\n',proposedSolvers{s}.getName);
                end
            end
        end
end
% If we reach here, all solvers failed. With a single proposed solver (a forced
% method name such as 'bound' or 'mva') the delegate's own diagnostic is the whole
% story, so surface it instead of the generic message that hides it.
if ~isempty(lastError)
    if numel(proposedSolvers) == 1
        rethrow(lastError);
    end
    line_error(mfilename, 'All candidate solvers failed to execute method ''%s''. Last error: %s', ...
        method, lastError.message);
end
line_error(mfilename, 'All candidate solvers failed to execute method ''%s''.', method);
end
