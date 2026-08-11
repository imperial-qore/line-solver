function varargout = delegate(self, method, nretout, varargin)
this_model = self.model;

switch class(this_model)
    case 'Network'
        % first try with chosen solver, if the method is not available
        %     or fails keep going with the other candidates        
        if length(self.candidates)>1
            chosenSolver = chooseSolver(self, method);
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
    case 'LayeredNetwork'
        % first try with chosen solver, if the method is not available
        %     or fails keep going with the other candidates
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
            catch
                if self.options.verbose
                    line_printf('Switching %s.\n',proposedSolvers{s}.getName);
                end
            end
        end
end
% If we reach here, all solvers failed
line_error(mfilename, 'All candidate solvers failed to execute method ''%s''.', method);
end
