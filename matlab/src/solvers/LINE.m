classdef LINE < SolverAuto
%{
%{
 % @class LINE
 % @brief Main entry point for the LINE solver.
 % @details Provides a unified interface to access various solvers (CTMC, MVA, SSA, JMT, Fluid, NC, MAM).
 % Can also be used to load models from files.
%}
%}
    % Alias for SolverAuto
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = LINE(model, varargin)
%{
%{
 % @brief Constructor for the LINE solver wrapper.
 % @param model The network model (Network or LayeredNetwork).
 % @param varargin Optional arguments passed to SolverAuto.
 % @return self Instance of the LINE solver.
%}
%}
            % SELF = LINE(MODEL, VARARGIN)
            self@SolverAuto(model, varargin{:});
        end
    end

    methods (Static)
        function result = load(varargin)
%{
%{
 % @brief Factory method to load models or instantiate solvers.
 % @details
 % Usage 1: Load a model from file.
 %    model = LINE.load(filename)
 %    model = LINE.load(filename, verbose)
 % Supported formats: .jsim, .jsimg, .jsimw (JMT), .jmva (JMVA), .xml/lqn/lqnx (LQN), .mat (MATLAB).
 %
 % Usage 2: Instantiate a solver with a specific method.
 %    solver = LINE.load(method, model, options...)
 %
 % @param varargin Arguments for loading a model or creating a solver.
 % @return result Loaded model object or Solver instance.
%}
%}
            % RESULT = LOAD(...)
            % Flexible loading method that handles both solver and model loading
            %
            % Usage 1: Load a model from file
            %   model = LINE.load(filename)
            %   model = LINE.load(filename, verbose)
            %
            % Usage 2: Load a solver with chosen method (legacy)
            %   solver = LINE.load(method, model, options...)

            % Check if this is model loading (1-2 args with first arg as string filename)
            if nargin >= 1 && nargin <= 2 && ischar(varargin{1}) && contains(varargin{1}, '.')
                % Check if this is actually a filename with a supported extension
                [~, ~, ext] = fileparts(varargin{1});
                supported_extensions = {'.jsim', '.jsimg', '.jsimw', '.jmva', '.xml', '.lqn', '.lqnx', '.mat'};

                if ismember(lower(ext), supported_extensions)
                    filename = varargin{1};
                    verbose = false;
                    if nargin == 2
                        verbose = varargin{2};
                    end

                    % Determine file format and load accordingly
                    switch lower(ext)
                    case {'.jsim', '.jsimg', '.jsimw'}
                        % Load JSIM/JMT model
                        result = JSIM2LINE(filename);
                    case {'.jmva'}
                        % Load JMVA model
                        result = JMVA2LINE(filename);
                    case {'.xml', '.lqn', '.lqnx'}
                        % Try to load as LQN model
                        try
                            result = LQN2MATLAB(filename);
                        catch
                            % If LQN fails, try as generic XML
                            error('LINE:load', 'Unsupported XML format. Only LQN XML files are supported.');
                        end
                    case '.mat'
                        % Load MATLAB saved model
                        data = load(filename);
                        if isfield(data, 'model')
                            result = data.model;
                        elseif isfield(data, 'network')
                            result = data.network;
                        else
                            % Find first Network object in the file
                            fields = fieldnames(data);
                            found = false;
                            for i = 1:length(fields)
                                if isa(data.(fields{i}), 'Network')
                                    result = data.(fields{i});
                                    found = true;
                                    break;
                                end
                            end
                            if ~found
                                error('LINE:load', 'No Network model found in MAT file.');
                            end
                        end

                        if verbose
                            fprintf('Loaded model "%s" from %s\n', result.name, filename);
                        end
                    otherwise
                        error('LINE:load', 'Unsupported file format: %s', ext);
                    end

                    return;
                end
                % If we get here, the string contains a dot but is not a supported file format
                % Continue to solver loading logic below
            end

            % Otherwise, this is solver loading (legacy behavior)
            if nargin < 2
                error('LINE:load', 'Solver loading requires at least method and model arguments.');
            end

            chosenmethod = varargin{1};
            model = varargin{2};
            options_args = varargin(3:end);

            options = Solver.parseOptions(options_args, Solver.defaultOptions);
            options.method = chosenmethod;
            switch options.method
                case {'default','auto'}
                    if strcmp(options.method,'auto'), options.method='default'; end
                    result = LINE(model, options);
                case {'ctmc','ctmc.gpu','gpu'}
                    if strcmp(options.method,'ctmc'), options.method='default'; end
                    options.method = erase(options.method,'ctmc.');
                    result = SolverCTMC(model, options);
                case {'mva','mva.exact','amva','mva.amva','qna','mva.qna','sqrt','mva.sqrt' ...
                        'amva.qd','mva.amva.qd', ...
                        'amva.bs','mva.amva.bs', ...
                        'amva.qli','mva.amva.qli', ...
                        'amva.fli','mva.amva.fli', ...
                        'amva.lin','mva.amva.lin', ...
                        'mm1','mmk','mg1','gm1','gig1','gim1','gig1.kingman', ...
                        'gigk', 'gigk.kingman_approx', ...
                        'gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal', ...
                        'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'gb.upper', 'gb.lower', 'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower'}
                        %'amva.aql','mva.amva.aql', ...
                        %'amva.qdaql','mva.amva.qdaql', ...
                    if strcmp(options.method,'mva'), options.method='default'; end
                    options.method = erase(options.method,'mva.');
                    result = SolverMVA(model, options);
                case {'nrm','ssa','ssa.serial','ssa.parallel','serial','parallel'}
                    if strcmp(options.method,'ssa'), options.method='default'; end
                    options.method = erase(options.method,'ssa.');
                    result = SolverSSA(model, options);
                case {'jmt','jsim','jmva','jmva.mva','jmva.recal','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin',... %,'jmva.ls','jmt.jmva.ls'
                        'jmt.jsim','jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin'}
                    if strcmp(options.method,'jmt'), options.method='default'; end
                    options.method = erase(options.method,'jmt.');
                    result = SolverJMT(model, options);
                case {'fluid','fluid.softmin','fluid.statedep','fluid.closing'}
                    if strcmp(options.method,'fluid'), options.method='default'; end
                    options.method = erase(options.method,'fluid.');
                    result = SolverFluid(model, options);
                case {'nc','nc.exact','nc.imci','nc.ls','comom','comomld','cub','ls','nc.le','le','mmint2','nc.panacea','nc.pana','nc.mmint2','nc.kt','nc.deterministic','nc.sampling','nc.propfair','nc.comom','nc.comomld','nc.mom','nc.cub','nc.brute','nc.rd', 'nc.nr.probit', 'nc.nr.logit','nc.gm'}
                    if strcmp(options.method,'nc'), options.method='default'; end
                    options.method = erase(options.method,'nc.');
                    result = SolverNC(model, options);
                case {'mam','mam.dec.source','mam.dec.mmap','mam.dec.poisson'}
                    if strcmp(options.method,'mam'), options.method='default'; end
                    options.method = erase(options.method,'mam.');
                    result = SolverMAM(model, options);
                otherwise
                    if strcmp(options.method,'auto'), options.method='default'; end
                    result = LINE(model, options);
            end
        end
    end
end
