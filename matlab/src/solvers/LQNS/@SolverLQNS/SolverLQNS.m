classdef SolverLQNS < LayeredNetworkSolver
    % A solver that interfaces the LQNS to LINE.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverLQNS(model, varargin)
            % SELF = SOLVERLQNS(MODEL, VARARGIN)
            self@LayeredNetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if ~SolverLQNS.isAvailable()
                line_error(mfilename,'SolverLQNS requires the lqns and lqsim commands to be available on the system path. Please visit: http://www.sce.carleton.ca/rads/lqns/');
            end
        end

        runtime = runAnalyzer(self, options);        
        [result, iterations] = parseXMLResults(self, filename);

        function varargout = getAvg(varargin)
            % [QN,UN,RN,TN] = GETAVG(SELF,~,~,~,~,USELQNSnaming)
            [varargout{1:nargout}] = getEnsembleAvg( varargin{:} );
        end
        
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self,~,~,~,~, useLQNSnaming);
    end
    
    methods (Static)
        
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','lqns','srvn','exactmva','srvn.exactmva','sim','lqsim','lqnsdefault'};
        end

        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            
            bool = true;
            if ispc % windows
                [~,ret] = dos('lqns -V -H');
                if containsstr(ret,'not recognized')
                    bool = false;
                end
                if containsstr(ret,'Version 5') || containsstr(ret,'Version 4') ...
                        || containsstr(ret,'Version 3') || containsstr(ret,'Version 2') ...
                        || containsstr(ret,'Version 1')
                    line_warning(mfilename,'Unsupported LQNS version. LINE requires Version 6.0 or greater.\n');
                    bool = true;
                end
            else % linux
                [~,ret] = unix('lqns -V -H');
                if containsstr(ret,'command not found')
                    bool = false;
                end
                if containsstr(ret,'Version 5') || containsstr(ret,'Version 4') ...
                        || containsstr(ret,'Version 3') || containsstr(ret,'Version 2') ...
                        || containsstr(ret,'Version 1')
                    line_warning(mfilename,'Unsupported LQNS version. LINE requires Version 6.0 or greater.\n');
                    bool = true;
                end
            end
        end
        
        function savedfname = plot(model)
            % FNAME = PLOT(model)
            %
            % Save a plot of an LQN model using the lqn2ps utility
            lastwd = pwd;
            cwd = [lineRootFolder,filesep,'workspace',filesep,'lqns'];
            cd(cwd);
            fname = tempname([lineRootFolder,filesep,'workspace',filesep,'lqns']);
            model.writeXML([fname,'.lqnx'],false);
            system(sprintf('lqn2ps %s',[fname,'.lqnx']));
            line_printf('Postscript file saved in: %s\n',[fname,'.ps'])
            cd(lastwd);
            savedfname = [fname,'.ps'];
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Queue',...
                'Coxian','Erlang','Exp','HyperExp',...
                'Buffer','Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_PS','SchedStrategy_FCFS','ClosedClass'});
            bool = true;
            for e=1:model.getNumberOfLayers()
                bool = bool && SolverFeatureSet.supports(featSupported, featUsed{e});
            end
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('LQNS');
        end

    end
end
