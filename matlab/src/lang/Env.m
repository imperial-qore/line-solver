classdef Env < Ensemble
    % An environment model defined by a collection of network sub-models
    % coupled with an environment transition rule that selects the active
    % sub-model.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        env;
        envGraph;
        proc; % Markovian representation of each stage transition
        holdTime; % holding times
        probEnv; % steady-stage probability of the environment
        probOrig; % probability that a request originated from phase
        resetFun; % function implementing the reset policy        
    end
    
    methods
        function self = Env(name)
            % SELF = ENVIRONMENT(MODELS, ENV)
            self@Ensemble({});
            self.name = name;
            self.envGraph = digraph();
            self.envGraph.Nodes.Model = cell(0);
            self.envGraph.Nodes.Type = cell(0);
        end
        
        function name = addStage(self, name, type, model)            
            wcfg = warning; % store warning configuration
            warning('off','MATLAB:table:RowsAddedExistingVars');
            self.envGraph = self.envGraph.addnode(name);
            warning(wcfg); % restore warning configuration
            self.envGraph.Nodes.Model{end} = model;
            self.envGraph.Nodes.Type{end} = type;
            E = height(self.envGraph.Nodes);
            if E>1
                if self.envGraph.Nodes.Model{1}.getNumberOfStatefulNodes ~= model.getNumberOfStatefulNodes
                    line_error(mfilename,'Unsupported feature. Random environment stages must map to networks with identical number of stateful nodes.');
                end
            end
            for e=E
                for h=1:E
                    self.env{e,h} = Disabled.getInstance();
                    self.env{h,e} = Disabled.getInstance();
                end
            end
            self.ensemble{E} = model;
        end
        
        function self = addTransition(self, fromName, toName, distrib, resetFun)            
           self.envGraph = self.envGraph.addedge(fromName, toName);            
           e = self.envGraph.findnode(fromName);
           h = self.envGraph.findnode(toName);
           self.env{e,h} = distrib;
           if nargin<5 %~exist('resetFun','var')
                self.resetFun{e,h} = @(q) q;
           else
                self.resetFun{e,h} = resetFun;
           end 
        end
        
        function init(self)
            E = height(self.envGraph.Nodes);
            T = height(self.envGraph.Edges);
            Pemb = zeros(E);

            for t=1:T
                e = self.envGraph.findnode(self.envGraph.Edges.EndNodes{t,1});
                h = self.envGraph.findnode(self.envGraph.Edges.EndNodes{t,2});
                self.envGraph.Edges.Distribution{t} = self.env{e,h};
            end
            
            % analyse holding times
            emmap = cell(E,1);
            for e=1:E
                emmap{e} = cell(1,E);
            end
            self.holdTime = {};
            for e=1:E
                for h=1:E
                    if isa(self.env{e,h},'Disabled')
                        emmap{e}{h} = {0,0}; % multiclass MMAP representation
                    else
                        emmap{e}{h} = self.env{e,h}.getRepres; % multiclass MMAP representation
                    end
                    for j = 1:E
                        if j == h
                            emmap{e}{h}{2+j} = emmap{e}{h}{2};
                        else
                            emmap{e}{h}{2+j} = 0 * emmap{e}{h}{2};
                        end
                    end
                end
                self.holdTime{e} = emmap{e}{e};
                for h=setdiff(1:E,e)
                    self.holdTime{e}{1} = krons(self.holdTime{e}{1},emmap{e}{h}{1});
                    for j = 2:(E+2)
                        self.holdTime{e}{j} = krons(self.holdTime{e}{j},emmap{e}{h}{j});
                        completion_rates = self.holdTime{e}{j}*ones(length(self.holdTime{e}{j}),1);
                        self.holdTime{e}{j} = 0*self.holdTime{e}{j};
                        self.holdTime{e}{j}(:,1) = completion_rates;
                    end
                    self.holdTime{e} = mmap_normalize(self.holdTime{e});
                end
                count_lambda = mmap_count_lambda(self.holdTime{e}); % completiom rates for the different transitions
                Pemb(e,:) = count_lambda/sum(count_lambda);
            end
            self.proc = emmap;
            
            %
            lambda = zeros(1,E);
            A = zeros(E); I=eye(E);
            for e=1:E
                lambda(e) = 1/map_mean(self.holdTime{e});
                for h=1:E
                    A(e,h) = -lambda(e)*(I(e,h)-Pemb(e,h));
                end
            end
            
            if all(lambda>0)
                penv = ctmc_solve(A);
                self.probEnv = penv;
                self.probOrig = zeros(E);
                for e = 1:E
                    for h = 1:E
                        self.probOrig(h,e) = penv(h) * lambda(h) * Pemb(h,e);
                    end
                    if penv(e) > 0
                        self.probOrig(:,e) = self.probOrig(:,e) / sum(self.probOrig(:,e));
                    end
                end
            else
                %T = A(lambda>0, lambda>0), % subgenerator
                %t = A(lambda>0, lambda==0), % exit vector
            end
        end
        
        function env = getEnv(self)
            env = self.env;
        end
        
        function self = setEnv(self, env)
            self.envGraph = env;
        end
        
        function self = setStageName(self, stageId, name)
            self.stageNames{stageId} = name;
        end
                
        function self = setStageType(self, stageId, stageCategory)                      
            if ischar(stageCategory)
                self.stageTypes(stageId) = categorical(stageCategory);
            elseif iscategorical(stageCategory)
                self.stageTypes(stageId) = stageCategory;
            else
                line_error(mfilename,'Stage type must be of type categorical, e.g., categorical("My Semantics").');
            end
        end
        
        function ET = getStageTable(self)
            E = height(self.envGraph.Nodes);
            Stage = [];
            HoldT = {};
            type = categorical([]);
            if isempty(self.probEnv)
                self.init;
            end
            for e=1:E
                Stage(e,1) = e;
                type(e,1) = self.envGraph.Nodes.Type{e};
            end
            Prob = self.probEnv(:);
            Name = categorical(self.envGraph.Nodes.Name(:));
            Model = self.ensemble(:);
            for e=1:E
                HoldT{e,1} = self.holdTime{e};
            end
            Type = type;
            ET = Table(Stage, Name, Type, Prob, HoldT, Model);
        end
    end
    
    methods (Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Ensemble(self);
            % Make a deep copy of each ensemble object
            clone.env = copy(self.envGraph);
        end
    end
end
