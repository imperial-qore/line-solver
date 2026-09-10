function ret = linemcr(varargin)
% LINEMCR is the main (wrapper) script of the LINECLI tool
% This function receives the model and solves it
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% arguments
%     NameValueArgs.F
%     NameValueArgs.S
%     NameValueArgs.A
% end
% ext=NameValueArgs.F;
% solver=NameValueArgs.S;
% analysis=NameValueArgs.A;

warning off;
%javaclasspath
%javaaddpath('/opt/matlabruntime/v99/java/jarext/matlab-websocket-1.6.jar');
%system('ls /opt/matlabruntime/v99/java/jarext/')
%system('cat /opt/matlabruntime/v99/toolbox/local/classpath.txt')
javaaddpath(which('matlab-websocket-1.6.jar'))
%javaclasspath

ret = [];
inputext = 'jsim';
solver = 'mva';
analysis = 'all';
outputext = 'json';
file = [];
verbosity = 'normal'; % default
randomSeed = 1+randi(1e5,1);
serverMode = false;
serverPort = 5463;
maxRequests = Inf;
for v=1:2:length(varargin)
    switch varargin{v}
        case {'-p','--port'}
            serverMode = true;
            serverPort = str2double(varargin{v+1});
        case {'-m','--maxreq'}
            maxRequests = str2double(varargin{v+1});
        case {'-s','--solver'}
            solver = varargin{v+1};
        case {'-a','--analysis'}
            analysis = varargin{v+1};
        case {'-f','--file'}
            file = varargin{v+1};
        case {'-v','--verbosity'}
            verbosity = varargin{v+1};
        case {'-i','--input'}
            inputext = varargin{v+1};
        case {'-o','--output'}
            outputext = varargin{v+1};
        case {'-d','--seed'}
            randomSeed = str2double(varargin{v+1});
        case {'-h','--help'}
            fprintf('--------------------------------------------------------------------\n');
            fprintf('LINE Solver - Command Line Interface\n');
            fprintf('Copyright (c) 2012-2026, QORE group, Imperial College London\n');
            fprintf(sprintf('Version %s. All rights reserved.\n',Model('').getVersion()));
            fprintf('--------------------------------------------------------------------\n');
            fprintf('-p, --port     : run in server mode on the specified port \n');
            fprintf('-i, --input    : input file format (jsim*, jsimg, jsimw, json, lqnx, xml, pnml) \n');
            fprintf('-o, --output   : output file format (json*, obj) \n');
            fprintf('-s, --solver   : solver (auto, ctmc, fluid, jmt, mva*, mam, nc, ssa, ldes, qns,\n');
            fprintf('                 env; ln, ln.mva, ln.nc, ln.comom, lqns for layered models) \n');
            fprintf('-d, --seed     : random number seed \n');
            fprintf('-m, --maxreq   : quit after processing the specified number of requests \n');
            fprintf('-a, --analysis : analysis type (all*, avg, sys, chain, node, nodechain,\n');
            fprintf('                 cache, item). The JAR, C++ and python CLIs serve the\n');
            fprintf('                 full set (prob, sample, reward, tran, cdf, sens, ...) \n');
            fprintf('-h, --help     : print help\n');
            fprintf('-v, --verbosity: set verbosity level \n');
            fprintf('    --version  : print version number\n');
            fprintf('\n');
            fprintf('Defaults marked by * \n');
            fprintf('\n');
            fprintf('EXAMPLE: cat myfile.jsimg | docker run -i --rm cli-ubuntu -i jsimg -o json -s mva -a sys\n');
            return
        case {'--version'}
            fprintf(sprintf('%s\n',Model('').getVersion()));
            return
    end
end

if serverMode
    lineserver(serverPort);  % start LINE server mode
    fprintf('--------------------------------------------------------------------\n');
    fprintf('LINE Solver - Command Line Interface\n');
    fprintf('Copyright (c) 2012-2026, QORE group, Imperial College London\n');
    fprintf(sprintf('Version %s. All rights reserved.\n',Model('').getVersion()));
    fprintf('--------------------------------------------------------------------\n');
    fprintf(sprintf('Running in server mode on port %d.\n',serverPort));
    fprintf('Enter Q to stop the server at any time.\n')
    while true
        cmd = input('','s');
        switch cmd
            case {'Q','q'}
                fprintf('Shutting down. Please hold on, it may take several seconds.\n')
                return
        end
    end
end

%persistent lineSplashScreenShown
%lineSplashScreenShown = false;

%setmcruserdata('ParallelProfile', 'lineClusterProfile.settings');

warning off
modelfile = [tempname,'.',inputext];
if isempty(file)
    fid = fopen(modelfile,'w+');
    filecontent = input('', 's');
    while ~isempty(filecontent)
        try
            filecontent = input('', 's');
            fprintf(fid,'%s',filecontent);
        catch
            break
        end
    end
    fclose(fid);
else
    %     switch inputext
    %         case 'zip'
    %             % recursively call solver on each zip content
    %             destFolder = tempname;
    %             unzip(file, destFolder);
    %             D = dir(destFolder);
    %             output = cell(length(D)-2,1);
    %             for d=3:length(D)
    %                 [~,~,dext]=fileparts(D(d).name);
    %                 output{d-2,1}=D(d).name;
    %                 output{d-2,2}=LINECLI('-i',[dext(2:end)],'-f',[destFolder,filesep,D(d).name],'-a',analysis,'-s',solver,'-v','silent','-o','bin');
    %             end
    %             outputMsg = jsonencode(output');
    %             switch outputext
    %                 case 'json'
    %                     if nargout>0
    %                         ret = outputMsg;
    %                     end
    %                     switch verbosity
    %                         case {'silent'}
    %                         case {'normal'}
    %                             fprintf(outputMsg);
    %                             fprintf('\n');
    %                         otherwise
    %                             error('Unknown verbosity level: %s',verbosity);
    %                     end
    %                 case 'bin'
    %                     ret = output;
    %             end
    %             return
    %         otherwise
    copyfile(file,modelfile,'f');
    %end
end
[fpath, name, fileext] = fileparts(modelfile);

%% choose analysis type
% AN UNKNOWN -a IS AN ERROR, NOT A SUBSTITUTION. The `otherwise` arm used to
% fall back to 'all', so `-a cache` or `-a tran-avg` -- both real analyses in
% the other CLIs -- silently produced the average and system tables instead,
% under the name of the analysis the caller asked for. A refusal that names the
% supported set is the only honest answer for an analysis this front end does
% not implement.
wantAvgChainTable = false;
wantAvgNodeTable = false;
wantAvgNodeChainTable = false;
wantAvgCacheTable = false;
wantAvgItemTable = false;
switch analysis
    case 'avg'
        wantAvgSysTable = false;
        wantAvgTable = true;
    case 'sys'
        wantAvgSysTable = true;
        wantAvgTable = false;
    case 'all'
        wantAvgSysTable = true;
        wantAvgTable = true;
    case 'chain'
        wantAvgSysTable = false;
        wantAvgTable = false;
        wantAvgChainTable = true;
    case 'node'
        wantAvgSysTable = false;
        wantAvgTable = false;
        wantAvgNodeTable = true;
    case 'nodechain'
        wantAvgSysTable = false;
        wantAvgTable = false;
        wantAvgNodeChainTable = true;
    case 'cache'
        wantAvgSysTable = false;
        wantAvgTable = false;
        wantAvgCacheTable = true;
    case 'item'
        wantAvgSysTable = false;
        wantAvgTable = false;
        wantAvgItemTable = true;
    otherwise
        error(['Unknown analysis type: %s. This front end serves avg, sys, all, ' ...
               'chain, node, nodechain, cache and item; the JAR (java -jar jline.jar), ' ...
               'the C++ line-cli and the python line CLI serve the full set.'], analysis);
end

%% choose solver
% `fluid` is MATLAB's spelling and `fld` the other CLIs'; `des` and `qnsolver`
% are the root wrapper's. Resolving them here is what lets one command line
% parse against every front end.
switch solver
    case 'fld',      solver = 'fluid';
    case 'des',      solver = 'ldes';
    case 'qnsolver', solver = 'qns';
end
switch fileext
    case {'.jsimg', '.jsimw', '.jsim', '.json'}
        % LINE's own portable model carries its own type, so ONE branch reads a
        % Network, a LayeredNetwork or an Environment and routes each to the
        % solver that can answer for it. `.json` was absent entirely, which made
        % this the only LINE CLI that could not read the format every codebase
        % interchanges through.
        if strcmp(fileext,'.json')
            model = linemodel_load(modelfile);
        else
            model = JMT2LINE(modelfile);
        end
        if isa(model,'LayeredNetwork')
            wantAvgSysTable = false;
            solverObj = linemcr_layered_solver(model, solver, randomSeed);
        elseif isa(model,'Environment')
            % An Environment accepts `env` and nothing else: every Network
            % solver would have to pick a stage and answer about a model the
            % document does not describe.
            if ~any(strcmp(solver,{'env','auto'}))
                error(['An Environment model is solved by -s env; %s solves a Network, ' ...
                       'and this document holds the coupling rather than any one stage.'], solver);
            end
            solverObj = SolverEnv(model, @(m) SolverFluid(m,'seed',randomSeed,'verbose',false), ...
                                  'seed',randomSeed,'verbose',false);
        else
            % EVERY BRANCH ENDS IN A REFUSAL. This switch had no `otherwise`, so
            % `-s mam`, `-s ldes` or `-s auto` left solverObj undefined and the
            % run died on `solverObj.getAvgTable` with "Undefined variable"
            % instead of saying which solvers the document accepts.
            switch solver
                case 'auto'
                    solverObj = SolverAUTO(model,'seed',randomSeed,'verbose',false);
                case 'ctmc'
                    solverObj = SolverCTMC(model,'seed',randomSeed,'verbose',false,'force',true);
                case 'fluid'
                    solverObj = SolverFluid(model,'seed',randomSeed,'verbose',false);
                case 'jmt'
                    solverObj = SolverJMT(model,'seed',randomSeed,'verbose',false);
                case 'mva'
                    solverObj = SolverMVA(model,'seed',randomSeed,'verbose',false);
                case 'mam'
                    solverObj = SolverMAM(model,'seed',randomSeed,'verbose',false);
                case 'nc'
                    solverObj = SolverNC(model,'seed',randomSeed,'verbose',false);
                case 'ssa'
                    solverObj = SolverSSA(model,'method','serial','seed',randomSeed,'verbose',false);
                case 'ldes'
                    solverObj = SolverLDES(model,'seed',randomSeed,'verbose',false);
                case 'qns'
                    solverObj = SolverQNS(model,'seed',randomSeed,'verbose',false);
                otherwise
                    error(['Solver %s cannot solve a %s document; use auto, ctmc, fluid, ' ...
                           'jmt, mva, mam, nc, ssa, ldes or qns.'], solver, fileext);
            end
        end
    case '.pnml'
        % A PNML document is a place/transition net (ISO/IEC 15909-2), so only
        % the solvers whose feature set declares Transition can answer for it;
        % the product-form solvers cannot, and naming one is an error rather
        % than a silent substitution.
        model = pnml_load(modelfile);
        switch solver
            case {'ctmc','auto'}
                solverObj = SolverCTMC(model,'seed',randomSeed,'verbose',false,'force',true);
            case 'jmt'
                solverObj = SolverJMT(model,'seed',randomSeed,'verbose',false);
            case 'ldes'
                solverObj = SolverLDES(model,'seed',randomSeed,'verbose',false);
            case 'ssa'
                solverObj = SolverSSA(model,'method','serial','seed',randomSeed,'verbose',false);
            otherwise
                error('Solver %s cannot solve a place/transition net; use ctmc, jmt, ldes or ssa.', solver);
        end
    case {'.lqnx','.xml'}
        %fprintf('Parsing LQN model: %s\n', modelfile);
        wantAvgSysTable = false;
        model = LQN2MATLAB(modelfile, name);
        solverObj = linemcr_layered_solver(model, solver, randomSeed);
    otherwise
        error(['Unknown model format: %s. This front end reads .json, .jsim/.jsimg/.jsimw, ' ...
               '.lqnx/.xml and .pnml.'], fileext);
end

%% run analyses
output = {};
try
    if wantAvgTable
        %fprintf(sprintf('Saving average performance metrics in:\n%s',[fpath,filesep,name,'_AvgTable.csv\n']));
        AvgTable = solverObj.getAvgTable;
        output{end+1} = AvgTable;
    end
catch ME
    getReport(ME,'basic')
end
try
    if wantAvgSysTable
        %fprintf(sprintf('Saving average performance metrics in:\n%s',[fpath,filesep,name,'_AvgSysTable.csv\n']));
        AvgSysTable = solverObj.getAvgSysTable;
        output{end+1} = AvgSysTable;
    end
catch ME
    getReport(ME,'basic')
end
% The views this front end gained with the -a vocabulary above. Each is one
% getter on the solver, and each was reachable from every other CLI.
tableGetters = { wantAvgChainTable,     @() solverObj.getAvgChainTable; ...
                 wantAvgNodeTable,      @() solverObj.getAvgNodeTable; ...
                 wantAvgNodeChainTable, @() solverObj.getAvgNodeChainTable; ...
                 wantAvgCacheTable,     @() solverObj.getAvgCacheTable; ...
                 wantAvgItemTable,      @() solverObj.getAvgItemTable };
for g = 1:size(tableGetters,1)
    if tableGetters{g,1}
        try
            output{end+1} = tableGetters{g,2}();
        catch ME
            getReport(ME,'basic')
        end
    end
end
outputMsg = jsonencode(output);
switch verbosity
    case {'silent'}
    case {'normal'}
        fprintf(outputMsg);
        fprintf('\n');
    otherwise
        error('Unknown verbosity level: %s',verbosity);
end
%% close
%fprintf('LINE Solver completed successfully.\n');
if nargout>0
    switch outputext
        case 'json'
            ret = outputMsg;
        case 'bin'
            ret = output;
    end
else
    ret = [];
end
return
end

function solverObj = linemcr_layered_solver(model, solver, randomSeed)
% LINEMCR_LAYERED_SOLVER  Layered solver for a token, shared by the .lqnx and
% line-model .json branches.
%
% THE `otherwise` ARM USED TO READ `error('Unknown solver name: %s',solverObj)`
% -- naming a variable that is by definition unassigned at that point, so the
% refusal itself raised "Undefined function or variable 'solverObj'" and the
% caller never saw which solvers an LQN accepts.
switch solver
    case 'lqns'
        solverObj = SolverLQNS(model,'seed',randomSeed,'verbose',false);
    case {'nc','ln.nc','ln.comom'}
        solverObj = SolverLN(model, @(m) SolverNC(m,'method','comom','seed',randomSeed,'verbose',false),'seed',randomSeed,'verbose',false);
    case {'mva','ln','ln.mva','auto'}
        % Bare `ln` is MVA layers, as it is in every other CLI.
        solverObj = SolverLN(model, @(m) SolverMVA(m,'seed',randomSeed,'verbose',false),'seed',randomSeed,'verbose',false);
    case 'ldes'
        solverObj = SolverLDES(model,'seed',randomSeed,'verbose',false);
    otherwise
        error(['Solver %s cannot solve a layered model; use ln, ln.mva, ln.nc, ' ...
               'ln.comom, lqns, mva, nc or ldes.'], solver);
end
end
