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
            fprintf('-i, --input    : input file format (jsim*, jsimg, jsimw, lqnx, xml) \n');
            fprintf('-o, --output   : output file format (json*, obj) \n');
            fprintf('-s, --solver   : solver (ctmc, fluid, jmt, mva*, nc, ssa, ln) \n');
            fprintf('-d, --seed     : random number seed \n');
            fprintf('-m, --maxreq   : quit after processing the specified number of requests \n');
            fprintf('-a, --analysis : analysis type (all*, sys, mva) \n');
            fprintf('-h, --help     : print help\n');
            fprintf('-v, --version  : print version number\n');
            fprintf('\n');
            fprintf('Defaults marked by * \n');
            fprintf('\n');
            fprintf('EXAMPLE: cat myfile.jsimg | docker run -i --rm cli-ubuntu -i jsimg -o json -s mva -a sys\n');
            return
        case {'-v','--version'}
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
    otherwise
        wantAvgSysTable = true;
        wantAvgTable = true;
end

%% choose solver
switch fileext
    case {'.jsimg', '.jsimw', '.jsim'}
        %fprintf(1,'Parsing JSIMgraph model: %s\n', modelfile);
        model = JMT2LINE(modelfile);
        %fprintf(1,'Model %s successfully parsed.\n', model.getName);
        switch solver
            case 'ctmc'
                M = model.getNumberOfStations;
                K = model.getNumberOfClasses;
                solverObj = SolverCTMC(model,'seed',randomSeed,'verbose',false,'force',true);
            case 'fluid'
                solverObj = SolverFluid(model,'seed',randomSeed,'verbose',false);
            case 'jmt'
                solverObj = SolverJMT(model,'seed',randomSeed,'verbose',false);
            case 'mva'
                solverObj = SolverMVA(model,'seed',randomSeed,'verbose',false);
            case 'nc'
                solverObj = SolverNC(model,'seed',randomSeed,'verbose',false);
            case 'ssa'
                solverObj = SolverSSA(model,'method','serial','seed',randomSeed,'verbose',false);
        end
    case {'.lqnx','.xml'}
        %fprintf('Parsing LQN model: %s\n', modelfile);
        wantAvgSysTable = false;
        model = LQN2MATLAB(modelfile, name);
        %fprintf('Model %s successfully parsed.\n', model.getName);
        switch solver
            case 'lqns'
                solverObj = SolverLQNS(model,'seed',randomSeed,'verbose',false);
            case {'nc','ln','ln.comom'}
                solverObj = SolverLN(model, @(m) SolverNC(m,'method','comom','seed',randomSeed,'verbose',false),'seed',randomSeed,'verbose',false);
            case {'mva','ln.mva'}
                solverObj = SolverLN(model, @(m) SolverMVA(m,'seed',randomSeed,'verbose',false),'seed',randomSeed,'verbose',false);
            otherwise
                error('Unknown solver name: %s',solverObj);
        end
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
