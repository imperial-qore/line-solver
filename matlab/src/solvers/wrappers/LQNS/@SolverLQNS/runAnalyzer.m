function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

tic;
if nargin<2
    options = self.getOptions;
end
line_ack('LQNS', options.verbose);

% Solver console: this wrapper leaves runAnalyzerChecks commented out (see
% below), so it opens its own run here. The guard must live until this
% function returns.
consoleGuard = LineConsole.beginRun(self, options); %#ok<NASGU>

line_debug(options, 'LQNS: starting (method=%s, multiserver=%s)', options.method, options.config.multiserver);

% see _kb/06-solver-catalog.md (Wrappers: two ways to reach an external binary)
useRemote = isfield(options.config, 'remote') && options.config.remote;

% Refuse a construct the binary cannot model BEFORE the file is written, so the
% answer is LINE's own named refusal rather than lqns reporting a syntax error
% for a modelling limit. See SolverLQNS.unsupportedLNConstructs.
SolverLQNS.assertSupported(self.model, options.method);

dirpath = lineTempName('lqns');
filename = [dirpath,filesep,'model.lqnx'];
LineConsole.step('writing the LQN model to %s', filename);
self.model.writeXML(filename);

%self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

% The binary's advisories and warnings are its own debug channel, not LINE's:
% they stay suppressed up to STD and are let through at DEBUG. Before the
% default verbosity became STD this test was `if options.verbose`, which at
% the old default (false) took the same branch.
if options.verbose == VerboseLevel.DEBUG
    verbose = '';
else
    verbose = '-a -w';
end

multiserver_praqma = '';
switch options.method
    case 'lqsim'
        %no-op
    otherwise
        switch options.config.multiserver
            case 'conway'
                multiserver_praqma='-Pmultiserver=conway';
            case 'rolia'
                multiserver_praqma='-Pmultiserver=rolia';
            case 'zhou'
                multiserver_praqma='-Pmultiserver=zhou';
            case 'suri'
                multiserver_praqma='-Pmultiserver=suri';
            case 'reiser'
                multiserver_praqma='-Pmultiserver=reiser';
            case 'schmidt'
                multiserver_praqma='-Pmultiserver=schmidt';
            case 'default'
                multiserver_praqma='-Pmultiserver=rolia';
        end
end

if isunix
    %                 switch options.method
    %                     case {'default','lqns'}
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' --iteration-limit=',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
    %                     case {'srvn'}
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' --iteration-limit=',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
    %                     case {'exactmva'}
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' --iteration-limit=',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
    %                     case {'srvnexact'}
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' --iteration-limit=',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
    %                     case {'sim','lqsim'}
    %                         cmd=['lqsim ',verbose,' ',multiserver_praqma,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
    %                     case {'lqnsdefault'}
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' -x ',filename]);
    %                     otherwise
    %                         cmd=['lqns ',verbose,' ',multiserver_praqma,' --iteration-limit=',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
    %                 end

    % --iteration-limit seems faulty as of 6.2.27
    if options.verbose
        switch options.method
            case {'default','lqns'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename];
            case {'srvn'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pstop-on-message-loss=false -x ',filename];
            case {'exactmva'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pmva=exact -Pstop-on-message-loss=false -x ',filename];
            case {'srvn.exactmva'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename];
            case {'sim','lqsim'}
                cmd=['lqsim ',verbose,' ',multiserver_praqma,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename];
            case {'lqnsdefault'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -x ',filename];
            otherwise
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename];
        end
    else
        switch options.method
            case {'default','lqns'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
            case {'srvn'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
            case {'exactmva'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pmva=exact -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
            case {'srvn.exactmva'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
            case {'sim','lqsim'}
                cmd=['lqsim ',verbose,' ',multiserver_praqma,' -A ',num2str(options.samples),' -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
            case {'lqnsdefault'}
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -x ',filename,' 2>&1'];
            otherwise
                cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename,' 2>&1'];
        end
    end
else
    switch options.method
        %         case {'default','lqns'}
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
        %         case {'srvn'}
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -i ',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
        %         case {'exactmva'}
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -i ',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
        %         case {'srvnexact'}
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -i ',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
        %         case {'sim','lqsim'}
        %             cmd=['lqsim ',verbose,' ',multiserver_praqma,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
        %         case {'lqnsdefault'}
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -x ',filename]);
        %         otherwise
        %             cmd=['lqns ',verbose,' ',multiserver_praqma,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
        case {'default','lqns'}
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename];
        case {'srvn'}
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pstop-on-message-loss=false -x ',filename];
        case {'exactmva'}
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pmva=exact -Pstop-on-message-loss=false -x ',filename];
        case {'srvn.exactmva'}
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename];
        case {'sim','lqsim'}
            %cmd=['lqsim ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename];
            cmd=['lqsim ',verbose,' ',multiserver_praqma,' -A ',num2str(options.samples),' -Pstop-on-message-loss=false -x ',filename];
        case {'lqnsdefault'}
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -x ',filename];
        otherwise
            cmd=['lqns ',verbose,' ',multiserver_praqma,' -Pstop-on-message-loss=false -x ',filename];
    end
end
if LineConsole.isActive()
    % the console already reports the command through the routed line_debug
    % below, so printing it again would duplicate the line
elseif options.verbose
%    line_printf('\nLQNS model: %s',filename);
    line_printf('\nLQNS command: %s\n',cmd);
end

% Check for remote execution
if useRemote
    line_debug(options, 'LQNS: using remote execution at %s', options.config.remote_url);
    if options.verbose
        line_printf('\nUsing remote LQNS at: %s\n', options.config.remote_url);
    end
    self.runRemoteLQNS(filename, options);
else
    line_debug(options, 'LQNS: using local execution, command: %s', cmd);
    % see _kb/06-solver-catalog.md (Wrappers: LQNS/lqsim LD_LIBRARY_PATH GLIBCXX strip)
    if isunix
        cmd = ['env -u LD_LIBRARY_PATH ', cmd];
    end
    LineConsole.step('running the lqns binary as a subprocess');
    system(cmd);
end
LineConsole.step('parsing the lqns XML results');
self.parseXMLResults(filename);

if ~options.keep
    rmdir(dirpath,'s');
end
runtime = toc;
end
