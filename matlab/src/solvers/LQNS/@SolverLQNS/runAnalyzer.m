function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

tic;
if nargin<2
    options = self.getOptions;
end
dirpath = lineTempName('lqns');
filename = [dirpath,filesep,'model.lqnx'];
self.model.writeXML(filename);

%self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

if options.verbose
    %verbose = '-v';
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
if options.verbose
%    line_printf('\nLQNS model: %s',filename);
    line_printf('\nLQNS command: %s\n',cmd);
end

% Check for remote execution
if isfield(options.config, 'remote') && options.config.remote
    % Remote execution via lqns-rest API
    if options.verbose
        line_printf('\nUsing remote LQNS at: %s\n', options.config.remote_url);
    end
    self.runRemoteLQNS(filename, options);
else
    % Local execution
    system(cmd);
end
self.parseXMLResults(filename);

if ~options.keep
    rmdir(dirpath,'s');
end
runtime = toc;
end
