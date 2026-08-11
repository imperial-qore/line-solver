function line_warning(caller, MSG, varargin)
% LINE_WARNING(CALLER, ERRMSG)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~coder.target('MATLAB')
    return;  % No-op in codegen mode
end

%global GlobalConstants.Verbose
persistent lastWarning;
persistent suppressedWarnings;
persistent suppressedWarningsTic;
persistent lastWarningTime;
persistent suppressedAnnouncement;

if GlobalConstants.Verbose == VerboseLevel.SILENT
    return
end

suppressedAnnouncement = false;
errmsg=sprintf(MSG, varargin{:});
w = warning('QUERY','ALL');
w(1).state = 'on'; % always print warnings by default
switch w(1).state
    case 'on'
        %warning('[%s] %s',caller,MSG);
        finalmsg = sprintf('Warning [%s.m]: %s',caller,errmsg);
        % Terminate the warning with a newline, as line_warning_always and the
        % Python/JAR loggers do. line_printf writes the string verbatim, so a
        % call site that does not embed a trailing '\n' leaves the next
        % printout glued to the warning (e.g. 'reference tasks only.AvgTable =').
        if isempty(finalmsg) || finalmsg(end) ~= sprintf('\n')
            finalmsg = [finalmsg, sprintf('\n')];
        end
        try
            % Time since suppression began. The previous form,
            % toc(suppressedWarningsTic)-toc(lastWarningTime), cancels the
            % current time and evaluates to lastWarningTime-suppressedWarningsTic,
            % a fixed gap between two past instants that is ~0 because both are
            % set in the same call. It therefore never exceeded 60 and an
            % identical message stayed suppressed for the whole session rather
            % than for a minute, diverging from the JAR and Python
            % implementations which both measure elapsed time.
            if ~strcmp(finalmsg, lastWarning) || toc(suppressedWarningsTic)>60
                line_printf(finalmsg);
                %warning(finalmsg);
                lastWarning = finalmsg;
                suppressedWarnings = false;                
                suppressedWarningsTic = tic;                
            else
                if ~suppressedWarnings && ~suppressedAnnouncement
                    %line_printf(finalmsg);
                    %warning(finalmsg);
                    finalmsg = sprintf('\nWarning [%s.m]: %s',caller,errmsg);
                    line_printf(sprintf('[%s.m] %s',caller,'Message casted more than once, repetitions will not be printed on screen for 60 seconds.\n'));
                    %warning(sprintf('[%s.m] %s',caller,'Message casted more than once, repetitions will not be printed on screen for 60 seconds.'));
                    suppressedAnnouncement = true;
                    suppressedWarnings = true;
                    suppressedWarningsTic = tic;
                end
            end
            lastWarningTime=tic;
        catch ME
            switch ME.identifier
                case 'MATLAB:toc:callTicFirstNoInputs'
                    %warning(finalmsg);
                    line_printf(finalmsg);
                    lastWarning = finalmsg;
                    suppressedWarnings = false;
                    suppressedWarningsTic = -1;
                    lastWarningTime=tic;
            end
        end
    case 'off'
        %line_printf(sprintf('Warning [%s.m]: %s',caller,errmsg));
end
end
