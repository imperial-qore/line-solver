function [loggerBefore,loggerAfter] = linkNetworkAndLog(self, nodes, classes, P, wantLogger, logPath)% obsolete - old name
% [LOGGERBEFORE,LOGGERAFTER] = LINKNETWORKANDLOG(NODES, CLASSES, P, WANTLOGGER, LOGPATH)% OBSOLETE - OLD NAME

[loggerBefore,loggerAfter] = linkAndLog(self, nodes, classes, P, wantLogger, logPath);
end