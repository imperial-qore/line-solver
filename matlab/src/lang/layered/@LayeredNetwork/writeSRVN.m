function writeSRVN(self, filename)
% WRITESRVN(FILENAME)
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
fid = fopen(filename,'w+');
fprintf(fid,'G\n');
fprintf(fid,['"',self.name,'"\n']);
fprintf(fid,[num2str(Solver.defaultOptions.iter_tol),'\n']);
fprintf(fid,[num2str(Solver.defaultOptions.iter_max),'\n']);
fprintf(fid,['-1\n\n']);

%%
fprintf(fid,'#Host block\n');
fprintf(fid,['P ',num2str(length(self.hosts)),'\n']);
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    switch SchedStrategy.toId(curProc.scheduling)
        case SchedStrategy.ID_INF
            fprintf(fid,['p ',curProc.name,' i','\n']);
        case SchedStrategy.ID_FCFS
            fprintf(fid,['p ',curProc.name,' f m ',num2str(curProc.multiplicity),'\n']);
        case SchedStrategy.ID_PS
            fprintf(fid,['p ',curProc.name,' s ',num2str(curProc.quantum),' m ',num2str(curProc.multiplicity),'\n']);
        case SchedStrategy.ID_SIRO
            fprintf(fid,['p ',curProc.name,' r m ',num2str(curProc.multiplicity),'\n']);
        case SchedStrategy.ID_HOL
            fprintf(fid,['p ',curProc.name,' h m ',num2str(curProc.multiplicity),'\n']);
        otherwise
            line_error(mfilename,'Unsupported scheduling policy.');
    end
end

%%
fprintf(fid,'-1\n\n#Task block\n');
numTasks = 0;
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    numTasks = numTasks + length(curProc.tasks);
end
fprintf(fid,['T ',num2str(numTasks),'\n']);
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            if e == 1
                entryNameList = curEntry.name;
            else
                entryNameList = [entryNameList,' ',curEntry.name];
            end
        end
        switch SchedStrategy.toId(curTask.scheduling)
            case SchedStrategy.ID_REF
                fprintf(fid,['t ',curTask.name,' r ',entryNameList,' -1 ',curProc.name,' z ',num2str(curTask.thinkTimeMean),' m ',num2str(curTask.multiplicity),'\n']);
            case SchedStrategy.ID_INF
                fprintf(fid,['t ',curTask.name,' i ',entryNameList,' -1 ',curProc.name,'\n']);
            case SchedStrategy.ID_FCFS
                fprintf(fid,['t ',curTask.name,' f ',entryNameList,' -1 ',curProc.name,' m ',num2str(curTask.multiplicity),'\n']);
            case SchedStrategy.ID_HOL
                fprintf(fid,['t ',curTask.name,' h ',entryNameList,' -1 ',curProc.name,' m ',num2str(curTask.multiplicity),'\n']);
            otherwise
                line_error(mfilename,'Unsupported scheduling policy.');
        end
    end
end

%%
fprintf(fid,'-1\n\n#Entry block\n');
numEntry = 0;
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        numEntry = numEntry + length(curTask.entries);
    end
end
fprintf(fid,['E ',num2str(numEntry),'\n']);
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            for a=1:length(curTask.activities)
                curAct = curTask.activities(a);
                if strcmp(curAct.boundToEntry,curEntry.name)
                    fprintf(fid,['A ',curEntry.name,' ',curAct.name,'\n']);
                    break
                end
            end
        end
    end
end

%%
fprintf(fid,'-1\n\n#Activity blocks\n');
for p=1:length(self.hosts)
    curProc = self.hosts(p);
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        fprintf(fid,['A ',curTask.name,'\n']);
        for a=1:length(curTask.activities)
            curAct = curTask.activities(a);
            fprintf(fid,['s ',curAct.name,' ',num2str(curAct.hostDemandMean),'\n']);
            fprintf(fid,['c ',curAct.name,' ',num2str(curAct.hostDemandSCV),'\n']);
            if strcmp(curAct.callOrder,'DETERMINISTIC')
                fprintf(fid,['f ',curAct.name,' 1\n']);
            else
                fprintf(fid,['f ',curAct.name,' 0\n']);
            end
            for d=1:numel(curAct.syncCallDests)
                fprintf(fid,['y ',curAct.name,' ',curAct.syncCallDests{d},' ',num2str(curAct.syncCallMeans),'\n']);
            end
            for d=1:numel(curAct.asyncCallDests)
                fprintf(fid,['z ',curAct.name,' ',curAct.asyncCallDests{d},' ',num2str(curAct.asyncCallMeans),'\n']);
            end
        end
        buffer = '';
        preReplyActNames = cell(0);
        for ap=1:length(curTask.precedences)
            curActPrec = curTask.precedences(ap);
            preActFields = curActPrec.preActs;
            for a=1:length(preActFields)
                for e=1:length(curTask.entries)
                    curEntry = curTask.entries(e);
                    if any(strcmp(preActFields{a},curEntry.replyActivity))
                        preReplyActNames{end+1,1} = preActFields{a};
                        preActFields{a} = [preActFields{a},'[',curEntry.name,']'];
                        break
                    end
                end
            end
            switch curActPrec.preType
                case ActivityPrecedenceType.PRE_SEQ
                    preActDelim = '';
                case ActivityPrecedenceType.PRE_AND
                    preActDelim = ' & ';
                case ActivityPrecedenceType.PRE_OR
                    preActDelim = ' + ';
                otherwise
                    line_error(mfilename,'Unsupported activity precedence.');
            end
            postActFields = curActPrec.postActs;
            switch curActPrec.postType
                case ActivityPrecedenceType.POST_SEQ
                    postActDelim = '';
                case ActivityPrecedenceType.POST_AND
                    postActDelim = ' & ';
                case ActivityPrecedenceType.POST_OR
                    for a=1:length(postActFields)
                        postActFields{a} = ['(',num2str(curActPrec.postParams(a)),') ',postActFields{a}];
                    end
                    postActDelim = ' + ';
                case ActivityPrecedenceType.POST_LOOP
                    for a=1:length(postActFields)-1
                        postActFields{a} = [num2str(curActPrec.postParams(a)),' * ',postActFields{a}];
                    end
                    postActDelim = ', ';
                otherwise
                    line_error(mfilename,'Unsupported activity precedence.');
            end
            prePrecSeg = join(preActFields,preActDelim);
            prePrecSeg = prePrecSeg{1};
            if strcmp(curActPrec.preType,ActivityPrecedenceType.PRE_AND) && ~isempty(curActPrec.preParams)
                prePrecSeg = [prePrecSeg,' (',num2str(curActPrec.preParams(1)),')'];
            end
            postPrecSeg = join(postActFields,postActDelim);
            postPrecSeg = postPrecSeg{1};
            buffer = [buffer,prePrecSeg,' -> ',postPrecSeg,';\n'];
        end
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            for r=1:length(curEntry.replyActivity)
                replyActName = curEntry.replyActivity{r};
                if ~any(strcmp(replyActName,preReplyActNames))
                    buffer = [buffer,replyActName,'[',curEntry.name,'];\n'];
                end
            end
        end
        if ~isempty(buffer)
            fprintf(fid,[':\n',buffer(1:end-3),'\n']);
        end
        fprintf(fid,['-1\n\n']);
    end
end

if fid~=1
    fclose(fid);
end
end
