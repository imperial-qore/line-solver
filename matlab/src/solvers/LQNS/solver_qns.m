function [QN,UN,RN,TN,CN,XN] = solver_qns(sn, options)
% [Q,U,R,T,C,X] = SOLVER_QNS(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

M = sn.nstations; % number of stations
K = sn.nclasses; % number of classes

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

filePath = lineTempName('qrs');
fileName = 'model';
fname = [fileName,'.jmva'];
logFileName = 'console';
logfname = [logFileName,'.out'];
outputFileName = [filePath,filesep,fname];
outputFileName = SolverJMT.writeJMVA(sn, outputFileName, options);
fileName = 'result';
ofname = [fileName,'.jmva'];
resultFileName = [filePath,filesep,ofname];
%stationames={sn.nodenames{sn.isstation}};

switch options.method
    case 'conway'
        options.config.multiserver='conway';
    case 'reiser'
        options.config.multiserver='reiser';
    case 'rolia'
        options.config.multiserver='rolia';
    case 'zhou'
        options.config.multiserver='zhou';
end
if any(sn.nservers>1 & sn.nservers<Inf)
    switch options.config.multiserver
        case {'default','conway'}
            if ispc
                cmd=['qnsolver -l ',outputFileName,' -mconway -o ',resultFileName,' > ',logfname];
            else
                cmd=['qnsolver -l ',outputFileName,' -mconway -o ',resultFileName,' > ',logfname,' 2>&1'];
            end
        case 'reiser'
            if ispc
                cmd=['qnsolver -l ',outputFileName,' -mreiser -o ',resultFileName,' > ',logfname];
            else
                cmd=['qnsolver -l ',outputFileName,' -mreiser -o ',resultFileName,' > ',logfname,' 2>&1'];
            end
        case 'rolia'
            if ispc
                cmd=['qnsolver -l ',outputFileName,' -mrolia -o ',resultFileName,' > ',logfname];
            else
                cmd=['qnsolver -l ',outputFileName,' -mrolia -o ',resultFileName,' > ',logfname,' 2>&1'];
            end
        case 'zhou'
            if ispc
                cmd=['qnsolver -l ',outputFileName,' -mzhou -o ',resultFileName,' > ',logfname];
            else
                cmd=['qnsolver -l ',outputFileName,' -mzhou -o ',resultFileName,' > ',logfname,' 2>&1'];
            end
    end
else
    if ispc
        cmd=['qnsolver -l ',outputFileName,' -o ',resultFileName,' > ',logfname];
    else
        cmd=['qnsolver -l ',outputFileName,' -o ',resultFileName,' > ',logfname,' 2>&1'];
    end
end

if GlobalConstants.Verbose == VerboseLevel.DEBUG
    line_printf('SolverQNS command:\n');
    disp(cmd)
end
system(cmd);
name = cell(sn.nstations*(sn.nchains+1),0);
Uchain = zeros(0,sn.nchains);
Qchain = zeros(0,sn.nchains);
Wchain = zeros(0,sn.nchains);
Tchain = zeros(0,sn.nchains);
[Lchain,STchain,Vchain,alpha,~,~,~] = snGetDemandsChain(sn);
%Lchain = zeros(sn.nstations,sn.nchains); % uncomment for starred output
i = 1;
statName = {};
try
    fid=fopen(resultFileName,'r');
    strline = 1;
    while strline>0
        strline = fgetl(fid);
        if sn.nclasses==1
            [Uchain, Qchain, Wchain, Tchain, statlabel] = parse_dollar_output_singleclass(strline, Uchain, Qchain, Wchain, Tchain);
        else
            [Uchain, Qchain, Wchain, Tchain, statlabel] = parse_dollar_output(strline, Uchain, Qchain, Wchain, Tchain);
        end
        if ~isempty(statlabel)
            statName{end+1} = statlabel;
        end
        %[Lchain, Uchain, Qchain, Wchain, Tchain] = parse_starred_output(strline, Lchain, Uchain, Qchain, Wchain, Tchain);
    end
    fclose(fid);
catch
    line_printf('Failed execution: cannot open the qnsolver output file at: ');
    disp(resultFileName)
    QN = nan(M,K);
    UN = nan(M,K);
    RN = nan(M,K);
    TN = nan(M,K);
    CN = nan(1,K);
    XN = nan(1,K);
    return
end

reorder = [];
for i=1:sn.nstations
    reorder(find(cellfun(@(x)strcmp(x,statName{i}),sn.nodenames)))=i;
end
Uchain = Uchain(reorder,:);
Qchain = Qchain(reorder,:);
Wchain = Wchain(reorder,:);
Tchain = Tchain(reorder,:);

ref= zeros(sn.nchains,1);
for c=1:sn.nchains
    chain = find(sn.chains(c,:));
    Xchain(c)=Tchain(sn.refstat(c),c);
end
%Rchain=Wchain./Vchain; % needs to be reinstated for starred output
Rchain=Wchain;
Rchain(isnan(Rchain))=0;
for i=1:sn.nstations
    if ~isinf(sn.nservers(i))
        Uchain(i,:) = Uchain(i,:) / sn.nservers(i);
    end
end

[QN,UN,RN,TN,CN,XN] = snDeaggregateChainResults(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
end

function [Uchain, Qchain, Wchain, Tchain, statName] = parse_dollar_output(strline, Uchain, Qchain, Wchain, Tchain)
statName = {};
if any(find(strline==','))
    if any(find(strline=='$'))
        % skip
    else
        R = size(Uchain,2);
        %strline
        str=strrep(strline,' ','');
        str=strsplit(str,',');
        Uchain(end+1,1:R) = 0;
        Qchain(end+1,1:R) = 0;
        Wchain(end+1,1:R) = 0;
        Tchain(end+1,1:R) = 0;
        ptr = 1;
        statName{end+1} = str{1};
        for r=1:R
            Qchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1 + R; % skip aggregate value
        for r=1:R
            Wchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1 + R; % skip aggregate value
        for r=1:R
            Uchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1 + R; % skip aggregate value
        for r=1:R
            Tchain(end,r)=str2num(str{ptr+r});
        end
    end
end
%Station, $Q(Chain01), $Q(Chain02), $Q, $R(Chain01), $R(Chain02), $R, $U(Chain01), $U(Chain02), $U, $X(Chain01), $X(Chain02), $X
%Delay1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.63033, 7.40631, 10.0366
%Queue1, 4, 4, 8, 1.52072, 0.54008, 0.797079, 1, 1, 2, 2.63033, 7.40631, 10.0366
end

function [Uchain, Qchain, Wchain, Tchain, statName] = parse_dollar_output_singleclass(strline, Uchain, Qchain, Wchain, Tchain)
statName = {};
if any(find(strline==','))
    if any(find(strline=='$'))
        % skip
    else
        R = size(Uchain,2);
        %strline
        str=strrep(strline,' ','');
        str=strsplit(str,',');
        Uchain(end+1,1:R) = 0;
        Qchain(end+1,1:R) = 0;
        Wchain(end+1,1:R) = 0;
        Tchain(end+1,1:R) = 0;
        ptr = 1;
        statName{end+1} = str{1};
        for r=1:R
            Qchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1;
        for r=1:R
            Wchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1;
        for r=1:R
            Uchain(end,r)=str2num(str{ptr+r});
        end
        ptr = ptr + 1;
        for r=1:R
            Tchain(end,r)=str2num(str{ptr+r});
        end
    end
end
%Station, $Q(Chain01), $Q(Chain02), $Q, $R(Chain01), $R(Chain02), $R, $U(Chain01), $U(Chain02), $U, $X(Chain01), $X(Chain02), $X
%Delay1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.63033, 7.40631, 10.0366
%Queue1, 4, 4, 8, 1.52072, 0.54008, 0.797079, 1, 1, 2, 2.63033, 7.40631, 10.0366
end

function [Lchain, Uchain, Qchain, Wchain, Tchain] = parse_starred_output(strline, Lchain, Uchain, Qchain, Wchain, Tchain)
if any(find(strline=='.'))
    str=strline(2:end);
    str=strrep(str,' ','');
    %[name,service,busypct,custnb,response,thruput]
    str=strsplit(str,'*');
    name = str{1};
    if name(1)=='('
        chainnum = strrep(name,'(Chain','');
        chainnum = str2num(strrep(chainnum,')',''));
        Lchain(i,chainnum) = str2num(str{2});
        Uchain(i,chainnum) = str2num(str{3});
        Qchain(i,chainnum) = str2num(str{4});
        Wchain(i,chainnum) = str2num(str{5});
        Tchain(i,chainnum) = str2num(str{6});
    else
        i=find(cellfun(@any,strfind(stationames,name)));
    end
end
end