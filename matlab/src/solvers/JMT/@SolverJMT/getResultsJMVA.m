function [result, parsed] = getResultsJMVA(self)
% [RESULT, PARSED] = GETRESULTSJMVA()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

try
    fileName = [getFileName(self),'.jmva-result.jmva'];
    filePath = [getFilePath(self),filesep,fileName];
    if exist(filePath,'file')
        Pref.Str2Num = 'always';
        parsed = xml_read(filePath,Pref);
    else
        line_error(mfilename,'JMT did not output a result file, the analysis has likely failed.');
    end
catch ME
    line_error(mfilename,'Unknown error upon parsing JMT result file. ');
end
self.result.('solver') = getName(self);
self.result.('model') = parsed.ATTRIBUTE;
self.result.('metric') = {};
self.result.('Prob') = struct();
try
    % older JMVA versions do not have the logValue field and will throw an
    % exception
    switch class(parsed.solutions.algorithm.normconst.ATTRIBUTE.logValue)
        case 'double'
            self.result.Prob.logNormConstAggr = parsed.solutions.algorithm.normconst.ATTRIBUTE.logValue;
        otherwise
            self.result.Prob.logNormConstAggr = NaN;
    end
catch
    self.result.Prob.logNormConstAggr = NaN;
end

sn = self.getStruct;
%%%
M = sn.nstations;    %number of stations
S = sn.nservers;
NK = sn.njobs';  % initial population per class
C = sn.nchains;
SCV = sn.scv;

% determine service times
ST = 1./sn.rates;
ST(isnan(sn.rates))=0;
SCV(isnan(SCV))=1;

[~,STchain,Vchain,alpha,~,~] = sn_get_demands_chain(sn);
%%%
statres = parsed.solutions.algorithm.stationresults;

for k=1:sn.nclasses
    for i=1:sn.nstations
        switch sn.nodetype(self.getStruct.stationToNode(i))
            case NodeType.Source
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = sn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = sn.rates(i,k);
                s.('maxSamples') = Inf;
                s.('meanValue') = sn.rates(i,k);
                s.('measureType') = MetricType.toText(MetricType.Tput);
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = sn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = sn.rates(i,k);
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = sn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = MetricType.toText(MetricType.QLen);
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = sn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = sn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = MetricType.toText(MetricType.RespT);
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = sn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
                
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = sn.classnames{k};
                s.('discardedSamples') = 0;
                s.('lowerLimit') = 0;
                s.('maxSamples') = Inf;
                s.('meanValue') = 0;
                s.('measureType') = MetricType.toText(MetricType.Util);
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = sn.nodenames{self.getStruct.stationToNode(i)};
                s.('successful') = 'true';
                s.('upperLimit') = 0;
                self.result.metric{end+1} = s;
        end
    end
end

%%
%Rchain
%Xchain
%Tchain
% for c=1:sn.nchains
%     inchain = sn.inchain{c};
%     for k=inchain(:)'
%         X(k) = Xchain(c) * alpha(sn.refstat(k),k);
%         for i=1:sn.nstations
%             if isinf(S(i))
%                 U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c)) * alpha(i,k);
%             else
%                 U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c)) * alpha(i,k) / S(i);
%             end
%             if Lchain(i,c) > 0
%                 Q(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
%                 T(i,k) = Tchain(i,c) * alpha(i,k);
%                 R(i,k) = Q(i,k) / T(i,k);
%             else
%                 T(i,k) = 0;
%                 R(i,k)=0;
%                 Q(i,k)=0;
%             end
%         end
%         C(k) = sn.njobs(k) / X(k);
%     end
% end
% Q=abs(Q); R=abs(R); X=abs(X); U=abs(U); T=abs(T); C=abs(C);
% T(~isfinite(T))=0; U(~isfinite(U))=0; Q(~isfinite(Q))=0; R(~isfinite(R))=0; X(~isfinite(X))=0; C(~isfinite(C))=0;
%%

for i=1:length(statres)
    classres = statres(i).classresults;
    for c=1:length(classres)
        inchain = sn.inchain{c};
        for m=1:length(classres(c).measure)
            for k=inchain(:)'
                s = struct();
                s.('alfa') = NaN;
                s.('analyzedSamples') = Inf;
                s.('class') = sn.classnames{k};
                s.('discardedSamples') = 0;
                s.('meanValue') = classres(c).measure(m).ATTRIBUTE.meanValue;
                if strcmp(s.meanValue,'NaN')
                    s.meanValue = NaN;
                end
                s.('maxSamples') = Inf;
                s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                switch classres(c).measure(m).ATTRIBUTE.measureType
                    case 'Utilization'
                        if isinf(sn.nservers(i))
                            s.meanValue = ST(i,k) * (s.meanValue / STchain(i,c)) * Vchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                        else
                            s.meanValue = ST(i,k) * (s.meanValue / STchain(i,c)) / Vchain(sn.refstat(k),c) * alpha(i,k) * min(sum(NK(isfinite(NK))), sn.nservers(i)) / sn.nservers(i);
                        end
                        s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                    case 'Throughput'
                        s.meanValue = s.meanValue * alpha(i,k);
                    case 'Number of Customers'
                        % Q(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                        s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                    case 'Residence time'
                        s.('measureType') = 'Response Time';
                        s.meanValue = s.meanValue / sn.visits{c}(i,k); % this is to convert from JMVA's residence into LINE's response time per visit
                        if isinf(sn.nservers(i))
                            s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                        else
                            s.meanValue = s.meanValue * ST(i,k) / STchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                        end
                    otherwise
                        s.('measureType') = classres(c).measure(m).ATTRIBUTE.measureType;
                end
                s.('lowerLimit') = s.meanValue;
                s.('upperLimit') = s.meanValue;
                s.('nodeType') = 'station';
                s.('precision') = Inf;
                s.('station') = statres(i).ATTRIBUTE.station;
                s.('successful') = classres(c).measure(m).ATTRIBUTE.successful;
                self.result.metric{end+1} = s;
            end
        end
    end
end
result = self.result;
end
