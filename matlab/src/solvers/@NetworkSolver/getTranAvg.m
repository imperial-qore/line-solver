function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self,Qt,Ut,Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF,QT,UT,TT)

% Return transient average station metrics
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if nargin == 1
    [Qt,Ut,Tt] = self.getTranHandles;
end
if nargin == 2
    handles = Qt;
    Qt=handles{1};
    Ut=handles{2};
    %Rt=handlers{3};
    Tt=handles{3};
end

QNclass_t={};
UNclass_t={};
%RNclass_t={};
TNclass_t={};

options = self.options;

sn = self.model.getStruct;
minrate = min(sn.rates(isfinite(sn.rates)));
if ~hasTranResults(self)
    if length(options.timespan)<2
        line_error(mfilename,'Timespan option requires to specify a two-dimensional vector, e.g., [0,1e3].\n');
    end
    if isinf(options.timespan(1)) && isinf(options.timespan(2))
        options.timespan = [0,30/minrate];
        line_warning(mfilename,'Timespan of transient analysis unspecified, setting the timespan option to [0, %d]. Use %s(model,''timespan'',[0,T]) to customize.\n',options.timespan(2),class(self));
    end
    if isinf(options.timespan(1))
            line_warning(mfilename,'Start time of transient analysis unspecified, setting the timespan option to [0,%d].\n',options.timespan(2));
        options.timespan(1) = 0;
    end
    if isinf(options.timespan(2))
        options.timespan(2) = 30/minrate;
        line_warning(mfilename,'End time of transient analysis unspecified, setting the timespan option to [%d,%d]. Use %s(model,''timespan'',[0,T]) to customize.\n',options.timespan(1),options.timespan(2),class(self));
    end
    runAnalyzer(self, options);
end

M = sn.nstations;
K = sn.nclasses;
if ~isempty(Qt)
    QNclass_t = cell(M,K);
    UNclass_t = cell(M,K);
    %RNclass_t = cell(M,K);
    TNclass_t = cell(M,K);
    for k=1:K
        for i=1:M
            %%
            if ~Qt{i,k}.disabled && ~isempty(self.result.Tran.Avg.Q)
                ret = self.result.Tran.Avg.Q{i,k};
                metricVal = struct();
                metricVal.handle = {self.model.stations{i}, self.model.classes{k}};
                metricVal.t = ret(:,2);
                metricVal.metric = ret(:,1);
                metricVal.isaggregate = true;
                QNclass_t{i,k} = metricVal;
            else
                ret = NaN;
            end
            
            %%
            if ~Ut{i,k}.disabled && ~isempty(self.result.Tran.Avg.U)
                ret = self.result.Tran.Avg.U{i,k};
                metricVal = struct();
                metricVal.handle = {self.model.stations{i}, self.model.classes{k}};
                metricVal.t = ret(:,2);
                metricVal.metric = ret(:,1);
                metricVal.isaggregate = true;
                UNclass_t{i,k} = metricVal;
            else
                ret = NaN;
            end
            
            %%
            if ~Tt{i,k}.disabled && ~isempty(self.result.Tran.Avg.T)
                ret = self.result.Tran.Avg.T{i,k};
                metricVal = struct();
                metricVal.handle = {self.model.stations{i}, self.model.classes{k}};
                metricVal.t = ret(:,2);
                metricVal.metric = ret(:,1);
                metricVal.isaggregate = true;
                TNclass_t{i,k} = metricVal;
            else
                ret = NaN;
            end
        end
    end
end
end
