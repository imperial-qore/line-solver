function [QNclass_t, UNclass_t, TNclass_t] = getTranAvg(self,Qt,Ut,Tt)
% [QNCLASS_T, UNCLASS_T, TNCLASS_T] = GETTRANAVG(SELF,QT,UT,TT)

% Return transient average station metrics
%
% Copyright (c) 2012-2026, Imperial College London
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

% Check if transient results exist
hasTran = isfield(self.result, 'Tran') && isfield(self.result.Tran, 'Avg');
hasQ = hasTran && isfield(self.result.Tran.Avg, 'Q') && ~isempty(self.result.Tran.Avg.Q);
hasU = hasTran && isfield(self.result.Tran.Avg, 'U') && ~isempty(self.result.Tran.Avg.U);
hasT = hasTran && isfield(self.result.Tran.Avg, 'T') && ~isempty(self.result.Tran.Avg.T);

if ~isempty(Qt)
    QNclass_t = cell(M,K);
    UNclass_t = cell(M,K);
    %RNclass_t = cell(M,K);
    TNclass_t = cell(M,K);
    for k=1:K
        for ist=1:M
            %%
            if ~Qt{ist,k}.disabled && hasQ
                ret = self.result.Tran.Avg.Q{ist,k};
                if ~isempty(ret) && size(ret,2) >= 2
                    metricVal = struct();
                    metricVal.handle = {self.model.stations{ist}, self.model.classes{k}};
                    metricVal.t = ret(:,2);
                    metricVal.metric = ret(:,1);
                    metricVal.isaggregate = true;
                    QNclass_t{ist,k} = metricVal;
                else
                    QNclass_t{ist,k} = NaN;
                end
            else
                QNclass_t{ist,k} = NaN;
            end

            %%
            if ~Ut{ist,k}.disabled && hasU
                ret = self.result.Tran.Avg.U{ist,k};
                if ~isempty(ret) && size(ret,2) >= 2
                    metricVal = struct();
                    metricVal.handle = {self.model.stations{ist}, self.model.classes{k}};
                    metricVal.t = ret(:,2);
                    metricVal.metric = ret(:,1);
                    metricVal.isaggregate = true;
                    UNclass_t{ist,k} = metricVal;
                else
                    UNclass_t{ist,k} = NaN;
                end
            else
                UNclass_t{ist,k} = NaN;
            end

            %%
            if ~Tt{ist,k}.disabled && hasT
                ret = self.result.Tran.Avg.T{ist,k};
                if ~isempty(ret) && size(ret,2) >= 2
                    metricVal = struct();
                    metricVal.handle = {self.model.stations{ist}, self.model.classes{k}};
                    metricVal.t = ret(:,2);
                    metricVal.metric = ret(:,1);
                    metricVal.isaggregate = true;
                    TNclass_t{ist,k} = metricVal;
                else
                    TNclass_t{ist,k} = NaN;
                end
            else
                TNclass_t{ist,k} = NaN;
            end
        end
    end
end
end
