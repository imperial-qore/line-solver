function [lst] = refreshLST(self,statSet,classSet)
% [LT] = REFRESHLAPLST(STATSET,CLASSSET)
% Refresh the Laplace-Stieltjes transforms in the NetworkStruct object

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
    lst = cell(M,1);
    for ist=1:M
        lst{ist,1} = cell(1,K);
    end
elseif nargin==2
    classSet = 1:K;
    lst = cell(M,1);
    for ist=1:M
        lst{ist,1} = cell(1,K);
    end
elseif nargin==3 && isfield(self.sn,'lt')
    % we are only updating selected stations and classes so use the
    % existing ones for the others
    lst = self.sn.lst;
else
    lst = cell(M,1);
    for ist=1:M
        lst{ist,1} = cell(1,K);
    end
end

source_i = self.getIndexSourceStation;
for ist=statSet
    for r=classSet
        if ist == source_i
            if  isa(self.stations{ist}.input.sourceClasses{r}{end},'Disabled')
                lst{ist}{r} = [];
            else
                lst{ist}{r} = @(s) self.stations{ist}.arrivalProcess{r}.evalLST(s);
            end
        else
            switch class(self.stations{ist})
                case {'Fork'}
                    lst{ist}{r} = [];
                case {'Join'}
                    lst{ist}{r} = [];
                otherwise
                    lst{ist}{r} = @(s) self.stations{ist}.serviceProcess{r}.evalLST(s);
            end
        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'mu')
    self.sn.lst = lst;
end
end