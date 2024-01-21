function [lst] = refreshLST(self,statSet,classSet)
% [LT] = REFRESHLAPLST(STATSET,CLASSSET)
% Refresh the Laplace-Stieltjes transforms in the NetworkStruct object

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
    lst = cell(M,1);
    for i=1:M
        lst{i,1} = cell(1,K);
    end
elseif nargin==2
    classSet = 1:K;
    lst = cell(M,1);
    for i=1:M
        lst{i,1} = cell(1,K);
    end
elseif nargin==3 && isfield(self.sn,'lt')
    % we are only updating selected stations and classes so use the
    % existing ones for the others
    lst = self.sn.lst;
else
    lst = cell(M,1);
    for i=1:M
        lst{i,1} = cell(1,K);
    end
end

source_i = self.getIndexSourceStation;
for i=statSet
    for r=classSet
        if i == source_i
            if  isa(self.stations{i}.input.sourceClasses{r}{end},'Disabled')
                lst{i}{r} = [];
            else
                lst{i}{r} = @(s) self.stations{i}.arrivalProcess{r}.evalLST(s);
            end
        else
            switch class(self.stations{i})
                case {'Fork'}
                    lst{i}{r} = [];
                case {'Join'}
                    lst{i}{r} = [];
                otherwise
                    lst{i}{r} = @(s) self.stations{i}.serviceProcess{r}.evalLST(s);
            end
        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'mu')
    self.sn.lst = lst;
end
end