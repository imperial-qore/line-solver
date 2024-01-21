function [ph, phases] = refreshServiceRepresentations(self)
% [PH, PHASES] = REFRESHSERVICEREPRESENTATIONS()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.


M = getNumberOfStations(self);
K = getNumberOfClasses(self);
ph = cell(M,1);
for i=1:M
    ph{i,1} = cell(1,K);
end
phases = zeros(M,K);
stations = self.stations;
for i=1:M
    if i == self.getIndexSourceStation
        ph_i = stations{i}.getSourceRates();
    else
        switch class(stations{i})
            case 'Fork'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
                ph_i = Coxian(mu_i,phi_i).getRepres;
            case 'Join'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                    ph_i{r} = Coxian(mu_i{r},phi_i{r}).getRepres;
                end
            otherwise
                ph_i = stations{i}.getServiceRates();
        end
    end
    for r=1:K
        ph{i}{r} = ph_i{r};
        if isempty(ph{i}{r}) % fluid fails otherwise
            phases(i,r) = 1;
        elseif any(isnan(ph{i}{r}{1}(:))) || any(isnan(ph{i}{r}{2}(:))) % disabled
            phases(i,r) = 0;
        else
            phases(i,r) = length(ph{i}{r}{1});
        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'mu')
	proc = ph;
	pie = cell(size(ph));
	for i=1:M
    	for r=1:K
        	map_ir = ph{i}{r};
        	if ~isempty(map_ir)
				%.proc{i}{r} = map_normalize(map_ir);
                proc{i}{r} = map_ir;
				pie{i}{r} = map_pie(map_ir);
    	    else
				pie{i}{r} = NaN;
    	    end
    	end
	end
	self.sn.proc = proc;
	self.sn.pie = pie;
	self.sn.phases = phases;
	self.sn.phasessz = max(self.sn.phases,ones(size(self.sn.phases)));
	self.sn.phasessz(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:)=phases(self.sn.nodeToStation(self.sn.nodetype == NodeType.Join),:);
	self.sn.phaseshift = [zeros(size(phases,1),1),cumsum(self.sn.phasessz,2)];	
end
end
