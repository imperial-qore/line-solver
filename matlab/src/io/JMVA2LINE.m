function model = JMVA2LINE(filename,modelName)
% MODEL = JMVA2LINE(FILENAME,MODELNAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
T0=tic;
% import model
Pref.Str2Num = 'always';
xDoc = xml_read(filename,Pref);
[~,fname]=fileparts(filename);
try
    xDoc = xDoc.sim;
end

% create network
if nargin<2
    modelName = fname;
end
model = Network(modelName);

%%
nodes = {};

% create stations
nInfStations = 0;
if isfield(xDoc.parameters.stations,'delaystation')
    nInfStations = length(xDoc.parameters.stations.delaystation);
    for i=1:nInfStations
        nodes{end+1} = Delay(model, xDoc.parameters.stations.delaystation(i).ATTRIBUTE.name);
    end
end

nLIStations = 0;
if isfield(xDoc.parameters.stations,'listation')
    nLIStations = length(xDoc.parameters.stations.listation);
    for i=1:nLIStations
        nodes{end+1} = Queue(model, xDoc.parameters.stations.listation(i).ATTRIBUTE.name, SchedStrategy.PS);
    end
end

nLDStations = 0;
if isfield(xDoc.parameters.stations,'ldstation')
    nLDStations = length(xDoc.parameters.stations.ldstation);
    for i=1:nLDStations
        nServers = 1;
        if isfield(xDoc.parameters.stations.ldstation(i).ATTRIBUTE, 'servers')
            nServers = xDoc.parameters.stations.ldstation(i).ATTRIBUTE.servers;
        end
        nodes{end+1} = Queue(model, xDoc.parameters.stations.ldstation(i).ATTRIBUTE.name, SchedStrategy.PS);
        nodes{end}.setNumberOfServers(nServers);
    end
end

%%
jobclass = {};
visited = false(0);
nOpenClasses = 0;
% create open classes
if isfield(xDoc.parameters.classes,'openclass')
    source = Source(model, 'Source');
    sink = Sink(model, 'Sink');
    nOpenClasses = length(xDoc.parameters.classes.openclass);
    for r=1:nOpenClasses
        name = xDoc.parameters.classes.openclass(r).ATTRIBUTE.name;
        jobclass{r} = OpenClass(model, name, 0);
        rate = xDoc.parameters.classes.openclass(r).ATTRIBUTE.rate;
        source.setArrival(jobclass{r}, Exp(rate));
    end
end

% create closed classes
nClosedClasses = 0;
if isfield(xDoc.parameters.classes,'closedclass')
    nClosedClasses = length(xDoc.parameters.classes.closedclass);
    for r=1:nClosedClasses
        name = xDoc.parameters.classes.closedclass(r).ATTRIBUTE.name;
        population = xDoc.parameters.classes.closedclass(r).ATTRIBUTE.population;
        if isfield(xDoc.parameters,'ReferenceStation')
            refstat = model.getStationByName(xDoc.parameters.ReferenceStation.Class(nOpenClasses+r).ATTRIBUTE.refStation);
        end
        jobclass{nOpenClasses + r} = ClosedClass(model, name, population, nodes{1}, 0);
    end
end

for i=1:nInfStations
    for r=1:(nOpenClasses+nClosedClasses)
        servTime = xDoc.parameters.stations.delaystation(i).servicetimes.servicetime(r).CONTENT;
        visits = xDoc.parameters.stations.delaystation(i).visits.visit(r).CONTENT;
        servClass = model.getClassByName(xDoc.parameters.stations.delaystation(i).servicetimes.servicetime(r).ATTRIBUTE.customerclass);
        demand = servTime*visits;
        if visits > 0
            nodes{i}.setService(servClass,Exp(1/demand));
            visited(i,servClass) = true;
        else
            nodes{i}.setService(servClass,Disabled.getInstance());
            visited(i,servClass) = false;
        end
    end
end
for i=1:nLIStations
    for r=1:(nOpenClasses+nClosedClasses)
        servTime = xDoc.parameters.stations.listation(i).servicetimes.servicetime(r).CONTENT;
        visits = xDoc.parameters.stations.listation(i).visits.visit(r).CONTENT;
        servClass = model.getClassByName(xDoc.parameters.stations.listation(i).servicetimes.servicetime(r).ATTRIBUTE.customerclass);
        demand = servTime*visits;
        if visits > 0
            nodes{nInfStations+i}.setService(servClass,Exp(1/demand));
            visited(nInfStations + i,servClass) = true;
        else
            nodes{nInfStations+i}.setService(servClass,Disabled.getInstance());
            visited(nInfStations + i,servClass) = false;
        end
    end
end


for i=1:nLDStations
    ldst = xDoc.parameters.stations.ldstation(i);
    nodeIdx = nInfStations + nLIStations + i;
    for r=1:(nOpenClasses+nClosedClasses)
        visits = ldst.visits.visit(r).CONTENT;
        servClass = model.getClassByName(ldst.servicetimes.servicetimes(r).ATTRIBUTE.customerclass);
        % Parse per-population service demands from LD station
        stData = ldst.servicetimes.servicetimes(r).servicetime;
        if ischar(stData) || isstring(stData)
            % Semicolon or comma-separated format: "d(0);d(1);...;d(N)"
            if contains(char(stData), ';')
                demands = str2double(strsplit(char(stData), ';'));
            else
                demands = str2double(strsplit(char(stData), ','));
            end
        elseif isnumeric(stData)
            demands = stData(:)';
        else
            % Array of structured entries
            demands = arrayfun(@(x) x.CONTENT, stData);
        end
        % Remove NaN entries
        demands = demands(~isnan(demands));
        % demands(1) is at population 0, demands(2) at population 1, etc.
        % Skip zero-population entry if present
        if length(demands) > 1 && demands(1) == 0
            demands = demands(2:end); % now demands(n) = demand at population n
        end
        if visits > 0 && ~isempty(demands) && demands(1) > 0
            baseDemand = demands(1) * visits;
            nodes{nodeIdx}.setService(servClass, Exp(1/baseDemand));
            visited(nodeIdx, servClass) = true;
            % Set load-dependent scaling: alpha(n) = baseDemand / (demands(n)*visits)
            if length(demands) > 1
                alpha = demands(1) ./ demands;
                alpha(demands <= 0) = 0;
                nodes{nodeIdx}.setLoadDependence(alpha);
            end
        else
            nodes{nodeIdx}.setService(servClass, Disabled.getInstance());
            visited(nodeIdx, servClass) = false;
        end
    end
end

P = model.initRoutingMatrix;
for r=1:nOpenClasses
    P{r,r} = Network.serialRouting({source,nodes{find(visited(:,r))},sink});
end
for r=1:nClosedClasses
    P{nOpenClasses+r,nOpenClasses+r} = Network.serialRouting(nodes{find(visited(:,nOpenClasses+r))});
end
model.link(P);
tot=toc(T0);
%line_printf(['JMT2LINE parsing time: ',num2str(Ttot),' s\n']);

end
