function ret = tget(AvgTable,station,class)
if GlobalConstants.DummyMode
    ret = [];
    return
end
if ~isstr(station) % inputs are objects
    if nargin==2
        if isa(station,'JobClass')
            class = station;
            station=[];
        else
            class=[];
        end
    end
    if isempty(station)
        ret = AvgTable(AvgTable.JobClass == class.name,:);
    elseif isempty(class)
        switch AvgTable.Properties.VariableNames{1}
            case 'Station'
                ret = AvgTable(AvgTable.Station == station.name,:);
            case 'Node'
                ret = AvgTable(AvgTable.Node == station.name,:);
        end
    else
        switch AvgTable.Properties.VariableNames{1}
            case 'Station'
                ret = AvgTable(AvgTable.Station == station.name & AvgTable.JobClass == class.name,:);
            case 'Node'
                ret = AvgTable(AvgTable.Node == station.name & AvgTable.JobClass == class.name,:);
        end
    end
else % inputs are strings
    inputstring = station;
    if nargin==2
        switch AvgTable.Properties.VariableNames{1}
            case 'Station'
                ret = AvgTable(AvgTable.Station == inputstring,:);
            case 'Node'
                ret = AvgTable(AvgTable.Node == inputstring,:);
        end
        if isempty(ret)
            ret = AvgTable( AvgTable.JobClass == inputstring,:);
        end
    else
        switch AvgTable.Properties.VariableNames{1}
            case 'Station'
                ret = AvgTable(AvgTable.Station == station & AvgTable.JobClass == class,:);
            case 'Node'
                ret = AvgTable(AvgTable.Node == station & AvgTable.JobClass == class,:);
        end
    end
end
end