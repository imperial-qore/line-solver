function ret = tget(AvgTable,station,class)
if GlobalConstants.DummyMode
    ret = [];
    return
end
if ~isstr(station) % inputs are objects
    if nargin==2
        if isa(station,'JobClass') || isa(station,'Chain')
            class = station;
            station=[];
        else
            class=[];
        end
    end
    if any(ismember(AvgTable.Properties.VariableNames,'JobClass'))
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
    else % Chain table
        if isempty(station)
            ret = AvgTable(AvgTable.Chain == class.name,:);
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
                    ret = AvgTable(AvgTable.Station == station.name & AvgTable.Chain == class.name,:);
                case 'Node'
                    ret = AvgTable(AvgTable.Node == station.name & AvgTable.Chain == class.name,:);
            end
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
            if any(ismember(AvgTable.Properties.VariableNames,'JobClass'))
                ret = AvgTable( AvgTable.JobClass == inputstring,:);
            else
                ret = AvgTable( AvgTable.Chain == inputstring,:);
            end
        end
    else
        if any(ismember(AvgTable.Properties.VariableNames,'JobClass'))
            switch AvgTable.Properties.VariableNames{1}
                case 'Station'
                    ret = AvgTable(AvgTable.Station == station & AvgTable.JobClass == class,:);
                case 'Node'
                    ret = AvgTable(AvgTable.Node == station & AvgTable.JobClass == class,:);
            end
        else % Chain table
            switch AvgTable.Properties.VariableNames{1}
                case 'Station'
                    ret = AvgTable(AvgTable.Station == station & AvgTable.Chain == class,:);
                case 'Node'
                    ret = AvgTable(AvgTable.Node == station & AvgTable.Chain == class,:);
            end
        end
    end
end
end