function [ algorithm ] = chooseMethod( data )

R = size(data,2)-1;

nbSamples = 0;
arrival = 1;
departure = 1;
responseTime = 1;
throughput = 1;
utilization = 1;

for i = 1:R
    if isempty(data{3,i})
        arrival = 0;
    else
        nbSamples = nbSamples + length(data{3,i});
    end
    if isempty(data{4,i})
        departure = 0;
    end
    if isempty(data{5,i})
        responseTime = 0;
    end
    if isempty(data{6,i})
        throughput = 0;
    end
end

if isempty(data{2,R+1})
    utilization = 0;
end

if arrival && departure && nbSamples < 10^4
    algorithm = 'ci';
    return
end

if arrival && departure && nbSamples > 10^4 && nbSamples < 10^5
    algorithm = 'gql';
    return
end

if arrival && departure && nbSamples > 10^5
    algorithm = 'erps';
    return
end

if throughput && utilization && responseTime
    algorithm = 'ubo';
    return
end

if throughput && utilization
    algorithm = 'ubr';
    return
end

end