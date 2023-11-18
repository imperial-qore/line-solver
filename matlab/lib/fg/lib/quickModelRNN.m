function model = quickModel(open, stations, classes, servers, jobs, routing)
%QUICKMODEL Generate simple closed queueing network based on given parameters
    numStations = size(stations,2);
    numClasses = size(classes, 1);
   
    if ~exist('servers','var')
        % servers parameter does not exist, so default it to no concurrency
        servers = ones(numStations);
    end

    if ~exist('jobs','var')
        % jobs parameter does not exist, so default it to 1
        jobs = ones(numClasses);
    end

    model = Network('quickModel');
    if open
        node{1} = Source(model, 'mySource');
        node{numStations+2} = Sink(model, 'mySink');
    end
    
    for i=1:numStations
        node{i} = Queue(model, strcat('QueueStation',num2str(i)), stations{i});
        node{i}.setNumberOfServers(servers(i));
    end

    for i=2:numStations
        node{i}.setState(randsample(1:max(jobs), 1));
    end
    
    P = model.initRoutingMatrix;
    for c=1:numClasses
        if ~open
            jobclass{c} = ClosedClass(model, strcat('Class', num2str(c)), jobs(c), node{1});

            if ~exist('routing','var')
                simpleRoutes = eye(numStations);
                P{jobclass{c}} = simpleRoutes([setdiff(1:size(simpleRoutes,1), [1]), [1]], :);
            else
                P{jobclass{c}} = routing{c};
            end
        else
            jobclass{c} = OpenClass(model, strcat('Class', num2str(c)));
        end
        
        for i=1:numStations
            node{i}.setService(jobclass{c}, Exp(classes(c,i)));
        end
    end

    if open
        P{c} = Network.serialRouting(node);
    end

    model.link(P);
end

