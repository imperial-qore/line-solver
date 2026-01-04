clear solver AvgTable;

n=5; m=[2,1]; alpha=1.0;
replStrat = {ReplacementStrategy.RR, ReplacementStrategy.FIFO, ReplacementStrategy.LRU};

AvgTable={};
solver={};
%fprintf(1,'replc-strategy, ctmc, mva, nc:\n');
for s = 1:length(replStrat)
    model = Network('model');
    source = Source(model, 'Source');
    cacheNode = Cache(model, 'Cache', n, m, replStrat{s});
    sink = Sink(model, 'Sink');

    jobClass = OpenClass(model, 'InitClass', 0);
    hitClass = OpenClass(model, 'HitClass', 0);
    missClass = OpenClass(model, 'MissClass', 0);

    source.setArrival(jobClass, Exp(2));

    pAccess = Zipf(alpha, n);  % Zipf-like item references
    cacheNode.setRead(jobClass, pAccess);

    cacheNode.setHitClass(jobClass, hitClass);
    cacheNode.setMissClass(jobClass, missClass);

    P = model.initRoutingMatrix;

    P{jobClass, jobClass}(source, cacheNode) =  1.0;
    P{hitClass, hitClass}(cacheNode, sink) =  1.0;
    P{missClass, missClass}(cacheNode, sink) =  1.0;

    model.link(P);

    solver{end+1} = CTMC(model,'keep',false,'cutoff',1,'verbose',false);
    AvgTable{end+1} = solver{end}.getAvgNodeTable;
    if ~isempty(cacheNode.getHitRatio)
        ctmcHitRatio(s) = cacheNode.getHitRatio;
    else
        ctmcHitRatio(s) = NaN;
    end

    if MVA.supports(model)
        solver{end+1} = MVA(model,'verbose',false);
        model.reset; AvgTable{end+1} = solver{end}.getAvgNodeTable;
        if ~isempty(cacheNode.getHitRatio)
            mvaHitRatio(s) = cacheNode.getHitRatio;
        else
            mvaHitRatio(s) = NaN;
        end
    end

    if NC.supports(model)
        solver{end+1} = NC(model,'verbose',false);
        model.reset; AvgTable{end+1} = solver{end}.getAvgNodeTable;
        if ~isempty(cacheNode.getHitRatio)
            ncHitRatio(s) = cacheNode.getHitRatio;
        else
            ncHitRatio(s) = NaN;
        end
    else
        ncHitRatio(s) = NaN;
    end

    fprintf(1,'%s: %.8f, %.8f, %.8f\n',ReplacementStrategy.toString(replStrat{s}),ctmcHitRatio(s),mvaHitRatio(s),ncHitRatio(s));
end
