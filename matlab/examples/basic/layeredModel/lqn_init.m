clear model solver options AvgTable RD

fprintf(1,'This example illustrates the initialization of LN using the output of LQNS.\n');

cwd = fileparts(which(mfilename));
model = LayeredNetwork.parseXML([cwd,filesep,'lqn_serial.xml']);

%% LQNS, whose solution is what an LN run can be initialized from
options = LQNS.defaultOptions;
options.keep = true; % keep the intermediate XML files of the translation

if SolverLQNS.isAvailable()
    solver{1} = LQNS(model, options);
    AvgTable{1} = solver{1}.getAvgTable();
    fprintf(1,'\nLQNS Results:\n');
    disp(AvgTable{1});
else
    fprintf(1,'\nLQNS solver not available - skipping.\n');
end

%% LN without initialization, for the elapsed time the initialization saves
fprintf(1,'\nSolve with LN without initialization:\n');
solver{2} = LN(model, @(x) MVA(x));
tic;
AvgTable{2} = solver{2}.getAvgTable();
timeElapsed = toc;
disp(AvgTable{2});
fprintf(1,'Time elapsed: %.3fs\n', timeElapsed);

%% CDF of response times, taken on the layer the calls terminate in
fprintf(1,'\nWe now obtain the CDF of response times:\n');
ensemble = model.getEnsemble();
if numel(ensemble) >= 3
    RD = FLD(ensemble{3}).getCdfRespT();
    fprintf(1,'RD (CDF of response times):\n');
    disp(RD);
else
    fprintf(1,'Model ensemble has fewer than three layers.\n');
end
