% Smoke test: exercise every model variant and every solver call used by
% poc_ansim_prio on the smallest configuration, so an API error surfaces in
% seconds rather than in the middle of the full sweep.
addpath('/data/gcasale/line-dev.git/experiments/scratch/poc_ansim');
run('/data/gcasale/line-dev.git/matlab/lineStart.m');

cfg = struct('D',[2,3],'Z',[4,6],'prio',[0,1],'Wdps',100, ...
             'maxIter',50,'tol',1e-6,'law','exp','N',[2,2]);

fprintf('\n--- model variants under exact CTMC ---\n');
for v = {'psprio','dps','hol','prs','fcfs'}
    try
        m = ansim_models(cfg, v{1});
        Q = SolverCTMC(m,'verbose',false,'force',true).getAvgQLen();
        fprintf('%-8s ok   %s\n', v{1}, mat2str(round(Q,4)));
    catch ME
        fprintf('%-8s FAIL %s\n', v{1}, ME.message);
    end
end

fprintf('\n--- competitor solvers ---\n');
calls = { ...
  'FLDmn',     'dps',    @(m) SolverFLD(m,'method','minnormal','verbose',false,'force',true); ...
  'MVAdps',    'dps',    @(m) SolverMVA(m,'verbose',false,'force',true); ...
  'NCdps',     'dps',    @(m) SolverNC(m,'verbose',false,'force',true); ...
  'MVAcl',     'hol',    @(m) SolverMVA(m,'config.np_priority','cl','verbose',false,'force',true); ...
  'MVAshadow', 'hol',    @(m) SolverMVA(m,'config.np_priority','shadow','verbose',false,'force',true); ...
  'MAMhol',    'hol',    @(m) SolverMAM(m,'verbose',false,'force',true); ...
  'MVAprs',    'prs',    @(m) SolverMVA(m,'verbose',false,'force',true); ...
  'MAMprs',    'prs',    @(m) SolverMAM(m,'verbose',false,'force',true); ...
  'MVAmarie',  'fcfs',   @(m) SolverMVA(m,'method','marie','verbose',false,'force',true); ...
  'LDES',      'psprio', @(m) SolverLDES(m,'seed',23000,'samples',1e4,'verbose',false,'force',true); ...
  'SSA',       'psprio', @(m) SolverSSA(m,'seed',23000,'samples',1e4,'verbose',false,'force',true); ...
};
for i = 1:size(calls,1)
    try
        m = ansim_models(cfg, calls{i,2});
        Q = calls{i,3}(m).getAvgQLen();
        fprintf('%-10s (%-6s) ok   %s\n', calls{i,1}, calls{i,2}, mat2str(round(Q,4)));
    catch ME
        fprintf('%-10s (%-6s) FAIL %s\n', calls{i,1}, calls{i,2}, strrep(ME.message,newline,' '));
    end
end

fprintf('\n--- ansim ---\n');
for md = {'ctmcP','live','mn'}
    [Q,t,it,ok,msg] = ansim_solve(cfg, md{1});
    if ok
        fprintf('ansim-%-6s ok   (%.2fs, %d it) %s\n', md{1}, t, it, mat2str(round(Q,4)));
    else
        fprintf('ansim-%-6s FAIL %s\n', md{1}, strrep(msg,newline,' '));
    end
end
fprintf('\nSMOKE DONE\n');
