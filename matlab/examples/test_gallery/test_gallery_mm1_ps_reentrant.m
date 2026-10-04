% Test gallery_mm1_ps_reentrant with MVA
model = gallery_mm1_ps_reentrant();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
