% Test gallery_mm1_reentrant with MVA
model = gallery_mm1_reentrant();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
