% Test gallery_hyperl1_reentrant with MVA
model = gallery_hyperl1_reentrant();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
