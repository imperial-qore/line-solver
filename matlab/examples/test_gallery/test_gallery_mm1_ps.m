% Test gallery_mm1_ps with MVA
model = gallery_mm1_ps();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
