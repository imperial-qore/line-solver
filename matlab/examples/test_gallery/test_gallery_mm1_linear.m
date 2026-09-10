% Test gallery_mm1_linear with MVA
model = gallery_mm1_linear();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
