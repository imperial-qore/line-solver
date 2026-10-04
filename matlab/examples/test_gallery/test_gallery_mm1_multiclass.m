% Test gallery_mm1_multiclass with MVA
model = gallery_mm1_multiclass();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
