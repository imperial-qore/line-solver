% Test gallery_mm1_ps_multiclass with MVA
model = gallery_mm1_ps_multiclass();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
