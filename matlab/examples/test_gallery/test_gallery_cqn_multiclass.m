% Test gallery_cqn_multiclass with MVA
model = gallery_cqn_multiclass();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
