% Test gallery_merlk with MVA
model = gallery_merlk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
