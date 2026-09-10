% Test gallery_hyperlk with MVA
model = gallery_hyperlk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
