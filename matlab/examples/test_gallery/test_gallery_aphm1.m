% Test gallery_aphm1 with MVA
model = gallery_aphm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
