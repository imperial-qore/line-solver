% Test gallery_erlm1 with MVA
model = gallery_erlm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
