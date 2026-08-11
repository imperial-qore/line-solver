% Test gallery_erlm1ps with MVA
model = gallery_erlm1ps();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
