% Test gallery_detm1 with MVA
model = gallery_detm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
