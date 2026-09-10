% Test gallery_erlerl1 with MVA
model = gallery_erlerl1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
