% Test gallery_mpar1 with MVA
model = gallery_mpar1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
