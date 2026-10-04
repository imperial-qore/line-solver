% Test gallery_merl1 with MVA
model = gallery_merl1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
