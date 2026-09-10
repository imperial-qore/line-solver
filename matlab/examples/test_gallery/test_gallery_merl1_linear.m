% Test gallery_merl1_linear with MVA
model = gallery_merl1_linear();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
