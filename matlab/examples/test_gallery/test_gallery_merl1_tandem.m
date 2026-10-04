% Test gallery_merl1_tandem with MVA
model = gallery_merl1_tandem();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
