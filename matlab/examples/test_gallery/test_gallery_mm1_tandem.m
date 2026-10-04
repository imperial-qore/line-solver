% Test gallery_mm1_tandem with MVA
model = gallery_mm1_tandem();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
