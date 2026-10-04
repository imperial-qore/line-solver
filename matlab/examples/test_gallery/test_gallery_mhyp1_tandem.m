% Test gallery_mhyp1_tandem with MVA
model = gallery_mhyp1_tandem();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
