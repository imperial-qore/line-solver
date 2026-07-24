% Test gallery_mapmk with MVA
model = gallery_mapmk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
