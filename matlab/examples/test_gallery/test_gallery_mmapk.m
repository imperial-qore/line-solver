% Test gallery_mmapk with MVA
model = gallery_mmapk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
