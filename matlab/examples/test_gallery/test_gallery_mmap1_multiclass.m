% Test gallery_mmap1_multiclass with MVA
model = gallery_mmap1_multiclass();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
