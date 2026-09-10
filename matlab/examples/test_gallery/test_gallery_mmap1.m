% Test gallery_mmap1 with MVA
model = gallery_mmap1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
