% Test gallery_mhypk with MVA
model = gallery_mhypk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
