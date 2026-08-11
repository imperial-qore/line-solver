% Test gallery_mdk with MVA
model = gallery_mdk();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
