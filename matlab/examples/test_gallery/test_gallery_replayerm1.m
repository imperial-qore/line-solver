% Test gallery_replayerm1 with MVA
model = gallery_replayerm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
