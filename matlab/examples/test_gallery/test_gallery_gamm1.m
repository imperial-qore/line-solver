% Test gallery_gamm1 with MVA
model = gallery_gamm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
