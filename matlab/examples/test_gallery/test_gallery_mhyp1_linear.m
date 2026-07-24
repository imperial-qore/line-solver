% Test gallery_mhyp1_linear with MVA
model = gallery_mhyp1_linear();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
