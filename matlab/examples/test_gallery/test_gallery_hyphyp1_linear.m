% Test gallery_hyphyp1_linear with MVA
model = gallery_hyphyp1_linear();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
