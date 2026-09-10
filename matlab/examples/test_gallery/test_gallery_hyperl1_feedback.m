% Test gallery_hyperl1_feedback with MVA
model = gallery_hyperl1_feedback();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
