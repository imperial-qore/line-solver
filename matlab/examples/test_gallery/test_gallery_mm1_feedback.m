% Test gallery_mm1_feedback with MVA
model = gallery_mm1_feedback();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
