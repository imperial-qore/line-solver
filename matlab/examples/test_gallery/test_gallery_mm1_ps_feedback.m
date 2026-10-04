% Test gallery_mm1_ps_feedback with MVA
model = gallery_mm1_ps_feedback();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
