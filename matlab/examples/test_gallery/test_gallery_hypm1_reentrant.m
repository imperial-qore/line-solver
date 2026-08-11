% Test gallery_hypm1_reentrant with MVA
model = gallery_hypm1_reentrant();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
