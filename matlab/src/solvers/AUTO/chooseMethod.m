function solver = chooseMethod(model)
% This function takes as input a QN model defined in LINE and returns
% a Solver object with the predicted method loaded
  
    dataVector = extractFeatures(model);    
    % Add derived features
    dataVector = [dataVector dataVector(:, 1:3) ./ sum(dataVector(:, 1:3), 2)]; % Percent FCFS, PS, Delay
    dataVector = [dataVector logical(dataVector(:, 4))]; % Has CS or not
    dataVector = [dataVector dataVector(:, 5) ./ sum(dataVector(:, 1:2), 2)]; % Avg svrs per Queue
    dataVector = [dataVector dataVector(:, 7) ./ dataVector(:, 6)]; % Num jobs per chain
    dataVector = [dataVector dataVector(:, 8:10) ./ sum(dataVector(:, 8:10), 2)]; % Percent distributions
    
    load('classifier.mat', 'classifier', 'methodNames', 'selected'); 
    if isa(classifier, 'cell')
        chosenMethod = predictEnsemble(classifier, dataVector(selected));
    else
        chosenMethod = predict(classifier, dataVector(selected));
    end
    
    solver = LINE.load(methodNames(chosenMethod), model);
end