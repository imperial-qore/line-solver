function [pred, scores] = predictEnsemble(learners, X)
% This function is a custom predict function for ensembles
% of decision trees trained using either umrbBagging or 
% smoteBagging. "learners" is a cell arry of MATLAB decision
% tree objects, and X is an N x M array, where N is the number
% of test instances, and M is the number of features used to
% train the learners.

    if ~iscolumn(learners)
        learners = learners';
    end   
    preCombinedScores = cell2mat(cellfun(@testPredict, learners, 'UniformOutput', false));
    numLearners = length(learners);
    numRows = size(X, 1);
    scores = zeros(numRows, size(preCombinedScores, 2));
    pred = zeros(numRows, 1);
    for i = 1 : numRows
        scores(i, :) = mean(preCombinedScores(i:numRows:i+(numLearners-1)*numRows, :));
        highestIdx = find(scores(i, :) == max(scores(i, :)));
        pred(i) = highestIdx(randi(length(highestIdx))); % Break ties using random choice
    end
    
    function scores = testPredict(learner)
        [~, scores] = learner.predict(X);
    end
end