function [trainedClassifier, validationAccuracy] = Classifier_k(CM,centroids)
% trainClassifier(trainingData)
%  returns a trained K means classifier
%
%   Input:
%       centroids: the kMeans centroids.
%
%   Output:
%       trainedClassifier: a struct containing the trained classifier.
%        The struct contains various fields with information about the
%        trained classifier.
%
%       trainedClassifier.predictFcn: a function to make predictions
%        on new data. It takes an input of the same form as this training
%        code (table or matrix) and returns predictions for the response.
%        If you supply a matrix, include only the predictors columns (or
%        rows).
%
%  Use the code to train the model with new data.
%  To retrain your classifier, call the function from the command line
%  with your original data or new data as the input argument trainingData.
%
%  For example, to retrain a classifier trained with the original data set
%  T, enter:
%    [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
%  To make predictions with the returned 'trainedClassifier' on new data T,
%  use
%    yfit = trainedClassifier.predictFcn(T)
%
%  To automate training the same classifier with new data, or to learn how
%  to programmatically train classifiers, examine the generated code.

centers=centroids{CM};
trainedClassifier.predictFcn = @(x) kmPredictFcn(x , centers);

% Add additional fields to the result struct
trainedClassifier.ColorMap = CM; %%%% added
trainedClassifier.Centroids = centers;

validationAccuracy=1;
end

function labels = kmPredictFcn(d,c)
    % 2 Clusters   
    d1=vecnorm(d-c(1,:),2,2).^2; %Distance of points to cluster 1
    d2=vecnorm(d-c(2,:),2,2).^2; %Distance of points to cluster 2
    labels=1+max((d1-d2)>0,0);   
end
