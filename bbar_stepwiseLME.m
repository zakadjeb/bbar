function [bestModel, selectionHistory] = bbar_stepwiseLME(data, responseVar, candidatePredictors, randomEffects, varargin)
% STEPWISELME Performs stepwise selection for Linear Mixed Effects models
%
% INPUTS:
%   data - table containing all variables
%   responseVar - string, name of response variable
%   candidatePredictors - cell array of predictor variable names
%   randomEffects - string, random effects specification (e.g., '(1|Subject)')
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'Direction' - 'forward', 'backward', or 'both' (default: 'both')
%   'Criterion' - 'AIC', 'BIC', or 'LRT' (default: 'AIC')
%   'Alpha' - significance level for LRT (default: 0.05)
%   'MaxIterations' - maximum number of iterations (default: 100)
%   'Verbose' - logical, display progress (default: true)
%   'StartingModel' - cell array of predictors to start with (default: {})
%
% OUTPUTS:
%   bestModel - final selected fitlme model
%   selectionHistory - structure containing selection process details
%
% Zakaria Djebbara, 08-2025

% Parse input arguments
p = inputParser;
addRequired(p, 'data', @istable);
addRequired(p, 'responseVar', @ischar);
addRequired(p, 'candidatePredictors', @iscell);
addRequired(p, 'randomEffects', @ischar);
addParameter(p, 'Direction', 'both', @(x) ismember(x, {'forward', 'backward', 'both'}));
addParameter(p, 'Criterion', 'AIC', @(x) ismember(x, {'AIC', 'BIC', 'LRT'}));
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && x > 0 && x < 1);
addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && x > 0);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'StartingModel', {}, @iscell);

parse(p, data, responseVar, candidatePredictors, randomEffects, varargin{:});

direction = p.Results.Direction;
criterion = lower(p.Results.Criterion);
alpha = p.Results.Alpha;
maxIter = p.Results.MaxIterations;
verbose = p.Results.Verbose;
startingModel = p.Results.StartingModel;

% Initialize
currentPredictors = startingModel;
remainingPredictors = setdiff(candidatePredictors, currentPredictors);
iteration = 0;
improved = true;
selectionHistory = struct();

if verbose
    fprintf('Starting stepwise LME selection...\n');
    fprintf('Direction: %s, Criterion: %s\n', direction, criterion);
    fprintf('Candidate predictors: %s\n', strjoin(candidatePredictors, ', '));
end

% Main stepwise loop
while improved && iteration < maxIter
    iteration = iteration + 1;
    improved = false;
    
    if verbose
        fprintf('\n--- Iteration %d ---\n', iteration);
        fprintf('Current model: %s\n', formatModelString(currentPredictors, randomEffects));
    end
    
    % Current model performance
    currentModel = fitCurrentModel(data, responseVar, currentPredictors, randomEffects);
    currentScore = getModelScore(currentModel, criterion);
    
    bestScore = currentScore;
    bestAction = 'none';
    bestPredictor = '';
    
    % Test all possible actions and find the best one
    candidateActions = [];
    
    % Forward step: try adding each remaining predictor
    if ismember(direction, {'forward', 'both'}) && ~isempty(remainingPredictors)
        if verbose
            fprintf('Testing forward additions...\n');
        end
        
        for i = 1:length(remainingPredictors)
            testPredictors = [currentPredictors, remainingPredictors{i}];
            
            try
                testModel = fitCurrentModel(data, responseVar, testPredictors, randomEffects);
                testScore = getModelScore(testModel, criterion);
                
                % Store candidate action
                candidateActions(end+1).action = 'add';
                candidateActions(end).predictor = remainingPredictors{i};
                candidateActions(end).score = testScore;
                candidateActions(end).model = testModel;
                
                if verbose
                    fprintf('  + %s: %s = %.4f\n', remainingPredictors{i}, criterion, testScore);
                end
                
            catch ME
                if verbose
                    fprintf('  + %s: FAILED (%s)\n', remainingPredictors{i}, ME.message);
                end
            end
        end
    end
    
    % Backward step: try removing each current predictor
    if ismember(direction, {'backward', 'both'}) && ~isempty(currentPredictors)
        if verbose
            fprintf('Testing backward removals...\n');
        end
        
        for i = 1:length(currentPredictors)
            testPredictors = currentPredictors;
            testPredictors(i) = [];
            
            try
                if isempty(testPredictors)
                    % Intercept-only model
                    testModel = fitCurrentModel(data, responseVar, {}, randomEffects);
                else
                    testModel = fitCurrentModel(data, responseVar, testPredictors, randomEffects);
                end
                testScore = getModelScore(testModel, criterion);
                
                % Store candidate action
                candidateActions(end+1).action = 'remove';
                candidateActions(end).predictor = currentPredictors{i};
                candidateActions(end).score = testScore;
                candidateActions(end).model = testModel;
                
                if verbose
                    fprintf('  - %s: %s = %.4f\n', currentPredictors{i}, criterion, testScore);
                end
                
            catch ME
                if verbose
                    fprintf('  - %s: FAILED (%s)\n', currentPredictors{i}, ME.message);
                end
            end
        end
    end
    
    % Find the best action among all candidates
    bestAction = 'none';
    bestPredictor = '';
    
    for i = 1:length(candidateActions)
        candidate = candidateActions(i);
        if isImprovement(candidate.score, bestScore, criterion, currentModel, candidate.model, alpha)
            bestScore = candidate.score;
            bestAction = candidate.action;
            bestPredictor = candidate.predictor;
            improved = true;
        end
    end
    
    % Apply best action
    if improved
        if strcmp(bestAction, 'add')
            currentPredictors{end+1} = bestPredictor;
            remainingPredictors = setdiff(remainingPredictors, {bestPredictor});
            if verbose
                fprintf('ADDED: %s (%s = %.4f)\n', bestPredictor, criterion, bestScore);
            end
        elseif strcmp(bestAction, 'remove')
            currentPredictors = setdiff(currentPredictors, {bestPredictor});
            remainingPredictors{end+1} = bestPredictor;
            if verbose
                fprintf('REMOVED: %s (%s = %.4f)\n', bestPredictor, criterion, bestScore);
            end
        end
        
        % Store history
        selectionHistory(iteration).iteration = iteration;
        selectionHistory(iteration).action = bestAction;
        selectionHistory(iteration).predictor = bestPredictor;
        selectionHistory(iteration).score = bestScore;
        selectionHistory(iteration).predictors = currentPredictors;
    else
        if verbose
            fprintf('No improvement found. Stopping.\n');
        end
    end
end

% Fit final model
bestModel = fitCurrentModel(data, responseVar, currentPredictors, randomEffects);

if verbose
    fprintf('\n=== FINAL MODEL ===\n');
    fprintf('Predictors: %s\n', strjoin(currentPredictors, ', '));
    fprintf('Formula: %s\n', formatModelString(currentPredictors, randomEffects));
    fprintf('%s = %.4f\n', criterion, getModelScore(bestModel, criterion));
    fprintf('Iterations: %d\n', iteration);
end

end

% Helper functions
function model = fitCurrentModel(data, responseVar, predictors, randomEffects)
    if isempty(predictors)
        formula = sprintf('%s ~ 1 + %s', responseVar, randomEffects);
    else
        formula = sprintf('%s ~ %s + %s', responseVar, strjoin(predictors, ' + '), randomEffects);
    end
    model = fitlme(data, formula);
end

function score = getModelScore(model, criterion)
    switch criterion
        case 'AIC'
            score = model.ModelCriterion.AIC;
        case 'BIC'
            score = model.ModelCriterion.BIC;
        case 'LRT'
            score = -2 * model.LogLikelihood; % Deviance
    end
end

function improved = isImprovement(newScore, currentScore, criterion, currentModel, newModel, alpha)
    switch criterion
        case {'AIC', 'BIC'}
            % Lower is better for AIC/BIC
            improved = newScore < currentScore;
        case 'LRT'
            % Use likelihood ratio test
            try
                % Compare nested models
                if currentModel.NumCoefficients > newModel.NumCoefficients
                    % Removing predictor: test if removal is justified
                    [~, pValue] = lratiotest(currentModel.LogLikelihood, newModel.LogLikelihood, ...
                        currentModel.NumCoefficients - newModel.NumCoefficients);
                    improved = pValue > alpha; % Non-significant means we can remove
                else
                    % Adding predictor: test if addition is justified
                    [~, pValue] = lratiotest(newModel.LogLikelihood, currentModel.LogLikelihood, ...
                        newModel.NumCoefficients - currentModel.NumCoefficients);
                    improved = pValue < alpha; % Significant means we should add
                end
            catch
                % If LRT fails, fall back to AIC
                improved = newScore < currentScore;
            end
    end
end

function str = formatModelString(predictors, randomEffects)
    if isempty(predictors)
        str = sprintf('1 + %s', randomEffects);
    else
        str = sprintf('%s + %s', strjoin(predictors, ' + '), randomEffects);
    end
end