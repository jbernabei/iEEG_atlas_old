function [optimal_threshold, best_accuracy] = get_optimal_threshold(non_resected_z_scores, resected_z_scores)
%GET_OPTIMAL_THRESHOLD returns the z-value that is most predictive of
%   resected vs. non-resected regions in a patient. Also returns the
%   corresponding accuracy for the z-value

% helper function
get_data = @(x) x(triu(true(size(x))));

min_z = min(min(non_resected_z_scores(:)),min(resected_z_scores(:)));
max_z = max(max(non_resected_z_scores(:)),max(resected_z_scores(:)));

% placeholders
optimal_threshold = min_z;
best_accuracy = 0;

% number of steps to take between the min and max z-score
steps = 100;

% get z-scores, discarding NaNs
non_res = get_data(non_resected_z_scores);
res = get_data(resected_z_scores);
    
non_res = non_res(~isnan(non_res));
res = res(~isnan(res));

% determine threshold values to test
for test_threshold = linspace(min_z,max_z,steps)
    
    correct_non_resected = non_res < test_threshold;
    correct_resected = res > test_threshold;
    
    % calculate proportion of correct guesses
    accuracy = mean([correct_non_resected; correct_resected]);
    
    % update optimal_threshold if the test threshold produces a new highest
    % accuracy
    if accuracy > best_accuracy
        optimal_threshold = test_threshold;
        best_accuracy = accuracy;
    end
    
end

