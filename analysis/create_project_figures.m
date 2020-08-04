%% create_project_figures.m
% In this script we create main project figures
% John Bernabei
% With assistance from Ian Ong
% Litt Laboratory
% Summer 2020
% https://github.com/jbernabei/iEEG_atlas

%% set up workspace

clear all

band_names = {'broadband','alpha-theta','beta','low-gamma','high-gamma'};

% suppress warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field, id_field,...
    implant_field, outcome_field, resect_field, roi_field, target_field,...
    therapy_field, region_list, region_names] = set_up_workspace(iEEG_atlas_path);

%% generate csv file with patient demographics

create_data_table(all_patients, hasData_field);

%% Figure 1A: construct adjacency matrix of all good outcome patients

set(0,'units','inches')
screen_dims = get(0,'ScreenSize');
figure_width = 10;

% initializing cell arrays
good_mean_conn = cell(1,5);
good_std_conn = cell(1,5);

% minimum sample size for an edge weight to be included in the atlas
test_threshold = 3;

for f = 1:5
    test_band = f;

    % run all good outcome patients in atlas
    cond = [hasData_field{:}] & strcmp(outcome_field,'good');
    [good_mean_conn{f}, good_std_conn{f}, num_conn] = create_atlas(conn_field(cond), roi_field(cond), ...
    resect_field(cond), region_list, test_band, test_threshold);

    % visualize adjacency matrices with labels added
    figure(f+1);clf;
    fig = gcf;
    set(fig,'defaultAxesTickLabelInterpreter','none'); 
    set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_width)/2, figure_width, figure_width-1])
    imagesc(good_mean_conn{f},'AlphaData',~isnan(good_mean_conn{f}))
    axis(gca,'equal');
    set(gca,'color',0*[1 1 1]);
    set(gca,'xtick',(1:90),'xticklabel',all_locs)
    xtickangle(45)
    set(gca,'ytick',(1:90),'yticklabel',all_locs)
    set(gca,'fontsize', 4)
    cb = colorbar;
    cb.FontSize = 8;
    title(sprintf('Connectivity atlas of non-resected regions in good outcome patients (band %d)',test_band),'fontsize',12)
    save_name = sprintf('output/supplemental_figures/non_resected_good_outcome_atlas_band_%d.png',test_band);
    fig.InvertHardcopy = 'off';
    saveas(fig,save_name) % save plot to output folder
end

save('output/good_outcome_pt_atlas.mat','good_mean_conn','good_std_conn')

%% Figure 1B: visualize number of patients with each connection
figure(1);clf
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_width)/2, figure_width, figure_width-1])
imagesc(num_conn,'AlphaData',~(num_conn==0))
axis(gca,'equal');
set(gca,'color',0*[1 1 1]);
set(gca,'xtick',(1:90),'xticklabel',all_locs)
xtickangle(45)
set(gca,'ytick',(1:90),'yticklabel',all_locs)
set(gca,'fontsize', 4)
cb = colorbar;
cb.FontSize = 8;
title(sprintf('Sample sizes for each edge in non-resected regions of good outcome patients'),'fontsize',12)
save_name = sprintf('output/supplemental_figures/non_resected_good_outcome_sample_sizes.png');
fig.InvertHardcopy = 'off';
saveas(fig,save_name) % save plot to output folder

%% Figure 2B: cross - validate out-of-bag predictions

test_threshold = 3;

num_good_patients = sum([hasData_field{:}] & strcmp(outcome_field,'good'));
num_poor_patients = sum([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good'));
poor_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_z_score_results = cell(num_good_patients,1);
good_resected_z_score_results = cell(num_good_patients,1);

threshold_results = NaN(num_good_patients,5);
threshold_accuracies = NaN(num_good_patients,5);

% to hold logistic regression results
mnr_results = cell(5,3);
mdl_results = cell(5,1);

set(0,'units','inches')
screen_dims = get(0,'ScreenSize');
figure_width = 12;

for test_band = 1:5
    fprintf('\n')
    % cross-validation of good-outcome patients
    for s = 1:length(good_patient_indices)
        test_patient = all_patients(good_patient_indices(s));
        cv_patients = all_patients(good_patient_indices);
        cv_patients(s) = [];

        line_length = fprintf('Testing %s...', test_patient.patientID);

        % get connectivity atlas of excluded patients
        [mean_conn, std_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

        % get connectivity atlas of test patient
        [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list, test_band);
        
        % get non-resected region labels of test patient
        [~, patient_roi, ~] = nifti_values(test_patient.coords(setdiff(1:length(test_patient.roi),test_patient.resect),:),'localization/AAL116_WM.nii');

        % test atlas
        good_z_score_results{s} = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi);

        % get resected region labels of test patient
        [~, patient_roi, ~] = nifti_values(test_patient.coords(test_patient.resect,:),'localization/AAL116_WM.nii');

        % test atlas
        good_resected_z_score_results{s} = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi);
        
        % place results into all_patients
        all_patients(good_patient_indices(s)).z_scores(test_band).data.non_resected = good_z_score_results{s};
        all_patients(good_patient_indices(s)).z_scores(test_band).data.resected = good_resected_z_score_results{s};
        
        if strcmp(test_patient.therapy,'Ablation')
            [threshold_results(s,test_band),threshold_accuracies(s,test_band)] = get_optimal_threshold(good_z_score_results{s},good_resected_z_score_results{s});
        end
        
        fprintf(repmat('\b',1,line_length))
        
    end

    good_z_score_mean = nanmean(cat(3,good_z_score_results{:}),3);
    good_resected_z_score_mean = nanmean(cat(3,good_resected_z_score_results{:}),3);

    % save results to output folder
    %save('output/figure_2B_good_data.mat','good_z_score_mean','good_resected_z_score_mean')

    % repeat cross-validation for poor outcome patients
    poor_z_score_results = cell(num_poor_patients,1);
    poor_resected_z_score_results = cell(num_poor_patients,1);

    % calculate atlas for good outcome patients only
    % this serves as the "model" for the cross-validation
    cv_patients = all_patients(good_patient_indices);
    [mean_conn, std_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

    % cross-validation of poor-outcome patients
    for s = 1:length(poor_patient_indices)
        test_patient = all_patients(poor_patient_indices(s));

        line_length = fprintf('Testing %s...', test_patient.patientID);

        % get connectivity atlas of test patient
        [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list, test_band);

        % get non-resected region labels of test patient
        [~, patient_roi, ~] = nifti_values(test_patient.coords(setdiff(1:length(test_patient.roi),test_patient.resect),:),'localization/AAL116_WM.nii');

        % test atlas
        poor_z_score_results{s} = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi);

        % get resected region labels of test patient
        [~, patient_roi, ~] = nifti_values(test_patient.coords(test_patient.resect,:),'localization/AAL116_WM.nii');

        % test atlas
        poor_resected_z_score_results{s} = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi);
        
        % place results into all_patients
        all_patients(poor_patient_indices(s)).z_scores(test_band).data.non_resected = poor_z_score_results{s};
        all_patients(poor_patient_indices(s)).z_scores(test_band).data.resected = poor_resected_z_score_results{s};
        
        fprintf(repmat('\b',1,line_length))
        
    end

    poor_z_score_mean = mean(cat(3,poor_z_score_results{:}),3,'omitnan');
    poor_resected_z_score_mean = mean(cat(3,poor_resected_z_score_results{:}),3,'omitnan');

    % save results to output folder
    %save('output/figure_2B_poor_data.mat','poor_z_score_mean','poor_resected_z_score_mean')

    figure
    set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, 2, figure_width, 6])
    
    % plot ALL z-score results for GOOD outcome patients
    subplot(1,2,1)
    get_data = @(x) x(triu(true(size(x))));
    bin_width = 0.2;
    % remove outliers and bottom triangle of data
    good_plot_data = cellfun(get_data,good_z_score_results,'UniformOutput',false);
    good_plot_data = cell2mat(good_plot_data);
    good_plot_data = good_plot_data(:);
    good_plot_data = rmoutliers(good_plot_data(~isnan(good_plot_data) & ~isinf(good_plot_data)));
    %good_plot_data = good_plot_data(~isnan(good_plot_data) & ~isinf(good_plot_data));

    good_resected_plot_data = cellfun(get_data,good_resected_z_score_results,'UniformOutput',false);
    good_resected_plot_data = cell2mat(good_resected_plot_data);
    good_resected_plot_data = good_resected_plot_data(:);
    good_resected_plot_data = rmoutliers(good_resected_plot_data(~isnan(good_resected_plot_data) & ~isinf(good_resected_plot_data)));
    %good_resected_plot_data = good_resected_plot_data(~isnan(good_resected_plot_data) & ~isinf(good_resected_plot_data));

    histogram(good_plot_data,'Normalization','probability','BinWidth',bin_width);
    hold on
    histogram(good_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
    title({sprintf('Z-scores of connectivity strengths in good outcome patients (band %d)',test_band),''})
    ylabel('Density')
    xlabel('Z-score')
    % draw lines representing the medians of both groups
    xline(mean(good_plot_data),'b','LineWidth',2);
    xline(mean(good_resected_plot_data),'r','LineWidth',2);
    legend('Non-resected regions','Resected regions','Non-resected mean','Resected mean')
    legend('Location','northeast','Box','off')
    % set(gca,'YScale','log')
    hold off

    % plot ALL z-score results for POOR outcome patients
    % remove outliers and bottom triangle of data
    subplot(1,2,2)
    poor_plot_data = cellfun(get_data,poor_z_score_results,'UniformOutput',false);
    poor_plot_data = cell2mat(poor_plot_data);
    poor_plot_data = poor_plot_data(:);
    poor_plot_data = rmoutliers(poor_plot_data(~isnan(poor_plot_data) & ~isinf(poor_plot_data)));
    %poor_plot_data = poor_plot_data(~isnan(poor_plot_data) & ~isinf(poor_plot_data));

    poor_resected_plot_data = cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false);
    poor_resected_plot_data = cell2mat(poor_resected_plot_data);
    poor_resected_plot_data = poor_resected_plot_data(:);
    poor_resected_plot_data = rmoutliers(poor_resected_plot_data(~isnan(poor_resected_plot_data) & ~isinf(poor_resected_plot_data)));
    %poor_resected_plot_data = poor_resected_plot_data(~isnan(poor_resected_plot_data) & ~isinf(poor_resected_plot_data));

    
    hold on
    histogram(poor_plot_data,'Normalization','probability','BinWidth',bin_width);
    histogram(poor_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
    title({sprintf('Z-scores of connectivity strengths in poor outcome patients (band %d)',test_band),''})
    ylabel('Density')
    xlabel('Z-score')
    % draw lines representing the medians of both groups
    xline(mean(poor_plot_data),'b','LineWidth',2);
    xline(mean(poor_resected_plot_data),'r','LineWidth',2);
    legend('Non-resected regions','Resected regions','Non-resected mean','Resected mean')
    legend('Location','northeast','Box','off')
    % set(gca,'YScale','log')
    save_name = sprintf('output/z_score_histogram_band_%d.png',test_band);
    saveas(gcf,save_name) % save plot to output folder
    hold off

    % Statistical testing on z-score distributions
    fprintf('\n\n===== TESTS ON BAND %d =====\n', test_band)
%     fprintf('\n=== Two-sample Kolmogorov-Smirnov test ===\n')
%     fprintf('H0: Z-scores in non-resected and resected regions come from the same distribution.\n')
%     [h,p] = kstest2([good_plot_data; poor_plot_data],[good_resected_plot_data; poor_resected_plot_data],'Alpha',0.05);
%     fprintf('p-value = %d',p)
%     if h, fprintf('*'); end
%     fprintf('\nH0: Z-scores in non-resected and resected regions of good outcome patients come from the same distribution.\n')
%     [h,p] = kstest2(good_plot_data,good_resected_plot_data,'Alpha',0.05);
%     fprintf('p-value = %d',p)
%     if h, fprintf('*'); end
%     fprintf('\nH0: Z-scores in non-resected and resected regions of poor outcome patients come from the same distribution.\n')
%     [h,p] = kstest2(poor_plot_data,poor_resected_plot_data,'Alpha',0.05);
%     fprintf('p-value = %d',p)
%     if h, fprintf('*'); end

    fprintf('\n=== Two-sample t-test ===\n')
    fprintf('H0: The mean z-scores of non-resected regions in good and poor outcome patients are equal.\n')
    [h,p] = ttest2(good_plot_data,poor_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-score of resected regions in good and poor outcome patients are equal.\n')
    [h,p] = ttest2(good_resected_plot_data,poor_resected_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-scores of non-resected regions and resected regions in good outcome patients are equal.\n')
    [h,p] = ttest2(good_resected_plot_data,good_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-scores of non-resected regions and resected regions in poor outcome patients are equal.\n')
    [h,p] = ttest2(poor_resected_plot_data,poor_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end

    fprintf('\n=== Welch''s t-test ===\n')
    fprintf('H0: The mean z-scores of non-resected regions in good and poor outcome patients are equal.\n')
    [h,p] = ttest2(good_plot_data,poor_plot_data,'Vartype','unequal','Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-score of resected regions in good and poor outcome patients are equal.\n')
    [h,p] = ttest2(good_resected_plot_data,poor_resected_plot_data,'Vartype','unequal','Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-scores of non-resected regions and resected regions in good outcome patients are equal.\n')
    [h,p] = ttest2(good_resected_plot_data,good_plot_data,'Vartype','unequal','Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The mean z-scores of non-resected regions and resected regions in poor outcome patients are equal.\n')
    [h,p] = ttest2(poor_resected_plot_data,poor_plot_data,'Vartype','unequal','Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    
    fprintf('\n=== Wilcoxon rank sum test ===\n')
    fprintf('H0: The z-scores of non-resected regions in good and poor outcome patients come from distributions with equal medians.\n')
    [p,h] = ranksum(good_plot_data,poor_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The z-scores of resected regions in good and poor outcome patients come from distributions with equal medians.\n')
    [p,h] = ranksum(good_resected_plot_data,poor_resected_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The median z-scores of non-resected regions and resected regions in good outcome patients are equal.\n')
    [p,h] = ranksum(good_resected_plot_data,good_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end
    fprintf('\nH0: The median z-scores of non-resected regions and resected regions in poor outcome patients are equal.\n')
    [p,h] = ranksum(poor_resected_plot_data,poor_plot_data,'Alpha',0.05);
    fprintf('p-value = %d',p)
    if h, fprintf('*'); end

    % logistic regression
    % predictors: distance between mean resected and mean non-resected 
    % z-scores, mean & variance of the distribution of all z scores together
    % response: good or poor outcome
    
    % helper function
    get_average_data = @(x) nanmean(x(triu(true(size(x)))));
    
    % get distance values for good outcome patients
    good_distances = cell2mat(cellfun(get_average_data,good_resected_z_score_results,'UniformOutput',false)) - cell2mat(cellfun(get_average_data,good_z_score_results,'UniformOutput',false));
    
    % get distance values for poor outcome patients and append them
    poor_distances = cell2mat(cellfun(get_average_data,poor_resected_z_score_results,'UniformOutput',false)) - cell2mat(cellfun(get_average_data,poor_z_score_results,'UniformOutput',false));
    
    % combine two arrays
    distances = [good_distances; poor_distances];
    
    % get mean z-scores
    good_means = nanmean([cell2mat(cellfun(get_data,good_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,good_z_score_results,'UniformOutput',false).')]);
    
    poor_means = nanmean([cell2mat(cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,poor_z_score_results,'UniformOutput',false).')]);
    
    z_score_means = [good_means, poor_means];
    
    % get variance of z-scores
    good_variances = nanvar([cell2mat(cellfun(get_data,good_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,good_z_score_results,'UniformOutput',false).')]);
    
    poor_variances = nanvar([cell2mat(cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,poor_z_score_results,'UniformOutput',false).')]);
    
    z_score_variances = [good_variances, poor_variances];
    
    % generate array of outcomes
    outcomes = categorical([repmat(["good"],1,length(good_patient_indices)), repmat(["poor"],1,length(poor_patient_indices))]).';

    predictors = [distances,z_score_means.',z_score_variances.'];
    
    % remove any patients with NaN predictor values
    good_rows = ~any(isnan(predictors),2);
    predictors = predictors(good_rows,:);
    outcomes = outcomes(good_rows,:);
    
    % perform logistic regression
    [B,dev,stats] = mnrfit(predictors,outcomes.','Model','hierarchical');
    mnr_results(test_band,:) = {B,dev,stats};
    
    mdl = fitglm(predictors,outcomes,'Distribution','binomial','Link','logit');
    mdl_results(test_band) = {mdl};
    
    close all
end
fprintf('\n')

% get average threshold values for each band
avg_thresholds = nanmean(threshold_results,1);

% remove some variables from memory
vars = {'B','get_average_data','good_distances','good_means', ...
    'good_plot_data','good_resected_plot_data','good_resected_z_score_mean', ...
    'good_rows','good_variances','good_z_score_mean','h','line_length' ...
    'mdl','outcomes','p','poor_distances','poor_means','poor_plot_data', ...
    'poor_resected_plot_data','poor_resected_z_score_mean','poor_variances' ...
    'poor_z_score_mean','predictors','s','save_name','screen_dims', ...
    'stats','test_band','test_threshold','z_score_means','z_score_variances'};
clear(vars{:})
clear('vars')

%% saving z-scores and z-score histograms on an individual basis
fprintf('\n')

% exclude RNS patients
cond = [hasData_field{:}] & (strcmp(outcome_field,'good') | strcmp(outcome_field,'poor'));

all_outcome_patients = all_patients(cond);

counter = 0;

for k = 1:length(all_outcome_patients)
    % save z-score data to a .mat file
    dat = all_outcome_patients(k).z_scores;
    patientID = all_outcome_patients(k).patientID;
    
    line_length = fprintf('Generating plots for %s...', patientID);
    
    % create folder for the patient if one does not exist
    if ~exist(sprintf('output/patient_specific/%s',patientID), 'dir')
       mkdir(sprintf('output/patient_specific/%s',patientID))
    end
    
    save(sprintf('output/patient_specific/%s/%s_z_scores.mat',patientID,patientID),'dat')
    
    % generate figure with a z-score histogram for each band
    fig = figure;
    fig.WindowState = 'maximized';
    get_triu_data = @(x) x(triu(true(size(x))));
    bin_width = 0.2;
    
    for test_band = 1:5
        
        subplot(2,3,test_band)
        
        non_resected_scores = get_triu_data(all_outcome_patients(k).z_scores(test_band).data.non_resected);
        resected_scores = get_triu_data(all_outcome_patients(k).z_scores(test_band).data.resected);

        % remove NaN values
        non_resected_scores = non_resected_scores(~isnan(non_resected_scores) & ~isinf(non_resected_scores));
        resected_scores = resected_scores(~isnan(resected_scores) & ~isinf(resected_scores));
        
        % default p-value to display
        p = NaN;
        % test fails by default
        h = 0;
        % get Mann-Whitney p-value
        if ~(isempty(non_resected_scores) || isempty(resected_scores)), [p,h] = ranksum(non_resected_scores,resected_scores,'Alpha',0.05); end
        
        legend_entries = {};
        if ~isempty(non_resected_scores)
            histogram(non_resected_scores,'Normalization','probability','BinWidth',bin_width);
            xline(mean(non_resected_scores),'b','LineWidth',2);
            legend_entries = {sprintf('Non-resected regions (n=%d)',length(non_resected_scores)),'Non-resected mean'};
        end
        if ~isempty(resected_scores)
            hold on
            histogram(resected_scores,'Normalization','probability','BinWidth',bin_width);
            xline(mean(resected_scores),'r','LineWidth',2);
            legend_entries = [legend_entries,{sprintf('Resected regions (n=%d)',length(resected_scores)),'Resected mean'}];
        end
        title({sprintf('%s, p = %.5f%s',band_names{test_band},p,char(42*h)),''})
        ylabel('Density')
        xlabel('Z-score')
        if ~isempty(legend_entries), legend(legend_entries); end
        legend('Location','northeast','Box','off')
        hold off
    end
    sgtitle({sprintf('%s: z-scores of connectivity strengths by band',patientID),sprintf('Surgical outcome: %s',all_outcome_patients(k).outcome),sprintf('Targeted region: %s',all_outcome_patients(k).target)},'FontSize',10)
    saveas(gcf,sprintf('output/patient_specific/%s/%s_z_score_histograms.png',patientID,patientID)) % save plot to output folder
    close all
    
    if isempty(non_resected_scores)
        fprintf(' (no non-resected z-scores)')
    end
    
    if isempty(resected_scores)
        fprintf(' (no resected z-scores)')
    end
    
    if ~(isempty(resected_scores) || isempty(non_resected_scores))
        fprintf(repmat('\b',1,line_length))
    else
        fprintf('\n')
        counter = counter + 1;
    end
end

fprintf('Total number of patients with missing data: %d\n',counter)

%% visualizing p-values by individual
fprintf('\n')

% exclude RNS patients
cond = [hasData_field{:}] & (strcmp(outcome_field,'good') | strcmp(outcome_field,'poor'));

all_outcome_patients = all_patients(cond);

% initialize matrices to hold all patient stats
patient_p_values = NaN(5,length(all_outcome_patients));
direction = false(5,length(all_outcome_patients));
distance = NaN(5,length(all_outcome_patients));

for k = 1:length(all_outcome_patients)
    patientID = all_outcome_patients(k).patientID;
    line_length = fprintf('Getting p-values for %s...', patientID);
    get_triu_data = @(x) x(triu(true(size(x))));
    
    for test_band = 1:5
        
        non_resected_scores = get_triu_data(all_outcome_patients(k).z_scores(test_band).data.non_resected);
        resected_scores = get_triu_data(all_outcome_patients(k).z_scores(test_band).data.resected);

        % remove NaN values
        non_resected_scores = non_resected_scores(~isnan(non_resected_scores) & ~isinf(non_resected_scores));
        resected_scores = resected_scores(~isnan(resected_scores) & ~isinf(resected_scores));
        
        % default p-value to display
        p = NaN;
        % test fails by default
        h = 0;
        % get Mann-Whitney p-value
        if ~(isempty(non_resected_scores) || isempty(resected_scores))
            [p,h] = ranksum(non_resected_scores,resected_scores,'Alpha',0.05);
            % also, calculate the difference between the means
            distance(test_band,k) = mean(resected_scores) - mean(non_resected_scores);
        end
        
        % add p-value to visualiztion matrix
        patient_p_values(test_band,k) = p;
        
        % add direction of distance between the resected and non-resected
        % means
        direction(test_band,k) = (median(resected_scores) - median(non_resected_scores) > 0);
        
    end
    fprintf(repmat('\b',1,line_length))
end

alpha = 0.05;

set(0,'units','inches')
screen_dims = get(0,'ScreenSize');
figure_width = 12;
figure_height = 4;

fig = gcf;
imagesc(patient_p_values,'AlphaData',~isnan(patient_p_values))
set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_height)/2, figure_width, figure_height])
set(gca,'color',0*[1 1 1]);
axis(gca,'equal');
set(gca,'xtick',1:length({all_outcome_patients(:).patientID}),'xticklabel',{all_outcome_patients(:).patientID})
set(gca,'ytick',(1:5),'yticklabel',band_names)
set(gca,'fontsize', 8)
xtickangle(90)
hcb = colorbar;
colormap(flip(parula(40)))
title(hcb,'Two-tailed p-value')
title(sprintf('Mann-Whitney U test p-values for each band of each patient'),'fontsize',12)

pos_dist_sig_indices = find((patient_p_values < alpha) & direction);
neg_dist_sig_indices = find((patient_p_values < alpha) & ~direction);
%pos_dist_indices = find((patient_p_values >= alpha) & direction);
%neg_dist_indices = find((patient_p_values >= alpha) & ~direction);
[X,Y] = ndgrid(1:1:length({all_outcome_patients(:).patientID}),1:1:5);
X = X.';
Y = Y.';

% add significance asterisks to plot
text(X(pos_dist_sig_indices),Y(pos_dist_sig_indices),char(42),'Color',[0 0 1],'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
text(X(neg_dist_sig_indices),Y(neg_dist_sig_indices),char(42),'Color',[1 0 0],'HorizontalAlignment', 'center','VerticalAlignment', 'middle')
%text(X(pos_dist_indices),Y(pos_dist_indices),'.','Color',[0 0 1],'HorizontalAlignment', 'center','VerticalAlignment', 'bottom')
%text(X(neg_dist_indices),Y(neg_dist_indices),'.','Color',[1 0 0],'HorizontalAlignment', 'center','VerticalAlignment', 'bottom')

% add characters signifying outcome to each column
text(find(strcmp({all_outcome_patients(:).outcome},'good')),6*ones(1,sum(strcmp({all_outcome_patients(:).outcome},'good'))),'G','Color','g','HorizontalAlignment','center','fontsize',10)
text(find(strcmp({all_outcome_patients(:).outcome},'poor')),6*ones(1,sum(strcmp({all_outcome_patients(:).outcome},'poor'))),'P','Color','r','HorizontalAlignment','center','fontsize',10)

hold on

h = zeros(2, 1);
h(1) = plot(NaN,NaN,'.b');
h(2) = plot(NaN,NaN,'.r');
lg = legend(h,'resected median > non-resected median','resected median < non-resected median');
set(lg,'color','w')

save_name = sprintf('output/supplemental_figures/patient_p_values.png');
fig.InvertHardcopy = 'off';
saveas(fig,save_name) % save plot to output folder

hold off

% plot the difference between resected mean and non-resected mean

figure;
fig = gcf;
imagesc(distance,'AlphaData',~isnan(patient_p_values))
set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_height)/2, figure_width, figure_height])
set(gca,'color',0*[1 1 1]);
axis(gca,'equal');
set(gca,'xtick',1:length({all_outcome_patients(:).patientID}),'xticklabel',{all_outcome_patients(:).patientID})
set(gca,'ytick',(1:5),'yticklabel',band_names)
set(gca,'fontsize', 8)
xtickangle(90)
hcb = colorbar;
colormap(flip(cool(40)))

bound = max(abs(distance),[],'all');
caxis([-3, 3]);

title(hcb,'Resected mean - non-resected mean')
title(sprintf('Difference between mean resected/non-resected z-score by patient'),'fontsize',12)

% add characters signifying outcome to each column
text(find(strcmp({all_outcome_patients(:).outcome},'good')),6*ones(1,sum(strcmp({all_outcome_patients(:).outcome},'good'))),'G','Color','g','HorizontalAlignment','center','fontsize',10)
text(find(strcmp({all_outcome_patients(:).outcome},'poor')),6*ones(1,sum(strcmp({all_outcome_patients(:).outcome},'poor'))),'P','Color','r','HorizontalAlignment','center','fontsize',10)

frmt = @(x) convertCharsToStrings(sprintf('%.1f',x));
formatted_distance = arrayfun(frmt,distance);
formatted_distance(strcmp(formatted_distance,'NaN')) = "";

% add distance values to plot
text(X(:),Y(:),formatted_distance(:),'HorizontalAlignment', 'center','VerticalAlignment', 'middle','fontsize',5,'Color',[0,0,0])

h=gca;
h.YAxis.TickLength = [0 0];

save_name = sprintf('output/supplemental_figures/patient_differences.png');
fig.InvertHardcopy = 'off';
saveas(fig,save_name) % save plot to output folder

%% histograms of individual z-score distances by band

% generate figure with a z-score histogram for each band
fig = figure;
fig.WindowState = 'maximized';
get_triu_data = @(x) x(triu(true(size(x))));
bin_width = 0.2;

for test_band = 1:5

    subplot(2,3,test_band)
    
    good_distances = distance(test_band,find(strcmp({all_outcome_patients.outcome},'good')));
    poor_distances = distance(test_band,find(strcmp({all_outcome_patients.outcome},'poor')));

    % remove NaN values
    good_distances = good_distances(~isnan(good_distances) & ~isinf(good_distances));
    poor_distances = poor_distances(~isnan(poor_distances) & ~isinf(poor_distances));

    [p,h] = ranksum(good_distances,poor_distances,'Alpha',0.05);

    histogram(good_distances,'Normalization','probability','BinWidth',bin_width);
    xline(mean(good_distances),'b','LineWidth',2);
    hold on
    histogram(poor_distances,'Normalization','probability','BinWidth',bin_width);
    xline(mean(poor_distances),'r','LineWidth',2);
    title({sprintf('%s, p = %.5f%s',band_names{test_band},p,char(42*h)),''})
    ylabel('Density')
    xlabel('Difference')
    legend({'Good outcome patients','Mean distance of good outcome patients','Poor outcome patients','Mean distance of poor outcome patients'});
    legend('Location','southoutside','Box','off')
    hold off
end
sgtitle({sprintf('Difference in mean z-scores (resected - non-resected) by band and outcome'),''},'FontSize',14)
saveas(gcf,sprintf('output/supplemental_figures/distance_and_outcome_comparison.png')) % save plot to output folder

%% plot logistic regression results

bound = 30;
interval = 5;
jitter = 1;
marker_size = 80;

fig = figure;
fig.WindowState = 'maximized';
hold on

% plot 3d volume
for test_band = 1:5
    coeffs = mnr_results{test_band,1};
    f = @(x,y,z) 1./(1+exp(-(coeffs(1)+(coeffs(2)*x)+(coeffs(3)*y)+(coeffs(4)*z))));
    [x,y,z] = ndgrid(-bound/2:interval:bound/2,-bound/2:interval:bound/2,interval:bound);
    x = x + jitter*(rand(size(x))-0.5);
    y = y + jitter*(rand(size(y))-0.5);
    z = z + jitter*(rand(size(z))-0.5);
    v = f(x,y,z);

%     p1 = patch(isosurface(x,y,z,v,0.25));
%     hold on
%     p2 = patch(isosurface(x,y,z,v,0.5));
%     p3 = patch(isosurface(x,y,z,v,0.75));
%     isonormals(x,y,z,v,p1);
%     isonormals(x,y,z,v,p2);
%     isonormals(x,y,z,v,p3);
%     p1.FaceColor = 'red';
%     p2.FaceColor = 'yellow';
%     p3.FaceColor = 'green';
%     p1.EdgeColor = 'none';
%     p2.EdgeColor = 'none';
%     p3.EdgeColor = 'green';
%     daspect([1,1,1])
%     view(3)
    
    x = reshape(x,[],1); 
    y = reshape(y,[],1);
    z = reshape(z,[],1);
    v = reshape(v,[],1);
    
    subplot(2,3,test_band)
    set(gca, 'LooseInset', get(gca,'TightInset'))
    scatter3(x,y,z,marker_size,v,'filled','o');
    alpha(2*interval/bound);
    view(3)
    axis tight
    
    cb = colorbar;
    cb.Label.String = 'Probability of good surgical outcome';
    
    colormap('parula')
    
    xlabel('Distance')
    ylabel('Mean z-score')
    zlabel('Variance of z-scores')
    
    title({sprintf('Logistic regression for patient outcome on band %d',test_band),''});
end

saveas(gcf,'output/3d_logistic_regression_results.png') % save plot to output folder

hold off

fig = figure;
fig.WindowState = 'maximized';
hold on

% plot 3d surface
for test_band = 1:5
    coeffs = mnr_results{test_band,1};
    f = @(x,y) 1./(1+exp(-(coeffs(1)+(coeffs(2)*x)+(coeffs(3)*y))));
    [x,y] = ndgrid(-bound/2:interval/2:bound/2);
    v = f(x,y);
    
    subplot(2,3,test_band)
    set(gca, 'LooseInset', get(gca,'TightInset'))
    surf(x,y,v);
    view(3)
    axis tight
    
    colormap('parula')
    
    xlabel('Distance')
    ylabel('Mean z-score')
    zlabel('Probability of good outcome')
    
    title({sprintf('Logistic regression for patient outcome on band %d',test_band),''});
end

saveas(gcf,'output/2d_logistic_regression_results.png') % save plot to output folder

hold off

fig = figure;
fig.WindowState = 'maximized';

format long

% plot ROC curves
for test_band = 1:5
    subplot(2,3,test_band)
    results = mdl_results{test_band};
    scores = results.Fitted.Probability;
    [X,Y,T,AUC] = perfcurve(outcomes,scores,'good');
    %fprintf('Band %d:\nArea under curve = %16.f\n',test_band,AUC);
    
    subplot(2,3,test_band)
    plot(X,Y)
    x = 0:0.1:1;
    y = x;
    hold on
    plot(x,y)
    xlabel('False positive rate') 
    ylabel('True positive rate')
    text(0.4,0.1,sprintf('Area under curve = %.5f',AUC));
    title({sprintf('ROC for Classification by Logistic Regression (band %d)',test_band),''})
end

saveas(gcf,'output/roc_results.png') % save plot to output folder

hold off

%% Clinical hypothesis testing

% part 1
% assemble set of bilateral (RNS) patients, and test whether they have
% higher inter-hemispheric z scores compared to unilateral good outcome
% patients
bilateral_pt_indices = find([hasData_field{:}] & strcmp(therapy_field,'RNS'));
num_bilateral_pts = length(bilateral_pt_indices);
bilateral_z_score_results = cell(num_bilateral_pts,1);

% define checkerboard matrices for extracting inter/intra hemispheric
% connections (odd-odd and even-even conns are intra, odd-even and even-odd conns are inter)
inter_hemisphere_conns = zeros(91*91,1); inter_hemisphere_conns(1:2:end) = 1;
inter_hemisphere_conns = reshape(inter_hemisphere_conns,[91,91]);
inter_hemisphere_conns(end,:) = []; inter_hemisphere_conns(:,end) = [];

intra_hemisphere_conns = zeros(91*91,1); intra_hemisphere_conns(2:2:end) = 1;
intra_hemisphere_conns = reshape(intra_hemisphere_conns,[91,91]);
intra_hemisphere_conns(end,:) = []; intra_hemisphere_conns(:,end) = [];

for test_band = 1:5
    for s = 1:length(bilateral_pt_indices)
        test_patient = all_patients(bilateral_pt_indices(s));

         % get connectivity atlas of test patient
         [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list, test_band);

         % get non-resected region labels of test patient
         [~, patient_roi, ~] = nifti_values(test_patient.coords,'localization/AAL116_WM.nii');

         % test atlas
         bilateral_z_score_results{s} = test_patient_conn(good_mean_conn{test_band}, good_std_conn{test_band}, region_list, patient_conn, patient_roi);
         
         interhemispheric_results(test_band,s) = nanmean(bilateral_z_score_results{s}(find(inter_hemisphere_conns)));
         intrahemispheric_results(test_band,s) = nanmean(bilateral_z_score_results{s}(find(intra_hemisphere_conns)));
    end
end

%% Clinical hypothesis testing
% here we will simulate how the atlas could be used in a prospective manner

patient_hypothesis_list = {'HUP116','HUP117','HUP130','HUP138','HUP139','HUP140',...
                           'HUP141','HUP157','HUP164','HUP165',...
                           'HUP171','HUP185','HUP188'};
                       
lobar_data = readtable('localization/lobes_aal.xlsx');

for s = 1:length(patient_hypothesis_list)
    % get test patient
    test_pt_ind = find(strcmp(id_field,patient_hypothesis_list{s}));
    test_patient = all_patients(test_pt_ind);

     % get connectivity atlas of test patient
     [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list, test_band);

     % get non-resected region labels of test patient
     [~, patient_roi, ~] = nifti_values(test_patient.coords,'localization/AAL116_WM.nii');

     % test atlas
     test_z_score_results{s} = test_patient_conn(good_mean_conn{test_band}, good_std_conn{test_band}, region_list, patient_conn, patient_roi);

    % step 2 extract z scores of ROIs within primary/secondary hypothesis
    % regions and scores from each region to all other regions
    primary_hypothesis = test_patient.hypothesis_1;
    secondary_hypothesis = test_patient.hypothesis_2;
    %true_resection = strcat(strcat(test_patient.laterality,'_'), test_patient.target);
    
    % find which ROI are contained in primary and secondary hypothesis
    ROI_primary = contains(lobar_data{:,3},primary_hypothesis);
    ROI_secondary = contains(lobar_data{:,3},secondary_hypothesis);
   
    % in-in primary
    in_in_1 = test_z_score_results{s}(ROI_primary,ROI_primary);
    
    % in-out primary
    in_out_1 = test_z_score_results{s}(ROI_primary,setdiff([1:90],ROI_primary));
    
    % in-in secondary
    in_in_2 = test_z_score_results{s}(ROI_secondary,ROI_secondary);
    
    % in-out secondary
    in_out_2 = test_z_score_results{s}(ROI_primary,setdiff([1:90],ROI_secondary));
    
    % out-out
    out_out = test_z_score_results{s}(setdiff([1:90],[ROI_primary,ROI_secondary]),setdiff([1:90],[ROI_primary,ROI_secondary]));
    
    % step 3 generate rendering of brain with primary/secondary regions labeled 
    % by colored nodes and all edge weights rendered
    
    [atlas_mni] = create_distance_matrix({test_patient.coords}, region_list);
    
    node_color = zeros(90,1);
    node_color(ROI_primary,1) = -1; % sets to blue
    node_color(ROI_secondary,1) = 1; % sets to red
    
    % create final electrode matrix
    final_elec_matrix = [atlas_mni, node_color, ones(90,1)];
    
    unsampled_roi = find(isnan(final_elec_matrix(:,1)));
    
    final_elec_matrix(unsampled_roi,:) = [];
    
    % get rid of unsampled regions
    test_z_score_results{s}(unsampled_roi,:) = [];
    test_z_score_results{s}(:,unsampled_roi) = [];
    
    % zero out any NaNs
    test_z_score_results{s}(find(isnan(test_z_score_results{s}))) = 0;
    
    z_score_results = test_z_score_results{s};
    
    dlmwrite('output/hypothesis_test.node',final_elec_matrix,'delimiter',' ','precision',5)
    save('output/hypothesis_test.edge','z_score_results','-ascii');
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/hypothesis_test.node','output/hypothesis_test.edge','output/hypothesis_test_options.mat',sprintf('output/hypothesis_test_%s.jpg',patient_hypothesis_list{s}))

    % step 4 plot boxplot of z scores for each of the four classes: in-in
    % primary, in-out primary, in-in secondary, in-out secondary
    
    boxplot_data = NaN.*zeros(90*90,5);
    plot_labels = {'IN-IN Region 1','IN-OUT Region 1','IN-IN Region 2','IN-OUT Region 2','OUT-OUT'};
    
    boxplot_data(1:length(in_in_1(:)),1) = in_in_1(:);
    boxplot_data(1:length(in_out_1(:)),2) = in_out_1(:);
    boxplot_data(1:length(in_in_2(:)),3) = in_in_2(:);
    boxplot_data(1:length(in_out_2(:)),4) = in_out_2(:); 
    boxplot_data(1:length(out_out(:)),5) = out_out(:);
    
    % create boxplot
    figure(s);clf;
    boxplot(boxplot_data,'labels',plot_labels)
    
    % calculate p values
    for h = 1:5
        for y = 1:5
            hypothesis_p_val(h,y) = ranksum(boxplot_data(:,h),boxplot_data(:,y));
        end
    end

end
%% algorithm for identifying targets from full patient connectivity
% across frequency bands. Use dice score for localization quality

% apply in good outcome patients

% apply in poor outcome patients

% quantify and compare results

%% testing create_atlas_by_edge
test_band = 1;
test_threshold = 3;
[mean_conn, std_conn, num_conn] = create_atlas_by_edge({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

figure;
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
%set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_width)/2, figure_width, figure_width-1])
imagesc(num_conn,'AlphaData',~(num_conn==0))
axis(gca,'equal');
set(gca,'color',0*[1 1 1]);
set(gca,'xtick',(1:90),'xticklabel',all_locs)
xtickangle(45)
set(gca,'ytick',(1:90),'yticklabel',all_locs)
set(gca,'fontsize', 4)
cb = colorbar;
cb.FontSize = 8;
title(sprintf('Electrode-level sample sizes for each edge in non-resected regions of good outcome patients'),'fontsize',12)
save_name = sprintf('output/supplemental_figures/non_resected_good_outcome_electrode_level_sample_sizes.png');
fig.InvertHardcopy = 'off';
saveas(fig,save_name) % save plot to output folder

for f = 1:5
    test_band = f;

    % run all good outcome patients in atlas
    [mean_conn, std_conn, num_conn] = create_atlas_by_edge({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

    % visualize adjacency matrices with labels added
    figure;
    fig = gcf;
    set(fig,'defaultAxesTickLabelInterpreter','none'); 
    %set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, (screen_dims(4)-figure_width)/2, figure_width, figure_width-1])
    imagesc(mean_conn,'AlphaData',~isnan(mean_conn))
    axis(gca,'equal');
    set(gca,'color',0*[1 1 1]);
    set(gca,'xtick',(1:90),'xticklabel',all_locs)
    xtickangle(45)
    set(gca,'ytick',(1:90),'yticklabel',all_locs)
    set(gca,'fontsize', 4)
    cb = colorbar;
    cb.FontSize = 8;
    title(sprintf('Electrode-level connectivity atlas of non-resected regions in good outcome patients (band %d)',test_band),'fontsize',12)
    save_name = sprintf('output/supplemental_figures/electrode_level_non_resected_good_outcome_atlas_band_%d.png',test_band);
    fig.InvertHardcopy = 'off';
    saveas(fig,save_name) % save plot to output folder
end