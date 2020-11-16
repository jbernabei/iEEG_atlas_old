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

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field, id_field,...
    implant_field, outcome_field, resect_field, roi_field, target_field,...
    therapy_field, region_list, region_names] = set_up_workspace(iEEG_atlas_path);

%% generate csv file with patient demographics

create_data_table(all_patients, hasData_field);

%% Figure 2B: cross - validate out-of-bag predictions

test_threshold = 3;

num_good_patients = sum([all_patients.hasData] & strcmp({all_patients.outcome},'good'));
num_poor_patients = sum([all_patients.hasData] & strcmp({all_patients.outcome},'poor'));

good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good'));
poor_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor'));
ablation_indices = find([all_patients.hasData] & strcmp({all_patients.therapy},'Ablation'));
resection_indices = find([all_patients.hasData] & strcmp({all_patients.therapy},'Resection'));
ecog_indices = find([all_patients.hasData] & strcmp({all_patients.implant},'ECoG'));
seeg_indices = find([all_patients.hasData] & strcmp({all_patients.implant},'SEEG'));

threshold_results = NaN(num_good_patients,5);
threshold_accuracies = NaN(num_good_patients,5);

% to hold rank-sum test results
validation_test_results = cell(5,9);
abl_res_test_results = cell(12,9);
implant_test_results = cell(12,9);

set(0,'units','inches')
screen_dims = get(0,'ScreenSize');
figure_width = 12;

for atlas_method = ["patient","edge"]
    for test_band = 5
        % cross-validation of good-outcome patients
        for s = 1:length(good_patient_indices)
            test_patient = all_patients(good_patient_indices(s));
            cv_patients = all_patients(good_patient_indices);
            cv_patients(s) = [];
            
            patient_id = test_patient.patientID;

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [out_out_scores ,in_in_scores, in_out_scores, all_scores] =...
                validate_patient(test_patient,cv_patients,region_list,test_band,test_threshold,"sem",atlas_method);

            % get all good outcome patient zscores
            good_patient_zscores(test_band).freq{s} = all_scores;
            
            % place results into all_patients
            all_patients(good_patient_indices(s)).z_scores(test_band).data.out_out = out_out_scores;
            all_patients(good_patient_indices(s)).z_scores(test_band).data.in_in = in_in_scores;
            all_patients(good_patient_indices(s)).z_scores(test_band).data.in_out = in_out_scores;

            if strcmp(test_patient.therapy,'Ablation')
                [threshold_results(s,test_band),threshold_accuracies(s,test_band)] =...
                    get_optimal_threshold(out_out_scores,in_in_scores);
            end

            fprintf(repmat('\b',1,line_length))

        end
        
        save(sprintf('data/%s/%s_all_scores_freq_%d.mat',patient_id,patient_id,test_band),'all_scores');

        % repeat cross-validation for poor outcome patients

        % calculate atlas for good outcome patients only
        % this serves as the "model" for the cross-validation
        cv_patients = all_patients(good_patient_indices);

        % cross-validation of poor-outcome patients
        for s = 1:length(poor_patient_indices)
            test_patient = all_patients(poor_patient_indices(s));

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [out_out_scores, in_in_scores, in_out_scores ,all_scores] = ...
                validate_patient(test_patient,cv_patients,region_list,test_band,test_threshold,"sem",atlas_method);

            % get all good outcome patient zscores
            poor_patient_zscores(test_band).freq{s} = all_scores;
            
            % place results into all_patients
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.out_out = out_out_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_in = in_in_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_out = in_out_scores;

            fprintf(repmat('\b',1,line_length))

        end

        % plot ALL z-score results for GOOD and POOR outcome patients

        figure
        set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, 2, figure_width, 6])

        bin_width = 1;

        good_scores = {all_patients(good_patient_indices).z_scores};
        poor_scores = {all_patients(poor_patient_indices).z_scores};
        sub_groups = {good_scores, poor_scores};

        edge_groups = {'out_out','in_out','in_in'};
        edge_pairs = {'out_out','in_in';'out_out','in_out';'in_out','in_in'};
        colors = {'#0072BD','#D95319','#EDB120'};
        titles = {{sprintf('Standardized scores of connectivity strengths\nin good outcome patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Standardized scores of connectivity strengths\nin poor outcome patients (%s)',band_names{test_band}),''}};

        figure(1);clf;
        for j = 1:length(sub_groups)

            subplot(1,2,j)
            for k = 1:length(edge_groups)

                get_score_data = @(x) x(test_band).data.(edge_groups{k})(triu(true(size(x(test_band).data.(edge_groups{k})))));
                % remove outliers and bottom triangle of data
                plot_data = cellfun(get_score_data,sub_groups{j},'UniformOutput',false);
                plot_data = reshape(cell2mat(plot_data),[],1);
                plot_data = plot_data(~isnan(plot_data) & ~isinf(plot_data));

                histogram(plot_data,'BinWidth',bin_width); % no rm_outliers
                hold on

                % plot median
                %xline(median(plot_data),'--','Color',colors{k},'LineWidth',1);

            end

            title(titles{j})
            ylabel('Count')
            xlabel('Score')

            legend('OUT-OUT','OUT-OUT median','IN-OUT','IN-OUT median','IN-IN','IN-IN median')
            legend('Location','northeast','Box','off')

            hold off
        end

        save_name = sprintf('output/supplemental_figures/%s_method/by_outcome/z_score_histogram_band_%d.png',atlas_method,test_band);
        saveas(gcf,save_name) % save plot to output folder

        % Statistical testing on z-score distributions
        % test edge type distributions within outcome group
        for j = 1:length(sub_groups)
            for k = 1:length(edge_pairs)
                edge_pair = edge_pairs(k,:);
                distribs = cell(1,2);
                for m = 1:2
                    get_score_data = @(x) x(test_band).data.(edge_pair{m})(triu(true(size(x(test_band).data.(edge_pair{m})))));
                    distrib = cellfun(get_score_data,sub_groups{j},'UniformOutput',false);
                    distrib = cell2mat(distrib);
                    distrib = distrib(:);
                    distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                    distribs{m} = distrib;
                end
                [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                validation_test_results{test_band,k+((j-1)*(length(edge_pairs)))} = ...
                    sprintf('%.5f%s',p,char(42*h));
            end
        end
        % test the same edge type distributions between outcome groups
        for j = 1:length(edge_groups)
            distribs = cell(1,2);
            for m = 1:length(sub_groups)
                get_score_data = @(x) x(test_band).data.(edge_groups{j})(triu(true(size(x(test_band).data.(edge_groups{j})))));
                distrib = cellfun(get_score_data,sub_groups{m},'UniformOutput',false);
                distrib = cell2mat(distrib);
                distrib = distrib(:);
                distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                distribs{m} = distrib;
            end
            [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                validation_test_results{test_band,j+(length(sub_groups)*(length(edge_pairs)))} = ...
                    sprintf('%.5f%s',p,char(42*h));
        end

        % plot results by therapy and outcome type

        figure
        set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, 1, figure_width, 8])

        bin_width = 1;

        good_scores_1 = {all_patients(intersect(good_patient_indices,ablation_indices)).z_scores};
        poor_scores_1 = {all_patients(intersect(poor_patient_indices,ablation_indices)).z_scores};
        good_scores_2 = {all_patients(intersect(good_patient_indices,resection_indices)).z_scores};
        poor_scores_2 = {all_patients(intersect(poor_patient_indices,resection_indices)).z_scores};
        sub_groups = {good_scores_1, poor_scores_1; good_scores_2, poor_scores_2};

        edge_groups = {'out_out','in_out','in_in'};
        edge_pairs = {'out_out','in_in';'out_out','in_out';'in_out','in_in'};
        colors = {'#0072BD','#D95319','#EDB120'};
        titles = {{sprintf('Good outcome ablation patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Poor outcome ablation patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Good outcome resection patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Poor outcome resection patients (%s)',band_names{test_band}),''}};

        figure(2);clf;
        for j = 1:length(sub_groups(:))

            subplot(2,2,j)
            for k = 1:length(edge_groups)

                get_score_data = @(x) x(test_band).data.(edge_groups{k})(triu(true(size(x(test_band).data.(edge_groups{k})))));
                % remove outliers and bottom triangle of data
                plot_data = cellfun(get_score_data,sub_groups{j},'UniformOutput',false);
                plot_data = reshape(cell2mat(plot_data),[],1);
                plot_data = plot_data(~isnan(plot_data) & ~isinf(plot_data));

                histogram(plot_data,'BinWidth',bin_width); % no rmoutliers
                hold on

                % plot mean
                %xline(median(plot_data),'--','Color',colors{k},'LineWidth',1);

            end

            title(titles{j})
            ylabel('Count')
            xlabel('Score')

            legend('OUT-OUT','OUT-OUT median','IN-OUT','IN-OUT median','IN-IN','IN-IN median')
            legend('Location','northeast','Box','off')

            hold off
        end

        save_name = sprintf('output/supplemental_figures/%s_method/by_therapy/abl_res_histogram_band_%d.png',atlas_method,test_band);
        saveas(gcf,save_name) % save plot to output folder

        % Statistical testing on z-score distributions
        % test edge type distributions within outcome group
        for row = 1:length(sub_groups)
            for j = 1:length(sub_groups)
                for k = 1:length(edge_pairs)
                    edge_pair = edge_pairs(k,:);
                    distribs = cell(1,2);
                    for m = 1:2
                        get_score_data = @(x) x(test_band).data.(edge_pair{m})(triu(true(size(x(test_band).data.(edge_pair{m})))));
                        distrib = cellfun(get_score_data,sub_groups{row,j},'UniformOutput',false);
                        distrib = cell2mat(distrib);
                        distrib = distrib(:);
                        distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                        distribs{m} = distrib;
                    end
                    [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                    abl_res_test_results{1+test_band+(row-1)*6,k+((j-1)*(length(edge_pairs)))} = ...
                        sprintf('%.5f%s',p,char(42*h));
                    %fprintf('Tested band %d, %s and %s. Row = %d, col = %d. Outputted to (%d,%d)\n',test_band,edge_pair{1},edge_pair{2},row,j,1+test_band+(row-1)*6,k+((j-1)*(length(edge_pairs))))
                end
            end
            % test the same edge type distributions between outcome groups
            for j = 1:length(edge_groups)
                distribs = cell(1,2);
                for m = 1:length(sub_groups)
                    get_score_data = @(x) x(test_band).data.(edge_groups{j})(triu(true(size(x(test_band).data.(edge_groups{j})))));
                    distrib = cellfun(get_score_data,sub_groups{row,m},'UniformOutput',false);
                    distrib = cell2mat(distrib);
                    distrib = distrib(:);
                    distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                    distribs{m} = distrib;
                end
                [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                    abl_res_test_results{1+test_band+(row-1)*6,j+(length(sub_groups)*(length(edge_pairs)))} = ...
                        sprintf('%.5f%s',p,char(42*h));
            end
        end
        
        % plot results by implant and outcome type

        figure
        set(gcf,'Units','inches','Position',[(screen_dims(3)-figure_width)/2, 1, figure_width, 8])

        bin_width = 1;

        good_scores_1 = {all_patients(intersect(good_patient_indices,ecog_indices)).z_scores};
        poor_scores_1 = {all_patients(intersect(poor_patient_indices,ecog_indices)).z_scores};
        good_scores_2 = {all_patients(intersect(good_patient_indices,seeg_indices)).z_scores};
        poor_scores_2 = {all_patients(intersect(poor_patient_indices,seeg_indices)).z_scores};
        sub_groups = {good_scores_1, poor_scores_1; good_scores_2, poor_scores_2};

        edge_groups = {'out_out','in_out','in_in'};
        edge_pairs = {'out_out','in_in';'out_out','in_out';'in_out','in_in'};
        colors = {'#0072BD','#D95319','#EDB120'};
        titles = {{sprintf('Good outcome ECoG patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Poor outcome ECoG patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Good outcome SEEG patients (%s)',band_names{test_band}),''}, ...
            {sprintf('Poor outcome SEEG patients (%s)',band_names{test_band}),''}};

        figure(3);clf;
        for j = 1:length(sub_groups(:))

            subplot(2,2,j)
            for k = 1:length(edge_groups)

                get_score_data = @(x) x(test_band).data.(edge_groups{k})(triu(true(size(x(test_band).data.(edge_groups{k})))));
                % remove outliers and bottom triangle of data
                plot_data = cellfun(get_score_data,sub_groups{j},'UniformOutput',false);
                plot_data = reshape(cell2mat(plot_data),[],1);
                plot_data = plot_data(~isnan(plot_data) & ~isinf(plot_data));

                histogram(plot_data,'BinWidth',bin_width); %no rmoutliers
                hold on

                % plot mean
                %xline(median(plot_data),'--','Color',colors{k},'LineWidth',1);

            end

            title(titles{j})
            ylabel('Count')
            xlabel('Score')

            legend('OUT-OUT','OUT-OUT median','IN-OUT','IN-OUT median','IN-IN','IN-IN median')
            legend('Location','northeast','Box','off')

            hold off
        end

        save_name = sprintf('output/supplemental_figures/%s_method/by_implant/implant_histogram_band_%d.png',atlas_method,test_band);
        saveas(gcf,save_name) % save plot to output folder

        % Statistical testing on z-score distributions
        % test edge type distributions within outcome group
        for row = 1:length(sub_groups)
            for j = 1:length(sub_groups)
                for k = 1:length(edge_pairs)
                    edge_pair = edge_pairs(k,:);
                    distribs = cell(1,2);
                    for m = 1:2
                        get_score_data = @(x) x(test_band).data.(edge_pair{m})(triu(true(size(x(test_band).data.(edge_pair{m})))));
                        distrib = cellfun(get_score_data,sub_groups{row,j},'UniformOutput',false);
                        distrib = cell2mat(distrib);
                        distrib = distrib(:);
                        distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                        distribs{m} = distrib;
                    end
                    [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                    implant_test_results{1+test_band+(row-1)*6,k+((j-1)*(length(edge_pairs)))} = ...
                        sprintf('%.5f%s',p,char(42*h));
                    %fprintf('Tested band %d, %s and %s. Row = %d, col = %d. Outputted to (%d,%d)\n',test_band,edge_pair{1},edge_pair{2},row,j,1+test_band+(row-1)*6,k+((j-1)*(length(edge_pairs))))
                end
            end
            % test the same edge type distributions between outcome groups
            for j = 1:length(edge_groups)
                distribs = cell(1,2);
                for m = 1:length(sub_groups)
                    get_score_data = @(x) x(test_band).data.(edge_groups{j})(triu(true(size(x(test_band).data.(edge_groups{j})))));
                    distrib = cellfun(get_score_data,sub_groups{row,m},'UniformOutput',false);
                    distrib = cell2mat(distrib);
                    distrib = distrib(:);
                    distrib = distrib(~isnan(distrib) & ~isinf(distrib));
                    distribs{m} = distrib;
                end
                [p,h] = ranksum(distribs{1},distribs{2},'Alpha',0.05);
                    implant_test_results{1+test_band+(row-1)*6,j+(length(sub_groups)*(length(edge_pairs)))} = ...
                        sprintf('%.5f%s',p,char(42*h));
            end
        end

        %close all
    end
% get average threshold values for each band
avg_thresholds = nanmean(threshold_results,1);

% validation_test_result_table = cell2table(validation_test_results, ...
%     'VariableNames',{'OUT-OUT to IN-IN (good)', 'OUT-OUT to IN-OUT (good)', ...
%     'IN-OUT to IN-IN (good)','OUT-OUT to IN-IN (poor)','OUT-OUT to IN-OUT (poor)', ...
%     'IN-OUT to IN-IN (poor)', 'OUT-OUT to OUT-OUT (between)', 'IN-OUT to IN-OUT (between)', ...
%     'IN-IN to IN-IN (between)'},'RowNames',band_names);
% 
% add_char = @(x) [x,':'];

% abl_res_test_result_table = cell2table(abl_res_test_results, ...
%     'VariableNames',{'OUT-OUT to IN-IN (good)', 'OUT-OUT to IN-OUT (good)', ...
%     'IN-OUT to IN-IN (good)','OUT-OUT to IN-IN (poor)','OUT-OUT to IN-OUT (poor)', ...
%     'IN-OUT to IN-IN (poor)', 'OUT-OUT to OUT-OUT (between)', 'IN-OUT to IN-OUT (between)', ...
%     'IN-IN to IN-IN (between)'},'RowNames',[{'Ablation'}, band_names(:)', {'Resection'}, cellfun(add_char,band_names(:)','UniformOutput',false)]);
% 
% implant_test_result_table = cell2table(implant_test_results, ...
%     'VariableNames',{'OUT-OUT to IN-IN (good)', 'OUT-OUT to IN-OUT (good)', ...
%     'IN-OUT to IN-IN (good)','OUT-OUT to IN-IN (poor)','OUT-OUT to IN-OUT (poor)', ...
%     'IN-OUT to IN-IN (poor)', 'OUT-OUT to OUT-OUT (between)', 'IN-OUT to IN-OUT (between)', ...
%     'IN-IN to IN-IN (between)'},'RowNames',[{'ECoG'}, band_names(:)', {'SEEG'}, cellfun(add_char,band_names(:)','UniformOutput',false)]);
% 
% writetable(validation_test_result_table,sprintf('output/supplemental_figures/%s_method/by_outcome/validation_test_results.xlsx',atlas_method),'WriteRowNames',true)
% writetable(abl_res_test_result_table,sprintf('output/supplemental_figures/%s_method/by_therapy/abl_res_test_results.xlsx',atlas_method),'WriteRowNames',true)
% writetable(implant_test_result_table,sprintf('output/supplemental_figures/%s_method/by_implant/implant_test_results.xlsx',atlas_method),'WriteRowNames',true)
% 
% fprintf('Statistical test results saved. (%s method)\n',atlas_method)

end

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
    
    [atlas_mni, distance_matrices] = create_distance_matrix({test_patient.coords}, region_list);
    
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

% (fix this stuff)


good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

poor_cond = [hasData_field{:}] & (strcmp(outcome_field,'poor'));
poor_outcome_patients = all_patients(poor_cond);

% apply in good outcome patients
[good_pt_dice] = localize_EZ_atlas(z_score_mat, {good_outcome_patients.roi}, {good_outcome_patients.resect}, {good_outcome_patients.coords});

% apply in poor outcome patients
[poor_pt_dice] = localize_EZ_atlas(z_score_mat, {good_outcome_patients.roi}, {good_outcome_patients.resect}, {good_outcome_patients.coords});

% quantify and compare results

%% Spatial extent of atlas -> correlation to outcome
% use modularity to quantify communities of 'abnormal' connectivity and
% find the spatial extent (mean interregional distance)

blue = [0, 0.4470, 0.7410];
beige = [254, 249, 213]./255;
red = [0.6350, 0.0780, 0.1840];

color_bar1 = [linspace(blue(1),beige(1),50)', linspace(blue(2),beige(2),50)', linspace(blue(3),beige(3),50)'];
color_bar2 = [linspace(beige(1),red(1),50)', linspace(beige(2),red(2),50)', linspace(beige(3),red(3),50)'];

color_bar = [color_bar1;color_bar2];

%save('color_bar.mat','color_bar')

color3 = [78 172 91]/255;
color4 = [103 55 155]/255;
beige = [254, 249, 213]./255;

color_bar3 = [linspace(color3(1),beige(1),50)', linspace(color3(2),beige(2),50)', linspace(color3(3),beige(3),50)'];
color_bar4 = [linspace(beige(1),color4(1),50)', linspace(beige(2),color4(2),50)', linspace(beige(3),color4(3),50)'];

color_bar_alt = [color_bar3;color_bar4];

figure(1);clf
imagesc(all_patients(51).z_scores(1).data.all,'AlphaData',~isnan(all_patients(51).z_scores(1).data.all))
colormap(color_bar)
colorbar
