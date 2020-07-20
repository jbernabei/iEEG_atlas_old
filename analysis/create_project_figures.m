%% create_project_figures.m
% John Bernabei
% With assistance from Ian Ong
% Litt Laboratory

%% Set up workspace

% suppress warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% load in data from excel spreadsheet
metadata = readtable("data/atlas_project_metadata.xlsx");

% enable warnings
warning ('on','all')

good_outcome_pts = {'HUP082','HUP086','HUP088','HUP094','HUP105','HUP106','HUP111','HUP116','HUP117',...
                    'HUP125','HUP130','HUP139','HUP140','HUP150','HUP151','HUP157','HUP163','HUP164',...
                    'HUP165','HUP173','HUP177','HUP179','HUP180','HUP181','HUP185'};

poor_outcome_pts = {'HUP060','HUP075','HUP078','HUP112','HUP133','HUP138','HUP141','HUP158','HUP170','HUP171','HUP172','HUP188'};

% place patients in a struct, extracting all relevant metadata
all_patients = struct('patientID',metadata.Patient, ...
'outcome', metadata.Outcome,'conn',cell(length(metadata.Patient),1), ...
'coords',cell(length(metadata.Patient),1), ...
'roi',cell(length(metadata.Patient),1), ...
'resect',cell(length(metadata.Patient),1), ...
'hasData',cell(length(metadata.Patient),1),...
'therapy',metadata.Therapy,'implant',metadata.Implant,...
'target',metadata.Target,'laterality',metadata.Laterality,...
'lesion_status',metadata.Lesion_status,'age_onset',metadata.Age_onset,...
'age_surgery',metadata.Age_surgery,'gender',metadata.Gender,...
'hypothesis_1',metadata.Hypothesis_1,'hypothesis_2',metadata.Hypothesis_2);

% Extract atlas indices and ROIs available from atlas (here AAL116 w/WM)
fileID = fopen('localization/AAL116_WM.txt');
atlas_info = textscan(fileID,'%s %s %d');
all_inds = [double(atlas_info{3})];
all_locs = [atlas_info{2}];

% set up arrays to store data
id_field = {all_patients.patientID};
conn_field = {all_patients.conn};
coords_field = {all_patients.coords};
roi_field = {all_patients.coords};
resect_field = {all_patients.resect};
outcome_field = {all_patients.outcome};
hasData_field = {all_patients.hasData};
therapy_field = {all_patients.therapy};
implant_field = {all_patients.implant};
target_field = {all_patients.target};

% if true, the script will automatically move problematic data to another
% directory
move_files = false;

% load in data from all patients
for k = 1:length(metadata.Patient)
    folderpath = sprintf('data/%s',id_field{k});
    datapath = sprintf('%s/patient_data.mat',folderpath);
    if isfile(datapath)
        fprintf('%s: ',datapath)
        d = load(datapath);
        conn_field{k} = d.II_conn;
        if sum(sum(~isnan(d.II_conn(1).data))) == 0
            hasData_field{k} = false;
            fprintf('(connectivity data is all NaNs!)\n')
            if move_files, movefile(folderpath,'data/exclude/no_conn_data'); end
            continue
        else
            hasData_field{k} = true;
        end
        coords_field{k} = d.mni_coords;
        
        try
            resect_field{k} = d.res_elec_inds;
        catch Error
            resect_field{k} = [];
        end
        % convert all electrode coordinates to region names
        try
            [~,electrode_regions,~] = nifti_values(coords_field{1,k},'localization/AAL116_WM.nii');
            roi_field{k} = electrode_regions;
            fprintf('loaded\n')
        catch ME
            fprintf('failed to load\n')
            warning('Problem converting MNI coordinates to region labels\n(%s)',datapath, ME.identifier)
            hasData_field{k} = false;
            if move_files, movefile(folderpath,'data/exclude/out_of_bound_electrodes'); end
        end
    else
        hasData_field{k} = false;
    end
end

% place data back into main struct
[all_patients.conn] = conn_field{:};
[all_patients.coords] = coords_field{:};
[all_patients.roi] = roi_field{:};
[all_patients.resect] = resect_field{:};
[all_patients.hasData] = hasData_field{:};

fprintf('\nAll patient data loaded.')

% load in region numbers
region_list = zeros(1,90); % for the 90 AAL regions we will be using
region_names = cell(1,90);
fi = fopen("localization/AAL116_WM.txt");
for j = 1:90
    label = split(fgetl(fi));
    region_list(j) = str2double(label{3});
    region_names{j} = label{2};
end

fprintf('\nRegion list loaded.\n')

%% generate csv file with patient demographics

% get patients who have data
study_patients = all_patients([hasData_field{:}]);
good_patients = study_patients(strcmp({study_patients.outcome},'good'));
poor_patients = study_patients(strcmp({study_patients.outcome},'poor'));

% create output table
row_names = {'Total number of subjects','Age at surgery',sprintf('Mean \x00B1 SD'), ...
    'Age at diagnosis',sprintf('Mean \x00B1 STD'),'Sex','Male','Female', ...
    'Electrode type','ECoG','SEEG','Resected/ablated region','MTL','Temporal', ...
    'MFL','Frontal','Parietal/FP','Insular','MRI','Lesional','Non-lesional', ...
    'Unknown','Therapy','Ablation','Resection'};

values_to_count = {'M','F','ECoG','SEEG','MTL','Temporal','MFL', ...
    'Frontal','Parietal','Insular','Lesional','Non-Lesional','', ...
    'Ablation','Resection'};

remove_rows = [1,2,3,4,5,6,9,12,19,23];
fill_rows = 1:length(row_names);
fill_rows(remove_rows) = [];

% dictionary for target lobe names
% target_map = containers.Map({'Frontal','Temporal','FP','Parietal', ...
%     'MTL','MFL','Insular'},{'FL/PL/FPL','TL','FL/PL/FPL','FL/PL/FPL', ...
%     'TL','FL/PL/FPL','FL/PL/FPL'});
% get_target = @(x) target_map(x);

% set up table

structs = {good_patients, poor_patients};
demographic_table = cell2table(cell(length(row_names),length(structs)));
for b = 1:length(structs)
    patients = structs{b};
    
    % prepare data to be read into table
    targets = {patients.target};
    targets{strcmp(targets,'FP')} = 'Parietal';
    [patients.target] = targets{:};
   
    patient_attributes = {patients.gender;patients.gender; ...
        patients.implant;patients.implant;patients.target; ...
        patients.target;patients.target;patients.target; ...
        patients.target;patients.target;patients.lesion_status; ...
        patients.lesion_status;patients.lesion_status; ...
        patients.therapy;patients.therapy}.';
%     targets = cellfun(get_target,{patients.target},'UniformOutput',false);
%     targets_and_lateralities = join([{patients.laterality}.',targets.'],2);
%     target_counts = countcats(categorical(targets_and_lateralities));
    data_column = cell(length(row_names),1);
    % Total number of subjects
    data_column{1} = length(patients);
    % Age at surgery
    data_column{3} = sprintf('%.1f \x00B1 %.1f',nanmean([patients.age_surgery]),nanstd([patients.age_surgery]));
    % Age at diagnosis
    data_column{5} = sprintf('%.1f \x00B1 %.1f',nanmean([patients.age_onset]),nanstd([patients.age_onset]));
    
    for w = 1:length(fill_rows)
        data_column{fill_rows(w)} = sum(strcmp(patient_attributes(:,w), values_to_count{w}));
    end
    
    demographic_table(:,b) = data_column;
end

demographic_table.Properties.RowNames = row_names;
demographic_table.Properties.VariableNames = {'Good surgical outcome','Poor surgical outcome'};

writetable(demographic_table,'output/patient_demographics.xlsx','WriteRowNames',true)

fprintf('Demographic table saved.\n')

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
    colorbar
    title(sprintf('Connectivity atlas of non-resected regions in good outcome patients (band %d)',test_band),'fontsize',12)
    save_name = sprintf('output/non_resected_good_outcome_atlas_band_%d.png',test_band);
    %saveas(fig,save_name) % save plot to output folder
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
colorbar
title(sprintf('Sample sizes for each edge in non-resected regions of good outcome patients'),'fontsize',12)
save_name = sprintf('output/non_resected_good_outcome_sample_sizes.png');
%saveas(fig,save_name) % save plot to output folder

%%
%dlmwrite('output/render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
%save('output/atlas.edge','good_mean_conn','-ascii');
%BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/render_elecs.node','output/atlas.edge','final_render.mat','output/elecs.jpg')
%% Figure 2A: cross - validate out-of-bag predictions on good outcome patients
% plot atlas of non-resected regions in good-outcome patients
% fig = figure;
% set(fig,'defaultAxesTickLabelInterpreter','none');  
% fig.WindowState = 'maximized';
% imagesc(good_mean_conn)
% title_text = sprintf('Atlas of non-resected regions in good outcome patients (band %d)',test_band);
% title(title_text)
% xticks((1:length(region_names)))
% yticks((1:length(region_names)))
% xlabel('Region')
% ylabel('Region')
% xticklabels(region_names)
% yticklabels(region_names)
% xtickangle(90)
% ax = gca;
% ax.XAxis.FontSize = 6;
% ax.YAxis.FontSize = 6;
% save_name = sprintf('output/good_non_resected_atlas_%d.png',test_band, test_threshold);
% saveas(gcf,save_name) % save plot to output folder

% plot matrices showing number of samples available for each edge
% samples = {all_samples,good_samples,poor_samples};
% title_suffixes = {'all patients','good outcome patients','poor outcome patients'};
% 
% mymap = colormap('hot');
% mymap = cat(1,[0 0 0],mymap);

% for a = 1:length(samples)
%     fig = figure;
%     set(fig,'defaultAxesTickLabelInterpreter','none');
%     fig.WindowState = 'maximized';
%     imagesc(samples{a})
%     colorbar(gca);
%     colormap(mymap);
%     title_text = sprintf('Sample sizes by edge (%s)',title_suffixes{a});
%     title(title_text)
%     xticks((1:length(region_names)))
%     yticks((1:length(region_names)))
%     xlabel('Region')
%     ylabel('Region')
%     xticklabels(region_names);
%     yticklabels(region_names);
%     xtickangle(90)
%     ax = gca;
%     ax.XAxis.FontSize = 6;
%     ax.YAxis.FontSize = 6;
%     save_name = sprintf('output/samples_available_%d.png',a);
%     saveas(gcf,save_name)
% end

%% Figure 2B: cross - validate out-of-bag predictions

test_threshold = 3;

num_good_patients = sum([hasData_field{:}] & strcmp(outcome_field,'good'));
num_poor_patients = sum([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good'));
poor_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_z_score_results = cell(num_good_patients,1);
good_resected_z_score_results = cell(num_good_patients,1);

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

        line_length = fprintf('Testing good outcome patient %d of %d...', s, length(good_patient_indices));

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

        line_length = fprintf('Testing poor outcome patient %d of %d...', s, length(poor_patient_indices));

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

    histogram(poor_plot_data,'Normalization','probability','BinWidth',bin_width);
    hold on
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
end
fprintf('\n')

%% plot logistic regression results

bound = 20;
interval = 2;
jitter = 0.5;
marker_size = 100;

fig = figure;
fig.WindowState = 'maximized';
hold on

% plot 3d volume
for test_band = 1:5
    coeffs = mnr_results{test_band,1};
    f = @(x,y,z) 1./(1+exp(-(coeffs(1)+(coeffs(2)*x)+(coeffs(3)*y)+(coeffs(4)*z))));
    [x,y,z] = ndgrid(-bound/2:interval:bound/2);
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
                           'HUP141','HUP150','HUP157','HUP164','HUP165',...
                           'HUP171','HUP185','HUP188'};
                       
lobar_data = readtable('localization/lobes_aal.xlsx');

for s = 1:length(patient_hypothesis_list)
    % get test patient
    test_pt_ind = find(strcmp(id_field,test_patient_hypothesis_list{s}));
    test_patient = all_patients(test_pt_ind);

     % get connectivity atlas of test patient
     [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list, test_band);

     % get non-resected region labels of test patient
     [~, patient_roi, ~] = nifti_values(test_patient.coords,'localization/AAL116_WM.nii');

     % test atlas
     bilateral_z_score_results{s} = test_patient_conn(good_mean_conn{test_band}, good_std_conn{test_band}, region_list, patient_conn, patient_roi);

    % step 2 extract z scores of ROIs within primary/secondary hypothesis
    % regions and scores from each region to all other regions
    primary_hypothesis = test_patient.hypothesis_1;
    secondary_hypothesis = test_patient.hypothesis_2;
    %true_resection = strcat(strcat(test_patient.laterality,'_'), test_patient.target);
    
    % find which ROI are contained in primary hypothesis
    ROI_primary = contains(lobar_data{:,3},primary_hypothesis)
    ROI_secondary = contains(lobar_data{:,3},secondary_hypothesis)
    
    % find which ROI are contained in secondary hypothesis

    % step 3 generate rendering of brain with primary/secondary regions labeled 
    % by colored nodes and all edge weights rendered

    % step 4 plot boxplot of z scores for each of the four classes: in-in
    % primary, in-out primary, in-in secondary, in-out secondary

end
