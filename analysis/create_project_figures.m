%% create_project_figures.m
% John Bernabei
% With assistance from Ian Ong
% Litt Laboratory

%% Set up workspace

% suppress warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% load in data from excel spreadsheet
metadata = readtable("data/atlas_project_metadata.xlsx");
% place patients in a struct
all_patients = struct("patientID",metadata.Patient, ...
"outcome", metadata.Outcome,"conn",cell(length(metadata.Patient),1), ...
"coords",cell(length(metadata.Patient),1), ...
"roi",cell(length(metadata.Patient),1), ...
"resect",cell(length(metadata.Patient),1), ...
"hasData",cell(length(metadata.Patient),1));

% Use AAL116WM to get white matter
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

% functions to help convert mni coordinates to regions of interest
f = @nifti_values;
g = @(x) f(x,"localization/AAL116_WM.nii");

% load in data from all patients
for k = 1:length(metadata.Patient)
    datapath = sprintf("data/%s/patient_data.mat",id_field{k});
    if isfile(datapath)
        d = load(datapath);
        conn_field{k} = d.II_conn;
        coords_field{k} = d.mni_coords;
        resect_field{k} = d.res_elec_inds;
        hasData_field{k} = true;
        % convert all electrode coordinates to region names
        [~,electrode_regions,~] = cellfun(g, coords_field(k), 'UniformOutput',false);
        roi_field{k} = electrode_regions{1};
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

fprintf("\nAll patient data loaded.")

% load in region numbers
region_list = zeros(1,90); % for the 90 AAL regions we will be using
region_names = cell(1,90);
fi = fopen("localization/AAL116_WM.txt");
for j = 1:90
    label = split(fgetl(fi));
    region_list(j) = str2double(label{3});
    region_names{j} = label{2};
end

fprintf("\nRegion list loaded.\n")

%% Figure 1A: visualize adjacency matrix (one panel in pipeline figure)

%% Figure 1B: visualize brain (one panel in pipeline figure)

%% Figure 2A: anatomical analysis

% all tests will be run on this band
test_band = 3;

% run all patients in atlas
cond = [hasData_field{:}];
[mean_conn, std_conn, all_samples] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list, test_band);

% run all good outcome patients in atlas
cond = [hasData_field{:}] & strcmp(outcome_field,'good');
[good_mean_conn, good_std_conn, good_samples] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list, test_band);

% run all poor outcome patients in atlas
cond = [hasData_field{:}] & strcmp(outcome_field,'poor');
[poor_mean_conn, poor_std_conn, poor_samples] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list, test_band);

% calculate difference between good and poor outcome means
difference_mean_conn = good_mean_conn - poor_mean_conn;

% save results to output folder
save('output/figure_2A_data.mat','mean_conn','std_conn','good_mean_conn', ...
'good_std_conn','poor_mean_conn', 'poor_std_conn', 'difference_mean_conn')

% plot atlas of non-resected regions in good-outcome patients
fig = figure;
set(fig,'defaultAxesTickLabelInterpreter','none');  
fig.WindowState = 'maximized';
imagesc(good_mean_conn)
title_text = sprintf('Atlas of non-resected regions in good outcome patients (band %d)',test_band);
title(title_text)
xticks((1:length(region_names)))
yticks((1:length(region_names)))
xlabel('Region')
ylabel('Region')
xticklabels(region_names)
yticklabels(region_names)
xtickangle(90)
ax = gca;
ax.XAxis.FontSize = 6;
ax.YAxis.FontSize = 6;
save_name = sprintf('output/good_non_resected_atlas_%d.png',test_band);
%saveas(gcf,save_name) % save plot to output folder

% plot matrices showing number of samples available for each edge
samples = {all_samples,good_samples,poor_samples};
title_suffixes = {'all patients','good outcome patients','poor outcome patients'};

mymap = colormap('hot');
mymap = cat(1,[0 0 0],mymap);

for a = 1:length(samples)
    fig = figure;
    set(fig,'defaultAxesTickLabelInterpreter','none');
    fig.WindowState = 'maximized';
    imagesc(samples{a})
    colorbar(gca);
    colormap(mymap);
    title_text = sprintf('Number of samples available for %s',title_suffixes{a});
    title(title_text)
    xticks((1:length(region_names)))
    yticks((1:length(region_names)))
    xlabel('Region')
    ylabel('Region')
    xticklabels(region_names);
    yticklabels(region_names);
    xtickangle(90)
    ax = gca;
    ax.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 6;
    save_name = sprintf('output/samples_available_%d.png',a);
    saveas(gcf,save_name)
end

%% Figure 2B: cross - validate out-of-bag predictions

num_good_patients = sum([hasData_field{:}] & strcmp(outcome_field,'good'));
num_poor_patients = sum([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good'));
poor_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_z_score_results = cell(num_good_patients,1);
good_resected_z_score_results = cell(num_good_patients,1);

% cross-validation of good-outcome patients
for s = 1:length(good_patient_indices)
    test_patient = all_patients(good_patient_indices(s));
    cv_patients = all_patients(good_patient_indices);
    cv_patients(s) = [];
    
    fprintf('\nTesting patient %d of %d:', s, length(good_patient_indices))
    
    % get connectivity atlas of excluded patients
    [mean_conn, std_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band);
    
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
end

good_z_score_mean = nanmean(cat(3,good_z_score_results{:}),3);
good_resected_z_score_mean = nanmean(cat(3,good_resected_z_score_results{:}),3);

% save results to output folder
save('output/figure_2B_good_data.mat','good_z_score_mean','good_resected_z_score_mean')

% repeat cross-validation for poor outcome patients
poor_z_score_results = cell(num_poor_patients,1);
poor_resected_z_score_results = cell(num_poor_patients,1);

% calculate atlas for good outcome patients only
% this serves as the "model" for the cross-validation
cv_patients = all_patients(good_patient_indices);
[mean_conn, std_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band);

% cross-validation of poor-outcome patients
for s = 1:length(poor_patient_indices)
    test_patient = all_patients(poor_patient_indices(s));
    
    fprintf('\nTesting patient %d of %d:', s, length(poor_patient_indices))
    
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
end

poor_z_score_mean = mean(cat(3,poor_z_score_results{:}),3,'omitnan');
poor_resected_z_score_mean = mean(cat(3,poor_resected_z_score_results{:}),3,'omitnan');

% save results to output folder
save('output/figure_2B_poor_data.mat','poor_z_score_mean','poor_resected_z_score_mean')

% plot AVERAGED z-score results for GOOD outcome patients
bin_width = 0.2;
figure
% remove outliers and bottom triangle of data
good_plot_data = rmoutliers(good_z_score_mean(triu(true(size(good_z_score_mean)))));
good_plot_data = good_plot_data(~isnan(good_plot_data) & ~isinf(good_plot_data));

good_resected_plot_data = rmoutliers(good_resected_z_score_mean);
good_resected_plot_data = good_resected_plot_data(~isnan(good_resected_plot_data) & ~isinf(good_resected_plot_data));

histogram(good_plot_data,'Normalization','probability','BinWidth',bin_width);
hold on
histogram(good_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
title({'Z-scores of averaged connectivity strengths','in good outcome patients'})
ylabel('Density')
xlabel('Z-score')
% draw lines representing the medians of both groups
xline(median(good_plot_data),'b','LineWidth',2);
xline(median(good_resected_plot_data),'r','LineWidth',2);
legend('Non-resected regions','Resected regions','Non-resected median','Resected median')
legend('Location','northeast','Box','off')
% set(gca,'YScale','log')
save_name = sprintf('output/avg_good_z_score_histogram.png');
saveas(gcf,save_name) % save plot to output folder
hold off

% plot AVERAGED z-score results for POOR outcome patients
figure
% remove outliers and bottom triangle of data
poor_plot_data = rmoutliers(poor_z_score_mean(triu(true(size(poor_z_score_mean)))));
poor_plot_data = poor_plot_data(~isnan(poor_plot_data) & ~isinf(poor_plot_data));

poor_resected_plot_data = rmoutliers(poor_resected_z_score_mean);
poor_resected_plot_data = poor_resected_plot_data(~isnan(poor_resected_plot_data) & ~isinf(poor_resected_plot_data));

histogram(poor_plot_data,'Normalization','probability','BinWidth',bin_width);
hold on
histogram(poor_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
title({'Z-scores of averaged connectivity strengths','in poor outcome patients'})
ylabel('Density')
xlabel('Z-score')
% draw lines representing the medians of both groups
xline(median(poor_plot_data),'b','LineWidth',2);
xline(median(poor_resected_plot_data),'r','LineWidth',2);
legend('Non-resected regions','Resected regions','Non-resected median','Resected median')
legend('Location','northeast','Box','off')
% set(gca,'YScale','log')
save_name = sprintf('output/avg_poor_z_score_histogram.png');
saveas(gcf,save_name) % save plot to output folder
hold off

% plot ALL z-score results for GOOD outcome patients
get_data = @(x) x(triu(true(size(x))));
bin_width = 0.2;
figure
% remove outliers and bottom triangle of data
good_plot_data = cellfun(get_data,good_z_score_results,'UniformOutput',false);
good_plot_data = cell2mat(good_plot_data);
good_plot_data = good_plot_data(:);
good_plot_data = rmoutliers(good_plot_data(~isnan(good_plot_data) & ~isinf(good_plot_data)));

good_resected_plot_data = cellfun(get_data,good_resected_z_score_results,'UniformOutput',false);
good_resected_plot_data = cell2mat(good_resected_plot_data);
good_resected_plot_data = good_resected_plot_data(:);
good_resected_plot_data = rmoutliers(good_resected_plot_data(~isnan(good_resected_plot_data) & ~isinf(good_resected_plot_data)));

histogram(good_plot_data,'Normalization','probability','BinWidth',bin_width);
hold on
histogram(good_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
title('Z-scores of connectivity strengths in good outcome patients')
ylabel('Density')
xlabel('Z-score')
% draw lines representing the medians of both groups
xline(median(good_plot_data),'b','LineWidth',2);
xline(median(good_resected_plot_data),'r','LineWidth',2);
legend('Non-resected regions','Resected regions','Non-resected median','Resected median')
legend('Location','northeast','Box','off')
% set(gca,'YScale','log')
save_name = sprintf('output/all_good_z_score_histogram.png');
saveas(gcf,save_name) % save plot to output folder
hold off

% plot ALL z-score results for POOR outcome patients
figure
% remove outliers and bottom triangle of data
poor_plot_data = cellfun(get_data,poor_z_score_results,'UniformOutput',false);
poor_plot_data = cell2mat(poor_plot_data);
poor_plot_data = poor_plot_data(:);
poor_plot_data = rmoutliers(poor_plot_data(~isnan(poor_plot_data) & ~isinf(poor_plot_data)));

poor_resected_plot_data = cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false);
poor_resected_plot_data = cell2mat(poor_resected_plot_data);
poor_resected_plot_data = poor_resected_plot_data(:);
poor_resected_plot_data = rmoutliers(poor_resected_plot_data(~isnan(poor_resected_plot_data) & ~isinf(poor_resected_plot_data)));

histogram(poor_plot_data,'Normalization','probability','BinWidth',bin_width);
hold on
histogram(poor_resected_plot_data,'Normalization','probability','BinWidth',bin_width); % specify data and number of bins
title('Z-scores of connectivity strengths in poor outcome patients')
ylabel('Density')
xlabel('Z-score')
% draw lines representing the medians of both groups
xline(median(poor_plot_data),'b','LineWidth',2);
xline(median(poor_resected_plot_data),'r','LineWidth',2);
legend('Non-resected regions','Resected regions','Non-resected median','Resected median')
legend('Location','northeast','Box','off')
% set(gca,'YScale','log')
save_name = sprintf('output/all_poor_z_score_histogram.png');
saveas(gcf,save_name) % save plot to output folder
hold off

%% Figure 2C: comparison against distance

%% Figure 3A: Distinguish between epilepsy type (mesial temporal vs neocortical)
% temporal vs extratemporal

%% Figure 3B: Relationship of connectivity to duration of epilepsy

%% Figure 3C: Relationship of connectivity to whether patient has generalized SZ or not

%% Figure 4A: localize in good outcome patients
% make figure based upon prior calculations

%% Supplement: choice of atlas / segmentation

%% Figure 4B: localize in poor outcome patients
for s = 1:length(all_poor_patients)
    
end