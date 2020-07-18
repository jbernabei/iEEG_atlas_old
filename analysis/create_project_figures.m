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
'age_surgery',metadata.Age_surgery,'gender',metadata.Gender);

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

% functions to help convert mni coordinates to regions of interest
%f = @nifti_values;
%g = @(x) f(x,'localization/AAL116_WM.nii');

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
            movefile(folderpath,'data/exclude/no_conn_data')
            continue
        else
            hasData_field{k} = true;
        end
        coords_field{k} = d.mni_coords;
        resect_field{k} = d.res_elec_inds;
        % convert all electrode coordinates to region names
        try
            %[~,electrode_regions,~] = cellfun(g, coords_field(k), 'UniformOutput',false);
            %roi_field{k} = electrode_regions{1};
            [~,electrode_regions,~] = nifti_values(coords_field{1,k},'localization/AAL116_WM.nii');
            roi_field{k} = electrode_regions;
            fprintf('loaded\n')
        catch ME
            fprintf('failed to load\n')
            warning('Problem converting MNI coordinates to region labels\n(%s)',datapath, ME.identifier)
            hasData_field{k} = false;
            movefile(folderpath,'data/exclude/out_of_bound_electrodes')
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

%% Figure 1A: Calculate atlas
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
    fig.WindowState = 'maximized';
    imagesc(good_mean_conn{f},'AlphaData',~isnan(good_mean_conn{f}))
    set(gca,'color',0*[1 1 1]);
    set(gca,'xtick',(1:90),'xticklabel',all_locs,'fontsize',6)
    xtickangle(45)
    set(gca,'ytick',(1:90),'yticklabel',all_locs,'fontsize',6)
    colorbar
    title(sprintf('Connectivity atlas of non-resected regions in good outcome patients (band %d)',test_band),'fontsize',12)
    save_name = sprintf('output/non_resected_good_outcome_atlas_band_%d.png',test_band);
    saveas(fig,save_name) % save plot to output folder
end


%% Figure 1B: visualize number of patients with each connection
figure(1);clf
fig = gcf;
set(fig,'defaultAxesTickLabelInterpreter','none'); 
fig.WindowState = 'maximized';
imagesc(num_conn,'AlphaData',~(num_conn==0))
set(gca,'color',0*[1 1 1]);
set(gca,'xtick',(1:90),'xticklabel',all_locs)
xtickangle(45)
set(gca,'ytick',(1:90),'yticklabel',all_locs)
set(gca,'fontsize', 6)
colorbar
title(sprintf('Sample sizes for each edge in non-resected regions of good outcome patients'),'fontsize',12)
save_name = sprintf('output/non_resected_good_outcome_sample_sizes.png');
saveas(fig,save_name) % save plot to output folder

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

end
fprintf('\n')

%% Clinical hypothesis testing

% part 1
% assemble set of bilateral (RNS) patients, and test whether they have
% higher inter-hemispheric z scores compared to unilateral good outcome
% patients

% part 2
% assemble set of good outcome patients that had temporal / extratemporal
% hypotheses and test whether connections abnormal within  / between these
% regions

% part 3
% assemble set of mesial / lateral temporal patients and test whether 

%%

%% algorithm

%% Figure 4B: localize in poor outcome patients
