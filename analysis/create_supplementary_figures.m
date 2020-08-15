%% create_supplementary_figures.m
% In this script we create tables and supplementary figures
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

%% In this section we derive distance - connectivity relationships for SEEG and ECoG separately

ECoG_indices = find([all_patients.hasData] & strcmp({all_patients.implant},'ECoG') & strcmp({all_patients.outcome},'good'));
SEEG_indices = find([all_patients.hasData] & strcmp({all_patients.implant},'SEEG') & strcmp({all_patients.outcome},'good'));

ecog_patients = all_patients(ECoG_indices);
seeg_patients = all_patients(SEEG_indices);

[curve_ECoG, dist_ecog, conn_ecog, mean_ecog] = compute_dist_reg({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect});
[curve_SEEG, dist_seeg, conn_seeg, mean_seeg] = compute_dist_reg({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect});

x_axis = [2:0.1:200];

for f = 1:5
y_ecog = (curve_ECoG(f).data.a.*x_axis.^curve_ECoG(f).data.b)+curve_ECoG(f).data.c;
y_seeg = (curve_SEEG(f).data.a.*x_axis.^curve_SEEG(f).data.b)+curve_SEEG(f).data.c;

figure(f);clf;
hold on
plot(dist_ecog(1:50:end), conn_ecog(1:50:end,f),'r.')
plot(dist_seeg(1:50:end), conn_seeg(1:50:end,f),'b.')
plot(x_axis, y_ecog,'r-')
plot(x_axis, y_seeg,'b-')
legend('ECoG','SEEG','ECoG','SEEG')
hold off

end

%% Create atlas distance matrices and save to patient specific folder

[atlas_mni, distance_matrices_ecog] = create_distance_matrix({ecog_patients.coords}, all_inds(1:90));
[atlas_mni, distance_matrices_seeg] = create_distance_matrix({seeg_patients.coords}, all_inds(1:90));

%% Find mean connectivity in each band in each patient
for f = 1:5
    figure(f);clf;
    subplot(2,1,1)
    plot(mean_ecog(:,f),'ko')
    subplot(2,1,2)
    plot(mean_seeg(:,f),'ko')
end

%% construct adjacency matrix of all good outcome patients

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
%% logistic regression stuff
        % TODO: logistic regression
        % predictors: distance between mean resected and mean non-resected 
        % z-scores, mean & variance of the distribution of all z scores together
        % response: good or poor outcome

    %     % helper function
    %     get_average_data = @(x) nanmean(x(triu(true(size(x)))));
    %     
    %     % get distance values for good outcome patients
    %     good_distances = cell2mat(cellfun(get_average_data,good_resected_z_score_results,'UniformOutput',false)) - cell2mat(cellfun(get_average_data,good_z_score_results,'UniformOutput',false));
    %     
    %     % get distance values for poor outcome patients and append them
    %     poor_distances = cell2mat(cellfun(get_average_data,poor_resected_z_score_results,'UniformOutput',false)) - cell2mat(cellfun(get_average_data,poor_z_score_results,'UniformOutput',false));
    %     
    %     % combine two arrays
    %     distances = [good_distances; poor_distances];
    %     
    %     % get mean z-scores
    %     good_means = nanmean([cell2mat(cellfun(get_data,good_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,good_z_score_results,'UniformOutput',false).')]);
    %     
    %     poor_means = nanmean([cell2mat(cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,poor_z_score_results,'UniformOutput',false).')]);
    %     
    %     z_score_means = [good_means, poor_means];
    %     
    %     % get variance of z-scores
    %     good_variances = nanvar([cell2mat(cellfun(get_data,good_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,good_z_score_results,'UniformOutput',false).')]);
    %     
    %     poor_variances = nanvar([cell2mat(cellfun(get_data,poor_resected_z_score_results,'UniformOutput',false).'); cell2mat(cellfun(get_data,poor_z_score_results,'UniformOutput',false).')]);
    %     
    %     z_score_variances = [good_variances, poor_variances];
    %     
    %     % generate array of outcomes
    %     outcomes = categorical([repmat(["good"],1,length(good_patient_indices)), repmat(["poor"],1,length(poor_patient_indices))]).';
    % 
    %     predictors = [distances,z_score_means.',z_score_variances.'];
    %     
    %     % remove any patients with NaN predictor values
    %     good_rows = ~any(isnan(predictors),2);
    %     predictors = predictors(good_rows,:);
    %     outcomes = outcomes(good_rows,:);
    %     
    %     % perform logistic regression
    %     [B,dev,stats] = mnrfit(predictors,outcomes.','Model','hierarchical');
    %     mnr_results(test_band,:) = {B,dev,stats};
    %     
    %     mdl = fitglm(predictors,outcomes,'Distribution','binomial','Link','logit');
    %     mdl_results(test_band) = {mdl};
    
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
%% do spatial extent of modules
good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

poor_cond = [hasData_field{:}] & (strcmp(outcome_field,'poor'));
poor_outcome_patients = all_patients(poor_cond);

[atlas_mni, distance_matrices_good] = create_distance_matrix({good_outcome_patients.coords}, all_inds(1:90));
[atlas_mni, distance_matrices_poor] = create_distance_matrix({poor_outcome_patients.coords}, all_inds(1:90));

[spatial_extent_good] = calculate_spatial_extent(good_patient_zscores(5).freq, distance_matrices_good);
[spatial_extent_poor] = calculate_spatial_extent(poor_patient_zscores(5).freq, distance_matrices_poor);

plotting_extent_good = [];
plotting_abnormal_good = [];

plotting_extent_poor = [];
plotting_abnormal_poor = [];

for i = 1:length(spatial_extent_good)
    plotting_abnormal_good  = [plotting_abnormal_good; spatial_extent_good(i).abnormality(:)];
    plotting_extent_good  = [plotting_extent_good; spatial_extent_good(i).distance(:)];
end

bad_distance = find(isnan(plotting_extent_good));
plotting_abnormal_good(bad_distance) = [];
plotting_extent_good(bad_distance) = [];

for i = 1:length(spatial_extent_poor)
    plotting_abnormal_poor  = [plotting_abnormal_poor; spatial_extent_poor(i).abnormality(:)];
    plotting_extent_poor  = [plotting_extent_poor; spatial_extent_poor(i).distance(:)];
end

bad_distance = find(isnan(plotting_extent_poor));
plotting_abnormal_poor(bad_distance) = [];
plotting_extent_poor(bad_distance) = [];

figure(1);clf;
hold on
plot(plotting_extent_poor,plotting_abnormal_poor,'r.')
plot(plotting_extent_good,plotting_abnormal_good,'bo')
