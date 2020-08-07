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

[curve_ECoG, dist_ecog, conn_ecog] = compute_dist_reg({ecog_patients.conn}, {ecog_patients.coords}, {ecog_patients.roi}, {ecog_patients.resect});
[curve_SEEG, dist_seeg, conn_seeg] = compute_dist_reg({seeg_patients.conn}, {seeg_patients.coords}, {seeg_patients.roi}, {seeg_patients.resect});

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