% john analysis

%% set up workspace

clear all

band_names = {'broadband','alpha-theta','beta','low-gamma','high-gamma'};

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field, id_field,...
    implant_field, outcome_field, resect_field, roi_field, target_field,...
    therapy_field, region_list, region_names] = set_up_workspace(iEEG_atlas_path);

%% find z scores

test_threshold = 7;

good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good'));
poor_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor'));

for atlas_method = ["patient","not edge"]
    for test_band = 1:5
        % cross-validation of good-outcome patients
        for s = 1:length(good_patient_indices)
            test_patient = all_patients(good_patient_indices(s));
            cv_patients = all_patients(good_patient_indices);
            cv_patients(s) = [];
            
            if strcmp(test_patient.target,'Temporal') || strcmp(test_patient.target,'MTL')
               target_good(s) = 1;
            elseif strcmp(test_patient.target,'Frontal') || strcmp(test_patient.target,'MFL')
                target_good(s) = 2;
            else
                target_good(s) = 0;
            end
            
            patient_id = test_patient.patientID;

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [out_out_scores ,in_in_scores, in_out_scores, all_scores] =...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,"std",atlas_method);

            % get all good outcome patient zscores
            good_patient_zscores(test_band).freq{s} = abs(all_scores);
            good_patient_spared_zscores(test_band).freq{s} = abs(out_out_scores);
            
            good_outcome_out_out(s,test_band) = nanmedian(nanmedian(abs(out_out_scores)));
            good_outcome_in_out(s,test_band) = nanmedian(nanmedian(abs(in_out_scores)));
            good_outcome_in_in(s,test_band) = nanmedian(nanmedian(abs(in_in_scores)));
            
            % place results into all_patients
            all_patients(good_patient_indices(s)).z_scores(test_band).data.all = all_scores;
            all_patients(good_patient_indices(s)).z_scores(test_band).data.out_out = out_out_scores;
            all_patients(good_patient_indices(s)).z_scores(test_band).data.in_in = in_in_scores;
            all_patients(good_patient_indices(s)).z_scores(test_band).data.in_out = in_out_scores;

            if strcmp(test_patient.therapy,'Ablation')
                [threshold_results(s,test_band),threshold_accuracies(s,test_band)] =...
                    get_optimal_threshold(out_out_scores,in_in_scores);
            end

            fprintf(repmat('\b',1,line_length))

        end
        
        cv_patients = all_patients(good_patient_indices);

        % cross-validation of poor-outcome patients
        for s = 1:length(poor_patient_indices)
            test_patient = all_patients(poor_patient_indices(s));
            
            if strcmp(test_patient.target,'Temporal') || strcmp(test_patient.target,'MTL')
               target_poor(s) = 1;
            elseif strcmp(test_patient.target,'Frontal') || strcmp(test_patient.target,'MFL')
                target_poor(s) = 2;
            else
                target_poor(s) = 0;
            end

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [out_out_scores, in_in_scores, in_out_scores ,all_scores] = ...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,"std",atlas_method);

            % get all good outcome patient zscores
            poor_patient_zscores(test_band).freq{s} = abs(all_scores);
            poor_patient_spared_zscores(test_band).freq{s} = abs(out_out_scores);
            
            poor_outcome_out_out(s,test_band) = nanmedian(nanmedian(abs(out_out_scores)));
            poor_outcome_in_out(s,test_band) = nanmedian(nanmedian(abs(in_out_scores)));
            poor_outcome_in_in(s,test_band) = nanmedian(nanmedian(abs(in_in_scores)));
            
            % place results into all_patients
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.all = all_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.out_out = out_out_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_in = in_in_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_out = in_out_scores;

            fprintf(repmat('\b',1,line_length))

        end
        
    end
end

%% plot z score results

for f = 1
   figure(f);clf;
   hold on
   subplot(1,3,1)
   boxplot([good_outcome_in_in(:,f),[poor_outcome_in_in(:,f);NaN*ones(10,1)]])
   ylim([0 8])
   subplot(1,3,2)
   boxplot([good_outcome_in_out(:,f),[poor_outcome_in_out(:,f);NaN*ones(10,1)]])
   ylim([0 8])
   subplot(1,3,3)
   boxplot([good_outcome_out_out(:,f),[poor_outcome_out_out(:,f);NaN*ones(10,1)]])
   ylim([0 8])
   
   in_in_p_vals(f) = ranksum(good_outcome_in_in(:,f),poor_outcome_in_in(:,f));
   in_out_p_vals(f) = ranksum(good_outcome_in_out(:,f),poor_outcome_in_out(:,f));
   out_out_p_vals(f) = ranksum(good_outcome_out_out(:,f),poor_outcome_out_out(:,f));
   
   in_in_vs_in_out_good(f) = signrank(good_outcome_in_in(:,f),good_outcome_in_out(:,f));
   in_in_vs_in_out_poor(f) = signrank(poor_outcome_in_in(:,f),poor_outcome_in_out(:,f));
   
   in_in_vs_out_out_good(f) = signrank(good_outcome_in_in(:,f),good_outcome_out_out(:,f));
   in_in_vs_out_out_poor(f) = signrank(poor_outcome_in_in(:,f),poor_outcome_out_out(:,f));
   
   in_out_vs_out_out_good(f) = signrank(good_outcome_in_out(:,f),good_outcome_out_out(:,f));
   in_out_vs_out_out_poor(f) = signrank(poor_outcome_in_out(:,f),poor_outcome_out_out(:,f));
   
end
%%

in_in_p_vals
in_out_p_vals
out_out_p_vals

in_in_vs_in_out_good
in_in_vs_in_out_poor

in_in_vs_out_out_good
in_in_vs_out_out_poor

in_out_vs_out_out_good
in_out_vs_out_out_poor

%% do localization based on z score

test_band = 1;

good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

poor_cond = [hasData_field{:}] & (strcmp(outcome_field,'poor'));
poor_outcome_patients = all_patients(poor_cond);

[all_pt_dice_good] = localize_EZ_atlas(good_patient_zscores(test_band).freq, {good_outcome_patients.resect}, {good_outcome_patients.coords});

%[all_pt_dice_poor] = localize_EZ_atlas(poor_patient_zscores(test_band).freq, {poor_outcome_patients.resect}, {poor_outcome_patients.coords});

%% do node abnormality study
% lets just pick broadband for now
test_band = 1;

edge_threshold = 1;
node_threshold = 0.1;

% compute node abnormality of entire network (pre-surgical)
[good_abn, good_frac_abn] = compute_node_abnormality(good_patient_zscores(test_band).freq, edge_threshold, node_threshold);
[poor_abn, poor_frac_abn] = compute_node_abnormality(poor_patient_zscores(test_band).freq, edge_threshold, node_threshold);

[good_spared_abn, good_spared_frac_abn] = compute_node_abnormality(good_patient_spared_zscores(test_band).freq, edge_threshold, node_threshold);
[poor_spared_abn, poor_spared_frac_abn] = compute_node_abnormality(poor_patient_spared_zscores(test_band).freq, edge_threshold, node_threshold);

figure(1);clf;
boxplot([good_frac_abn', [poor_frac_abn';NaN*ones(10,1)], good_spared_frac_abn', [poor_spared_frac_abn';NaN*ones(10,1)]])
ylim([-0.1 1])
ylabel('Fraction of abnormal nodes')

ranksum(good_frac_abn,poor_frac_abn)
ranksum(good_spared_frac_abn,poor_spared_frac_abn)

signrank(good_frac_abn,good_spared_frac_abn)
signrank(poor_frac_abn,poor_spared_frac_abn)

%% spatial extent

% we use broadband CC
test_band = 1;

% get matrix of distances from MNI coordinates
[good_distances] = create_pt_distance({good_outcome_patients.coords});
[poor_distances] = create_pt_distance({poor_outcome_patients.coords});

% calculating spatial extent, which is paired values of mean internodal
% distances within module, and mean abnormality values for those modules
[good_spatial_extent] = calculate_spatial_extent(good_patient_zscores(test_band).freq, good_distances);
[poor_spatial_extent] = calculate_spatial_extent(poor_patient_zscores(test_band).freq, poor_distances);

% assemble into plottable results
num_good_pts = length(good_distances);
good_abn_all = [];
good_dist_all = [];
for pt = 1:num_good_pts
    [Y, I] = max(good_spatial_extent(pt).abnormality);
    good_abn_all = [good_abn_all, Y];
    good_dist_all = [good_dist_all, good_spatial_extent(pt).distance(I)];
end

zero_regions = find(good_dist_all==0);
good_abn_all(zero_regions) = [];
good_dist_all(zero_regions) = [];

num_poor_pts = length(poor_distances);
poor_abn_all = [];
poor_dist_all = [];
for pt = 1:num_poor_pts
    [Y, I] = max(poor_spatial_extent(pt).abnormality);
    poor_abn_all = [poor_abn_all, Y];
    poor_dist_all = [poor_dist_all, poor_spatial_extent(pt).distance(I)];
end

zero_regions = find(poor_dist_all==0);
poor_abn_all(zero_regions) = [];
poor_dist_all(zero_regions) = [];

% one regression for both good and poor outcome
[r1, m1, b1] = regression([good_dist_all,poor_dist_all] , [good_abn_all,poor_abn_all])
x1 = [min([good_dist_all,poor_dist_all]), max([good_dist_all,poor_dist_all])];

figure(1);clf;
hold on
plot(good_dist_all, good_abn_all,'bo')
plot(poor_dist_all, poor_abn_all,'ro')
plot(x1,[m1.*x1+b1],'k-')
%plot(x2,[m2.*x2+b2],'r-')
legend('good outcome','poor outcome')
xlabel('Mean internodal distance')
ylabel('Mean absolute z score')
txt = sprintf('rho = %d',r1);
t = text(50,5,txt)
txt2 = 'p-value = 1.25e-05';
t2 = text(50,4.75,txt2)

% do separate analysis of patients with known lesions versus patients
% without known lesions

%% Clinical hypothesis testing
% here we will simulate how the atlas could be used in a prospective manner

patient_hypothesis_list = {'HUP116','HUP117','HUP130','HUP138','HUP140',...
                           'HUP141','HUP157','HUP164','HUP165',...
                           'HUP171','HUP185','HUP188'};
                       
lobar_data = readtable('localization/lobes_aal.xlsx');

for s = 1:length(patient_hypothesis_list)
    % get test patient
    test_pt_ind = find(strcmp(id_field,patient_hypothesis_list{s}));
    test_patient = all_patients(test_pt_ind);

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


%% regress patient length of epilepsy with global abnormalities
test_band = 1;

all_cond = find([hasData_field{:}]&(strcmp(outcome_field,'good')|strcmp(outcome_field,'poor')));
all_data_patients = all_patients(all_cond);
all_age_onset = good_outcome_patients(1).age_onset(all_cond);
all_age_surg = good_outcome_patients(1).age_surgery(all_cond);

all_epilepsy_duration = all_age_surg-all_age_onset;

for s = 1:length(all_data_patients)
    abs_mean_abn(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.all)));
end

nan_inds = find(isnan(all_epilepsy_duration));
all_epilepsy_duration(nan_inds) = [];
abs_mean_abn(nan_inds) = [];

[r1, m1, b1] = regression(all_epilepsy_duration', abs_mean_abn)
x1 = [min(all_epilepsy_duration), max(all_epilepsy_duration)];

figure(1);clf;
hold on
plot(all_epilepsy_duration,abs_mean_abn,'k.')
plot(x1,[m1.*x1+b1],'k-')
xlabel('Duration of epilepsy')
ylabel('Mean atlas Z score')
txt = sprintf('rho = %d',r1);
t = text(40,1.6,txt)
txt2 = 'p-value = 0.001';
t2 = text(40,1.5,txt2)
hold off
%% predict patient outcome 
