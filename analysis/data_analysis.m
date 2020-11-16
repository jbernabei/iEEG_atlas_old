% john analysis

%% set up workspace

clear all

band_names = {'broadband','alpha-theta','beta','low-gamma','high-gamma'};

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, coords_field, hasData_field, id_field,...
    implant_field, outcome_field, resect_field, roi_field, target_field,...
    therapy_field, region_list, region_name, lesion_field,...
    sz_field] = set_up_workspace(iEEG_atlas_path);

load color_bar
load color_bar_alt
%% get basic atlas info
test_band = 1;

test_threshold = 3;
good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good'));
cv_patients = all_patients(good_patient_indices);
[mean_conn, std_conn, num_samples, sem_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

num_samples(num_samples==0) = NaN;
median_connectivity(num_samples==0) = NaN;

% find which roi belong to which lobes 
lobe_table = readtable('lobes_aal.xlsx')

figure(1);clf;
subplot(1,2,2)
imagesc(mean_conn)
title('median connectivity')
colormap(color_bar)
colorbar
subplot(1,2,1)
imagesc(num_samples)
title('number of samples')
colorbar

%% find z scores

test_threshold = 7; % main results for = 7

good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good'));
poor_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor'));

for atlas_method = ["patient","not edge"]
    for test_band = 1:5
        % cross-validation of good-outcome patients
        for s = 1:length(good_patient_indices)
            test_patient = all_patients(good_patient_indices(s));
            cv_patients = all_patients(good_patient_indices);
            cv_patients(s) = [];
            
            if strcmp(test_patient.therapy,'Resection')
               target_good(s) = 1;
            elseif strcmp(test_patient.therapy,'Ablation')
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
            
            good_outcome_out_out(s,test_band) = nanmedian(nanmedian((out_out_scores)));
            good_outcome_in_out(s,test_band) = nanmedian(nanmedian((in_out_scores)));
            good_outcome_in_in(s,test_band) = nanmedian(nanmedian((in_in_scores)));
            
            good_outcome_out_out_2(s,test_band) = nanvar(out_out_scores(:));
            good_outcome_in_out_2(s,test_band) = nanvar(in_out_scores(:));
            good_outcome_in_in_2(s,test_band) = nanvar(in_in_scores(:));
            
            good_outcome_all_var(s,test_band) = nanvar(all_scores(:));
            
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
            
            if strcmp(test_patient.therapy,'Resection') 
               target_poor(s) = 1;
            elseif strcmp(test_patient.therapy,'Ablation') 
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
            
            poor_outcome_out_out(s,test_band) = nanmedian(nanmedian((out_out_scores)));
            poor_outcome_in_out(s,test_band) = nanmedian(nanmedian((in_out_scores)));
            poor_outcome_in_in(s,test_band) = nanmedian(nanmedian((in_in_scores)));
            
            poor_outcome_out_out_2(s,test_band) = nanvar(out_out_scores(:));
            poor_outcome_in_out_2(s,test_band) = nanvar(in_out_scores(:));
            poor_outcome_in_in_2(s,test_band) = nanvar(in_in_scores(:));
            
            poor_outcome_all_var(s,test_band) = nanvar(all_scores(:));
            
            % place results into all_patients
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.all = all_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.out_out = out_out_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_in = in_in_scores;
            all_patients(poor_patient_indices(s)).z_scores(test_band).data.in_out = in_out_scores;

            fprintf(repmat('\b',1,line_length))

        end
        
    end
end

%% set up indices
%ablation / resection
g_abl_placehold_inds = find([all_patients.hasData] & strcmp({all_patients.outcome},'good') & strcmp({all_patients.therapy},'Ablation'));
gpai = ismember(good_patient_indices,g_abl_placehold_inds);

g_res_placehold_inds = find([all_patients.hasData] & strcmp({all_patients.outcome},'good') & strcmp({all_patients.therapy},'Resection'));
gpri = ismember(good_patient_indices,g_res_placehold_inds);

p_abl_placehold_inds = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor') & strcmp({all_patients.therapy},'Ablation'));
ppai= ismember(poor_patient_indices,p_abl_placehold_inds);

p_res_placehold_inds = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor') & strcmp({all_patients.therapy},'Resection'));
ppri = ismember(poor_patient_indices,p_res_placehold_inds);

%% plot z score results all

color1 = [0, 0.4470, 0.7410];
colormid = [254, 249, 213]./255;
color2 = [0.6350, 0.0780, 0.1840];

abs_good_outcome_in_in = abs(good_outcome_in_in);
abs_good_outcome_in_out = abs(good_outcome_in_out);
abs_good_outcome_out_out = abs(good_outcome_out_out);
abs_poor_outcome_in_in = abs(poor_outcome_in_in);
abs_poor_outcome_in_out = abs(poor_outcome_in_out);
abs_poor_outcome_out_out = abs(poor_outcome_out_out);

for f = 1
   figure(f);clf;
   
   subplot(1,2,1)
   plot_data = abs([[good_outcome_in_in(:,f),[poor_outcome_in_in(:,f);NaN*(ones(size(good_outcome_in_in,1)-size(poor_outcome_in_in,1),1))]],...
       [good_outcome_in_out(:,f),[poor_outcome_in_out(:,f);NaN*(ones(size(good_outcome_in_in,1)-size(poor_outcome_in_in,1),1))]],...
       [good_outcome_out_out(:,f),[poor_outcome_out_out(:,f);NaN*(ones(size(good_outcome_in_in,1)-size(poor_outcome_in_in,1),1))]]])


   boxplot(plot_data)
   
   h = findobj(gca,'Tag','Box');
   colors = [color2; color1; color2; color1; color2; color1];
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    
    legend('Poor outcome','Good outcome')
    
    subplot(1,2,2)
   plot_data = [[good_outcome_all_var(:,f),[poor_outcome_all_var(:,f);NaN*(ones(size(good_outcome_in_in,1)-size(poor_outcome_in_in,1),1))]]]

    boxplot(plot_data)
    
    h = findobj(gca,'Tag','Box');
   colors = [color2; color1; color2; color1; color2; color1];
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
   
   in_in_p_vals(f) = ranksum(abs_good_outcome_in_in(:,f),abs_poor_outcome_in_in(:,f));
   in_out_p_vals(f) = ranksum(abs_good_outcome_in_out(:,f),abs_poor_outcome_in_out(:,f));
   out_out_p_vals(f) = ranksum(abs_good_outcome_out_out(:,f),abs_poor_outcome_out_out(:,f));
   
   in_in_vs_in_out_good(f) = signrank(abs_good_outcome_in_in(:,f),abs_good_outcome_in_out(:,f));
   in_in_vs_in_out_poor(f) = signrank(abs_poor_outcome_in_in(:,f),abs_poor_outcome_in_out(:,f));
   
   in_in_vs_out_out_good(f) = signrank(abs_good_outcome_in_in(:,f),abs_good_outcome_out_out(:,f));
   in_in_vs_out_out_poor(f) = signrank(abs_poor_outcome_in_in(:,f),abs_poor_outcome_out_out(:,f));
   
   in_out_vs_out_out_good(f) = signrank(abs_good_outcome_in_out(:,f),abs_good_outcome_out_out(:,f));
   in_out_vs_out_out_poor(f) = signrank(abs_poor_outcome_in_out(:,f),abs_poor_outcome_out_out(:,f));
   
    pvals_2(f) = ranksum(good_outcome_all_var(:,f),poor_outcome_all_var(:,f));
   
end

% want to also split up ECoG vs SEEG for this analysis, and resection vs
% ablation


%% extract p values for basic z score results

in_in_p_vals
in_out_p_vals
out_out_p_vals

in_in_vs_in_out_good
in_in_vs_in_out_poor

in_in_vs_out_out_good
in_in_vs_out_out_poor

in_out_vs_out_out_good
in_out_vs_out_out_poor

pvals_2


%% logistic regression predict outcome based on mean & variance of atlas Z scores
% use LOOCV
feat_matrix = []
all_patient_indices = find([all_patients.hasData]);
for pt = 1:length(all_patient_indices)
    patient_index = all_patient_indices(pt);
    
    %feat1 = nanmean(all_patients(patient_index).z_scores(1).data.in_in(:));
    feat1 = nanvar(all_patients(patient_index).z_scores(1).data.in_in(:));
    feat2 = nanvar(all_patients(patient_index).z_scores(1).data.in_out(:));
    %feat5 = nanmean(all_patients(patient_index).z_scores(1).data.out_out(:));
    feat3 = nanvar(all_patients(patient_index).z_scores(1).data.out_out(:));
    
    feat_matrix(pt,:) = [feat1, feat2, feat3];
    if strcmp(outcome_field{pt},'poor')
        label_matrix(pt,:) = 0;
    else
        label_matrix(pt,:) = 1;
    end
end

label_matrix = label_matrix+1;

for i = 1:64
    X_train = feat_matrix;
    Y_train = label_matrix;
    
    X_train(i,:) = [];
    Y_train(i,:) = [];
    
    X_test = feat_matrix(i,:);
    Y_test = label_matrix(i,:);
    
    B = mnrfit(X_train,Y_train);
    pihat = mnrval(B,X_test);
    pred_scores(i,:) = pihat;
end

labels = pt_outcomes2;
posclass = 2;
[X_AUC,Y_AUC,T,AUC] = perfcurve(labels,pred_scores(:,1),posclass);

figure(1);clf;
hold on
plot(X_AUC,Y_AUC,'b-','LineWidth',2)
plot([0,1],[0,1],'k-.')
%legend(sprintf('Pre-ictal epochs, AUC = %.2f',AUC_pre),sprintf('Ictal epochs, AUC = %.2f',AUC_ict),sprintf('Both epochs, AUC = %.2f',AUC_both),'Chance','Location','southeast')
xlabel('False positive rate') 
ylabel('True positive rate')
title(sprintf('AUC = %.2f',AUC))

%% do localization based on z score

% @@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@ doesn't work @@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@

test_band = 1;

good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

poor_cond = [hasData_field{:}] & (strcmp(outcome_field,'poor'));
poor_outcome_patients = all_patients(poor_cond);

[pt_dice_good,all_pt_dice_good] = localize_EZ_atlas({good_outcome_patients.z_scores}, {good_outcome_patients.resect}, {good_outcome_patients.coords});

[pt_dice_poor,all_pt_dice_poor] = localize_EZ_atlas({poor_outcome_patients.z_scores}, {poor_outcome_patients.resect}, {poor_outcome_patients.coords});

figure(1);clf;
boxplot([all_pt_dice_good',[all_pt_dice_poor';NaN*ones((length(good_outcome_patients)-length(poor_outcome_patients)),1)]])

%% do node abnormality study
% lets just pick broadband for now
test_band = 1;

color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];

edge_threshold = 1;
node_threshold = 0.1;

% compute node abnormality of entire network (pre-surgical)
[good_abn, good_frac_abn] = compute_node_abnormality(good_patient_zscores(test_band).freq, edge_threshold, node_threshold);
[poor_abn, poor_frac_abn] = compute_node_abnormality(poor_patient_zscores(test_band).freq, edge_threshold, node_threshold);

[good_spared_abn, good_spared_frac_abn] = compute_node_abnormality(good_patient_spared_zscores(test_band).freq, edge_threshold, node_threshold);
[poor_spared_abn, poor_spared_frac_abn] = compute_node_abnormality(poor_patient_spared_zscores(test_band).freq, edge_threshold, node_threshold);

figure(1);clf;
boxplot([good_frac_abn', [poor_frac_abn';NaN*ones((length(good_frac_abn)-length(poor_frac_abn)),1)], good_spared_frac_abn', [poor_spared_frac_abn';NaN*ones((length(good_frac_abn)-length(poor_frac_abn)),1)]])
ylim([-0.1 1])
ylabel('Fraction of abnormal nodes')
h = findobj(gca,'Tag','Box');
colors = [color2; color1; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end


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

% zero_regions = find(good_dist_all==0);
% good_abn_all(zero_regions) = [];
% good_dist_all(zero_regions) = [];

num_poor_pts = length(poor_distances);
poor_abn_all = [];
poor_dist_all = [];
for pt = 1:num_poor_pts
    [Y, I] = max(poor_spatial_extent(pt).abnormality);
    poor_abn_all = [poor_abn_all, Y];
    poor_dist_all = [poor_dist_all, poor_spatial_extent(pt).distance(I)];
end

% zero_regions = find(poor_dist_all==0);
% poor_abn_all(zero_regions) = [];
% poor_dist_all(zero_regions) = [];

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
txt = 'R = -0.62'%sprintf('rho = %d',r1);
t = text(50,5,txt)
txt2 = 'p-value = 1.25e-05';
t2 = text(50,4.75,txt2)

% should find whether this most abnormal module localized to resection
% zone...

%% do separate analysis of patients with known lesions versus patients
% without known lesions

% find indices of good outcome patients that are lesional / nonlesional
first_cond = find([hasData_field{:}] & (strcmp(outcome_field,'good')));
second_cond = find([hasData_field{:}] & (strcmp(outcome_field,'good')) & strcmp(lesion_field,'Lesional'));
lesional_good = find(ismember(first_cond, second_cond));
non_lesional_good = find(ismember(first_cond, second_cond)==0);

% find indices of poor outcome patients that are lesional / nonlesional
first_cond = find([hasData_field{:}] & (strcmp(outcome_field,'poor')));
second_cond = find([hasData_field{:}] & (strcmp(outcome_field,'poor')) & strcmp(lesion_field,'Lesional'));
lesional_poor = find(ismember(first_cond, second_cond));
non_lesional_poor = find(ismember(first_cond, second_cond)==0);


all_lesional_dist = [good_dist_all(lesional_good), poor_dist_all(lesional_poor)];
all_lesional_abn = [good_abn_all(lesional_good), poor_abn_all(lesional_poor)];

all_non_lesional_dist = [good_dist_all(non_lesional_good), poor_dist_all(non_lesional_poor)];
all_non_lesional_abn = [good_abn_all(non_lesional_good), poor_abn_all(non_lesional_poor)];

figure(1);clf;
hold on
plot(all_lesional_dist, all_lesional_abn,'bo')
plot(all_non_lesional_dist, all_non_lesional_abn,'ro')
plot(x1,[m1.*x1+b1],'k-')
%plot(x2,[m2.*x2+b2],'r-')
legend('Lesional','Non-Lesional')
xlabel('Mean internodal distance')
ylabel('Mean absolute z score')
txt = 'R = -0.62'%sprintf('rho = %d',r1);
t = text(50,5,txt)
txt2 = 'p-value = 1.25e-05';
t2 = text(50,4.75,txt2)
hold off

NL_HS_MCD = [1 2 3 5 6 7 14 17];
NL_GL = [4 15 18 25 27 28 29];

%% regress patient length of epilepsy with 
% a. all global abnormalities
% b. in-in abnormalities
% c. in-out abnormalities
% d. out-out abnormalities

good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

test_band = 1;

all_cond = find([hasData_field{:}]&(strcmp(outcome_field,'good')|strcmp(outcome_field,'poor')));
all_data_patients = all_patients(all_cond);
all_age_onset = good_outcome_patients(1).age_onset(all_cond);
all_age_surg = good_outcome_patients(1).age_surgery(all_cond);

all_epilepsy_duration = all_age_surg-all_age_onset;

for s = 1:length(all_data_patients)
    abs_mean_abn_all(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.all)));
    abs_mean_abn_in_in(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.in_in)));
    abs_mean_abn_in_out(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.in_out)));
    abs_mean_abn_out_out(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.out_out)));
end

nan_inds = find(isnan(all_epilepsy_duration));
all_epilepsy_duration(nan_inds) = [];
all_age_surg(nan_inds) = [];

abs_mean_abn_all(nan_inds) = [];
abs_mean_abn_in_in(nan_inds) = [];
abs_mean_abn_in_out(nan_inds) = [];
abs_mean_abn_out_out(nan_inds) = [];

features = [all_epilepsy_duration]%,all_age_surg];

f1 = fitlm(features, abs_mean_abn_all');
p1 = f1.Coefficients{2,4};
b1 = f1.Coefficients{1,1}; m1a = f1.Coefficients{2,1}; %m1b = f1.Coefficients{3,3};

f2 = fitlm(features, abs_mean_abn_in_in');
p2 = f2.Coefficients{2,4};
b2 = f2.Coefficients{1,1}; m2a = f2.Coefficients{2,1}; %m2b = f2.Coefficients{3,3};

f3 = fitlm(features, abs_mean_abn_in_out');
p3 = f3.Coefficients{2,4};
b3 = f3.Coefficients{1,1}; m3a = f3.Coefficients{2,1}; %m3b = f3.Coefficients{3,3};

f4 = fitlm(features, abs_mean_abn_out_out');
p4 = f4.Coefficients{2,4};
b4 = f4.Coefficients{1,1}; m4a = f4.Coefficients{2,1}; %m4b = f4.Coefficients{3,3};

x1 = [min(all_epilepsy_duration), max(all_epilepsy_duration)];
%x2 = [min(all_age_surg), max(all_age_surg)];

%
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
purple = [103 55 155]/255;

figure(1);clf;
subplot(1,3,1)
hold on
plot(all_epilepsy_duration,abs_mean_abn_in_in,'Color',red,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m2a.*x1+b2],'k-','LineWidth',2)
txt2 = sprintf('p = %d',p2)
%t2 = text(40,1.5,txt2)
title('IN-IN connections')
xlabel('Duration between first seizure and surgery (years)')
ylabel('Mean atlas Z score')
hold off
subplot(1,3,2)
hold on
plot(all_epilepsy_duration,abs_mean_abn_in_out,'Color',purple,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m3a.*x1+b3],'k-','LineWidth',2)
txt3 = sprintf('p = %d',p3)
%t3 = text(40,1.5,txt3)
title('IN-OUT connections')
xlabel('Duration between first seizure and surgery (years)')
ylabel('Mean atlas Z score')
hold off
subplot(1,3,3)
hold on
plot(all_epilepsy_duration,abs_mean_abn_out_out,'Color',blue,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m4a.*x1+b4],'k-','LineWidth',2)
txt4 = sprintf('p = %d',p4)
%t4 = text(40,1.5,txt4)
title('OUT-OUT connections')
xlabel('Duration between first seizure and surgery (years)')
ylabel('Mean atlas Z score')
hold off

% figure 2: all 3 on same axes
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
purple = [103 55 155]/255;

% figure(2);clf;
% hold on
% plot(all_epilepsy_duration,abs_mean_abn_in_in,'Color',[red+(1-red)*0.5],'Marker','.','LineStyle','none','MarkerSize',24)
% 
% plot(all_epilepsy_duration,abs_mean_abn_in_out,'Color',[purple+(1-purple)*0.5],'Marker','.','LineStyle','none','MarkerSize',24)
% 
% plot(all_epilepsy_duration,abs_mean_abn_out_out,'Color',[blue+(1-blue)*0.5],'Marker','.','LineStyle','none','MarkerSize',24)
% 
% plot(x1,[m2.*x1+b2],'Color',red,'LineWidth',2)
% plot(x1,[m3.*x1+b3],'Color',purple,'LineWidth',2)
% plot(x1,[m4.*x1+b4],'Color',blue,'LineWidth',2)
% xlabel('Duration between first seizure and surgery (years)')
% ylabel('Mean atlas Z score')
% legend('IN-IN','IN-OUT','OUT-OUT')
% hold off

%% test whether different global z score in patients with/without generalized seizures

%%%%%%% how to assess if there are generalized sz? %%%%%%%%
%%%% were recorded in emu? %%%% ever had a generalized sz? %%%% had before
%%%% AED?

test_band = 5;

yes_cond = find([hasData_field{:}]&(strcmp(outcome_field,'good')|strcmp(outcome_field,'poor'))& strcmp(sz_field,'yes'));
no_cond = find([hasData_field{:}]&(strcmp(outcome_field,'good')|strcmp(outcome_field,'poor'))& strcmp(sz_field,'no'));

yes_patients = all_patients(yes_cond);
no_patients = all_patients(no_cond);

for s = 1:length(yes_patients)
    yes_mean_abn(s) = nanmean(nanmean(((all_patients(yes_cond(s)).z_scores(test_band).data.in_in))));
end

for s = 1:length(no_patients)
    no_mean_abn(s) = nanmean(nanmean(((all_patients(no_cond(s)).z_scores(test_band).data.in_in))));
end

figure(1);clf;
boxplot([yes_mean_abn', [no_mean_abn';NaN*ones((length(yes_mean_abn)-length(no_mean_abn)),1)]])

ranksum(yes_mean_abn,no_mean_abn)
%% Clinical hypothesis testing
% here we will simulate how the atlas could be used in a prospective manner

% @@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@ doesn't work @@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@

patient_hypothesis_list = {'HUP116','HUP117','HUP138','HUP140',...
                           'HUP141','HUP157','HUP164','HUP165',...
                           'HUP171','HUP173','HUP181','HUP185','HUP188'};
                       
true_targets = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1];
                       
lobar_data = readtable('localization/lobes_aal.xlsx');

for s = [1:length(patient_hypothesis_list)]
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
    
    zscores = test_patient.z_scores(test_band).data.all;
    
    [in_in_1, in_in_2] = patient_hypothesis_test(zscores, test_patient.roi, region_list, ROI_primary, ROI_secondary);
   
%     % in-in primary
%     in_in_1 = test_z_score_results{s}(ROI_primary,ROI_primary);
%     
%     % in-out primary
%     in_out_1 = test_z_score_results{s}(ROI_primary,setdiff([1:90],ROI_primary));
%     
%     % in-in secondary
%     in_in_2 = test_z_score_results{s}(ROI_secondary,ROI_secondary);
%     
%     % in-out secondary
%     in_out_2 = test_z_score_results{s}(ROI_primary,setdiff([1:90],ROI_secondary));
%     
%     % out-out
%     out_out = test_z_score_results{s}(setdiff([1:90],[ROI_primary,ROI_secondary]),setdiff([1:90],[ROI_primary,ROI_secondary]));
%     
%     % step 3 generate rendering of brain with primary/secondary regions labeled 
%     % by colored nodes and all edge weights rendered
%     
%     [atlas_mni, distance_matrices] = create_distance_matrix({test_patient.coords}, region_list);
%     
%     node_color = zeros(90,1);
%     node_color(ROI_primary,1) = -1; % sets to blue
%     node_color(ROI_secondary,1) = 1; % sets to red
%     
%     % create final electrode matrix
%     final_elec_matrix = [atlas_mni, node_color, ones(90,1)];
%     
%     unsampled_roi = find(isnan(final_elec_matrix(:,1)));
%     
%     final_elec_matrix(unsampled_roi,:) = [];
%     
%     % get rid of unsampled regions
%     test_z_score_results{s}(unsampled_roi,:) = [];
%     test_z_score_results{s}(:,unsampled_roi) = [];
%     
%     % zero out any NaNs
%     test_z_score_results{s}(find(isnan(test_z_score_results{s}))) = 0;
%     
%     z_score_results = test_z_score_results{s};
%     
%     dlmwrite('output/hypothesis_test.node',final_elec_matrix,'delimiter',' ','precision',5)
%     save('output/hypothesis_test.edge','z_score_results','-ascii');
%     BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/hypothesis_test.node','output/hypothesis_test.edge','output/hypothesis_test_options.mat',sprintf('output/hypothesis_test_%s.jpg',patient_hypothesis_list{s}))

    % step 4 plot boxplot of z scores for each of the four classes: in-in
    % primary, in-out primary, in-in secondary, in-out secondary
    
    boxplot_data = NaN.*zeros(90*90,2);
    %plot_labels = {'IN-IN Region 1','IN-OUT Region 1','IN-IN Region 2','IN-OUT Region 2','OUT-OUT'};
    
    boxplot_data(1:length(in_in_1(:)),1) = in_in_1(:);
    boxplot_data(1:length(in_in_2(:)),2) = in_in_2(:);
    %boxplot_data(1:length(in_in_2(:)),3) = in_in_2(:);
    %boxplot_data(1:length(in_out_2(:)),4) = in_out_2(:); 
    %boxplot_data(1:length(out_out(:)),5) = out_out(:);
    
    % create boxplot
    figure(s);clf;
    boxplot(boxplot_data)%,'labels',plot_labels)
    
    pval(s) = ranksum(boxplot_data(:,1),boxplot_data(:,2))
    
    if (nanmedian(abs(in_in_1(:))) > nanmedian(abs(in_in_2(:)))) %&& (pval(s)<0.05)
        hyp_pred(s) = 1;
    elseif (nanmedian(abs(in_in_1(:))) < nanmedian(abs(in_in_2(:)))) %&& (pval(s)<0.05)
        hyp_pred(s) = 2;
    else
        hyp_pred(s) = NaN;
    end
    
    % calculate p values
%     for h = 1:5
%         for y = 1:5
%             hypothesis_p_val(h,y) = ranksum(boxplot_data(:,h),boxplot_data(:,y));
%         end
%     end

end

pred_acc = mean((hyp_pred==true_targets),'omitnan')
