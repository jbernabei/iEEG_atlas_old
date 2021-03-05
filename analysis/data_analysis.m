%% set up workspace
clear all

iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';
metadata = readtable("data/atlas_project_metadata.xlsx");

all_engel_scores = metadata{:,8:10}; % extracted 3 columns of engel scores (1.2 = 1B, 2.1 = 2A, etc)

[all_patients_raw, all_inds, all_locs, coords_field, hasData_field, hasVar_field, id_field,...
    implant_field, outcome_field, target_field,...
    therapy_field, region_list, region_name, lesion_field] = set_up_workspace(iEEG_atlas_path);

% Fix HUP140
% Fix HUP172

% extract outcome
for s = 1:length(all_patients_raw)
    engel_scores = metadata{s,8:10};
    early_outcome(s) = floor(engel_scores(1));
    late_outcome(s) = floor(engel_scores(end));
end

% need to modify so we remove any un-used patients.
remove_patients = ~[all_patients_raw.hasData];

age_onset = all_patients_raw(1).age_onset;
age_surgery = all_patients_raw(1).age_surgery;

all_patients_raw(remove_patients) = [];
early_outcome(remove_patients) = [];
late_outcome(remove_patients) = [];
age_onset(remove_patients) = [];
age_surgery(remove_patients) = [];
therapy_field(remove_patients) = [];
lesion_field(remove_patients) = [];
target_field(remove_patients) = [];

all_engel_scores(remove_patients,:) = [];

metadata(remove_patients,:) = [];

load color_bar
load color_bar_alt

all_patients = remove_all_wm(all_patients_raw);

%% Do all patient localization renderings
for s = [15,17,18]
    this_pt = all_patients(s).patientID
    all_coords = all_patients(s).coords;
    % set up .node file
    final_elec_matrix = [all_coords, ones(size(all_coords,1),1), ones(size(all_coords,1),1)];  
    dlmwrite('output/all_render.node',final_elec_matrix,'delimiter',' ','precision',5)
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/all_render.node','output/all_render.mat',sprintf('output/all_render/%s_all_node.jpg',this_pt))

end

%% Do rendering for all coordinates together
all_coords = [];
all_pt_ind = [];
for s = 2:2:length(all_patients)
    this_coords = all_patients(s).coords;
    all_coords = [all_coords;this_coords];
    all_pt_ind = [all_pt_ind;s*ones(size(this_coords,1),1)];
end

% all_coords = all_coords([1:3:end],:);
% all_pt_ind = all_pt_ind([1:3:end])

final_elec_matrix = [all_coords, all_pt_ind, ones(size(all_coords,1),1)];  
dlmwrite('output/all_render.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/all_render.node','output/all_pt_render.mat','output/all_render/all_pt_node.jpg')

%% get outcome data & modify based on current measures

figure(1);clf;
hold on
title('patient outcome scores')
ylabel('patient number')
imagesc(floor(all_engel_scores),'AlphaData',~isnan(all_engel_scores))
xticks([1:3])
xticklabels({'6 month','12 month','24 month'})
colorbar % will need to change color axis stuff here but otherwise its fine.
hold off

%% generate clinical metadata table for main text
% define outcome by 6 month
good_patient_indices = find(early_outcome==1);
poor_patient_indices = find(early_outcome>1);

% extract good and poor outcome metadata
good_metadata = metadata(good_patient_indices,:);
poor_metadata = metadata(poor_patient_indices,:);

% get patient numbers
patient_nums = [length(good_patient_indices), length(poor_patient_indices)];

% number of male/female
num_female = [sum(strcmp(good_metadata{:,18},'F')),sum(strcmp(poor_metadata{:,18},'F'))];
num_male = [sum(strcmp(good_metadata{:,18},'M')),sum(strcmp(poor_metadata{:,18},'M'))];

p_gender = chi2Tests([num_female;num_male]);

% number of right/left
num_right = [sum(strcmp(good_metadata{:,14},'R')),sum(strcmp(poor_metadata{:,14},'R'))];
num_left = [sum(strcmp(good_metadata{:,14},'L')),sum(strcmp(poor_metadata{:,14},'L'))];

p_laterality = chi2Tests([num_right;num_left]);

% number of lesional / non-lesional
num_lesional = [sum(strcmp(good_metadata{:,15},'Lesional')),sum(strcmp(poor_metadata{:,15},'Lesional'))];
num_NL = [sum(strcmp(good_metadata{:,15},'Non-Lesional')),sum(strcmp(poor_metadata{:,15},'Non-Lesional'))];

p_lesional = chi2Tests([num_lesional;num_NL]);

% age surgery & age onset
age_onset_mean_std = [mean(good_metadata{:,16}),std(good_metadata{:,16}),mean(poor_metadata{:,16}),std(poor_metadata{:,16})];
age_surgery_mean_std = [mean(good_metadata{:,17}),std(good_metadata{:,17}),mean(poor_metadata{:,17}),std(poor_metadata{:,17})];

p_age_onset = ranksum(good_metadata{:,16},poor_metadata{:,16});
p_age_surgery = ranksum(good_metadata{:,17},poor_metadata{:,17});

% target (temporal, frontal/FP, insular)
num_temporal = [(sum(strcmp(good_metadata{:,13},'Temporal'))+sum(strcmp(good_metadata{:,13},'MTL'))),(sum(strcmp(poor_metadata{:,13},'Temporal'))+sum(strcmp(poor_metadata{:,13},'MTL')))];
num_frontal = [(sum(strcmp(good_metadata{:,13},'Frontal'))+sum(strcmp(good_metadata{:,13},'MFL'))+sum(strcmp(good_metadata{:,13},'FP'))),...
               (sum(strcmp(poor_metadata{:,13},'Frontal'))+sum(strcmp(poor_metadata{:,13},'MFL'))+sum(strcmp(poor_metadata{:,13},'FP')))];
num_insular = [sum(strcmp(good_metadata{:,13},'Insular')),sum(strcmp(poor_metadata{:,13},'Insular'))];

p_target = chi2Tests([num_temporal;num_frontal;num_insular]);

% type of implant (ECoG / SEEG)
num_ECoG = [sum(strcmp(good_metadata{:,12},'ECoG')),sum(strcmp(poor_metadata{:,12},'ECoG'))];
num_SEEG = [sum(strcmp(good_metadata{:,12},'SEEG')),sum(strcmp(poor_metadata{:,12},'SEEG'))];

p_implant = chi2Tests([num_ECoG;num_SEEG]);

% type of surgery (resection/ablation)
num_resection = [sum(strcmp(good_metadata{:,11},'Resection')),sum(strcmp(poor_metadata{:,11},'Resection'))];
num_ablation = [sum(strcmp(good_metadata{:,11},'Ablation')),sum(strcmp(poor_metadata{:,11},'Ablation'))];

p_surgery = chi2Tests([num_resection;num_ablation]);

% number of total non-wm nodes & nodes targeted by surgery
for s = 1:length(all_patients)
    num_non_wm_all(s) = length(all_patients(s).coords);
    num_resect_all(s) = length(all_patients(s).resect);
end

num_nodes_mean_std = [mean(num_non_wm_all(good_patient_indices)),std(num_non_wm_all(good_patient_indices)),mean(num_non_wm_all(poor_patient_indices)),std(num_non_wm_all(poor_patient_indices))];
num_resect_mean_std = [mean(num_resect_all(good_patient_indices)),std(num_resect_all(good_patient_indices)),mean(num_resect_all(poor_patient_indices)),std(num_resect_all(poor_patient_indices))];

p_nodes = ranksum(num_non_wm_all(good_patient_indices),num_non_wm_all(poor_patient_indices)); % rank-sum test
p_resect = ranksum(num_resect_all(good_patient_indices),num_resect_all(poor_patient_indices)); % rank-sum test

% number that relapsed to engel 2+ by year 2
num_relapse = [sum([early_outcome==1].*[late_outcome>1]),NaN];

%% get basic atlas info
test_band = 1;

for test_threshold = 1:10

data_patient_indices = find([all_patients.hasData]);

atlas_patients = all_patients(data_patient_indices);
[mean_conn, std_conn, num_samples, sem_conn] = create_atlas({atlas_patients.conn}, {atlas_patients.roi}, {atlas_patients.resect}, region_list, test_band, test_threshold);
[mean_var, std_var, num_samples, sem_var] = create_atlas({atlas_patients.var}, {atlas_patients.roi}, {atlas_patients.resect}, region_list, test_band, test_threshold);

num_samples(num_samples==0) = NaN;
mean_connectivity(num_samples==0) = NaN;

% figure(1);clf;
% subplot(1,3,1)
% imagesc(num_samples)
% title('number of samples')
% colorbar
% subplot(1,3,2)
% imagesc(mean_conn,'AlphaData',~isnan(mean_conn))
% title('mean connectivity')
% colormap(color_bar)
% colorbar
% subplot(1,3,3)
% imagesc(std_conn,'AlphaData',~isnan(std_conn))
% title('std connectivity')
% colormap(color_bar)
% colorbar

std_conn_all(test_threshold) =  nanmedian(std_conn(:))


all_mean_conn = mean_conn(:);
all_std_conn = std_conn(:);
all_pt_count = num_samples(:);
figure(test_threshold);clf;
hold on
for i = 1:5:length(all_mean_conn)
    try plot(all_mean_conn(i),all_std_conn(i),'ko','MarkerSize',(all_pt_count(i)+1))
    catch ME
    end
end
xlabel('Median connectivity across patients')
ylabel('Standard deviation of connectivity across patients')
title('Atlas edge variability across patients')
end
% color by part of brain
% size by how many patients have that connection

figure(3);clf;
subplot(1,2,1)
title('edges > 1 sample')
imagesc(num_samples>1)
subplot(1,2,2)
title('edges > 3 samples')
imagesc(num_samples>3)

%% do distance regression -> check each frequency band and each roi separately 

% do distance regression in good outcome patients -> non resected (at first)
% need to make distance matrix
for s = 1:length(all_patients)
    dist_mat = [];
    this_pt = s;
    coord_mat = all_patients(this_pt).coords;
    num_coords = size(coord_mat,1);
    for e1 = 1:num_coords
        for e2 = 1:num_coords
            dist_mat(e1,e2) = sqrt(sum((coord_mat(e1,:) - coord_mat(e2,:)).^2));
        end
    end
     all_patients(this_pt).dist_mat = dist_mat;
end


%% do same regression in left and right
for r = 1:45
    
    % loop through frequency bands and make the base structure 
    intra_roi_conn(r).dist = [];
    for f = 1
        intra_roi_conn(r).freq(f).data = [];    
        intra_roi_var(r).freq(f).data = [];  
    end
    
    for s = 1:length(all_patients)
        this_pt = s;
        this_pt_roi = all_patients(this_pt).roi;
        
        % get right and left region numerical codes
        left_region = find(this_pt_roi==all_inds(2*r));
        right_region = find(this_pt_roi==all_inds(2*r-1));
        
        if length(left_region)>1

            dist_matrix = triu(all_patients(this_pt).dist_mat);
            this_roi_dist = dist_matrix(left_region,left_region);
            dist_vec = this_roi_dist(:);
            intra_roi_conn(r).dist = [intra_roi_conn(r).dist;dist_vec];
            intra_roi_conn(r).dist(intra_roi_conn(r).dist==0) = [];

            fprintf('more than 1 node on left\n')
            for f = 1
                triu_matrix = triu(all_patients(this_pt).conn(f).data);
                triu_matrix_var = triu(all_patients(this_pt).var(f).data);

                this_roi_var = triu_matrix_var(left_region,left_region);
                this_roi_conn = triu_matrix(left_region,left_region);

                roi_conn_vec = this_roi_conn(:);
                roi_var_vec = this_roi_var(:);

                intra_roi_conn(r).freq(f).data = [intra_roi_conn(r).freq(f).data;roi_conn_vec];
                intra_roi_conn(r).freq(f).data(intra_roi_conn(r).freq(f).data==0) = [];
                
                intra_roi_var(r).freq(f).data = [intra_roi_var(r).freq(f).data;roi_var_vec];
                intra_roi_var(r).freq(f).data(intra_roi_var(r).freq(f).data==0) = [];


            end
        end      

        if length(right_region)>1
            fprintf('more than 1 node on right\n')

            dist_matrix = triu(all_patients(this_pt).dist_mat);
            this_roi_dist = dist_matrix(right_region,right_region);
            dist_vec = this_roi_dist(:);
            intra_roi_conn(r).dist = [intra_roi_conn(r).dist;dist_vec];
            intra_roi_conn(r).dist(intra_roi_conn(r).dist==0) = [];

            for f = 1
                triu_matrix = triu(all_patients(this_pt).conn(f).data);
                triu_matrix_var = triu(all_patients(this_pt).var(f).data);

                this_roi_var = triu_matrix_var(right_region,right_region);
                this_roi_conn = triu_matrix(right_region,right_region);

                roi_conn_vec = this_roi_conn(:);
                roi_var_vec = this_roi_var(:);

                intra_roi_conn(r).freq(f).data = [intra_roi_conn(r).freq(f).data;roi_conn_vec];
                intra_roi_conn(r).freq(f).data(intra_roi_conn(r).freq(f).data==0) = [];
                
                intra_roi_var(r).freq(f).data = [intra_roi_var(r).freq(f).data;roi_var_vec];
                intra_roi_var(r).freq(f).data(intra_roi_var(r).freq(f).data==0) = [];
            end
        end
        
    end
end
%%
figure(1);clf;
for f = 1
    all_conn = [];
    all_dist = [];
    all_var = [];
for r = [1:45]
    if length(intra_roi_conn(r).freq(f).data)~=length(intra_roi_conn(r).dist)
    else
    all_conn = [all_conn;intra_roi_conn(r).freq(f).data];
    all_var = [all_var;intra_roi_var(r).freq(f).data];
    all_dist = [all_dist;intra_roi_conn(r).dist];
    end
end
    mdl=fit(all_dist,all_conn,'rat11');
    
    xvals = linspace(min(all_dist+1),max(all_dist),100);
    ypred = (mdl.p1.*xvals+mdl.p2)./(xvals+mdl.q1);
    figure(1)
    subplot(1,3,1)
    hold on
    plot(all_dist,all_conn,'ko')
    plot(xvals,ypred,'r-','LineWidth',2)
    xlabel('Internodal distance (mm)')
    ylabel('Broadband cross-correlation')
    legend('Atlas edge weights','non-linear polynomial regression')
    title('Intra-regional internodal distance regression')
    hold off
    
    intra_roi_curve(f).mdl = mdl;

    subplot(1,3,2)
    ypred = (mdl.p1.*all_dist+mdl.p2)./(all_dist+mdl.q1);
    residuals = all_conn-ypred;
    plot(all_dist,residuals,'ko')
    xlabel('Internodal distance (mm)')
    ylabel('Broadband cross-correlation')
    title('Residual to best-fit polynomial')

    num_bins = 20;
    for j = 1:num_bins
        upper_bound = (80/num_bins)*j;
        lower_bound = (80/num_bins)*(j-1)+1;
        which_inds = find([all_dist>lower_bound].*[all_dist<upper_bound]);
        std_val(j) = std(abs(residuals(which_inds))) 

    end
    subplot(1,3,3)
    plot([1:4:80],std_val,'ko')
    hold on
    xlabel('Internodal distance (mm)')
    ylabel('Correlation variance')
    title('Variance of Broadband cross-correlation residual')
    std_mdl = fitlm([1:4:80],std_val)
    plot([1,80],[0.127-[1,80].*0.001355],'r-','LineWidth',2)
    
end

%% find z scores
clear results_struct
for s = 1:length(all_patients)
    for f = 1
    results_struct(s).freq(f).all_scores = [];
    results_struct(s).freq(f).virtual_resect = [];
    results_struct(s).freq(f).corr_val = [];
    
    var_results_struct(s).freq(f).all_scores = [];
    var_results_struct(s).freq(f).virtual_resect = [];
    var_results_struct(s).freq(f).corr_val = [];
    end
end

conn_type = 1; % mean/median
test_band = 1;

good_patient_indices = find(early_outcome==1);
poor_patient_indices = find(early_outcome>1);

for atlas_method = ["native"]
    for test_threshold = 3
        % cross-validation of good-outcome patients
        for s = 1:length(good_patient_indices)
            test_patient = all_patients(good_patient_indices(s));
            cv_patients = all_patients(good_patient_indices);
            cv_patients(s) = [];           
            
            patient_id = test_patient.patientID;

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [all_scores, corr_val, virtual_resect] =...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,atlas_method,conn_type);
                
            results_struct(good_patient_indices(s)).freq(test_band).all_scores = all_scores;
            results_struct(good_patient_indices(s)).freq(test_band).corr_val = corr_val;
            results_struct(good_patient_indices(s)).freq(test_band).virtual_resect = virtual_resect;
            
            corr_val_thresh(good_patient_indices(s),test_threshold) = corr_val;
            
            fprintf(repmat('\b',1,line_length))

        end
        
        cv_patients = all_patients(good_patient_indices);

        % cross-validation of poor-outcome patients
        for s = 1:length(poor_patient_indices)
            test_patient = all_patients(poor_patient_indices(s));

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [all_scores, corr_val, virtual_resect] = ...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,atlas_method,conn_type);

            results_struct(poor_patient_indices(s)).freq(test_band).all_scores = all_scores;
            results_struct(poor_patient_indices(s)).freq(test_band).corr_val = corr_val;
            results_struct(poor_patient_indices(s)).freq(test_band).virtual_resect = virtual_resect;

            fprintf(repmat('\b',1,line_length)  )
            
            corr_val_thresh(poor_patient_indices(s),test_threshold) = corr_val;

        end
        
    end
end

%% Analyze results
% for this paper focus on in-out and out-out connections
% if in-in we need to restrict resection & ablation separately 
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color5 = [103 55 155]/255;
color6 = [78 172 91]/255;

test_band = 1;

% get different patients
stable_good = find([late_outcome==1]);%.*[strcmp(target_field,'MTL')+strcmp(target_field,'Temporal')]);
relapse = find([early_outcome==1].*[late_outcome>1]);%.*[strcmp(target_field,'MTL')+strcmp(target_field,'Temporal')]);
stable_poor = find([early_outcome>1].*[late_outcome>1]);%.*[strcmp(target_field,'MTL')+strcmp(target_field,'Temporal')]);
all_poor = [relapse, stable_poor];

% divide by type of surgery
stable_good_res = find([late_outcome==1].*strcmp(therapy_field,'Resection'));
stable_good_abl = find([late_outcome==1].*strcmp(therapy_field,'Ablation'));
relapse_res = find([early_outcome==1].*[late_outcome>1].*strcmp(therapy_field,'Resection'));
relapse_abl = find([early_outcome==1].*[late_outcome>1].*strcmp(therapy_field,'Ablation'));
stable_poor_res = find([early_outcome>1].*[late_outcome>1].*strcmp(therapy_field,'Resection'));
stable_poor_abl = find([early_outcome>1].*[late_outcome>1].*strcmp(therapy_field,'Ablation'));
all_poor_res = [relapse_res, stable_poor_res];
all_poor_abl = [relapse_abl, stable_poor_abl];

% loop through all patients
for s = 1:length(all_patients)
    % extract virtual resection scores & resected electrodes
    pt_scores = results_struct(s).freq(test_band).all_scores;
    num_elecs = size(results_struct(s).freq(test_band).all_scores,1);
    resect_elecs = all_patients(s).resect;
    %resect_elecs = find(sum(all_patients(s).roi(all_patients(s).resect)'==region_list)>0);
    non_res_elecs = setdiff([1:num_elecs],resect_elecs);
    
    % out-out scores -> get rid of resected rows and columns
    out_out_scores = pt_scores;
    out_out_scores(resect_elecs,:) = NaN;
    out_out_scores(:,resect_elecs) = NaN;

    % in-in scores -> get rid of non-resected rows and columns
    in_in_scores = pt_scores;
    in_in_scores(non_res_elecs,:) = NaN;
    in_in_scores(:,non_res_elecs) = NaN;

    % in-out scores -> need to separately get rid of non-resected rows and
    % columns then add, then get rid of in-in resected region
    in_out_scores = pt_scores;
    in_out_scores(resect_elecs,resect_elecs) = NaN;
    in_out_scores(non_res_elecs,non_res_elecs) = NaN;
    
    in_in_all(s) = nanmedian(in_in_scores(:));
    in_out_all(s) = nanmedian(in_out_scores(:));
    out_out_all(s) = nanmedian(out_out_scores(:));
end

p1 = ranksum(out_out_all(stable_good),out_out_all(stable_poor));
p2 = ranksum(out_out_all(stable_good),out_out_all(all_poor));
p3 = ranksum(out_out_all(stable_good),out_out_all(relapse));
p4 = ranksum(out_out_all(stable_poor),out_out_all(relapse));

[p1 p2 p3 p4]

p5 = ranksum(in_out_all(stable_good),in_out_all(stable_poor));
p6 = ranksum(in_out_all(stable_good),in_out_all(all_poor));
p7 = ranksum(in_out_all(stable_good),in_out_all(relapse));
p8 = ranksum(in_out_all(stable_poor),in_out_all(relapse));

[p5 p6 p7 p8]

p9 = ranksum(rmmissing(in_in_all(stable_good)),rmmissing(in_in_all(stable_poor)));
p10 = ranksum(rmmissing(in_in_all(stable_good)),rmmissing(in_in_all(all_poor)));
p11 = ranksum(rmmissing(in_in_all(stable_good)),rmmissing(in_in_all(relapse)));
p12 = ranksum(rmmissing(in_in_all(stable_poor)),rmmissing(in_in_all(relapse)));

[p9 p10 p11 p12]

p1 = signrank(in_in_all(stable_good),in_out_all(stable_good));
p2 = signrank(in_out_all(stable_good),out_out_all(stable_good));
p3 = signrank(in_in_all(stable_good),out_out_all(stable_good));

[p1 p2 p3]

p4 = signrank(in_in_all(stable_poor),in_out_all(stable_poor));
p5 = signrank(in_out_all(stable_poor),out_out_all(stable_poor));
p6 = signrank(in_in_all(stable_poor),out_out_all(stable_poor));

[p4 p5 p6]

p7 = signrank(in_in_all(relapse),in_out_all(relapse));
p8 = signrank(in_out_all(relapse),out_out_all(relapse));
p9 = signrank(in_in_all(relapse),out_out_all(relapse));

[p7 p8 p9]

p10 = signrank(in_in_all(all_poor),in_out_all(all_poor));
p11 = signrank(in_out_all(all_poor),out_out_all(all_poor));
p12 = signrank(in_in_all(all_poor),out_out_all(all_poor));

[p10 p11 p12]

plotting_structure = NaN*zeros(length(all_patients),6);
plotting_structure(1:length(in_in_all(stable_good)),1) = in_in_all(stable_good);
plotting_structure(1:length(in_out_all(stable_good)),3) = in_out_all(stable_good);
plotting_structure(1:length(out_out_all(stable_good)),5) = out_out_all(stable_good);
plotting_structure(1:length(in_in_all(all_poor)),2) = in_in_all(all_poor);
plotting_structure(1:length(in_out_all(all_poor)),4) = in_out_all(all_poor);
plotting_structure(1:length(out_out_all(all_poor)),6) = out_out_all(all_poor);

figure(1);clf;
which_color = color2;
hold on
for r = 1:6
    if mod(r,2)==1
        which_color = color1;
    else
        which_color = color2;
    end
    scatter(r.*ones(size(plotting_structure,1),1),plotting_structure(:,r),'MarkerEdgeColor',which_color,'MarkerFaceColor',which_color,'jitter','on')
    plot([r-0.25 r+0.25],[nanmedian(plotting_structure(:,r)) nanmedian(plotting_structure(:,r))],'k-','LineWidth',2)
end

plot([1,5], [11,11], '-k', 'LineWidth',1)
plot(3.5, 11.5, 'kd')
plot([3,5], [8,8], '-k', 'LineWidth',1)
plot(4, 8.5, 'kd')

hold off

% plotting_structure_2 = NaN*zeros(length(all_patients),6);
% plotting_structure_2(1:length(in_in_all(stable_good_res)),1) = in_in_all(stable_good_res);
% plotting_structure_2(1:length(in_out_all(stable_good_res)),3) = in_out_all(stable_good_res);
% plotting_structure_2(1:length(out_out_all(stable_good_res)),5) = out_out_all(stable_good_res);
% plotting_structure_2(1:length(in_in_all(all_poor_res)),2) = in_in_all(all_poor_res);
% plotting_structure_2(1:length(in_out_all(all_poor_res)),4) = in_out_all(all_poor_res);
% plotting_structure_2(1:length(out_out_all(all_poor_res)),6) = out_out_all(all_poor_res);
% 
% figure(1);clf;
% which_color = color2;
% hold on
% for r = 1:6
%     if mod(r,3)==1
%         which_color = color1;
%     elseif mod(r,3)==2
%         which_color = color5;
%     else
%         which_color = color2;
%     end
%     scatter(r.*ones(size(plotting_structure_2,1),1),plotting_structure_2(:,r),'MarkerEdgeColor',which_color,'MarkerFaceColor',which_color,'jitter','on')
%     plot([r-0.25 r+0.25],[nanmedian(plotting_structure_2(:,r)) nanmedian(plotting_structure_2(:,r))],'k-','LineWidth',2)
% end
% hold off
%%
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];

p13 = ranksum(in_in_all(strcmp(lesion_field,'Lesional')),in_in_all(strcmp(lesion_field,'Non-Lesional')))
p14 = ranksum(in_out_all(strcmp(lesion_field,'Lesional')),in_out_all(strcmp(lesion_field,'Non-Lesional')))
p15 = ranksum(out_out_all(strcmp(lesion_field,'Lesional')),out_out_all(strcmp(lesion_field,'Non-Lesional')))

figure(2);clf;
hold on
a1 = in_in_all(strcmp(lesion_field,'Lesional'))';a2 = in_in_all(strcmp(lesion_field,'Non-Lesional'))';
a2([14,25]) = [];
%boxplot([a1,a2])
scatter(ones(length(a1),1),a1,'MarkerEdgeColor',color5,'MarkerFaceColor',color5,'jitter','on')
plot([0.75 1.25],[nanmedian(a1) nanmedian(a1)],'k-','LineWidth',2)
scatter(2*ones(length(a2),1),a2,'MarkerEdgeColor',color6,'MarkerFaceColor',color6,'jitter','on')
plot([1.75 2.25],[nanmedian(a2) nanmedian(a2)],'k-','LineWidth',2)
xticks([1:2])
ylabel('Mean abnormality within resection zone')
xticklabels({'Lesional','Non-Lesional'})
plot([1,2], [11,11], '-k', 'LineWidth',1)
plot(1.5, 11.5, '*k')

%% regress patient length of epilepsy with 
% b. in-in abnormalities
% c. in-out abnormalities
% d. out-out abnormalities

test_band = 1;

all_epilepsy_duration = age_surgery-age_onset;

features = [all_epilepsy_duration];


f2 = fitlm(features, in_in_all','RobustOpts','on')
p2 = f2.Coefficients{2,4};
b2 = f2.Coefficients{1,1}; m2a = f2.Coefficients{2,1}; %m2b = f2.Coefficients{3,1};

f3 = fitlm(features, in_out_all','RobustOpts','on')
p3 = f3.Coefficients{2,4};
b3 = f3.Coefficients{1,1}; m3a = f3.Coefficients{2,1}; %m3b = f3.Coefficients{3,1};

f4 = fitlm(features, out_out_all','RobustOpts','on')
p4 = f4.Coefficients{2,4};
b4 = f4.Coefficients{1,1}; m4a = f4.Coefficients{2,1}; %m4b = f4.Coefficients{3,1};

x1 = [min(all_epilepsy_duration), max(all_epilepsy_duration)];

%
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
purple = [103 55 155]/255;

figure(1);clf;
subplot(1,3,1)
hold on
plot(all_epilepsy_duration,[in_in_all],'Color',red,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m2a.*x1+b2],'k-','LineWidth',2)
txt2 = sprintf('p = %d',p2)
%t2 = text(40,1.5,txt2)
title(sprintf('Resected edges, p = %2f',p2))
xlabel('Duration between first seizure and implant (years)')
ylabel('Mean atlas Z score')
hold off
subplot(1,3,2)
hold on
plot(all_epilepsy_duration,[in_out_all],'Color',purple,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m3a.*x1+b3],'k-','LineWidth',2)
txt3 = sprintf('p = %d',p3)
%t3 = text(40,1.5,txt3)
title(sprintf('Non-resected edges, p = %2f',p3))
xlabel('Duration between first seizure and surgery (years)')
ylabel('Mean atlas Z score')
hold off
subplot(1,3,3)
hold on
plot(all_epilepsy_duration,[out_out_all],'Color',blue,'Marker','.','LineStyle','none','MarkerSize',24)
plot(x1,[m4a.*x1+b4],'k-','LineWidth',2)
txt4 = sprintf('p = %d',p4)
%t4 = text(40,1.5,txt4)
title(sprintf('Non-resected edges, p = %2f',p4))
xlabel('Duration between first seizure and implant (years)')
ylabel('Mean atlas Z score')
hold off

% figure 2: all 3 on same axes
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
purple = [103 55 155]/255;

%% second example patient analysis
% -> find one where mapping agrees & good outcome

% not seizure free initially post-laser, but now free after ATL
% initial MRI non-lesional -> show post-surgical ATL

% set patient
test_band = 1;
which_pt = 31; % HUP133
secondary_inds = region_list([81:89])
tertiary_inds = region_list([47,55])

% define first region electrodes
resect_elecs = all_patients(which_pt).resect;

% define second region electrodes
elecs_2 = find(sum(all_patients(which_pt).roi==secondary_inds')>0)';
elecs_3 = find(sum(all_patients(which_pt).roi==tertiary_inds')>0)';
elecs_3 = setdiff(elecs_3,resect_elecs);
which_roi = all_patients(which_pt).roi(elecs_2);

% get in-in res and in-in 2
pt_scores = results_struct(which_pt).freq(test_band).all_scores;
abn_matrix = pt_scores;
in_in_res = pt_scores(resect_elecs,resect_elecs);
in_in_2 = pt_scores(elecs_2,elecs_2);
in_in_3 = pt_scores(elecs_3,elecs_3);


figure(1);clf
hold on
scatter(ones(length(in_in_res(:)),1),in_in_res(:),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'jitter','on')
plot([0.75 1.25],[nanmedian(in_in_res(:)),nanmedian(in_in_res(:))],'k-','LineWidth',2)
scatter(2*ones(length(in_in_2(:)),1),in_in_2(:),'MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'jitter','on')
plot([1.75 2.25],[nanmedian(in_in_2(:)),nanmedian(in_in_2(:))],'k-','LineWidth',2)
scatter(3*ones(length(in_in_3(:)),1),in_in_3(:),'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'jitter','on')
plot([2.75 3.25],[nanmedian(in_in_3(:)),nanmedian(in_in_3(:))],'k-','LineWidth',2)
xlim([0.5 3.5])
xticks([1:3])
xticklabels({'target','secondary','tertiary'})
hold off

p1 = ranksum(in_in_res(:),in_in_2(:))
p2 = ranksum(in_in_res(:),in_in_3(:))
p3 = ranksum(in_in_2(:),in_in_3(:))

% get the abnormality matrix
figure(2);clf
hold on
imagesc(abn_matrix,'AlphaData',~isnan(abn_matrix))
caxis([-5 10])
colorbar
hold off

% now do the rendering w/ brain net viewer

% get coordinates
all_coords = all_patients(which_pt).coords;
all_coords(:,3) = all_coords(:,3);

% set up .edge file
abn_matrix(isnan(abn_matrix)) = 0;
save('output/hypothesis_test.edge','abn_matrix','-ascii');

% set up .node file
final_elec_matrix = [all_coords, 4*ones(size(all_coords,1),1), ones(size(all_coords,1),1)];  
final_elec_matrix(resect_elecs,4) = 1;
final_elec_matrix(elecs_2,4) = 2;
final_elec_matrix(elecs_3,4) = 3; 
figure(3);clf;
dlmwrite('output/hypothesis_test.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/hypothesis_test.node','output/hypothesis_test.edge','output/fig_3_option.mat','output/hypothesis_test_HUP133.jpg')

%% third example patient analysis
% -> find one where mapping doesnt work

% not seizure free initially post-laser, but now free after ATL
% initial MRI non-lesional -> show post-surgical ATL

% set patient
test_band = 1;
which_pt = 47; % HUP177
secondary_inds = region_list([58:2:70]); % R parietal (cephalgic aura)
tertiary_inds = [3001 3002]; % insula (oral automatisms)

% define first region electrodes
resect_elecs = all_patients(which_pt).resect;

% define second region electrodes
elecs_2 = find(sum(all_patients(which_pt).roi==secondary_inds')>0)';
elecs_3 = find(sum(all_patients(which_pt).roi==tertiary_inds')>0)';
which_roi = all_patients(which_pt).roi(elecs_2);

% get in-in res and in-in 2
pt_scores = results_struct(which_pt).freq(test_band).all_scores;
abn_matrix = pt_scores;
in_in_res = pt_scores(resect_elecs,resect_elecs);
in_in_2 = pt_scores(elecs_2,elecs_2);
in_in_3 = pt_scores(elecs_3,elecs_3);

all_potential_elecs = [resect_elecs; elecs_2];
pt_scores(all_potential_elecs,:) = NaN;
pt_scores(:,all_potential_elecs) = NaN;
out_out_scores = pt_scores;

figure(1);clf
hold on
scatter(ones(length(in_in_res(:)),1),in_in_res(:),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'jitter','on')
plot([0.75 1.25],[nanmedian(in_in_res(:)),nanmedian(in_in_res(:))],'k-','LineWidth',2)
scatter(2*ones(length(in_in_2(:)),1),in_in_2(:),'MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'jitter','on')
plot([1.75 2.25],[nanmedian(in_in_2(:)),nanmedian(in_in_2(:))],'k-','LineWidth',2)
scatter(3*ones(length(in_in_3(:)),1),in_in_3(:),'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'jitter','on')
plot([2.75 3.25],[nanmedian(in_in_3(:)),nanmedian(in_in_3(:))],'k-','LineWidth',2)
xlim([0.5 3.5])
ylim([-5 10])
xticks([1:3])
xticklabels({'target','secondary','tertiary'})
hold off

p1 = ranksum(in_in_res(:),in_in_2(:))
p2 = ranksum(in_in_res(:),in_in_3(:))

% get the abnormality matrix
figure(2);clf
hold on
imagesc(abn_matrix,'AlphaData',~isnan(abn_matrix))
caxis([-5 10])
colorbar
hold off

% now do the rendering w/ brain net viewer

% get coordinates
all_coords = all_patients(which_pt).coords;
all_coords(:,3) = all_coords(:,3);

% set up .edge file
abn_matrix(isnan(abn_matrix)) = 0;
save('output/hypothesis_test.edge','abn_matrix','-ascii');

% set up .node file
final_elec_matrix = [all_coords, 4*ones(size(all_coords,1),1), ones(size(all_coords,1),1)];  
final_elec_matrix(resect_elecs,4) = 1;
final_elec_matrix(elecs_2,4) = 2;
final_elec_matrix(elecs_3,4) = 3;
dlmwrite('output/hypothesis_test.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/hypothesis_test.node','output/hypothesis_test.edge','output/fig_3_option.mat','output/hypothesis_test_HUP177.jpg')

%% example patient analysis
% -> find one where mapping agrees & good outcome

% set patient
test_band = 1;
which_pt = 39; % HUP148
secondary_inds = region_list([37,39,41]); % L MTL
tertiary_inds = region_list([81:2:90]); % L temp neocortical

% define first region electrodes
resect_elecs = all_patients(which_pt).resect;

% define second region electrodes
elecs_2 = find(sum(all_patients(which_pt).roi==secondary_inds')>0)';
elecs_3 = find(sum(all_patients(which_pt).roi==tertiary_inds')>0)';
which_roi = all_patients(which_pt).roi(elecs_2);

% get in-in res and in-in 2
pt_scores = results_struct(which_pt).freq(test_band).all_scores;
abn_matrix = pt_scores;
in_in_res = pt_scores(resect_elecs,resect_elecs);
in_in_2 = pt_scores(elecs_2,elecs_2);
in_in_3 = pt_scores(elecs_3,elecs_3);

all_potential_elecs = [resect_elecs; elecs_2];
pt_scores(all_potential_elecs,:) = NaN;
pt_scores(:,all_potential_elecs) = NaN;
out_out_scores = pt_scores;

figure(1);clf
hold on
scatter(ones(length(in_in_res(:)),1),in_in_res(:),'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],'jitter','on')
plot([0.75 1.25],[nanmedian(in_in_res(:)),nanmedian(in_in_res(:))],'k-','LineWidth',2)
scatter(2*ones(length(in_in_2(:)),1),in_in_2(:),'MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'jitter','on')
plot([1.75 2.25],[nanmedian(in_in_2(:)),nanmedian(in_in_2(:))],'k-','LineWidth',2)
scatter(3*ones(length(in_in_3(:)),1),in_in_3(:),'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410],'jitter','on')
plot([2.75 3.25],[nanmedian(in_in_3(:)),nanmedian(in_in_3(:))],'k-','LineWidth',2)
xlim([0.5 3.5])
ylim([-5 10])
xticks([1:3])
xticklabels({'target','secondary','tertiary'})
hold off

p1 = ranksum(in_in_res(:),in_in_2(:))
p2 = ranksum(in_in_res(:),in_in_3(:))

% get the abnormality matrix
figure(2);clf
hold on
imagesc(abn_matrix,'AlphaData',~isnan(abn_matrix))
caxis([-5 10])
colorbar
hold off

% now do the rendering w/ brain net viewer

% get coordinates
all_coords = all_patients(which_pt).coords;
all_coords(:,3) = all_coords(:,3);

% set up .edge file
abn_matrix(isnan(abn_matrix)) = 0;
save('output/hypothesis_test.edge','abn_matrix','-ascii');

% set up .node file
final_elec_matrix = [all_coords, 4*ones(size(all_coords,1),1), ones(size(all_coords,1),1)];  
final_elec_matrix(resect_elecs,4) = 1;
final_elec_matrix(elecs_2,4) = 2;
final_elec_matrix(elecs_3,4) = 3; 
dlmwrite('output/hypothesis_test.node',final_elec_matrix,'delimiter',' ','precision',5)
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/hypothesis_test.node','output/hypothesis_test.edge','output/fig_3_option.mat','output/hypothesis_test_HUP144.jpg')




