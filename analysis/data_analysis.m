% john analysis

%% set up workspace

clear all

band_names = {'broadband','alpha-theta','beta','low-gamma','high-gamma'};

base_path = '/Users/jbernabei/Documents/PhD_Research/ecog_seeg/ecog_vs_seeg';
iEEG_atlas_path = '/Users/jbernabei/Documents/PhD_Research/atlas_project/iEEG_atlas';

[all_patients, all_inds, all_locs, conn_field, var_field, coords_field, hasData_field, hasVar_field, id_field,...
    implant_field, outcome_field, resect_field, roi_field, target_field,...
    therapy_field, region_list, region_name, lesion_field,...
    sz_field] = set_up_workspace(iEEG_atlas_path);

load color_bar
load color_bar_alt

% identify 3 categories:
% 1) durable good outcome
% 2) relapse
% 3) poor outcome

%% remove WM from all structures
data_patient_indices = find([all_patients.hasData] & hasVar_field & ~strcmp({all_patients.patientID},'HUP111')& ~strcmp({all_patients.patientID},'HUP117'));
for s = 1:length(data_patient_indices)
    
    patient_data = all_patients(data_patient_indices(s));
    
    num_elecs = length(patient_data.roi);
    resect_bool = zeros(num_elecs,1);
    resect_bool(patient_data.resect) = 1;
    
    wm_elecs = find(patient_data.roi==9171); % get which electrodes are in WM
    resect_bool(wm_elecs) = [];
    
    final_res_elecs = find(resect_bool);
    
    patient_data.resect = final_res_elecs;
    
    for f = 1:5
    base_adj = patient_data.conn(f).data;
    final_adj = base_adj;
    final_adj(:,wm_elecs) = [];
    final_adj(wm_elecs,:) = [];
    
    patient_data.conn(f).data = final_adj;
    
    base_var = patient_data.var(f).data;
    final_var = base_var;
    final_var(:,wm_elecs) = [];
    final_var(wm_elecs,:) = [];
    
    patient_data.var(f).data = final_var;
    end
    
    final_roi = patient_data.roi;
    final_roi(wm_elecs) = [];
    
    patient_data.roi = final_roi;
    
    final_coords = patient_data.coords;
    final_coords(wm_elecs,:) = [];
    patient_data.coords = final_coords;
    
    all_patients(data_patient_indices(s)) = patient_data;
    
    size(all_patients(data_patient_indices(s)).conn(f).data)
    size(all_patients(data_patient_indices(s)).coords)
    size(all_patients(data_patient_indices(s)).roi)
    
end

%save('output/node_str_results/node_str_D_rs.mat','D_rs_all','pt_name','pt_outcome')
%% get basic atlas info
test_band = 1;

test_threshold = 5;
data_patient_indices = find([all_patients.hasData]);

atlas_patients = all_patients(data_patient_indices);
[mean_conn, std_conn, num_samples, sem_conn] = create_atlas({atlas_patients.conn}, {atlas_patients.roi}, {atlas_patients.resect}, region_list, test_band, test_threshold);
%[mean_var, std_var, num_samples, sem_var] = create_atlas({cv_patients.var}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

num_samples(num_samples==0) = NaN;
mean_connectivity(num_samples==0) = NaN;

figure(1);clf;
subplot(1,3,1)
imagesc(num_samples)
title('number of samples')
colorbar
subplot(1,3,2)
imagesc(mean_conn)
title('mean connectivity')
colormap(color_bar)
colorbar
subplot(1,3,3)
imagesc(std_conn)
title('std connectivity')
colormap(color_bar)
colorbar

nanmedian(std_conn(:))
figure(2);clf;
hold on
plot(mean_conn(:),std_conn(:),'ko')
xlabel('Median connectivity across patients')
ylabel('Standard deviation of connectivity across patients')
title('Atlas edge variability across patients')

inter_reg_std_mdl = fitlm(mean_conn(:),std_conn(:))

figure(3);clf;
subplot(1,2,1)
title('edges > 1 sample')
imagesc(num_samples>1)
subplot(1,2,2)
title('edges > 5 samples')
imagesc(num_samples>5)

%% do distance regression -> check each frequency band and each roi separately 
good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good') & ~strcmp({all_patients.patientID},'HUP111')& ~strcmp({all_patients.patientID},'HUP117'));

% do distance regression in good outcome patients -> non resected (at first)
% need to make distance matrix
for s = 1:length(good_patient_indices)
    dist_mat = [];
    this_pt = good_patient_indices(s);
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
    for f = 1:5
        intra_roi_conn(r).freq(f).data = [];    
    end
    
    for s = 1:length(good_patient_indices)
        this_pt = good_patient_indices(s);
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
            for f = 1:5
                triu_matrix = triu(all_patients(this_pt).conn(f).data);

                this_roi_conn = triu_matrix(left_region,left_region);

                roi_conn_vec = this_roi_conn(:);

                intra_roi_conn(r).freq(f).data = [intra_roi_conn(r).freq(f).data;roi_conn_vec];
                intra_roi_conn(r).freq(f).data(intra_roi_conn(r).freq(f).data==0) = [];


            end
        end      

        if length(right_region)>1
            fprintf('more than 1 node on right\n')

            dist_matrix = triu(all_patients(this_pt).dist_mat);
            this_roi_dist = dist_matrix(right_region,right_region);
            dist_vec = this_roi_dist(:);
            intra_roi_conn(r).dist = [intra_roi_conn(r).dist;dist_vec];
            intra_roi_conn(r).dist(intra_roi_conn(r).dist==0) = [];

            for f = 1:5
                triu_matrix = triu(all_patients(this_pt).conn(f).data);

                this_roi_conn = triu_matrix(right_region,right_region);

                roi_conn_vec = this_roi_conn(:);

                intra_roi_conn(r).freq(f).data = [intra_roi_conn(r).freq(f).data;roi_conn_vec];
                intra_roi_conn(r).freq(f).data(intra_roi_conn(r).freq(f).data==0) = [];


            end
        end
        
    end
end
%% intra - roi regression
figure(1);clf;
for f = 1:5
    all_conn = [];
    all_dist = [];
for r = [1:45]
    if length(intra_roi_conn(r).freq(f).data)~=length(intra_roi_conn(r).dist)
    else
    all_conn = [all_conn;intra_roi_conn(r).freq(f).data];
    all_dist = [all_dist;intra_roi_conn(r).dist];
    end
end
    mdl=fit(all_dist,all_conn,'rat11');
    
    xvals = linspace(min(all_dist),max(all_dist),100);
    ypred = (mdl.p1.*xvals+mdl.p2)./(xvals+mdl.q1);
    figure(1)
    subplot(1,5,f)
    hold on
    plot(all_dist,all_conn,'ko')
    plot(xvals,ypred,'r-','LineWidth',2)
    hold off
    
    intra_roi_curve(f).mdl = mdl;
    
    if f==2
        figure(2);clf;
        subplot(1,2,1)
        ypred = (mdl.p1.*all_dist+mdl.p2)./(all_dist+mdl.q1);
        residuals = all_conn-ypred;
        plot(all_dist,residuals,'ko')
        
        for j = 1:8
            upper_bound = 10*j;
            lower_bound = 10*(j-1)+1;
            which_inds = find([all_dist>lower_bound].*[all_dist<upper_bound]);
            std_val(j) = std(abs(residuals(which_inds))) 
            
        end
        std_val(7) = [];
        subplot(1,2,2)
        plot(std_val,'ko')
        std_mdl = fitlm([5:10:55,75],std_val)
    end
end

%% inter - roi regression

for f = 1
    all_inter_roi_conn(f).data = [];
end
all_inter_roi_dist = [];

for s = 1:length(good_patient_indices)
    
    s
    
    this_pt = good_patient_indices(s);
    this_pt_roi = all_patients(this_pt).roi;
    
    

    for e1 = 1:length(this_pt_roi)
        roi1 = this_pt_roi(e1);
        for e2 = 1:length(this_pt_roi)
            roi2 = this_pt_roi(e2);
            if roi1==roi2
            else
                inter_elec_dist = all_patients(this_pt).dist_mat(e1,e2);
                all_inter_roi_dist = [all_inter_roi_dist;inter_elec_dist];
                
                for f = 1
                    freq_conn = all_patients(this_pt).conn(f).data(e1,e2);
                    
                    all_inter_roi_conn(f).data = [all_inter_roi_conn(f).data;freq_conn];
                end
            end
        end
    end
    
end

%% find z scores

test_threshold = 5; % main results for = 7

good_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'good') & ~strcmp({all_patients.patientID},'HUP111')& ~strcmp({all_patients.patientID},'HUP117'));
poor_patient_indices = find([all_patients.hasData] & strcmp({all_patients.outcome},'poor'));

for atlas_method = ["patient","not edge"]
    for test_band = 3
        % cross-validation of good-outcome patients
        for s = 1:length(good_patient_indices)
            test_patient = all_patients(good_patient_indices(s));
            cv_patients = all_patients(good_patient_indices);
            cv_patients(s) = [];           
            
            patient_id = test_patient.patientID;

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [all_scores, corr_val, virtual_resect] =...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,"std",atlas_method,'conn');
                
            all_corr_vals(s,test_band) = corr_val;
            
            fprintf(repmat('\b',1,line_length))

        end
        
        cv_patients = all_patients(good_patient_indices);

        % cross-validation of poor-outcome patients
        for s = 1:length(poor_patient_indices)
            test_patient = all_patients(poor_patient_indices(s));

            line_length = fprintf('Testing %s...', test_patient.patientID);

            [all_scores, corr_val, virtual_resect] = ...
                get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,"std",atlas_method,'conn');

            all_corr_vals(s+length(good_patient_indices),test_band) = corr_val;

            fprintf(repmat('\b',1,line_length))

        end
        
    end
end


%% do node abnormality study
% lets just pick broadband for now
test_band = 3;

color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];

edge_threshold = 2.5;
node_threshold = 0.1;

% compute node abnormality of entire network (pre-surgical)
[good_abn, good_frac_abn] = compute_node_abnormality(good_patient_zscores(test_band).freq, edge_threshold, node_threshold);
[poor_abn, poor_frac_abn] = compute_node_abnormality(poor_patient_zscores(test_band).freq, edge_threshold, node_threshold);

for s = 1:length(good_abn)
    
    D_rs = compute_D_rs(good_abn{s}',all_patients(good_patient_indices(s)).resect);
    
    D_rs_abn(s,1) = D_rs;
end

for s = 1:length(poor_abn)
    
    D_rs = compute_D_rs(poor_abn{s}',all_patients(poor_patient_indices(s)).resect);
    
    D_rs_abn((s+length(poor_abn)),1) = D_rs;
end

mean(D_rs_abn(1:length(good_abn)))
nanmean(D_rs_abn((length(good_abn)+1):end))
mean(D_rs_all(pt_outcome==1))
mean(D_rs_all(pt_outcome==0))

ranksum(D_rs_abn,D_rs_all)
%% get all node abnormality for surgical network and then do D_RS for node
% abnormality between surgical and spared network
for s = 1:length(data_patient_indices)
    % get resected electrodes
    res_elecs = all_patients(s).resect;
    
    % get node strengths
    pt_node_str = [];
    for f = 1:5
        pt_node_str = nanmean(all_patients(s).z_scores(1).data.all)';

        res_result = pt_node_str(res_elecs);
        non_res_result = pt_node_str;
        non_res_result(res_elecs) = [];

        mwu_result = mwwtest(res_result',non_res_result');

        D_rs_val = max(mwu_result.U(2))./(length(res_result).*length(non_res_result));
        D_rs_all_abn(s,f) = D_rs_val;
    end
    
    pt_name_abn{s,1} = all_patients(data_patient_indices(s)).patientID;
    pt_outcome_abn{s,1} = all_patients(data_patient_indices(s)).outcome;
end
%% regress patient length of epilepsy with 
% a. all global abnormalities
% b. in-in abnormalities
% c. in-out abnormalities
% d. out-out abnormalities

good_cond = [hasData_field{:}] & (strcmp(outcome_field,'good'));
good_outcome_patients = all_patients(good_cond);

test_band = 3;

all_cond = find([hasData_field{:}]&(strcmp(outcome_field,'good')|strcmp(outcome_field,'poor')));
all_data_patients = all_patients(all_cond);
all_age_onset = good_outcome_patients(1).age_onset(all_cond);
all_age_surg = good_outcome_patients(1).age_surgery(all_cond);

all_epilepsy_duration = all_age_surg-all_age_onset;

for s = 1:length(all_data_patients)
    %abs_mean_abn_all(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.all)));
    abs_mean_abn_in_in(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.in_in)));
    abs_mean_abn_in_out(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.in_out)));
    abs_mean_abn_out_out(s) = nanmean(nanmean((all_patients(all_cond(s)).z_scores(test_band).data.out_out)));
end

nan_inds = find(isnan(all_epilepsy_duration));
all_epilepsy_duration(nan_inds) = [];
all_age_surg(nan_inds) = [];

%abs_mean_abn_all(nan_inds) = [];
abs_mean_abn_in_in(nan_inds) = [];
abs_mean_abn_in_out(nan_inds) = [];
abs_mean_abn_out_out(nan_inds) = [];

features = [all_epilepsy_duration]%,all_age_surg];


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

