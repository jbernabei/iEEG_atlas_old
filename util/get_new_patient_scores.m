function [out_out_scores, in_in_scores, in_out_scores, all_scores, corr_val] = get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,score_denominator,atlas_method,mapping_type)

% score_denominator is a string that can be either "std" or "sem"
% atlas_method is a string that can be either "patient" or "edge"

% sets default denominator to "std" if none is given
if ~exist('score_denominator','var'), score_denominator = "std"; end
% sets default method to "patient" if none is given
if ~exist('atlas_method','var'), atlas_method = "patient"; end

% get connectivity atlas of excluded patients
[mean_conn, std_conn, ~, sem_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);
%[mean_var, std_var, ~, sem_var] = create_atlas({cv_patients.var}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

% set up dummy corr value
corr_val = [];

if strcmp(mapping_type,'conn')
if score_denominator=='sem'
    [score_matrix, corr_val] = test_native_adj(test_patient.conn, test_patient.coords, test_patient.roi, mean_conn, nanmean(nanmean(sem_conn)), region_list, test_band);
else
    [score_matrix, corr_val] = test_native_adj(test_patient.conn, test_patient.coords, test_patient.roi, mean_conn, nanmean(nanmean(std_conn)), region_list, test_band);
end

else
if score_denominator=='sem'
    [score_matrix, corr_val] = test_native_adj(test_patient.var, test_patient.coords, test_patient.roi, mean_var, nanmean(nanmean(sem_var)), region_list, test_band);
else
    [score_matrix, corr_val] = test_native_adj(test_patient.var, test_patient.coords, test_patient.roi, mean_var, nanmean(nanmean(std_var)), region_list, test_band);
end
end

num_elecs = size(test_patient.conn(1).data,1);
non_res_elecs = setdiff([1:num_elecs],test_patient.resect);

% output matrices
all_scores = score_matrix; % added output of all scores

% out-out scores -> get rid of resected rows and columns
out_out_scores = score_matrix;
out_out_scores(test_patient.resect,:) = NaN;
out_out_scores(:,test_patient.resect) = NaN;

% in-in scores -> get rid of non-resected rows and columns
in_in_scores = score_matrix;
in_in_scores(non_res_elecs,:) = NaN;
in_in_scores(:,non_res_elecs) = NaN;

% in-out scores -> need to separately get rid of non-resected rows and
% columns then add, then get rid of in-in resected region
in_out_scores = score_matrix;
in_out_scores(test_patient.resect,test_patient.resect) = NaN;
in_out_scores(non_res_elecs,non_res_elecs) = NaN;

end

