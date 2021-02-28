function [all_scores, corr_val, virtual_resect] = get_new_patient_scores(test_patient,cv_patients,region_list,test_band,test_threshold,atlas_method,conn_type)

% atlas_method is a string that can be either "native" or "atlas"

% get connectivity atlas of excluded patients
if conn_type==1
    [mean_conn, std_conn, ~, sem_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);
else
    [mean_var, std_var, ~, sem_var] = create_atlas({cv_patients.var}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);
end

% set up dummy corr value
corr_val = [];

if strcmp(atlas_method,'atlas')
    if conn_type==1
        [pt_atlas, pt_std, ~, sem_conn] = create_atlas({test_patient.conn}, {test_patient.roi}, {[]}, region_list, test_band, 1);
    else
        [pt_atlas, pt_std, ~, sem_conn] = create_atlas({test_patient.var}, {test_patient.roi}, {[]}, region_list, test_band, 1);
    end
    % need to fix the resected electrodes
    [score_matrix, corr_val, virtual_resect] = test_native_adj(pt_atlas, test_patient.coords, region_list, mean_conn, nanmean(nanmean(sem_conn)), region_list, test_band);
    % need to do it without coordinates
elseif strcmp(atlas_method,'native')
    if conn_type==1
        [score_matrix, corr_val, virtual_resect] = test_native_adj(test_patient.conn, test_patient.coords, test_patient.roi, mean_conn, nanmean(nanmean(std_conn)), region_list, test_band);
    else
        [score_matrix, corr_val, virtual_resect] = test_native_adj(test_patient.var, test_patient.coords, test_patient.roi, mean_conn, nanmean(nanmean(std_conn)), region_list, test_band);
    end
end


num_elecs = size(test_patient.conn(1).data,1);


if strcmp(atlas_method,'atlas')
    res_elecs = find(sum(test_patient.roi(test_patient.resect)'==region_list)>0);
else
    res_elecs = test_patient.resect;
end
non_res_elecs = setdiff([1:num_elecs],res_elecs);

% output matrices
all_scores = score_matrix; % added output of all scores

% out-out scores -> get rid of resected rows and columns
out_out_scores = score_matrix;
out_out_scores(res_elecs,:) = NaN;
out_out_scores(:,res_elecs) = NaN;

% in-in scores -> get rid of non-resected rows and columns
in_in_scores = score_matrix;
in_in_scores(non_res_elecs,:) = NaN;
in_in_scores(:,non_res_elecs) = NaN;

% in-out scores -> need to separately get rid of non-resected rows and
% columns then add, then get rid of in-in resected region
in_out_scores = score_matrix;
in_out_scores(res_elecs,res_elecs) = NaN;
in_out_scores(non_res_elecs,non_res_elecs) = NaN;

end

