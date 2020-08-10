function [out_out_scores, in_in_scores, in_out_scores, all_scores] = validate_patient(test_patient,cv_patients,region_list,test_band,test_threshold,score_denominator,atlas_method)

% score_denominator is a string that can be either "std" or "sem"
% atlas_method is a string that can be either "patient" or "edge"

% sets default denominator to "std" if none is given
if ~exist('score_denominator','var'), score_denominator = "std"; end
% sets default method to "patient" if none is given
if ~exist('atlas_method','var'), atlas_method = "patient"; end

% get connectivity atlas of excluded patients
[mean_conn, std_conn, ~, sem_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list, test_band, test_threshold);

% get connectivity atlas of test patient
if atlas_method == "edge"
    patient_conn = create_atlas_by_edge({test_patient.conn}, {test_patient.roi}, {[]}, region_list, test_band);
else
    patient_conn = create_atlas({test_patient.conn}, {test_patient.roi}, {[]}, region_list, test_band);
end
    
% test patient
if score_denominator == "sem"
    score_matrix = test_patient_conn(mean_conn, sem_conn, region_list, patient_conn);
else
    score_matrix = test_patient_conn(mean_conn, std_conn, region_list, patient_conn);
end

% get labels of regions which contain non-resected electrodes
regions_with_non_res_elec = test_patient.roi(setdiff(1:length(test_patient.roi),test_patient.resect));

% get labels of regions which contain resected electrodes
regions_with_res_elec = test_patient.roi(test_patient.resect);

% get regions that only contain non-resected electrodes
non_res_regions = setdiff(regions_with_non_res_elec,regions_with_res_elec);

% get regions that only contain resected electrodes
res_regions = setdiff(regions_with_res_elec,regions_with_non_res_elec);

% convert region labels to indices on score_matrix
[~,non_res_indices] = ismember(non_res_regions,region_list);
[~,res_indices] = ismember(res_regions,region_list);
non_res_indices(non_res_indices == 0) = [];
res_indices(res_indices == 0) = [];

% calculate matrix to determine which edge is of which type
type_matrix = zeros(size(score_matrix));
type_matrix(non_res_indices,:) = type_matrix(non_res_indices,:)+1;
type_matrix(:,non_res_indices) = type_matrix(:,non_res_indices)+1;
type_matrix(res_indices,:) = type_matrix(res_indices,:)-1;
type_matrix(:,res_indices) = type_matrix(:,res_indices)-1;
type_matrix = floor(type_matrix./2);

% output matrices
all_scores = score_matrix; % added output of all scores

out_out_scores = score_matrix;
out_out_scores(~logical(logical(type_matrix)+type_matrix)) = NaN;

in_in_scores = score_matrix;
in_in_scores(~logical(-logical(type_matrix)+type_matrix)) = NaN;

in_out_scores = score_matrix;
in_out_scores(logical(type_matrix)) = NaN;

end

