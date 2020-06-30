function [z_score_mat] = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi)

% calculate z-score of all edges
z_score_mat = (patient_conn - mean_conn)./std_conn;

% extract z-scores involving relevant regions
roi_boolean = ismember(region_list,patient_roi);
z_score_mat(~roi_boolean) = NaN;

end