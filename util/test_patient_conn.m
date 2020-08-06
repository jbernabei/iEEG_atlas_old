function [score_mat, corr_val, residuals] = test_patient_conn(mean_conn, spread_conn, region_list, patient_conn, patient_roi)

% output whole z-score matrix if no ROIs are specified
if ~exist('patient_roi','var'), patient_roi = region_list; end

% calculate standardized score of all edges
score_mat = (patient_conn - mean_conn)./spread_conn;

% calculate residuals of all edges
residuals = patient_conn - mean_conn;

% extract atlas and patient edges
patient_edges = patient_conn(:);
atlas_edges = mean_conn(:);

% find NaNs
patient_NaN = find(isnan(patient_edges));
atlas_NaN = find(isnan(atlas_edges));

% get rid of NaNs
all_NaN = unique([patient_NaN;atlas_NaN]);
patient_edges(all_NaN) = [];
atlas_edges(all_NaN) = [];

% get correlation value
size(atlas_edges);
size(patient_edges);
if isempty(patient_edges)
    corr_val = NaN;
else
    corr_val = corr(patient_edges,atlas_edges);
end

% extract z-scores involving relevant regions
roi_boolean = ismember(region_list,patient_roi);
score_mat(:,~roi_boolean) = NaN;
score_mat(~roi_boolean,:) = NaN;

residuals(:,~roi_boolean) = NaN;
residuals(~roi_boolean,:) = NaN;

end