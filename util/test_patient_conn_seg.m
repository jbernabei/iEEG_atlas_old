function [score_mat, corr_val,corr_val_out,corr_val_in,corr_val_intra,corr_val_inter,roi_mat] = test_patient_conn_seg(mean_conn, spread_conn, curve_model, patient_coords, patient_conn, patient_roi, res_elecs)

num_elecs = size(patient_roi,1);

score_mat = NaN.*zeros(size(patient_conn));
pred_mat = NaN.*zeros(size(patient_conn));
intra_roi_mat = NaN.*zeros(size(patient_conn));
inter_roi_mat = NaN.*zeros(size(patient_conn));

roi_mat = NaN.*zeros(3,3,length(patient_conn));

% also look at frontal, temporoparietal, paralimbic
frontal_inds = [1:5,19:23];
temporoparietal_inds = [6:13,24:31];
paralimbic_inds = [14:18,32:36];
a = 0;
% loop through all electrodes
for e1 = 1:num_elecs
    which_roi1 = patient_roi(e1);
    
    if ~isempty(intersect(frontal_inds,which_roi1))
        c_row = 1;
    elseif ~isempty(intersect(temporoparietal_inds,which_roi1))
        c_row = 2;
    else
        c_row = 3;
    end
    
    for e2 = 1:num_elecs
        which_roi2 = patient_roi(e2);
        
        a = a+1;
        
        % check ROI for labeling
        if ~isempty(intersect(frontal_inds,which_roi2))
            c_col = 1;
        elseif ~isempty(intersect(temporoparietal_inds,which_roi2))
            c_col = 2;
        else
            c_col = 3;
        end
        
        
        if which_roi1~=which_roi2
            % the default comparison against the mean/spread
            atlas_conn_val = mean_conn(which_roi1,which_roi2);
            score_mat(e1,e2) = (patient_conn(e1,e2)-atlas_conn_val)./spread_conn;
            pred_mat(e1,e2) = atlas_conn_val;
            inter_roi_mat(e1,e2) = 1;
            
        else
            % comparison against within-roi distance regression
            if e1==e2
            else
            % calculate inter-elec distance
            elec1_coords = patient_coords(e1,:);
            elec2_coords = patient_coords(e2,:);
            inter_elec_dist = sqrt(sum((elec1_coords - elec2_coords).^2));
            
            % calculate expected connectivity
            exp_conn = (curve_model.p1.*inter_elec_dist+curve_model.p2)./(inter_elec_dist+curve_model.q1);
            
            % calculate expected variance
            exp_std = 0.1-0.0012*inter_elec_dist;
            
            score_mat(e1,e2) = (patient_conn(e1,e2)-exp_conn)./exp_std;
            pred_mat(e1,e2) = exp_conn;
            intra_roi_mat(e1,e2) = 1;
            end
        end
        roi_mat(c_row,c_col,a) = score_mat(e1,e2);
    end
end

% compute correlation accuracy of predicted and actual 
% *** actually should stipulate outside of resection zone ***
intra_mat = pred_mat;
intra_mat(isnan(intra_roi_mat)) = NaN;
non_nan_preds = find(~isnan(intra_mat));
[corr_val_intra, pval] = corr(patient_conn(non_nan_preds),pred_mat(non_nan_preds));

inter_mat = pred_mat;
inter_mat(inter_roi_mat~=1) = NaN;
non_nan_preds = find(~isnan(inter_mat));
[corr_val_inter, pval] = corr(patient_conn(non_nan_preds),pred_mat(non_nan_preds))

outside_rz_mat = pred_mat;
outside_rz_mat(res_elecs,:) = NaN;
outside_rz_mat(:,res_elecs) = NaN;
non_nan_preds = find(~isnan(outside_rz_mat));
[corr_val_out, pval] = corr(patient_conn(non_nan_preds),pred_mat(non_nan_preds));

inside_rz_mat = pred_mat;
non_res_elecs = setdiff([1:num_elecs],res_elecs);
inside_rz_mat(non_res_elecs,non_res_elecs) = NaN;
non_nan_preds = find(~isnan(inside_rz_mat));
try corr_val_in = corr(patient_conn(non_nan_preds),pred_mat(non_nan_preds));
catch anyerror
    corr_val_in = NaN;
end

non_nan_preds = find(~isnan(pred_mat));
corr_val = corr(patient_conn(non_nan_preds),pred_mat(non_nan_preds));

roi_mat(roi_mat==0) = NaN;

end