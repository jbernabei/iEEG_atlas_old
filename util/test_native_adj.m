function [native_adj_scores, corr_val, virtual_resect_mat] = test_native_adj(all_patient_conn, pt_mni, patient_roi, atlas_conn, atlas_var, all_inds, test_band)
    
    patient_conn = all_patient_conn(test_band).data;

    num_elecs = size(patient_conn,1);
    
    native_adj_scores = zeros(num_elecs);
    
    pred_mat = zeros(num_elecs);
    
    % double loop through electrodes
    for i = 1:num_elecs
        for j = 1:num_elecs
            % get regions from patient roi
            region1 = patient_roi(i);
            region2 = patient_roi(j);
            
            if region1==region2
                if i==j
                else
                    % calculate inter-elec distance
                    elec1_coords = pt_mni(i,:);
                    elec2_coords = pt_mni(j,:);
                    inter_elec_dist = sqrt(sum((elec1_coords - elec2_coords).^2));
            
                    % expected conn
                    exp_conn = (0.0334.*inter_elec_dist+3.769)./(inter_elec_dist+2.139); % code this better
                    
                    % expected std
                    exp_std = 0.1-0.0012*inter_elec_dist; % code this better
            
                    native_adj_scores(i,j) = (patient_conn(i,j)-exp_conn)./exp_std;
                    
                    pred_mat(i,j) = exp_conn;
                end
            else
                %atlas region1 
                atlas_row = find(all_inds==region1);
                atlas_col = find(all_inds==region2);

                atlas_edge_val = atlas_conn(atlas_row, atlas_col);
                
                atlas_var = 0.73.*atlas_edge_val-0.146;

                new_score = (patient_conn(i,j)-atlas_edge_val)./atlas_var;

                if isempty(new_score)
                    native_adj_scores(i,j) = NaN;
                    pred_mat(i,j) = NaN;
                else
                    native_adj_scores(i,j) = new_score;
                    pred_mat(i,j) = atlas_edge_val;
                end
            end
        end
    end
    
    % find predictions which are not NaN
    non_nan_preds = find(~isnan(pred_mat));
    
    % set up real vectors of real connectivity and predicted connectivity
    real_conn = patient_conn(non_nan_preds);
    pred_conn = pred_mat(non_nan_preds);
    
    
    
    % get a global correlation
    corr_val = corr(real_conn, pred_conn);
    
    % get number of edges
    num_real_edges = length(real_conn);
    
    % loop through edges and do 'virtual edge resection'
    for e = 1:num_real_edges
        real_conn_res = real_conn;
        real_conn_res(e) = [];
        
        pred_conn_res = pred_conn;
        pred_conn_res(e) = [];
        
        % see whether correlation goes up (abnormal) or down (normal)
        virtual_resect_edge(e) = (corr(real_conn_res, pred_conn_res)-corr_val)./corr_val;
    end
    
    virtual_resect_mat = NaN*zeros(num_elecs);
    virtual_resect_mat(non_nan_preds) = virtual_resect_edge;
end