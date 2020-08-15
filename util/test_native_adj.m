function [native_adj_scores] = test_native_adj(all_patient_conn, patient_roi, atlas_conn, atlas_var, all_inds, test_band)
    
    patient_conn = all_patient_conn(test_band).data;

    num_elecs = size(patient_conn,1);
    
    native_adj_scores = zeros(num_elecs);
    
    % double loop through electrodes
    for i = 1:num_elecs
        for j = 1:num_elecs
            % get regions from patient roi
            region1 = patient_roi(i);
            region2 = patient_roi(j);
            
            %atlas region1 
            atlas_row = find(all_inds==region1);
            atlas_col = find(all_inds==region2);
            
            atlas_edge_val = atlas_conn(atlas_row, atlas_col);
            
            new_score = (patient_conn(i,j)-atlas_edge_val)./atlas_var;
            
            if isempty(new_score)
                native_adj_scores(i,j) = NaN;
            else
                native_adj_scores(i,j) = new_score;
            end
        end
    end
    
end