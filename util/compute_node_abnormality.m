function [node_abnormality, percent_abnormal] = compute_node_abnormality(atlas_zscores, edge_threshold, node_threshold)
    
    % get number of patients
    num_patients = length(atlas_zscores);
   
    % loop through patients
    for pt = 1:num_patients
        % extract patient's zscore matrix
        pt_zscore = atlas_zscores{pt};
        
        num_nodes = size(pt_zscore,1);
        
        % find abnormal edges
        abn_edge = pt_zscore>edge_threshold;
        
        node_scores = mean(abn_edge,'omitnan');
        
        abn_node = find(node_scores>node_threshold);
        
        node_abnormality{pt} = abn_node;
        
        percent_abnormal(pt) = length(abn_node)./num_nodes;
    end
    
    
end
