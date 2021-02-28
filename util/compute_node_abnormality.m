function [node_abnormality, percent_abnormal] = compute_node_abnormality(atlas_zscores, edge_threshold, node_threshold)
    
    % get number of patients
    num_patients = length(atlas_zscores);
   
    % loop through patients
    for pt = 1:num_patients
        % extract patient's zscore matrix
        pt_zscore = atlas_zscores{pt};
        
        num_nodes = size(pt_zscore,1);
        
        % find abnormal edges
        abn_edge = abs(pt_zscore)>edge_threshold;
        
        node_scores = sum(abn_edge,'omitnan');
        
        abn_node = find(node_scores>1);
        
        node_abnormality{pt} = node_scores;
        
        percent_abnormal(pt) = length(abn_node)./num_nodes;
    end
    
    
end
