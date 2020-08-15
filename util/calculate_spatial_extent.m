function [spatial_extent] = calculate_spatial_extent(z_score_matrices, all_dists)
    
    % get number of patients
    num_patients = length(z_score_matrices)
    
    for pt = 1:num_patients
        
        clear mean_in_in_abnormality
        clear mean_distance
        
        %[~, atlas_distance] = create_distance_matrix(mni_coords{pt}, atlas_inds);
        atlas_distance = all_dists{pt};
        
        % extract patients z score matrices
        pt_zscores = z_score_matrices{pt};
        
        % find which connections are totally unsampled
        unsampled_conns = find(isnan(nanmean(pt_zscores)));
        
        % reduce samples
        reduced_pt_zscores = pt_zscores;
        reduced_pt_zscores(unsampled_conns,:) = [];
        reduced_pt_zscores(:,unsampled_conns) = [];
        
        % zero out NaNs
        reduced_pt_zscores(find(isnan(reduced_pt_zscores))) = 0;
        
        % compute modularity
        [modules, Q] = modularity_und(reduced_pt_zscores,1);
        
        % get number of modules
        num_modules = length(unique(modules));
        
        % loop through modules
        for m = 1:num_modules
           % find this modules nodes
           module_nodes = find(modules==m);
           
           mean_in_in_abnormality(m) = nanmean(nanmean(reduced_pt_zscores(module_nodes,module_nodes)));
           
           mean_distance(m) = nanmean(nanmean(atlas_distance(module_nodes,module_nodes)));
           
           
        end
        
        % extract spatial extent of abnormality modules and level of abnormality
        spatial_extent(pt).abnormality = mean_in_in_abnormality;
        spatial_extent(pt).distance = mean_distance;
        
    end
    
end