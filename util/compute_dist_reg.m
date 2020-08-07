function [curve, dist_data, conn_data] = compute_dist_reg(adj_matrices, mni_coordinates, patient_roi, resected_elecs)
    
    all_conn = []; % initialize elec x 5 freq
    all_dist = [];

    % get number of patients
    num_patients = length(adj_matrices);
    
    % loop through patients
    for pt = 1:num_patients
        
        clear pt_conn
        clear pt_dist
        
        pt_roi = patient_roi{pt};
        
        pt_res_elecs = resected_elecs{pt};
        
        % get number of electrodes
        num_elecs = size(mni_coordinates{pt},1);
        
        % assign into mni coordinates
        mni_coords = mni_coordinates{pt};
        
        % double loop through electrodes and calculate distance matrix
        for i = 1:num_elecs
            for j = 1:num_elecs
                % have to add in resecte elecs here
                pt_dist(i,j) = sqrt(sum((mni_coords(i,:)-mni_coords(j,:)).^2));
            end
        end
        
        bad_inds = [find(pt_roi>9000), find(pt_roi==0), pt_res_elecs(:)'];
        pt_dist(bad_inds,:) = [];
        pt_dist(:,bad_inds) = [];
        
        all_pt_dist = pt_dist(:);
        
        for f = 1:5
            adj = adj_matrices{pt}(f).data;
            adj(bad_inds,:) = [];
            adj(:,bad_inds) = [];
            pt_conn(:,f) = adj(:);
        end
        
        all_conn = [all_conn;pt_conn];
        all_dist = [all_dist;all_pt_dist];
        
    end
    
    % extract self connections from data
    self_conn = find(all_dist == 0);
    all_dist(self_conn) = [];
    all_conn(self_conn,:) = [];
    
    dist_data = all_dist;
    conn_data = all_conn;
    
    % do the curve fitting
    for f = 1:5
        curve(f).data = fit(all_dist,all_conn(:,f),'power2');
        %[fitresult, gof] = powerFit(all_dist,all_conn);
        %curve(f).data = gof;
    end

end