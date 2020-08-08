function [all_pt_dice] = localize_EZ_atlas(z_score_mat, pt_roi, resected_elecs, mni_coordinates)
% [mean_conn, std_conn, num_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
% takes in an array of conectivity structs, an array of 3D mni coordinate
% arrays, an array of resected electrode vectors, and a vector containing
% all mni labels, and outputs two matrices representing the average and
% standard deviation of connections between all regions.
%
% Input:
%   z_score_mat (double): (i,j) a ROI x ROI matrix of z scores for all of that patient's
%   connectivity when compared to the atlas
%   threshold (float): average z score for each node to be localized to
%   epileptogenic zone
%   patient_ROI (cell): a cell array of the ROI for each electrode index
%
% Output:
%   EZ_ROI (cell): a list of ROI which are determined to likely comprise
%   the epileptogenic zone
%   EZ_inds (double): a list of electrodes to target with resection based
%   on EZ_ROI and patient_ROI
%   frac_abnl_edge (double): a list of the fraction of edges for each node
%   that are deemed to be abnormal by the algorithm
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 7/6/2020

% get number of patients
num_patients = length(z_score_mat);

% loop through number of patients 
for pt = 1:num_patients
    
    clear dice_scores
    clear metric_resect
    clear reduced_adj
    clear node_group
    
    % extract z score matrix for that patient
    this_patient_mat = z_score_mat{pt};
    
    % find fraction of abnormal edges connected to each node
    %frac_abnl_edge = nanmean(this_patient_mat>edge_thresh);
    
    % find regions which have fractions of abnormal edges above the nodal
    % threshold
    %EZ_roi = find(frac_abnl_edge>node_thresh);
    
    % get num_roi
    patient_roi = pt_roi{pt};
    
    mni_coords = mni_coordinates{pt};
    
    % need to get centroids of each ROI
     % find unique brain regions
    unique_roi = unique(patient_roi{pt});

    a = 0;

    % loop through unique brain regions
    for i = 1:length(unique_roi)

        % find which nodes are in these
        nodes_1 = find(patient_roi{pt}==unique_roi(i));

        if size(resected_elecs{pt},2)==1
            resected_elecs{pt} = resected_elecs{pt}';
        end

        % find centroid of these regions
        centroid_1 = mean(mni_coords{pt}(nodes_1,:),1);
        pt_centroid{pt}.data(i,:) = centroid_1;

        % old code
         % check if its a resected region
        if sum(sum(nodes_1'==resected_elecs{pt}))>0
            a = a+1;
            resect_region(pt).data(a) = i;
        end

        % loop through regions again
        for j = 1:length(unique_roi)

            % find second set of nodes
            nodes_2 = find(patient_roi{pt}==unique_roi(j));

            % find second set of centroidss
            centroid_2 = mean(mni_coords{pt}(nodes_2,:),1);
            if i==j
            else

            end
        end

    end
    
    % once we have centroids of ROI loop through all in iterative algorithm
    
    % get number of regions to use in iterative algorithm
    for r = 1:length(unique_roi)
        % get number of regions to use
        
        
        % now we want to find spatially constrained number of regions with
        % best dice -> assess by intra-ROI connectivity -> find inflection
        % point probably
        
        % loop through nodes
        for i = 1:length(unique_roi)

            % get coords from centroid
            this_roi_coords = pt_centroid{pt}.data(i,:);

            % find all dists
            for j = 1:num_nodes
                dist(pt).data(j)= sqrt(sum((this_node_coords-pt_centroid{pt}.data(j,:)).^2));
            end

            % sort by distance
            [B,I] = sort(dist(pt).data,'ascend');

            % we'll try an increasing number of targets
            num_targets = r;

            % select group of nodes
            node_group(i).data = I(1:(num_targets));

            % get the metric required
            metric_resect(i) = nanmean(nanmean(z_score_mat{pt}(node_group(i).data,node_group(i).data)));

        end
        
     % find resection w/ greatest total node strength resected
    [y, which_resection] = max(metric_resect);

    targets_1 = node_group(which_resection).data;
            
    %targets_1
    %resected_elecs{pt}
    % compute dice
    dice_scores(r) = dice(targets_1,resect_region(pt).data);
            
    end
    
    all_pt_dice(pt).data = dice_scores;
   

end

end