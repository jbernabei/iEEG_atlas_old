function [data_structure_new] = process_roi_seg(data_structure_old,all_roi_list)

num_patients = length(data_structure_old);

for i = 1:num_patients
    pt_roi_list = data_structure_old(i).roi_list;
    
    
    
    num_elecs = length(pt_roi_list);
    
    pt_new_roi = zeros(num_elecs,1);
    
    laterality = contains(pt_roi_list,'Right');
    
    num_unique_roi = length(all_roi_list);
    % find new roi list
    for j = 1:(num_unique_roi)
 
        which_roi = find(contains(pt_roi_list,all_roi_list(j)));
        
        pt_new_roi(which_roi) = j;       
    end
    
    clear distance_matrix
    for j = 1:num_elecs
        
            elec1_coords = data_structure_old(i).coords(j,:);
            for m = 1:num_elecs
                elec2_coords = data_structure_old(i).coords(m,:);
                
                elec_dist = sqrt(sum((elec1_coords-elec2_coords).^2));
                
                distance_matrix(m,j) = elec_dist;
            end
    
            if laterality(j)==1 && pt_new_roi(j,1)~=0
                pt_new_roi(j,1) = pt_new_roi(j,1)+length(all_roi_list);
            end
    end
    
    % remove everything that doesn't have a new ROI
    remove_elecs = find(pt_new_roi==0);
    for q = 1:5
        data_structure_old(i).conn(q).data(:,remove_elecs) = [];
        data_structure_old(i).conn(q).data(remove_elecs,:) = [];
        
        data_structure_old(i).var(q).data(:,remove_elecs) = [];
        data_structure_old(i).var(q).data(remove_elecs,:) = [];
    end

    data_structure_old(i).roi_list(pt_new_roi==0) = [];
    data_structure_old(i).coords(pt_new_roi==0,:) = [];
    
    distance_matrix(pt_new_roi==0,:) = [];
    distance_matrix(:,pt_new_roi==0) = [];
    
    resect_bool = zeros(1,length(pt_new_roi));
    resect_bool(data_structure_old(i).resect) = 1;
    resect_bool(pt_new_roi==0) = [];
    data_structure_old(i).resect = find(resect_bool);
    
    pt_new_roi(pt_new_roi==0) = [];
    data_structure_old(i).final_roi = pt_new_roi; 
    data_structure_old(i).dist_mat = distance_matrix;
    
    size(data_structure_old(i).var(q).data)
    size(data_structure_old(i).dist_mat)
    size(pt_new_roi)
    
end

    data_structure_new = data_structure_old;

end