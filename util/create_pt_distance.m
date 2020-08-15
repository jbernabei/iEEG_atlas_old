function [pt_distances] = create_pt_distance(all_mni_coords)
    num_pt = length(all_mni_coords);
    
    for pt = 1:num_pt
        this_pt_coords = all_mni_coords{pt};
        num_elecs = size(this_pt_coords,1);
        
        clear this_pt_distance
        
        for i = 1:num_elecs
            for j = 1:num_elecs
                
                this_pt_distance(i,j) = sqrt(sum((this_pt_coords(i,:)-this_pt_coords(j,:)).^2));
                
            end
        end
        
        pt_distances{pt} = this_pt_distance;
        
    end
end