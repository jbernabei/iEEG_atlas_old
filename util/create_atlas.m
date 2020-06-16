function [mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)


% get number of patients
num_patients = length(all_conn);

% get number of regions
num_regions = length(region_list);

mean_conn = zeros(num_regions);
std_conn = zeros(num_regions);

% double loop through patients and regions
for i = 1:num_regions
    for j = 1:num_regions
        for s = 1:num_patients
            
            % extract non-resected connections
            
            % find connections from region i to region j
            
            % find average
            
            % find stdev
            
        end
    end
end



end