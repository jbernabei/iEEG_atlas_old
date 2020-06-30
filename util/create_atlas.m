function [mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
% [mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
% takes in an array of conectivity structs, an array of 3D mni coordinate
% arrays, an array of resected electrode vectors, and a vector containing
% all mni labels, and outputs two matrices representing the average and
% standard deviation of connections between all regions.
%
% Input:
%   all_conn (cell): cell array containing patient connectivity structs in
%   order
%   all_roi (cell): cell array containing regions of interest corresponding
%   to each electrode for each patient in order
%   all_resect (cell): cell array containing patient resected electrode
%   arrays in order
%   region_list (double): array containing all region labels
%
% Output:
%   mean_conn (double): (i,j) matrix of mean connectivity strengths between 
%   region_list(i) and region_list(j)
%   std_conn (double): (i,j) matrix of standard deviations of connectivity 
%   strengths between region_list(i) and region_list(j)
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 6/27/2020

% specify which band to use
band = 1;

% get number of patients
num_patients = length(all_conn);

% get number of regions
num_regions = length(region_list);

% initialize output arrays
mean_conn = zeros(num_regions);
std_conn = zeros(num_regions);

fprintf("\nCalculating connections from region        ")

% double loop through patients and regions
for i = 1:num_regions % first region
    
    % print progress to console
    fprintf("\b\b\b\b\b\b\b%d...",region_list(i))
    
    for j = i:num_regions % second region
        
        % initialize list of connection strengths between region i and j
        reg_conn_strengths = zeros(1, num_patients);
        
        for p = 1:num_patients
            
            % get electrode regions for the patient
            patient_electrode_regions = all_roi{p};
            % get resected electrodes for the patient
            res_elec_inds = all_resect{p};
            % orient vector vertically
            res_elec_size = size(res_elec_inds);
            if res_elec_size(1) == 1
                res_elec_inds = res_elec_inds.';
            end
            
            % calculate logical with resected indices
            resect_boolean = cat(2,accumarray(res_elec_inds,1).',...
            zeros(1,length(patient_electrode_regions)-max(res_elec_inds)));
            % get electrodes contained within first region
            first_reg_elec = find(patient_electrode_regions == region_list(i) & ~resect_boolean);
            % get electrodes contained within second region
            second_reg_elec = find(patient_electrode_regions == region_list(j) & ~resect_boolean);
            
            % get number of electrodes in each region
            num_first_reg_elec = length(first_reg_elec);
            num_second_reg_elec = length(second_reg_elec);
            
            % initialize list of connection strengths corresponding to the
            % two regions in this patient
            patient_conn_strengths = zeros(1, num_first_reg_elec*num_second_reg_elec);
            
            % get the patient's adjacency matrices
            adj_matrix = all_conn{p};
            
            % find strengths of connections from region i to region j
            for b = 1:num_first_reg_elec
                for c = 1:num_second_reg_elec
                    
                    % get electrode numbers
                    elec1 = first_reg_elec(b);
                    elec2 = second_reg_elec(c);
                    
                    % get current index
                    current_index = ((b-1)*num_second_reg_elec)+c;
                    
                    if elec1 ~= elec2
                    
                        % get connection strength from the specified band from 
                        % the connectivity matrices
                        patient_conn_strengths(current_index) = adj_matrix(band).data(elec1,elec2);
                        
                    else
                        
                        % mark for deletion if both electrodes are the same
                        patient_conn_strengths(current_index) = NaN;
                        
                    end
                end
            end
                      
            % remove NaN values
            patient_conn_strengths = patient_conn_strengths(~isnan(patient_conn_strengths));
            
            % store the average connection strength of the current patient
            reg_conn_strengths(p) = mean(patient_conn_strengths);
 
        end
        
        % remove NaN values
        reg_conn_strengths = reg_conn_strengths(~isnan(reg_conn_strengths));
        
        % find the average connection strength across all patients
        avg_reg_conn_strength = mean(reg_conn_strengths);
            
        % find stdev
        stdev_reg_conn_strength = std(reg_conn_strengths);
            
        % place avg and stdev in output arrays (this fills the upper
        % triangle of both arrays)
        mean_conn(i,j) = avg_reg_conn_strength;
        std_conn(i,j) = stdev_reg_conn_strength;
        
    end
end

% symmetrize both output matrices
mean_conn = triu(mean_conn) + tril(mean_conn.',-1);
std_conn = triu(std_conn) + tril(std_conn.',-1);

fprintf("\nSuccessfully generated atlas for band %d.\n", band)

end