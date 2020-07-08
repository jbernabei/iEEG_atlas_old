function [mean_conn, std_conn, num_conn] = create_atlas(all_conn, all_roi, all_resect, region_list, band)
% [mean_conn, std_conn, num_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
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
%   band (int): frequency band to be used
%
% Output:
%   mean_conn (double): (i,j) matrix of mean connectivity strengths between 
%   region_list(i) and region_list(j)
%   std_conn (double): (i,j) matrix of standard deviations of connectivity 
%   strengths between region_list(i) and region_list(j)
%   num_conn (double): (i,j) matrix of number of patients with each edge
%   measured (off diagonal) and each node measured (on diagonal)
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 7/5/2020

% get number of patients
num_patients = length(all_conn);

% get number of regions
num_regions = length(region_list);

% initialize output array to store connection strengths
mean_conn = NaN(90,90,num_patients);

fprintf("\nCalculating connections  ")

% for aesthetics
spinner = ['|','/','-','\'];

for p = 1:num_patients
    
    % display spinner
    fprintf("\b%c",spinner(floor(mod(p/8,4)+1)))
    
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

    % get connection strengths between each pair of electrodes
    adj_matrix = all_conn{p};
    % extract band data
    band_matrix = adj_matrix(band).data;

    % double loop through regions
    for i = 1:90 % first region
        
        % get electrodes contained within first region
        first_reg_elec = (patient_electrode_regions == region_list(i) & ~resect_boolean);
        
        % extract rows corresponding to electrodes in the first region
        first_reg_strengths = band_matrix(first_reg_elec,:);
        
        % calculate connections within the region
        patient_strengths = first_reg_strengths(:,first_reg_elec);
        patient_strength = mean(patient_strengths(triu(true(size(patient_strengths)),1)));
        
        % add to output array if the region contains any electrodes
        if ~isnan(patient_strength)
            mean_conn(i,i,p) = patient_strength;
        else
            % skip calculations for this region pair if the first region
            % does not contain any electrodes
            continue
        end

        for j = (i+1):90 % second region
            
            % get electrodes contained within second region
            second_reg_elec = (patient_electrode_regions == region_list(j) & ~resect_boolean);

            % extract connection strengths between the two regions
            patient_strengths = first_reg_strengths(:,second_reg_elec);
            
            % average connection strengths between the two regions
            patient_strength = sum(sum(patient_strengths))/numel(patient_strengths);
            
            % add to output array
            mean_conn(i,j,p) = patient_strength;
        end
    end
end

% calculate num conn
for i = 1:90
    for j = 1:90
        num_conn(i,j) = sum(~isnan(mean_conn(i,j,:)));
    end
end

% take the standard deviation element-wise to get the std matrix
std_conn = std(mean_conn,0,3,'omitnan');

% divide out the number of patients element-wise to get the mean matrix
mean_conn = mean(mean_conn,3,'omitnan');

% symmetrize both output matrices
mean_conn = triu(mean_conn) + tril(mean_conn.',-1);
std_conn = triu(std_conn) + tril(std_conn.',-1);
num_conn = triu(num_conn) + tril(num_conn.',-1);

fprintf("\b\b...\nSuccessfully generated atlas for band %d.\n", band)

end