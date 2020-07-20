function [atlas_mni] = create_distance_matrix(all_mni, region_list)
% add function info
%
% Input:
%   all_mni (cell): cell array containing patient mni coordinates in order
%   all_roi (cell): cell array containing regions of interest corresponding
%   to each electrode for each patient in order
%   all_resect (cell): cell array containing patient resected electrode
%   arrays in order
%   region_list (double): array containing all region labels
%
% Output:
%   atlas_mni (double): (n x 3) matrix of ROI centroids
%   mean_mat (double): (i,j) matrix of mean distances between 
%   region_list(i) and region_list(j)
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 7/5/2020

num_pts = length(all_mni);

all_elecs = [];

for s = 1:num_pts
    all_elecs = [all_elecs; all_mni{s}];
    
end

[~, all_roi, ~] = nifti_values(all_elecs,'localization/AAL116_WM.nii');

for r = 1:90
    this_roi = region_list(r);
    
    this_roi_inds = find(all_roi==this_roi);
    
    atlas_mni(r,1:3) = mean(all_elecs(this_roi_inds,:));
end

end