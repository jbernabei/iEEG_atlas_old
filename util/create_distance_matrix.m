function [dist_mat] = create_distance_matrix(all_mni, all_roi, all_resect, region_list, band)
% add function info
%
% Input:
%   all_mni (cell): cell array containing patient mni coordinates in order
%   all_roi (cell): cell array containing regions of interest corresponding
%   to each electrode for each patient in order
%   all_resect (cell): cell array containing patient resected electrode
%   arrays in order
%   region_list (double): array containing all region labels
%   band (int): frequency band to be used
%
% Output:
%   mean_mat (double): (i,j) matrix of mean distances between 
%   region_list(i) and region_list(j)
%
% John Bernabei and Ian Ong
% johnbe@seas.upenn.edu
% ianzyong@seas.upenn.edu
% 7/5/2020

end