function [EZ_roi, EZ_inds, frac_abnl_edge] = localize_EZ_atlas(z_score_mat, edge_thresh, node_thresh, patient_ROI)
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

frac_abnl_edge = nanmean(z_score_mat>edge_thresh);
EZ_roi = find(frac_abnl_edge>node_thresh);

end