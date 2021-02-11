function [node_str_all, final_res_elecs] = node_str_no_wm(patient_data)
    % this function computes distinguishability between resected and
    % non-resected regions for 
    
    num_elecs = length(patient_data.roi);
    resect_bool = zeros(num_elecs,1);
    resect_bool(patient_data.resect) = 1;
    
    wm_elecs = find(patient_data.roi==9171); % get which electrodes are in WM
    resect_bool(wm_elecs) = [];
    
    final_res_elecs = find(resect_bool);
    
    base_adj = patient_data.conn(1).data;
    final_adj = base_adj;
    final_adj(:,wm_elecs) = [];
    final_adj(wm_elecs,:) = [];
    
    node_str_all = sum(final_adj)'; % sum the adjacency matrix (w/o WM) to get final node strengths
    

end