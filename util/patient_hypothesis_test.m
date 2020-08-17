function [in_in_1_conns, in_in_2_conns] = patient_hypothesis_test(zscores, patient_ROI, region_list, region1, region2)
    % find electrode contacts in ROI 1
    region_1_nums = region_list(region1);
    
    % find electroded contacts in ROI 2
    region_2_nums = region_list(region2);
    
    for i = 1:length(patient_ROI)
        if sum(patient_ROI(i) == region_1_nums)>0
            this_roi(i) = 1;
        elseif sum(patient_ROI(i) == region_2_nums)>0
            this_roi(i) = 2;
        else
            this_roi(i) = 0;
        end
    end
    
    zscores = zscores-zscores.*eye(size(zscores));
    
    in_in_1_conns = zscores((this_roi==1),:);
    in_in_2_conns = zscores((this_roi==2),:);
    
end