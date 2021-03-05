function all_patients = remove_all_wm(all_patients_raw)

all_patients = all_patients_raw;
data_patient_indices = find([all_patients.hasData]);
for s = 1:length(data_patient_indices)
    
    patient_data = all_patients(data_patient_indices(s));
    
    num_elecs = length(patient_data.roi);
    resect_bool = zeros(num_elecs,1);
    resect_bool(patient_data.resect) = 1;
    
    wm_elecs = [find(patient_data.roi==9171),find(patient_data.roi==0),find(patient_data.in_brain==0)']; % get which electrodes are in WM
    resect_bool(wm_elecs) = [];
    
    final_res_elecs = find(resect_bool);
    
    patient_data.resect = final_res_elecs;
    
    for f = 1:5
    base_adj = patient_data.conn(f).data;
    final_adj = base_adj;
    final_adj(:,wm_elecs) = [];
    final_adj(wm_elecs,:) = [];
    
    patient_data.conn(f).data = final_adj;
    
    base_var = patient_data.var(f).data;
    final_var = base_var;
    final_var(:,wm_elecs) = [];
    final_var(wm_elecs,:) = [];
    
    patient_data.var(f).data = final_var;
    end
    
    final_roi = patient_data.roi;
    final_roi(wm_elecs) = [];
    
    patient_data.roi = final_roi;
    
    final_coords = patient_data.coords;
    final_coords(wm_elecs,:) = [];
    patient_data.coords = final_coords;
    
    all_patients(data_patient_indices(s)) = patient_data;
    
    size(all_patients(data_patient_indices(s)).conn(f).data)
    size(all_patients(data_patient_indices(s)).coords)
    size(all_patients(data_patient_indices(s)).roi)
    
end