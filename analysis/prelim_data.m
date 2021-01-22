% prelim_data 
% john bernabei
clear all

metadata_table = readtable('iEEG_atlas/data/atlas_project_metadata.xlsx');

whether_res = strcmp(metadata_table{:,5},'x');
whether_conn = strcmp(metadata_table{:,6},'x');
whether_seg = strcmp(metadata_table{:,8},'x');

patients_to_use = whether_res.*whether_conn.*whether_seg;

patients = metadata_table{find(patients_to_use),1:2};

%% load all data
for i = 1:size(patients,1)
    patient_path = sprintf('data/%s/patient_data.mat',patients{i,1});
    
    load(patient_path)
    
    data_struct(i).II_conn = II_conn;
    data_struct(i).II_var = II_var;
    data_struct(i).II_var = II_var;
end
%%