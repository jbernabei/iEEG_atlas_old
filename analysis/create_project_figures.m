%% create_project_figures.m
% John Bernabei
% With assistance from Ian Ong
% Litt Laboratory

%% Set up workspace

% suppress warning
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

% load in data from excel spreadsheet
metadata = readtable("data/atlas_project_metadata.xlsx");
% place patients in a struct
all_patients = struct("patientID",metadata.Patient, ...
"outcome", metadata.Outcome,"conn",cell(length(metadata.Patient),1), ...
"roi",cell(length(metadata.Patient),1), ...
"resect",cell(length(metadata.Patient),1), ...
"hasData",cell(length(metadata.Patient),1));

% Use AAL116WM to get white matter
fileID = fopen('localization/AAL116_WM.txt');
atlas_info = textscan(fileID,'%s %s %d');
all_inds = [double(atlas_info{3})];
all_locs = [atlas_info{2}];

%% Figure 1A: visualize adjacency matrix (one panel in pipeline figure)

%% Figure 1B: visualize brain (one panel in pipeline figure)

%% Figure 2A: anatomical analysis

% set up arrays to store data
id_field = {all_patients.patientID};
conn_field = {all_patients.conn};
roi_field = {all_patients.roi};
resect_field = {all_patients.resect};
outcome_field = {all_patients.outcome};
hasData_field = {all_patients.hasData};

% load in data from all poor outcome patients
for k = 1:length(metadata.Patient)
    datapath = sprintf("data/%s/patient_data.mat",id_field{k});
    if isfile(datapath)
        d = load(datapath);
        conn_field{k} = d.II_conn;
        roi_field{k} = d.mni_coords;
        resect_field{k} = d.res_elec_inds;
        hasData_field{k} = true;
    else
        hasData_field{k} = false;
    end
end

% place data back into main struct
[all_patients.conn] = conn_field{:};
[all_patients.roi] = roi_field{:};
[all_patients.resect] = resect_field{:};
[all_patients.hasData] = hasData_field{:};

fprintf("\nAll patient data loaded.")

% load in region numbers
region_list = zeros(1,117); % for the 117 AAL regions
fi = fopen("localization/AAL116_WM.txt");
for j = 1:117
    label = split(fgetl(fi));
    region_list(j) = str2double(label{3});
end

fprintf("\nRegion list loaded.\n")

% run all patients in atlas
cond = [hasData_field{:}];
[mean_conn, std_conn] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list);

% run all good outcome patients in atlas
cond = [hasData_field{:}] & strcmp(outcome_field,'good');
[good_mean_conn, good_std_conn] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list);

% run all poor outcome patients in atlas
cond = [hasData_field{:}] & strcmp(outcome_field,'poor');
[poor_mean_conn, poor_std_conn] = create_atlas(conn_field(cond), roi_field(cond), ...
resect_field(cond), region_list);

% calculate difference between good and poor outcome means
difference_mean_conn = good_mean_conn - poor_mean_conn;

%% Figure 2B: cross - validate out-of-bag predictions

for s = 1:length(all_good_patients)
    test_patient = all_good_patients(s);
    cv_patients = all_good_patients;
    cv_patients(s) = [];
    
    % call neuroimaging atlas
    [mni_coords, mni_labels, NN_flag] = nifti_values(mni_input_coords,'AAL116_WM.nii')
    
    % get connectivity atlas
    [mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
    
    % test atlas
    %[mean_conn, std_conn] = create_atlas(all_conn, all_roi, all_resect, region_list)
end

%% Figure 2C: comparison against distance

%% Figure 3A: Distinguish between epilepsy type (mesial temporal vs neocortical)
% temporal vs extratemporal

%% Figure 3B: Relationship of connectivity to duration of epilepsy

%% Figure 3C: Relationship of connectivity to whether patient has generalized SZ or not

%% Figure 4A: localize in good outcome patients
% make figure based upon prior calculations

%% Supplement: choice of atlas / segmentation

%% Figure 4B: localize in poor outcome patients
for s = 1:length(all_poor_patients)
    
end