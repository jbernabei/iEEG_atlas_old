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
"coords",cell(length(metadata.Patient),1), ...
"roi",cell(length(metadata.Patient),1), ...
"resect",cell(length(metadata.Patient),1), ...
"hasData",cell(length(metadata.Patient),1));

% Use AAL116WM to get white matter
fileID = fopen('localization/AAL116_WM.txt');
atlas_info = textscan(fileID,'%s %s %d');
all_inds = [double(atlas_info{3})];
all_locs = [atlas_info{2}];

% set up arrays to store data
id_field = {all_patients.patientID};
conn_field = {all_patients.conn};
coords_field = {all_patients.coords};
roi_field = {all_patients.coords};
resect_field = {all_patients.resect};
outcome_field = {all_patients.outcome};
hasData_field = {all_patients.hasData};

% functions to help convert mni coordinates to regions of interest
f = @nifti_values;
g = @(x) f(x,"localization/AAL116_WM.nii");

% load in data from all poor outcome patients
for k = 1:length(metadata.Patient)
    datapath = sprintf("data/%s/patient_data.mat",id_field{k});
    if isfile(datapath)
        d = load(datapath);
        conn_field{k} = d.II_conn;
        coords_field{k} = d.mni_coords;
        resect_field{k} = d.res_elec_inds;
        hasData_field{k} = true;
        % convert all electrode coordinates to region names
        [~,electrode_regions,~] = cellfun(g, coords_field(k), 'UniformOutput',false);
        roi_field{k} = electrode_regions{1};
    else
        hasData_field{k} = false;
    end
end

% place data back into main struct
[all_patients.conn] = conn_field{:};
[all_patients.coords] = coords_field{:};
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

%% Figure 1A: visualize adjacency matrix (one panel in pipeline figure)

%% Figure 1B: visualize brain (one panel in pipeline figure)

%% Figure 2A: anatomical analysis

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

num_good_patients = sum([hasData_field{:}] & strcmp(outcome_field,'good'));
num_poor_patients = sum([hasData_field{:}] & strcmp(outcome_field,'poor'));

good_patient_indices = find([hasData_field{:}] & strcmp(outcome_field,'good'));

z_score_results = cell(num_good_patients,1);

% cross validate the non-resected region of good outcome patients
for s = 1:length(good_patient_indices)
    test_patient = all_patients(good_patient_indices(s));
    cv_patients = all_patients(good_patient_indices);
    cv_patients(s) = [];
    
    fprintf('\nTesting patient %d of %d', s, length(good_patient_indices))
    
    % get connectivity atlas of excluded patients
    [mean_conn, std_conn] = create_atlas({cv_patients.conn}, {cv_patients.roi}, {cv_patients.resect}, region_list);
    
    % get connectivity atlas of test patient
    [patient_conn, patient_std] = create_atlas({test_patient.conn}, {test_patient.roi}, {test_patient.resect}, region_list);
    
    % get non-resected region labels of test patient
    [~, patient_roi, ~] = nifti_values(test_patient.coords(setdiff(1:length(test_patient.roi),test_patient.resect),:),'localization/AAL116_WM.nii');
    
    % test atlas
    z_score_results{s} = test_patient_conn(mean_conn, std_conn, region_list, patient_conn, patient_roi);
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