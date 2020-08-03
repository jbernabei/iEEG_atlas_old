function [] = create_data_table(all_patients, hasData_field)

% get patients who have data
study_patients = all_patients([hasData_field{:}]);
good_patients = study_patients(strcmp({study_patients.outcome},'good'));
poor_patients = study_patients(strcmp({study_patients.outcome},'poor'));

% create output table
row_names = {'Total number of subjects','Age at surgery',sprintf('Mean \x00B1 SD'), ...
    'Age at diagnosis',sprintf('Mean \x00B1 STD'),'Sex','Male','Female', ...
    'Electrode type','ECoG','SEEG','Resected/ablated region','MTL','Temporal', ...
    'MFL','Frontal','Parietal/FP','Insular','MRI','Lesional','Non-lesional', ...
    'Unknown','Therapy','Ablation','Resection'};

values_to_count = {'M','F','ECoG','SEEG','MTL','Temporal','MFL', ...
    'Frontal','Parietal','Insular','Lesional','Non-Lesional','', ...
    'Ablation','Resection'};

remove_rows = [1,2,3,4,5,6,9,12,19,23];
fill_rows = 1:length(row_names);
fill_rows(remove_rows) = [];

% dictionary for target lobe names
% target_map = containers.Map({'Frontal','Temporal','FP','Parietal', ...
%     'MTL','MFL','Insular'},{'FL/PL/FPL','TL','FL/PL/FPL','FL/PL/FPL', ...
%     'TL','FL/PL/FPL','FL/PL/FPL'});
% get_target = @(x) target_map(x);

% set up table

structs = {good_patients, poor_patients};
demographic_table = cell2table(cell(length(row_names),length(structs)));
for b = 1:length(structs)
    patients = structs{b};
    
    % prepare data to be read into table
    targets = {patients.target};
    targets{strcmp(targets,'FP')} = 'Parietal';
    [patients.target] = targets{:};
   
    patient_attributes = {patients.gender;patients.gender; ...
        patients.implant;patients.implant;patients.target; ...
        patients.target;patients.target;patients.target; ...
        patients.target;patients.target;patients.lesion_status; ...
        patients.lesion_status;patients.lesion_status; ...
        patients.therapy;patients.therapy}.';
%     targets = cellfun(get_target,{patients.target},'UniformOutput',false);
%     targets_and_lateralities = join([{patients.laterality}.',targets.'],2);
%     target_counts = countcats(categorical(targets_and_lateralities));
    data_column = cell(length(row_names),1);
    % Total number of subjects
    data_column{1} = length(patients);
    % Age at surgery
    data_column{3} = sprintf('%.1f \x00B1 %.1f',nanmean([patients.age_surgery]),nanstd([patients.age_surgery]));
    % Age at diagnosis
    data_column{5} = sprintf('%.1f \x00B1 %.1f',nanmean([patients.age_onset]),nanstd([patients.age_onset]));
    
    for w = 1:length(fill_rows)
        data_column{fill_rows(w)} = sum(strcmp(patient_attributes(:,w), values_to_count{w}));
    end
    
    demographic_table(:,b) = data_column;
end

demographic_table.Properties.RowNames = row_names;
demographic_table.Properties.VariableNames = {'Good_surgical_outcome','Poor_surgical_outcome'};

writetable(demographic_table,'output/patient_demographics.xlsx','WriteRowNames',true)

fprintf('Demographic table saved.\n')

end