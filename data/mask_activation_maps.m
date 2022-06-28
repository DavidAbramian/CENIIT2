
% Parameters
z = 5;
c = 0;

% Find subject folders
dataPath = '/flush/davab27/CENIIT/data';
dirSubj = dir(dataPath);
names = string({dirSubj(:).name}');
ISubjects = cellfun(@(x)~isempty(x), regexp(names, '[0-9]{5}'));  % Check if the file/folder name consists of 5 numbers
dirSubj = dirSubj(ISubjects);

nSubjects = length(dirSubj);

tasks = {'motor'; 'verb'; 'word'};
FWHMs = [4, 6];
contrasts = {["motor_finger"; "motor_foot"; "motor_lips"]; "verb_generation"; "word_repetition" };
motion_params = 1:2;

% a = -ones(20,6);
% b = -ones(20,6);
% c = -ones(20,6);
% d = -ones(20,6);

% for s = 1:nSubjects
for s = 6
    
    subject_id = dirSubj(s).name;
    
    % Find tumor voxel indices
%     tumorFile = fullfile(dataPath, subject_id, 'anat', 'tumor_mask_iso.nii.gz');
    tumorFile = fullfile(dataPath, subject_id, 'anat', 'tumor_mask_new.nii.gz');
    [h,v] = ml_load_nifti(tumorFile);
    I_tumor = find(v);
    
    % Add margin to tumor mask
    SE = strel('sphere', 2);    
    v2mm = convn(v, SE.Neighborhood, 'same') > 0;
    I_tumor_2mm = find(v2mm);
    
    SE = strel('sphere', 4);    
    v4mm = convn(v, SE.Neighborhood, 'same') > 0;
    I_tumor_4mm = find(v4mm);

    dirActivationMaps = fullfile(dataPath, subject_id, 'fmri_stats', ['z', num2str(z), '_c', num2str(c)]);
%     files = dir(fullfile(dirActivationMaps, 'verb*0.nii.gz'));
    files = dir(fullfile(dirActivationMaps, '*mm.nii.gz'));
    
    %% Mask activation maps with brain tumor mask
    
    for i = 1:length(files)
        
        fprintf('%s, %s \n', subject_id, files(i).name)
        
        % Load activation map
        fileActMap = fullfile(dirActivationMaps, files(i).name);
        [h,v] = ml_load_nifti(fileActMap);
        
%         fprintf('%i \n', nnz(v(I_tumor)))
        fprintf('%i \n', nnz(v))
        
%         a(i,s) = nnz(v);
%         b(i,s) = nnz(v(I_tumor));
%         c(i,s) = nnz(v(I_tumor_2mm));
%         d(i,s) = nnz(v(I_tumor_4mm));

        % Mask out tumor voxels
        v(I_tumor) = 0;
        
        % Save masked activation map
        h.fname = [fileActMap(1:end-7), '_m.nii.gz'];
        h.dt(1) = 2;
        spm_write_vol(h,v);
        
        % Mask out tumor voxels
        v(I_tumor_2mm) = 0;
        
        % Save masked activation map
        h.fname = [fileActMap(1:end-7), '_m2.nii.gz'];
        h.dt(1) = 2;
        spm_write_vol(h,v);
        
        % Mask out tumor voxels
        v(I_tumor_4mm) = 0;
        
        % Save masked activation map
        h.fname = [fileActMap(1:end-7), '_m4.nii.gz'];
        h.dt(1) = 2;
        spm_write_vol(h,v);
        
    end
    
end
