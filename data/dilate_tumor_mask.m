
% Parameters
z = 5;
c = 0;

margin = 2;

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

a = -ones(20,6);
b = -ones(20,6);

% for s = 1:nSubjects
for s = 1
    
    subject_id = dirSubj(s).name;
    
    % Find tumor voxel indices
    tumorFile = fullfile(dataPath, subject_id, 'anat', 'tumor_mask_new.nii.gz');
    [h,v] = ml_load_nifti(tumorFile);
    dim = h.dim;
    I_tumor = find(v);
    
    SE = strel('sphere', margin);
    
    v2 = convn(v, SE.Neighborhood, 'same') > 0;
    
%     dirActivationMaps = fullfile(dataPath, subject_id, 'fmri_stats', ['z', num2str(z), '_c', num2str(c)]);
%     files = dir(fullfile(dirActivationMaps, '*0.nii.gz'));
%     
%     %% Mask activation maps with brain tumor mask
%     
%     for i = 1:length(files)
%         
%         fprintf('%s, %s \n', subject_id, files(i).name)
%         
%         % Load activation map
%         fileActMap = fullfile(dirActivationMaps, files(i).name);
%         [h,v] = ml_load_nifti(fileActMap);
%         
% %         fprintf('%i \n', nnz(v(I_tumor)))
%         fprintf('%i \n', nnz(v))
%         
%         a(i,s) = nnz(v);
%         b(i,s) = nnz(v(I_tumor));
% 
% %         % Mask out tumor voxels
% %         v(I_tumor) = 0;
% %         
% %         % Save masked activation map
% %         h.fname = [fileActMap(1:end-7), '_masked_new.nii.gz'];
% %         h.dt(1) = 2;
% %         spm_write_vol(h,v);
%         
%     end
    
end
