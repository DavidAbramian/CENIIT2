% For every subject: transform T1.nii to 1mm isotropic. Transform z-maps  
% into same 1mm isotropic space. Conserves correct matrices so it is 
% possible to overlay activation maps on top of structural images.


% Find subject folders
dataPath = '/flush/davab27/CENIIT/data';
dirSubj = dir(dataPath);
names = string({dirSubj(:).name}');
ISubjects = cellfun(@(x)length(x)==5, regexp(names, '[0-9]'));  % Check if the file/folder name consists of 5 numbers
dirSubj = dirSubj(ISubjects);

nSubjects = length(dirSubj);

for s = 1:nSubjects
    
    subject_id = dirSubj(s).name;
    
    %% Convert T1 into 1mm isotropic
    
    fprintf('subj: %s, T1\n', subject_id)

    fileT1 = fullfile(dataPath, subject_id, 'anat', 'T1.nii.gz');
    [h,v] = ml_load_nifti(fileT1);
    
    % Resize T1 to 1mm isotropic
    v2 = imresize3(v, [256,203,256], 'linear');
    
    % Insert resized T1 volume into 256 cube volume
    vOut = zeros(256,256,256);
    vOut(:,1:203,:) = v2;
    
    % Adjust matrix to new voxel size
    mat = h.mat(1:3,1:3);
    [U,S,V] = svd(mat);
    mat2 = U*eye(3)*V';
    
    % Define new header
    hOut = h;
    hOut.mat(1:3,1:3) = mat2;
    hOut.dt(1) = 16;
    hOut.dim = size(vOut);
    fileOut = fullfile(dataPath, subject_id, 'anat', 'T1_iso.nii');
    hOut.fname = fileOut;
    
    % Save new T1 volume and gzip it
    spm_write_vol(hOut, vOut);
    system(['pigz -f ', fileOut]);
    
    %% Convert T1_brain into 1mm isotropic
    
    fprintf('subj: %s, T1_brain\n', subject_id)

    fileT1 = fullfile(dataPath, subject_id, 'anat', 'T1_brain.nii.gz');
    [h,v] = ml_load_nifti(fileT1);
    
    % Resize T1 to 1mm isotropic
    v2 = imresize3(v, [256,203,256], 'linear');
    
    % Insert resized T1 volume into 256 cube volume
    vOut = zeros(256,256,256);
    vOut(:,1:203,:) = v2;
    
    % Adjust matrix to new voxel size
    mat = h.mat(1:3,1:3);
    [U,S,V] = svd(mat);
    mat2 = U*eye(3)*V';
    
    % Define new header
    hOut = h;
    hOut.mat(1:3,1:3) = mat2;
    hOut.dt(1) = 16;
    hOut.dim = size(vOut);
    fileOut = fullfile(dataPath, subject_id, 'anat', 'T1_brain_iso.nii');
    hOut.fname = fileOut;
    
    % Save new T1 volume and gzip it
    spm_write_vol(hOut, vOut);
    system(['pigz -f ', fileOut]);
    
    %% Convert tumor mask into 1mm isotropic

    fprintf('subj: %s, tumor mask\n', subject_id)

    fileTumor = fullfile(dataPath, subject_id, 'anat', 'tumor_mask.nii.gz');
    [h,v] = ml_load_nifti(fileTumor);
    
    % Resize T1 to 1mm isotropic
    v2 = imresize3(v, [256,203,256], 'linear');
    
    % Insert resized T1 volume into 256 cube volume
    vOut = zeros(256,256,256);
    vOut(:,1:203,:) = v2;
    
    % Adjust matrix to new voxel size
    mat = h.mat(1:3,1:3);
    [U,S,V] = svd(mat);
    mat2 = U*eye(3)*V';
    
    % Define new header
    hOut = h;
    hOut.mat(1:3,1:3) = mat2;
    hOut.dt(1) = 16;
    hOut.dim = size(vOut);
    fileOut = fullfile(dataPath, subject_id, 'anat', 'tumor_mask_iso.nii');
    hOut.fname = fileOut;
    
    % Save new T1 volume and gzip it
    spm_write_vol(hOut, vOut);
    system(['pigz -f ', fileOut]);
    
end
