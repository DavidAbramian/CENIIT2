
% Parameters
% voxel_threshold = 3.1;
% clustersize_threshold = 10;

% voxel_threshold = 4.8917;  % tinv(1 - 0.05/nVox, dof)
% clustersize_threshold = 10;

voxel_threshold = 5;  % tinv(1 - 0.05/nVox, dof)
clustersize_threshold = 0;

clustersize_threshold_T1 = clustersize_threshold * 4 * 4 * 4; % T1 volume is 1 x 1 x 1 mm, fMRI volume is 4 x 4 x 4 mm
thresh_str = ['z', num2str(voxel_threshold),'_c', num2str(clustersize_threshold)];

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
motion_params = 1;

% for s = 1:nSubjects
for s = 6
    
    subject_id = dirSubj(s).name;
    
    dirOut = fullfile(dataPath, subject_id, 'fmri_stats', thresh_str);
    if ~exist(dirOut, 'dir')
        mkdir(dirOut);
    end
    
    %% Save activation maps in cube volumes
    
    % Load brain activity maps and tumor mask from NifTI files
    for t = 1:length(tasks)
        
        task = tasks{t};
        
        % Skip missing tasks
        if ~exist(fullfile(dataPath, subject_id, 'fmri', [task, '.nii.gz']), 'file')
            continue
        end
        
        for fwhm = FWHMs
        
            for m = motion_params

                switch m
                    case 1
                        mot_str = 'std-mot';
                    case 2
                        mot_str = 'ext-mot';
                end

                for c = 1:length(contrasts{t})

                    fprintf('subj: %s, task: %s, fwhm: %i, mot: %s, cont: %i\n', subject_id, task, fwhm, mot_str, c)

                    contrast = char(contrasts{t}(c));
                    fNameContrast = [contrast, '_', num2str(fwhm), 'mm_', mot_str, '.nii.gz'];
                    fileContrast = fullfile(dataPath, subject_id, 'fmri_stats', 'zmaps', fNameContrast);
                    [h,v] = ml_load_nifti(fileContrast);

                    % Z-value threshold
                    v = double(v > voxel_threshold);

                    % Cluster size threshold, remove small clusters
                    cc = bwconncomp(v);
                    for comp = 1:cc.NumObjects
                        if length(cc.PixelIdxList{comp}) < clustersize_threshold_T1
                            v(cc.PixelIdxList{comp}) = 0;
                        end
                    end

                    % Save as new nifti file
                    fNameOut = [replace(contrast, '_', '-'), '_', num2str(fwhm), 'mm.nii'];
                    fileOut = fullfile(dirOut, fNameOut);
                    h.fname = fileOut;
                    h.dt(1) = 2;
                    spm_write_vol(h,v);
                    system(['pigz -f ', fileOut]);

                end

            end

        end
        
    end

    
end
