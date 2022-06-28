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

% fMRI parameters
tasks = {'motor'; 'verb'; 'word'};
FWHMs = [4, 6];
contrasts = {["motor_finger"; "motor_foot"; "motor_lips"]; "verb_generation"; "word_repetition" };
motion_params = 1:2;

for s = 1:nSubjects
    
    subject_id = dirSubj(s).name;
        
    %% Convert fMRI to same space as T1
    
    % Get reference matrix from T1 volume
    fileT1 = fullfile(dataPath, subject_id, 'anat', 'T1_iso.nii.gz');
    hRef = ml_load_nifti(fileT1);
    matRef = hRef.mat;
    dimRef = hRef.dim;
    [x, y, z] = ndgrid(1:dimRef(1),1:dimRef(2),1:dimRef(3));
    
    % Create output directory
    dirOut = fullfile(dataPath, subject_id, 'fmri_stats', 'zmaps');
    if ~exist(dirOut, 'dir')
        mkdir(dirOut);
    end
        
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

                dirNameFMRI = [task, '_', num2str(fwhm), 'mm_', mot_str, '.feat'];
                dirFMRI = fullfile(dataPath, subject_id, 'fmri_out', dirNameFMRI, 'stats');

                for c = 1:length(contrasts{t})

                    fprintf('subj: %s, task: %s, fwhm: %i, mot: %s, cont: %i\n', subject_id, task, fwhm, mot_str, c)

                    fileFMRI = fullfile(dirFMRI, ['zstat', num2str(c), '.nii.gz']);
                    [h,v] = ml_load_nifti(fileFMRI);

                    % Resample z-maps
                    mat = h.mat;
                    A = affine3d((mat\matRef)');

                    [xTrans,yTrans,zTrans] = transformPointsForward(A, x, y, z);
                    vOut = interpn(v, xTrans, yTrans, zTrans, 'cubic', 0);

                    % Define new header
                    hOut = hRef;
                    hOut.dt(1) = 16;
                    fileOut = fullfile(dirOut, [char(contrasts{t}(c)),'_',num2str(fwhm),'mm_',mot_str,'.nii']);
                    hOut.fname = fileOut;

                    % Save new T1 volume and gzip it
                    spm_write_vol(hOut, vOut);
                    system(['pigz -f ', fileOut]);

                end

            end

        end
    end
    
end
