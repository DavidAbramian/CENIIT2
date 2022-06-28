
% Find subject folders
dataPath = '/flush/davab27/CENIIT/data';
dirSubj = dir(dataPath);
names = string({dirSubj(:).name}');
ISubjects = cellfun(@(x)length(x)==5, regexp(names, '[0-9]'));  % Check if the file/folder name consists of 5 numbers
dirSubj = dirSubj(ISubjects);

nSubjects = length(dirSubj);

% Load DICOM header
fileDICOM = 'template.dcm';
infoDICOM = dicominfo(fileDICOM);

for s = 1:nSubjects
% for s = 1
    
    subject_id = dirSubj(s).name;
    
    %% Convert T1 to DICOM
    
    % Create output folder for DICOM files
    dirOutDicom = fullfile(dataPath, subject_id, 'dicom', [subject_id, '_dicom']);
    if ~exist(dirOutDicom, 'dir')
        mkdir(dirOutDicom)
    end
        
    % Load T1w brain volume from NifTI file
    fileT1 = fullfile(dataPath, subject_id, 'anat', 'T1_iso.nii.gz');
    vol = niftiread(fileT1);
    
    % Split single 3D volume into one DICOM file per slice
    nSlices = size(vol,3);
    for i = 1:nSlices
        
        progresss(i, nSlices, 'Saving DICOM... ')
        
        % Take single slice. Flip for correct patient orientation
        %         slice = uint16(squeeze(vol(:,i,:))');
        
        slice = uint16(rot90(vol(:,:,i)));
        
        %imagesc(slice); colormap gray
        %pause(0.1)
        
        % Create DICOM header
        md = infoDICOM;
        md.ImagePositionPatient(1) = 0;
        md.ImagePositionPatient(2) = 0;
        md.ImagePositionPatient(3) = 0;
        md.Rows = size(slice,1);
        md.Columns = size(slice,2);
        
        % Not necessary?
        %     md.SpecificCharacterSet = 'ISO_IR 100';
        %     md.LargestImagePixelValue = max(slice(:));
        
        % These fields will be visible when loading data in GammaPlan
%         md.PatientName = [subject_id '_Test']; % Update when creating new files
%         md.PatientID = [subject_id '_Test']; % Update when creating new files
        md.PatientName = [subject_id]; % Update when creating new files
        md.PatientID = [subject_id]; % Update when creating new files
        %     md.SeriesNumber = 1002;  % To have several MR images for the same patient, not necessary ?
        md.SeriesDescription = 'T1w';
        md.MRAcquisitionType = '2D';
        
        md.InstanceNumber = i;  % Slice number
        
        % From standard: "specifies the x, y, and z coordinates of the upper left hand corner of the image"
        md.ImagePositionPatient(3) = md.ImagePositionPatient(3)+i;
        
        % Replacing field from the 3D DICOM to the 2D version
        %md.SliceLocation = md.SliceLocationVector(i);
        md.SliceLocation = i;
        md = rmfield(md, 'SliceLocationVector');
        
        md.PatientPosition = 'HFS';  % Head-first position
        
        md.ImageOrientationPatient(5) = 1;  % This is necessary for GammaPlan
        
        % Save slice as DICOM file
        fOut = fullfile(dirOutDicom, ['slice', num2str(i, '%.3i'), '.dcm']);
        x = dicomwrite(slice, fOut, md, 'ObjectType', 'MR Image Storage');
        
    end
    
end
