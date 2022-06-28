#!/bin/bash

maxJobs=6 # Max number of parallel jobs
check_jobs () {
	while [[ $(jobs -r | wc -l) -ge $maxJobs ]]; do
		sleep 1
	done
}
wait_finish () {
	while [[ $(jobs -r | wc -l) -gt 0 ]]; do
	  sleep 1
	done
}

# Function that carries out all the processing for a single subject
process_subj () {

  local subj=$1

  echo $subj

  dirSubj="sub-$subj"
  dirAnat=$dirSubj/anat
  dirFunc=$dirSubj/func

  # Create forlders for data
  mkdir -p $dirSubj/{anat,func}

  # Extract and rename T1
  local fileT1=$(ls $subj/*COR_3D*.zip)
  unzip -q $fileT1 -d $dirAnat
  fileT1=$dirAnat/sub-${subj}_T1w.nii.gz
  gzip -c $dirAnat/*COR_3D*/anon_*.nii > $fileT1
  rm -r $dirAnat/*COR_3D*

  # Fix T1 orientation
  fslswapdim $fileT1 x -z y $fileT1

  # Extract and rename T2
  local fileT2=$(ls $subj/*Axial_T2*.zip)
  unzip -q $fileT2 -d $dirAnat
  fileT2=$dirAnat/sub-${subj}_T2w.nii.gz
  gzip -c $dirAnat/*Axial_T2*/anon_*.nii > $fileT2
  rm -r $dirAnat/*Axial_T2*

  # Extract and rename tumor segmentation
  local fileTumor=$(ls $subj/tissue_classes.zip)
  unzip -q $fileTumor -d $dirAnat
  fileTumor=$dirAnat/sub-${subj}_T1w_label-lesion_roi.nii.gz
  gzip -c $dirAnat/c3anon_*.nii > $fileTumor
  rm $dirAnat/*anon*.nii

  # Fix tumor orientation
  fslswapdim $fileTumor x -z y $fileTumor

  # # Create brainmask
  # bet $subj/anat/T1.nii.gz $subj/anat/T1_brain.nii.gz -R -f 0.5

  # Extract and process motor fMRI
  local fileMotor=$(ls $subj/*finger*.zip 2> /dev/null)

  if [[ -n $fileMotor ]]; then
    unzip -q $fileMotor -d $dirFunc
    fileMotor=${dirFunc}/sub-${subj}_task-motor_bold.nii.gz
    fslmerge -tr $fileMotor $dirFunc/*finger*/*.nii 2.5
    fslroi $fileMotor $fileMotor 4 180
    rm -r $dirFunc/*finger*/
  fi

  # Extract and process verb generation fMRI
  local fileVerb=$(ls $subj/*silent*.zip 2> /dev/null)

  if [[ -n $fileVerb ]]; then
    unzip -q $fileVerb -d $dirFunc
    fileVerb=${dirFunc}/sub-${subj}_task-verb_bold.nii.gz
    fslmerge -tr $fileVerb $dirFunc/*silent*/*.nii 2.5
    fslroi $fileVerb $fileVerb 4 169
    rm -r $dirFunc/*silent*/
  fi

  # Extract and process word repetition fMRI
  local fileWord=$(ls $subj/*word*.zip 2> /dev/null)

  if [[ -n $fileWord ]]; then
    unzip -q $fileWord -d $dirFunc
    fileWord=${dirFunc}/sub-${subj}_task-word_bold.nii.gz
    fslmerge -tr $fileWord $dirFunc/*word*/*.nii 2.5
    fslroi $fileWord $fileWord 4 72
    rm -r $dirFunc/*word*/
  fi

  # Delete zip files
  rm -r $subj

}


# Extract subject data
# echo "untaring files"
# for file in orig_files/*.tar; do
#
#   tar -xf $file
#
# done

# echo ""

dir_fmriprep=./fmriprep
dir_fmriprep_data=$dir_fmriprep/data
dir_fmriprep_out=$dir_fmriprep/out-T1w-SyN-noFS

dir_out=./data

if [[ ! -d $dir_out ]]; then
  mkdir $dir_out
fi

for dir_fmriprep_data_subj in ${dir_fmriprep_data}/sub-[0-9][0-9][0-9][0-9][0-9] ; do

  subj=${dir_fmriprep_data_subj##*-}
  echo $subj

  dir_fmriprep_out_subj=$dir_fmriprep_out/sub-$subj
  dir_out_subj=${dir_out}/$subj

  # Create output folders
  mkdir -p $dir_out_subj/{anat,fmri}

  # Copy T1w volume
  cp $dir_fmriprep_out_subj/anat/sub-*_desc-preproc_T1w.nii.gz ${dir_out_subj}/anat/T1.nii.gz

  # Copy tumor mask
  # TODO: This is wrong mask. Replace with clean mask.
  cp $dir_fmriprep_data_subj/anat/sub*_T1w_label-lesion_roi.nii.gz ${dir_out_subj}/anat/tumor_mask.nii.gz

  # Copy brain mask
  cp $dir_fmriprep_out_subj/anat/sub*_desc-brain_mask.nii.gz ${dir_out_subj}/anat/brain_mask.nii.gz

  # Apply brain mask to T1w volume
  fslmaths ${dir_out_subj}/anat/T1.nii.gz -mas ${dir_out_subj}/anat/brain_mask.nii.gz ${dir_out_subj}/anat/T1_brain.nii.gz

  # Copy fMRI data
  for task in {motor,verb,word}; do
    cp $dir_fmriprep_out_subj/func/sub-*_task-${task}_space-T1w_desc-preproc_bold.nii.gz ${dir_out_subj}/fmri/${task}.nii.gz  2> /dev/null
    cp $dir_fmriprep_out_subj/func/sub-*_task-${task}_space-T1w_desc-brain_mask.nii.gz ${dir_out_subj}/fmri/${task}_mask.nii.gz  2> /dev/null
    cp $dir_fmriprep_out_subj/func/sub-*_task-${task}_desc-confounds_timeseries.tsv ${dir_out_subj}/fmri/${task}_confounds.tsv  2> /dev/null

    if [[ -e ${dir_out_subj}/fmri/${task}_confounds.tsv ]]; then

      # Find where relevant confounds begin
      i=$(awk -v b='trans_x' '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}' ${dir_out_subj}/fmri/${task}_confounds.tsv)

      # Motion parameters
      awk -v i=$i '{ print $i, $(i+4), $(i+8), $(i+12), $(i+16), $(i+20)}' ${dir_out_subj}/fmri/${task}_confounds.tsv > ${dir_out_subj}/fmri/${task}_mot_par.txt

      # Motion parameters + derivatives
      awk -v i=$i '{ print $i, $(i+1), $(i+4), $(i+5), $(i+8), $(i+9), $(i+12), $(i+13), $(i+16), $(i+17), $(i+20), $(i+21)}' ${dir_out_subj}/fmri/${task}_confounds.tsv > ${dir_out_subj}/fmri/${task}_mot_par_der.txt
    fi

  done

  # Start a new job as soon as another one finishes
  # check_jobs

done

# Wait for remaining jobs to finish
wait_finish
