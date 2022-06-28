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

	local dir_fmriprep_data_subj=${dir_fmriprep_data}/sub-$subj
	local dir_fmriprep_out_subj=$dir_fmriprep_out/sub-$subj
  local dir_out_subj=${dir_out}/$subj

  # Create output folders
  mkdir -p $dir_out_subj/{anat,fmri}

  # Copy T1w volume
  cp $dir_fmriprep_out_subj/anat/sub-${subj}_desc-preproc_T1w.nii.gz ${dir_out_subj}/anat/T1.nii.gz
	fslorient -swaporient ${dir_out_subj}/anat/T1.nii.gz
	fslswapdim ${dir_out_subj}/anat/T1.nii.gz -x y z ${dir_out_subj}/anat/T1.nii.gz > /dev/null

  # Copy tumor mask
  # TODO: This is wrong mask. Replace with clean mask.
  # cp $dir_fmriprep_data_subj/anat/sub*_T1w_label-lesion_roi.nii.gz ${dir_out_subj}/anat/tumor_mask.nii.gz
	# fslswapdim ${dir_out_subj}/anat/tumor_mask.nii.gz -x y z ${dir_out_subj}/anat/tumor_mask.nii.gz > /dev/null

	# Correct mask
	cp _other/Edinburgh/$subj/tumor_mask_clean.nii.gz ${dir_out_subj}/anat/tumor_mask.nii.gz
	fslswapdim ${dir_out_subj}/anat/tumor_mask.nii.gz x -z y ${dir_out_subj}/anat/tumor_mask.nii.gz > /dev/null

  # Copy brain mask
  cp $dir_fmriprep_out_subj/anat/sub*_desc-brain_mask.nii.gz ${dir_out_subj}/anat/brain_mask.nii.gz
	fslorient -swaporient ${dir_out_subj}/anat/brain_mask.nii.gz
	fslswapdim ${dir_out_subj}/anat/brain_mask.nii.gz -x y z ${dir_out_subj}/anat/brain_mask.nii.gz > /dev/null

  # Apply brain mask to T1w volume
  fslmaths ${dir_out_subj}/anat/T1.nii.gz -mas ${dir_out_subj}/anat/brain_mask.nii.gz ${dir_out_subj}/anat/T1_brain.nii.gz

  # Copy fMRI data
  for task in {motor,verb,word}; do

		if [[ -e $dir_fmriprep_out_subj/func/sub-${subj}_task-${task}_space-T1w_desc-preproc_bold.nii.gz ]]; then

			# Copy task data
			cp $dir_fmriprep_out_subj/func/sub-${subj}_task-${task}_space-T1w_desc-preproc_bold.nii.gz ${dir_out_subj}/fmri/${task}.nii.gz
			fslorient -swaporient ${dir_out_subj}/fmri/${task}.nii.gz
			fslswapdim ${dir_out_subj}/fmri/${task}.nii.gz -x y z ${dir_out_subj}/fmri/${task}.nii.gz > /dev/null

			# Copy brain mask
	    cp $dir_fmriprep_out_subj/func/sub-${subj}_task-${task}_space-T1w_desc-brain_mask.nii.gz ${dir_out_subj}/fmri/${task}_mask.nii.gz
			fslorient -swaporient ${dir_out_subj}/fmri/${task}_mask.nii.gz
			fslswapdim ${dir_out_subj}/fmri/${task}_mask.nii.gz -x y z ${dir_out_subj}/fmri/${task}_mask.nii.gz > /dev/null

			# Copy confounds table
	    cp $dir_fmriprep_out_subj/func/sub-${subj}_task-${task}_desc-confounds_timeseries.tsv ${dir_out_subj}/fmri/${task}_confounds.tsv

			# Extract motion parameters
      # Find where relevant confounds begin
      local i=$(awk -v b='trans_x' '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}' ${dir_out_subj}/fmri/${task}_confounds.tsv)

      # Motion parameters
      awk -v i=$i '{ print $i, $(i+4), $(i+8), $(i+12), $(i+16), $(i+20)}' ${dir_out_subj}/fmri/${task}_confounds.tsv > ${dir_out_subj}/fmri/${task}_mot_par.txt

      # Motion parameters + derivatives
      awk -v i=$i '{ print $i, $(i+1), $(i+4), $(i+5), $(i+8), $(i+9), $(i+12), $(i+13), $(i+16), $(i+17), $(i+20), $(i+21)}' ${dir_out_subj}/fmri/${task}_confounds.tsv > ${dir_out_subj}/fmri/${task}_mot_par_der.txt

		fi

  done

}


dir_fmriprep=./fmriprep
dir_fmriprep_data=$dir_fmriprep/data
dir_fmriprep_out=$dir_fmriprep/out-T1w-SyN-noFS

dir_out=./data

if [[ ! -d $dir_out ]]; then
  mkdir $dir_out
fi

for dir_fmriprep_data_subj in ${dir_fmriprep_data}/sub-[0-9][0-9][0-9][0-9][0-9] ; do

  subj=${dir_fmriprep_data_subj##*-}

	# Process subject
	process_subj $subj &

  # Start a new job as soon as another one finishes
  check_jobs

done

# Wait for remaining jobs to finish
wait_finish
