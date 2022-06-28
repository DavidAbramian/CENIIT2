#!/bin/bash

# Work order

# 1. Analyze fMRI data (analyze_all_subjects.sh)
# 2. Re-scale anatomical T1w volume to 1 mm isotropic, transform activity maps to this new T1 volume (transform_motor.sh)
# 3. Fit nifti volumes into 256 cubes, and then convert T1w brain volume to DICOM format (one DICOM file per slice) (nifti_to_dicom.m)
# 4. Run this script to Create RT struct files for tumor mask and brain activity maps

z=$1  # t-threshold
c=$2  # cluster extent threshold
f=$3  # smoothign FWHM
m=$4  # mot

# for subject in [0-9][0-9][0-9][0-9][0-9] ; do
for subject in 17904 ; do

	echo $subject

	dirRT=$subject/dicom/${subject}_dicom_RTSTRUCT/
	if [[ ! -d $dirRT ]] ; then
		mkdir $dirRT
	fi

	dirT1DICOM=$subject/dicom/${subject}_dicom/
	dirAnat=$subject/anat/
	dirfMRI=$subject/fmri_stats/z${z}_c${c}/

	# Tumor mask
	if [[ ! -e $dirRT/${subject}_tumor_mask_new_RTSTRUCT.dcm ]]; then
		echo ${subject}_tumor_mask_new
		python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask_new.nii.gz $dirT1DICOM $dirRT ${subject}_tumor_mask_new Red > /dev/null
	fi

	# # Tumor mask
	# if [[ ! -e $dirRT/${subject}_tumor_mask_RTSTRUCT.dcm ]]; then
	# 	echo ${subject}_tumor_mask
	# 	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask_iso.nii.gz $dirT1DICOM $dirRT ${subject}_tumor_mask Red > /dev/null
	# fi
	#
	# if [[ -e $subject/fmri/motor.nii.gz ]]; then
	#
	# 	# Finger
	# 	echo ${subject}_motor_finger_${f}mm_${m}-mot_z${z}_c${c}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/motor_finger_${f}mm_${m}-mot_z${z}_c${c}.nii.gz $dirT1DICOM $dirRT ${subject}_motor_finger_${f}mm_${m}-mot_z${z}_c${c} Green > /dev/null
	#
	# 	# Foot
	# 	echo ${subject}_motor_foot_${f}mm_${m}-mot_z${z}_c${c}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/motor_foot_${f}mm_${m}-mot_z${z}_c${c}.nii.gz $dirT1DICOM $dirRT ${subject}_motor_foot_${f}mm_${m}-mot_z${z}_c${c} Blue > /dev/null
	#
	# 	# Lips
	# 	echo ${subject}_motor_lips_${f}mm_${m}-mot_z${z}_c${c}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/motor_lips_${f}mm_${m}-mot_z${z}_c${c}.nii.gz $dirT1DICOM $dirRT ${subject}_motor_lips_${f}mm_${m}-mot_z${z}_c${c} Yellow > /dev/null
	#
	# fi
	#
	# if [[ -e $subject/fmri/verb.nii.gz ]]; then
	#
	# 	# Verb generation
	# 	echo ${subject}_verb_generation_${f}mm_${m}-mot_z${z}_c${c}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/verb_generation_${f}mm_${m}-mot_z${z}_c${c}.nii.gz $dirT1DICOM $dirRT ${subject}_verb_generation_${f}mm_${m}-mot_z${z}_c${c} Cyan > /dev/null
	#
	# fi
	#
	# if [[ -e $subject/fmri/word.nii.gz ]]; then
	#
	# 	# Word repetition
	# 	echo ${subject}_word_repetition_${f}mm_${m}-mot_z${z}_c${c}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/word_repetition_${f}mm_${m}-mot_z${z}_c${c}.nii.gz $dirT1DICOM $dirRT ${subject}_word_repetition_${f}mm_${m}-mot_z${z}_c${c} Magenta > /dev/null
	#
	# fi

done
