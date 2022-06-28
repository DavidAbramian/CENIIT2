#!/bin/bash

z=$1  # t-threshold
c=$2  # cluster extent threshold
f=$3  # smoothign FWHM
m=$4

# for subject in [0-9][0-9][0-9][0-9][0-9] ; do
# for subject in 17904 ; do
# for subject in 18582 ; do
# for subject in 18975 ; do
# for subject in 19015 ; do
for subject in 19849 ; do

	echo $subject

	dirRT=$subject/dicom/${subject}_dicom_RTSTRUCT/
	if [[ ! -d $dirRT ]] ; then
		mkdir $dirRT
	fi

	dirT1DICOM=$subject/dicom/${subject}_dicom/
	dirAnat=$subject/anat/
	dirfMRI=$subject/fmri_stats/z${z}_c${c}/

	# Tumor mask
	if [[ ! -e $dirRT/${subject}_tumor_mask_RTSTRUCT.dcm ]]; then
		echo ${subject}_tumor_mask
		python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask_new.nii.gz $dirT1DICOM $dirRT ${subject}_tumor_mask red > /dev/null
	fi

	# # Tumor mask
	# if [[ ! -e $dirRT/${subject}_tumor_mask_RTSTRUCT.dcm ]]; then
	# 	echo ${subject}_tumor_mask
	# 	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask_iso.nii.gz $dirT1DICOM $dirRT ${subject}_tumor_mask Red > /dev/null
	# fi

	if [[ -e $subject/fmri/motor.nii.gz ]]; then

		# Finger
		echo ${subject}_motor-finger_${f}mm${m}
		python3.7 convert_to_rtstruct.py $dirfMRI/motor-finger_${f}mm${m}.nii.gz $dirT1DICOM $dirRT ${subject}_motor-finger_${f}mm${m} lime > /dev/null

		# Foot
		echo ${subject}_motor-foot_${f}mm${m}
		python3.7 convert_to_rtstruct.py $dirfMRI/motor-foot_${f}mm${m}.nii.gz $dirT1DICOM $dirRT ${subject}_motor-foot_${f}mm${m} orange > /dev/null

		# Lips
		echo ${subject}_motor-lips_${f}mm${m}
		python3.7 convert_to_rtstruct.py $dirfMRI/motor-lips_${f}mm${m}.nii.gz $dirT1DICOM $dirRT ${subject}_motor-lips_${f}mm${m} cyan > /dev/null

	fi

	if [[ -e $subject/fmri/verb.nii.gz ]]; then

		# Verb generation
		echo ${subject}_verb-generation_${f}mm${m}
		python3.7 convert_to_rtstruct.py $dirfMRI/verb-generation_${f}mm${m}.nii.gz $dirT1DICOM $dirRT ${subject}_verb-generation_${f}mm${m} magenta > /dev/null

	fi

	# if [[ -e $subject/fmri/word.nii.gz ]]; then

	# 	# Word repetition
	# 	echo ${subject}_word-repetition_${f}mm${m}
	# 	python3.7 convert_to_rtstruct.py $dirfMRI/word-repetition_${f}mm${m}.nii.gz $dirT1DICOM $dirRT ${subject}_word-repetition_${f}mm${m} dark-violet > /dev/null

	# fi

done
