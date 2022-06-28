#!/bin/bash

z=$1  # t-threshold
c=$2  # cluster extent threshold
f=$3  # smoothign FWHM
m=$4

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
	echo red
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT red red > /dev/null

	echo tomato
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT tomato tomato > /dev/null

	echo orange
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT orange orange > /dev/null

	echo yellow-green
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT yellow-green yellow-green > /dev/null

	echo lime
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT lime lime > /dev/null

	echo cyan
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT cyan cyan > /dev/null

	echo turquoise
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT turquoise turquoise > /dev/null

	echo deep-blue-sky
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT deep-blue-sky deep-blue-sky > /dev/null

	echo blue
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT blue blue > /dev/null

	echo dark-violet
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT dark-violet dark-violet > /dev/null

	echo magenta
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT magenta magenta > /dev/null

	echo brown
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT brown brown > /dev/null

	echo tan
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT tan tan > /dev/null

	echo pink
	python3.7 convert_to_rtstruct.py $dirAnat/tumor_mask.nii.gz $dirT1DICOM $dirRT pink pink > /dev/null

done
