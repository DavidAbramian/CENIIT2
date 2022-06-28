#!/bin/bash

# Work order

# 1. Analyze fMRI data by running this script
# 2. Re-scale anatomical T1w volume to 1 mm isotropic, transform activity maps to this new T1 volume (transform_motor.sh)
# 3. Fit nifti volumes into 256 cubes, and then convert T1w brain volume to DICOM format (one DICOM file per slice) (nifti_to_dicom.m)
# 4. Create RT struct files for tumor mask and activity maps (create_RTSTRUCTS.sh)


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

dirData="/flush/davab27/CENIIT/data"

for subject in [0-9][0-9][0-9][0-9][0-9] ; do
# for subject in 17904 ; do

	mkdir -p $subject/fmri_out

	for task in {motor,verb,word}; do

		if [[ -e $subject/fmri/$task.nii.gz ]]; then

			dirTemplate="$dirData/fmri/$task"

			for motion in 1 2 ; do

				case $motion in
					1)
						motion_str=std-mot
						confounds_file=${task}_mot_par.txt
						;;
					2)
						motion_str=ext-mot
						confounds_file=${task}_mot_par_der.txt
						;;
				esac

				for fwhm in 4 6 ; do

					echo $subject, $task, $motion_str, $fwhm

					# Copy template design
					design=$subject/fmri_out/${task}_f${fwhm}_${motion_str}.fsf
					cp $dirTemplate/${task}_template.fsf $design

					# Change subject name
					sed -i "s/<subject>/$subject/g" $design

					# Change smoothing output
					sed -i "s/<fwhm>/$fwhm/g" $design

					# Change motion output
					sed -i "s/<motion_str>/$motion_str/g" $design
					# sed -i "s/<motion>/$motion/g" $design
					sed -i "s/<confounds_file>/$confounds_file/g" $design

					# Run analyses in parallel
					feat $design > /dev/null &

					check_jobs

				done

			done

		fi

	done

	# Remove design files (they are available inside each feat folder)
	rm $subject/fmri_out/*.fsf

done

# Wait for remaining jobs to finish
wait_finish
