#!/bin/bash

systemctl start docker

docker run -ti --rm -e TZ=Europe/Stockholm \
	-v /flush/davab27/br-tum/br-tum-bids/:/data:ro \
	-v /flush/davab27/br-tum/br-tum-out-syn-MNI-lesion:/out \
	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
	-v /flush/davab27/br-tum/work/br-tum-work-syn-MNI-lesion:/scratch \
	poldracklab/fmriprep:latest /data /out participant \
	--fs-no-reconall \
	--fs-license-file /license.txt \
	--nthreads 12 \
	--omp-nthreads 12 \
	--force-syn \
	-w /scratch \
	|& tee /flush/davab27/br-tum/fmriprep-log-syn-lesion.txt

#docker run -ti --rm -e TZ=Europe/Stockholm \
#	-v /flush/davab27/br-tum-bids/:/data:ro \
#	-v /flush/davab27/br-tum-out-syn:/out \
#	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
#	-v /flush/davab27/br-tum-bids-syn-work:/scratch \
#	poldracklab/fmriprep:latest /data /out participant \
#	--fs-no-reconall \
#	--fs-license-file /license.txt \
#	--nthreads 12 \
#	--omp-nthreads 12 \
#	--use-syn-sdc \
#	--output-space /data/sub-01/anat/sub-01_T1w.nii.gz \
#	-w /scratch \
#	|& tee /flush/davab27/fmriprep-log-syn.txt
