#!/bin/bash

systemctl start docker

docker run -ti --rm -e TZ=Europe/Stockholm \
	-v /flush/davab27/br-tum/br-tum-bids/:/data:ro \
	-v /flush/davab27/br-tum/br-tum-out-lesion:/out \
	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
	-v /flush/davab27/br-tum/work/br-tum-work-lesion:/scratch \
	poldracklab/fmriprep:latest /data /out participant \
	--fs-no-reconall \
	--fs-license-file /license.txt \
	--nthreads 12 \
	--omp-nthreads 12 \
	-w /scratch \
	|& tee /flush/davab27/br-tum/fmriprep-log-lesion.txt

#docker run -ti --rm -e TZ=Europe/Stockholm \
#	-v /flush/davab27/br-tum/br-tum-bids/:/data:ro \
#	-v /flush/davab27/br-tum/br-tum-out-lesion:/out \
#	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
#	-v /flush/davab27/br-tum/work/br-tum-work-lesion:/scratch \
#	poldracklab/fmriprep:latest /data /out participant \
#	--fs-no-reconall \
#	--fs-license-file /license.txt \
#	--nthreads 12 \
#	--omp-nthreads 12 \
#	--participant-label 11 \
#	--anat-only \
#	-w /scratch \
#	|& tee /flush/davab27/br-tum/fmriprep-log-lesion.txt
