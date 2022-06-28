#!/bin/bash

systemctl start docker

docker run -ti --rm -e TZ=Europe/Stockholm \
	-v /flush/davab27/br-tum-bids/:/data:ro \
	-v /flush/davab27/br-tum-out-fs:/out \
	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
	poldracklab/fmriprep:latest /data /out participant \
	--nthreads 12 \
	--omp-nthreads 12 \
	--fs-license-file /license.txt \
	|& tee /flush/davab27/fmriprep-log-fs.txt
