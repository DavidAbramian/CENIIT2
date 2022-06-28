#!/bin/bash

systemctl start docker

docker run -ti --rm -e TZ=Europe/Stockholm \
	-v /flush/davab27/HCP_1200/bids/data/:/data:ro \
	-v /flush/davab27/HCP_1200/bids/out:/out \
	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
	-v /flush/davab27/HCP_1200//bids/work/:/scratch \
	poldracklab/fmriprep:latest /data /out participant \
	--fs-no-reconall \
	--fs-license-file /license.txt \
	--nthreads 8 \
	--omp-nthreads 8 \
	--output-space T1w \
	-w /scratch \
	|& tee /flush/davab27/HCP_1200/bids/fmriprep-log.txt
