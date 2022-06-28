#!/bin/bash

systemctl start docker

docker run -ti --rm -e TZ=Europe/Stockholm \
	-v /flush/davab27/br-tum-bids/:/data:ro \
	-v /flush/davab27/br-tum-out-fs-syn:/out \
	-v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
	poldracklab/fmriprep:latest /data /out participant \
	--fs-license-file /license.txt \
	--use-syn-sdc \
	|& tee /flush/davab27/fmriprep-log-fs-syn.txt
