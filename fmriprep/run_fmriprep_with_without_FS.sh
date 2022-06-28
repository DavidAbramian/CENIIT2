
# T1w output space, SyN-SDC, no FreeSurfer
docker run --rm -it -v /flush/davab27/CENIIT/fmriprep/data:/data:ro \
  -v /flush/davab27/CENIIT/fmriprep/out-T1w-SyN-noFS:/out \
  -v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
  nipreps/fmriprep:21.0.0 \
  --random-seed 100 \
  --skull-strip-t1w auto \
  --output-spaces T1w \
  --use-syn-sdc \
  --fs-no-reconall \
  --fs-license-file /license.txt \
  /data /out participant | tee out-T1w-SyN-noFS.txt

# T1w output space + SyN-SDC
docker run --rm -it -v /flush/davab27/CENIIT/fmriprep/data:/data:ro \
  -v /flush/davab27/CENIIT/fmriprep/out-T1w-SyN:/out \
  -v /home/davab27/.licenses/freesurfer-license.txt:/license.txt \
  nipreps/fmriprep:21.0.0 \
  --random-seed 100 \
  --skull-strip-t1w auto \
  --output-spaces T1w \
  --use-syn-sdc \
  --fs-license-file /license.txt \
  /data /out participant | tee out-T1w-SyN.txt
