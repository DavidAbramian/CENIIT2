[h,v] = ml_load_nifti('17904/fmri_stats/z5_c0/motor-foot_6mm.nii.gz');

v = logical(v);
CC = bwconncomp(v);

%%

nVox = cellfun(@length, CC.PixelIdxList)';

vOut = false(h.dim);

% vOut(CC.PixelIdxList{3}) = true;
% vOut(CC.PixelIdxList{5}) = true;
vOut(CC.PixelIdxList{2}) = true;

h.fname = 'test3.nii';
h.dt(1) = 2;

spm_write_vol(h,vOut);