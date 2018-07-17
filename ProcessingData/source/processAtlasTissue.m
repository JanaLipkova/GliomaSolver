%===================================
%
%  Create brain mask from atlas tissue
% -----------------------------------
% 1) combine tissue segmenatation to get binary brain mask
% 2) use Gaussian filter to convert the brain mask to soft mask -> used for
% pff
% 3) save both masks to nii
%===================================


 function processAtlasTissue

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%path to data
path2atlas = '../../../DataForPaper/Atlas/';

wm  = MRIread([path2atlas,'wm_high.nii']);
gm  = MRIread([path2atlas,'gm_high.nii']);
csf = MRIread([path2atlas,'csf_high.nii']);

%threshold to get binary mask
mask = wm.vol + gm.vol + csf.vol;
mask(mask(:)>0.1) = 1;
mask(mask(:)<0.1) = 0;

% Apply gaussian smootheing filter to get soft mask -> used as pff
% Gaussin filter with param. 1 lead to smoothening around 5 gp, with values
% close to 0.5 at the border
pff = imgaussfilt(mask,1.0);

%save to nii
wm.vol = mask;
gm.vol = pff;

MRIwrite(wm,[path2atlas,'mask.nii']);
MRIwrite(gm,[path2atlas,'pff.nii']);
