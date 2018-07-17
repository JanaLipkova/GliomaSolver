%===================================
%
%   Convert folder with nii to dat
% -----------------------------------
%  INPUT:  path dat folder
%  OuTPUT: folder with nii
%
%===================================

% function convert_folder_content_nii2dat

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

inputPath =  '../../../DataForPaper/Atlas/';
outputPath = inputPath;

if( exist(outputPath,'dir') == 0 )
    mkdir(outputPath)
end;

bRotate = 0;

%4) i)  get names of files in the input folder
%   ii) convert them to nii
%   iii) save nii to output folder
files = dir(inputPath);
filesNames = {files.name};
filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));

for i = 1:length(filesNames)
    
    inFilename  = filesNames{i}
    outFilename = [inFilename(1:end-3),'dat'];
    
    nii_data = MRIread([inputPath,inFilename]);
    data = nii_data.vol;    
    
    if(bRotate)
        data = rotate90_3D(data,1);
        data = rotate90_3D(data,1);
        data = rotate90_3D(data,1);
    end;
   
    typeId = 1; % 0 double, 1 float
    writeMatrix(data,[outputPath,outFilename],typeId);
end;

