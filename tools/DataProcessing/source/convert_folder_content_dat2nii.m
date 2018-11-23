%===================================
%
%   Convert folder with dat to nii
% -----------------------------------
%  INPUT:  path dat folder
%  OuTPUT: folder with nii
%
%===================================

% function convert_folder_content_dat2nii

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

inputPath  = '../../18/Results/Propagation/';
outputPath = [inputPath(1:end-1),'_nii/'];


bResize = 1;
bRotate = 0;


if( exist(outputPath,'dir') == 0 )
    mkdir(outputPath)
end;

%3) Read in any nifty volume to learn nifty header
%  nii_data = MRIread('/Volumes/FileStorage/GLIOMA/DataForPaper/01/Results/ResultsNii/P01_T1w_out.nii');
 nii_data = MRIread('/Volumes/FileStorage/GLIOMA/DataForPaper/18/RawData/MRI/FLAIR.nii');


%4) i)  get names of files in the input folder
%   ii) convert them to nii
%   iii) save nii to output folder

files = dir(inputPath);
filesNames = {files.name};
filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));

for i = 1:length(filesNames)
    
    inFilename  = filesNames{i};
    outFilename = [inFilename(1:end-3),'nii'];
    
    if( inFilename(end-3:end) == '.dat')
        inFilename
        datVolume = loadMatrix([inputPath,inFilename]);
        
        if(bRotate)
            datVolume = rotate90_3D(datVolume,1);
%            datVolume = rotate90_3D(datVolume,1);
%            datVolume = rotate90_3D(datVolume,1);
        end;
        
        if(bResize)
            [Nx,Ny,Nz] = size(nii_data.vol);
            datVolume = resize_data(datVolume,Nx,Ny,Nz);
        end;
        
         nii_data.vol = datVolume;
        MRIwrite(nii_data, [outputPath,outFilename]);
    end
end;
