%===================================
%
%   Convert folder with nii to dat
% -----------------------------------
%  INPUT:  path dat folder
%  OuTPUT: folder with nii
%
%===================================

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

% FILL IN HERE: Provide full path to data or paht relative to this fodler
inputPath =  '../../Examples/PatientSpecificInference/InputData/';
outputPath = [inputPath(1:end-1),'_dat/'];

bRotate = 1;

if( exist(outputPath,'dir') == 0 )
    mkdir(outputPath)
end;

if( exist(outputPath,'dir') == 0 )
    mkdir(outputPath)
end;


%4) i)  get names of files in the input folder
%   ii) convert them to nii
%   iii) save nii to output folder
files = dir(inputPath);
filesNames = {files.name}
filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));

fprintf('\n Converting files: \n------------------\n');

for i = 1:length(filesNames)
    
    inFilename  = filesNames{i};
    fprintf('%s \n', inFilename);
    
    if(~isempty(strfind(inFilename,'nii')))
            
        if( inFilename(end-3:end) == '.nii')
            outFilename = [inFilename(1:end-3),'dat'];
        else if( inFilename(end-6:end) == '.nii.gz')
                outFilename = [inFilename(1:end-6),'dat'];
            end;
        end;
        
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
end;

