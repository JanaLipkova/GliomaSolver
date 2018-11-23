%===========================================================
%  
%                     Convert nii to dat 
% ---------------------------------------------------------
%  INPUT:  Path to folder or file with nifty data
%  OUTPUT: Folder or file converted to dat -> used as MRAG input
%  
%   NOTE: 
%   - if input is folder, output is folder, at the same location with same +_dat/
%   - if input is file, output is file with .dat format, at the same
%   location as input
%
%  Created by Lipkova on 22/11/18.
%  Copyright (c) 2018 Lipkova
%===================================


function nii2dat(InputPath)

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

if(isdir(InputPath))    

    if(~strcmp(InputPath(end),'/'))
        InputPath=[InputPath,'/'];
    end;
    
    outputPath = [InputPath(1:end-1),'_dat/'];
    
    if( exist(outputPath,'dir') == 0 )
        mkdir(outputPath)
    end;
    
    % i) get names of files in the input folder
    % ii) convert them to nii
    % iii) save nii to output folder
    files = dir(InputPath);
    filesNames = {files.name};
    filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));
    
    % check if there are nifty data on the path
    if(isempty(strfind(filesNames,'nii')))
        fprintf('Aborting, since there are no nifty data in the input path=%s \n',InputPath);
        return
    end; 
    
    fprintf('\n Converting files: \n------------------\n');
    
    for i = 1:length(filesNames)
        
        inFilename  = filesNames{i};
        
        if(~isempty(strfind(inFilename,'nii')))
            
            fprintf('%s \n', inFilename);
            
            if( strcmp(inFilename(end-3:end),'.nii'))
                outFilename = [inFilename(1:end-3),'dat'];
            else if( strcmp(inFilename(end-6:end), '.nii.gz'))
                    outFilename = [inFilename(1:end-6),'dat'];
                end;
            end;
            
            nii_data = MRIread([correctEmptySpaceInPathName(InputPath),inFilename]);
            data = nii_data.vol;
            
            typeId = 1; % 0 double, 1 float
            writeMatrix(data,fullfile([outputPath,outFilename]),typeId);
        end;
    end;
    
    
else
    if(~isempty(strfind(InputPath,'nii')))
        
        if( strcmp(InputPath(end-3:end), '.nii'))
            outFilename = [InputPath(1:end-3),'dat'];
        else if( strcmp(InputPath(end-6:end), '.nii.gz'))
                outFilename = [InputPath(1:end-6),'dat'];
            end;
        end
        
        nii_data = MRIread(correctEmptySpaceInPathName(InputPath));
        data = nii_data.vol;
        
        typeId = 1; % 0 double, 1 float
        writeMatrix(data,outFilename,typeId)
    else
        fprintf('Your input file: %s \n does not contain nifty data \n aborting',inputPath);
    end;
    
end;
