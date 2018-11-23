%===========================================================
%  
%                     Convert dat to nii 
% ---------------------------------------------------------
%  INPUT:  Path to folder or file with .dat data
%  OUTPUT: Folder or file converted to .nii 
%  
%   NOTE: 
%   - if input is folder, output is folder, at the same location with same +_dat2nii/
%   - if input is file, output is file with .nii format, at the same
%   location as input
%
%  Created by Lipkova on 22/11/18.
%  Copyright (c) 2018 Lipkova
%===================================


function dat2nii(InputPath, InputNiiPath)

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/NIfTI_20140122/')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')


%3) Read in any nifty volume to learn nifty header
 nii_data = MRIread(['correctEmptySpaceInPathName(InputNiiPath),CSF.nii.gz']);

if(isdir(InputPath))    

    if(~strcmp(InputPath(end),'/'))
        InputPath=[InputPath,'/'];
    end;
    
    outputPath = [InputPath(1:end-1),'_dat2nii/'];
    
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
        
        if(~isempty(strfind(inFilename,'dat')))
            fprintf('%s \n', inFilename);
            
            if( strcmp(inFilename(end-3:end),'.dat'))
                outFilename = [inFilename(1:end-3),'nii.gz'];
            end;
            
            datVolume = loadMatrix([inputPath,inFilename]);
            
            if(bResize)
                [Nx,Ny,Nz] = size(nii_data.vol);
                datVolume = resize_data(datVolume,Nx,Ny,Nz);
            end;
            
            nii_data.vol = datVolume;
            MRIwrite(nii_data, [outputPath,outFilename]);
        end;
    end;
    
    
else
    if(~isempty(strfind(InputPath,'dat')))
        
        if( strcmp(InputPath(end-3:end), '.dat'))
            outFilename = [InputPath(1:end-3),'nii.gz'];
        end
        
        datVolume = loadMatrix(correctEmptySpaceInPathName(InputPath));
        
        if(bResize)
            [Nx,Ny,Nz] = size(nii_data.vol);
            datVolume = resize_data(datVolume,Nx,Ny,Nz);
        end;
        
        nii_data.vol = datVolume;
        MRIwrite(nii_data, outFilename);
    else
        fprintf('Your input file: %s \n does not contain nifty data \n aborting',inputPath);
    end;
    
end;
