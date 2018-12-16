%===========================================================
%  
%                     Convert dat to nii 
% ---------------------------------------------------------
%  INPUT:  
%  InputDatPath  = path to input .dat data (folder or file)
%  InputNiiPath  = path to nifty data to which .dat should be mapped
%                = of not need to map data, provide path to any nifty
%                volume so i can learn its header
%  Optional parameters:
%  bRotate       = rotate data back - use if nii2dat used bRotate
%  bResize       = if used, resize data to resolution of input nifty volume
%
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


function dat2nii(varargin)

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

bRotate=1;
bResize=0;

switch nargin
    case 0
        error('You have to provide at least 2 input arguments.');
    case 1
        error('You have to provide at least 2 input arguments.');
    case 2
        InputDatPath=fullfile(varargin{1});
        InputNiiPath=fullfile(varargin{2}); 
    case 3
        InputDatPath=fullfile(varargin{1});
        InputNiiPath=fullfile(varargin{2}); 
        bRotate=varargin{3}; 
    case 4
        InputDatPath=fullfile(varargin{1});
        InputNiiPath=fullfile(varargin{2}); 
        bRotate=varargin{3}; 
        bResize=varargin{4}; 
    otherwise
        error('Too many input arguments. Expects min 2 max 4.')
end


%3) Read in any nifty volume to learn nifty header
if(isempty(strfind(InputNiiPath,'nii')))
    error('Aborting, input nifty path does not contain nifty data');
end;

if( strcmp(InputNiiPath(end-6:end), '.nii.gz'))
    InputNiiPath=correctEmptySpaceInPathName(InputNiiPath);
end;

nii_data = MRIread(InputNiiPath);

if(bRotate)
    % MRIread rotate data, so rotate to original reference space
    nii_data.vol = rotate90_3D(nii_data.vol,1);
    nii_data.vol = rotate90_3D(nii_data.vol,1);
    nii_data.vol = rotate90_3D(nii_data.vol,1);
end;


if(isdir(InputDatPath))    

    outputPath = [InputDatPath(1:end-5),'_dat2nii/'];
    if( exist(outputPath,'dir') == 0 )
        mkdir(outputPath)
    end;
    
    % i) get names of files in the input folder
    % ii) convert them to nii
    % iii) save nii to output folder
    files = dir(InputDatPath);
    filesNames = {files.name};
    filesNames = filesNames(~ismember(filesNames,{'.','..','.DS_Store'}));
    
    % check if there are dat data on the path
    if(isempty(strfind(filesNames,'dat')))
        error('Aborting, there are no .dat data in the provided input path');
    end; 
    
    fprintf('\n Converting files: \n------------------\n');
    
    for i = 1:length(filesNames)
        
        inFilename  = filesNames{i};
        
        if(~isempty(strfind(inFilename,'dat')))
            fprintf('%s \n', inFilename);
            
            if( strcmp(inFilename(end-3:end),'.dat'))
                outFilename = [inFilename(1:end-3),'nii'];
            end;
            
            datVolume = loadMatrix([InputDatPath,inFilename]);
            
            if(bResize)
                [Nx,Ny,Nz] = size(nii_data.vol);
                datVolume = resize_data(datVolume,Nx,Ny,Nz);
            end;
            
            if(bRotate)
                datVolume = rotate90_3D(datVolume,1);
            end;
            
            nii_data.vol = datVolume;
            MRIwrite(nii_data,[outputPath,outFilename]);
        end;
    end;
    
    
else
    if(~isempty(strfind(InputDatPath,'dat')))
        
        if( strcmp(InputDatPath(end-3:end), '.dat'))
            outFilename = [InputDatPath(1:end-3),'nii'];
        end
        
        datVolume = loadMatrix(InputDatPath);
        
        if(bResize)
            [Nx,Ny,Nz] = size(nii_data.vol);
            datVolume = resize_data(datVolume,Nx,Ny,Nz);
        end;
        
        if(bRotate)
            datVolume = rotate90_3D(datVolume,1);
        end;
        
        nii_data.vol = datVolume;
        MRIwrite(nii_data, outFilename);
        
    else
        fprintf('Your input file: %s \n does not contain nifty data \n aborting',InputDatPath);
    end;
    
end;
