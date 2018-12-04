%===========================================================
%
%                     Convert nii to dat
% ---------------------------------------------------------
%  INPUT:
%         - Path to nifty data (folder or file)
%  Optional parameters:
%         - bRotate = rotate data to avoid issue with MRIread rotation
%         - bResize = resize data to 256 resolution
%  OUTPUT: Folder or file converted to dat -> used as MRAG input
%
% IMPORTANT: Output is padd with zeros to 256^3 resolution to ensure
% compatibility with MRAG grid and avoid interpolation erros
%
%   NOTE:
%   - if input is folder, output is folder, at the same location with same +_dat/
%   - if input is file, output is file with .dat format, at the same
%   location as input
%
%  Created by Lipkova on 22/11/18.
%  Copyright (c) 2018 Lipkova
%===================================


function nii2dat(varargin)

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/Matlab2C/matrixMatlab2Cpp/matlab/')
addpath('../lib/')

% Optional parameters
bRotate=0;
bResize=0;

switch nargin
    case 0
        error('You have to provide at least one input argument.');
    case 1
        InputPath=fullfile(varargin{1});
    case 2
        InputPath=fullfile(varargin{1});
        bRotate=varargin{2};
    case 3
        InputPath=fullfile(varargin{1});
        bRotate=varargin{2};
        bResize=varargin{3};
    otherwise
        error('Too many input arguments. Expect min 1, max 3 arguments.')
end


if(isdir(InputPath))

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
        error('There are no nifty data in the provided input path');
    end;
    
    fprintf('\n Converting files: \n------------------\n');
    
    for i = 1:length(filesNames)
        
        inFilename  = filesNames{i};
        
        if(~isempty(strfind(inFilename,'nii')))
            
            fprintf('%s \n', inFilename);
            InputNiiFile = [InputPath,inFilename];
            
            if( strcmp(inFilename(end-3:end),'.nii'))
                outFilename  = [inFilename(1:end-3),'dat'];
            else if( strcmp(inFilename(end-6:end), '.nii.gz'))
                    outFilename = [inFilename(1:end-6),'dat'];
                    InputNiiFile = correctEmptySpaceInPathName(InputNiiFile);
                end;
            end;
            
            nii_data = MRIread(InputNiiFile);
            data = nii_data.vol;
            
            if(bRotate)
                % MRIread rotate data, so rotate to original reference space
                data = rotate90_3D(data,1);
                data = rotate90_3D(data,1);
                data = rotate90_3D(data,1);
            end;
            
            if(bResize)
                [Nx,Ny,Nz] = size(data);
                if( max([Nx,Ny,Nz]) > 256)
                    error('Data resolution is too big to fit in 256^3 grid');
                else
                    fprintf('Resizing the input data to 256^3 resolution by padding with zeros.\n');
                    pX = 256 - Nx;
                    pY = 256 - Ny;
                    pZ = 256 - Nz;
                    data = padarray(data,[pX,pY,pZ],0,'post');
                end;
            end;
            
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
                InputPath = correctEmptySpaceInPathName(InputPath);
            end;
        end
        
        nii_data = MRIread(InputPath);
        data = nii_data.vol;
        
        if(bRotate)
            % MRIread rotate data, so rotate to original reference space
            data = rotate90_3D(data,1);
            data = rotate90_3D(data,1);
            data = rotate90_3D(data,1);
        end;
        
        if(bResize)
            [Nx,Ny,Nz] = size(data);
            if(max(Nx,Ny,Nz) > 256)
                error('Data resolution is too big to fit in 256^3 grid');
            else
                fprintf('Resizing the input data to 256^3 resolution by padding with zeros.');
                pX = 256 - Nx;
                pY = 256 - Ny;
                pZ = 256 - Nz;
                data = padarray(data,[pX,pY,pZ],0,'post');
            end;
        end;
        
        typeId = 1; % 0 double, 1 float
        writeMatrix(data,outFilename,typeId)
    else
        fprintf('Your input file: %s \n does not contain nifty data. Aborting \n',InputPath);
    end;
    
end;

fprintf('------------------\n');

