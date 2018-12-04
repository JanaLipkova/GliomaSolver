%=======================================
%
% Visualise multimodal medical data
%
% @author: Jana Lipkova
% Created on 01.06.2018
%
%---------------------------------------
%
% InputFolder = folder with the input modalities to be visualised
% OutputPath  = where the final visualisaiton will be stored
%=======================================

function plotAllNiftyData(varargin)

addpath('../lib/toolbox_matlab_nifti')
addpath('../lib/vi');
addpath('../lib/')

switch nargin
    case 0
        error('You have to provide at least one input argument.');
    case 1
        vis.InputFolder  = varargin{1};
        vis.OutputFolder = varargin{1};
        vis.OutputName   = 'Overview';
    case 2
        vis.InputFolder = varargin{1};
        vis.OutputFolder  = varargin{2};
        vis.OutputName  = 'Overview';
    otherwise
        error('Too many input arguments. Expect min 1, max 2 arguments.')
end

% 1) Create a list of inpur nifty files and their color maps (hot for
% PET, jet for MAP, gray else)
files = dir(vis.InputFolder);
filesNames = {files.name};

j=1;
MAPindex=1;
[~,idx] = sort([files.datenum]);

for i = idx
    if(~isempty(strfind(filesNames{i},'nii')))
        vis.Modalities(j)= {filesNames{i}} ;
        vis.Colormap(j) = {'gray'};        
        j=j+1;
    end;
end;

% reset color map for selected cases
for i = 1:length(vis.Modalities)
    if(~isempty(strfind(vis.Modalities{i},'MAP')))
        vis.Colormap(i) = {'jet'};
        MAPindex=i;
    else
        if(~isempty(strfind(vis.Modalities{i},'FET')))
            vis.Colormap(i) = {'hot'};
        else
            if(~isempty(strfind(vis.Modalities{i},'GT')))
                vis.Colormap(i) = {'jet'};
                MAPindex=i;
            end;
        end;
    end;
end;


% 2) Compute number of subplots (+1 to show MAP + GM)
Nplots = length(vis.Modalities);
vis.nCols = ceil(sqrt(Nplots));
vis.nRows = ceil( Nplots/vis.nCols);

% 3) Compute slice that should be visualise based on MAP.nii (center of z-direction)
inMAP=fullfile(vis.InputFolder,filesNames{MAPindex});

if( strcmp(inMAP(end-6:end), '.nii.gz'))
    inMAP=correctEmptySpaceInPathName(inMAP);
end;

data = MRIread(inMAP);
[i1, i2, i3] = ind2sub(size(data.vol), find(data.vol));  % indecies of non-zero elements in each dimension
vis.Slice = i3(1) + floor((i3(end) - i3(1)) / 2);

% 4) Other visualisation set up
% vis.FigPosition = [0.1,0.1,0.4,0.89];
% UBx = min(1, 0.1*vis.nCols);
% UBy = min(1, 0.25*vis.nRows);

UBx =  min(0.9, 0.1*vis.nCols);
UBy =  min(0.9, 0.22*vis.nRows);

vis.FigPosition = [0.1,0.1,UBx,UBy];
vis.bPositive  = 1; % keep only positive values

visualise_multimodal_data(vis)


