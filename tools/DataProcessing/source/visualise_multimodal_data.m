%=======================================
%
% Visualise multimodal medical data
%
% @author: Jana Lipkova
% Created on 01.06.2018
%
%---------------------------------------
% Visualise in one joint plot all input modalities and their name
%
% INPUT:
% vis = visualisation structure with following atributes:
%
% InputFolder = folder with the input modalities (assume all modalities are in the folder)
% OuputFolder = folder where the final visualisaiton will be stored
% OutputName  = name of the output image that will be stored in the OutputFolder
% Modalities  = list of file names of modalities that shoule be visualised
% Colormap    = color map to be used for each modality (order same as in Modalities) (see matlab colormap options)
% nCols       = number of columns of the final image
% nRows       = number of rows in the final image
% Slice       = slice of the data to be visualised
% bPositive   = set negative values to zero if 1, else nothing
% FigPosition = figure scaling and positioning on screen, - set w.r. to
% final image, i.e.  [0.1,0.1,0.4,0.95] ir [0.1,0.1,0.4,0.4]
%
% OUTPUT:
% figure with all modalities
%=======================================

function visualise_multimodal_data(vis)

% addpath('../lib/toolbox_matlab_nifti')
% addpath('../lib/')


% Create figure, set position/size to almost full screen.
figure
set( gcf, 'Units', 'normalized', 'Position', vis.FigPosition ) ;

% Create grid of axes to fit subplots nicely
LB = 0.;
UB = 1;
[blx, bly] = meshgrid( LB: UB/vis.nCols: UB, LB: UB/vis.nRows:UB ) ;
hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, UB*UB/vis.nCols, UB*UB/vis.nRows] ), blx, bly, 'UniformOutput', false ) ;


row=1;
col=1;

for i = 1:length(vis.Modalities)
    
    % get data
    
    inputPath = fullfile(vis.InputFolder,vis.Modalities{i});
    if( strcmp(inputPath(end-6:end), '.nii.gz'))
        inputPath=correctEmptySpaceInPathName(inputPath);
    end;
    
    data      = MRIread(inputPath);
    
    if(vis.bPositive)
        data.vol(data.vol(:)<0) = 0;
    end;
    
    baseName  = strsplit(vis.Modalities{i},'.');
    modName   = baseName{1};
    
    % plot on given position
    axes( hAxes{vis.nRows - row + 1,col} )
    pcolor(data.vol(:,:,vis.Slice));
    daspect([1 1 1])
    shading flat;
    colormap(gca,vis.Colormap{i})
    
    % title
    v = axis;
    posX = v(1) + 0.1* (v(2) - v(1));
    posY = v(3) + 0.95* (v(4) - v(3));
    text( posX, posY, modName, 'FontSize', 15, 'FontWeight','bold','Color', [0.99,0.99,0.99],'Interpreter','none' );
    
    % turn off axis numbering
    set( gca, 'Visible', 'off' ) ;
    
    col=col+1;
    
    if(mod(i,vis.nCols)==0)
        col=1;
        row=row+1;
    end;    
end;

% Save output
if( exist(vis.OutputFolder,'dir') == 0 )
    sprintf('Output folder does not exist, creating it in: \n %s', vis.OutputFolder)
    mkdir(vis.OutputFolder)
end;

print(gcf,fullfile(vis.OutputFolder,vis.OutputName),'-djpeg')






