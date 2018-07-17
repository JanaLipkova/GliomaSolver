function NIFTIexport(volume, filename, opt_compress, pixdim, units)
%% function NIFTIexport(volume, filename, opt_compress)
%
% Export data to NIFTI format
% volume ... the data that should be exported
% filename ... the filename for the exported data (optional)
% opt_compress ... (optional) argument whether the nifti should be gzipped
% pixdim ... (optional) using the original pixel dimension
% units ... (optional) "m", "mm", or "um" (all other inputs ignored) 
%
% also necessary to download the NIFTI Package from:
% http://www.mathworks.com/matlabcentral/fileexchange/8797
%
% added on 2008-11-20 by eva.dittrich@student.tuwien.ac.at
% enhanced 2011       by andreas.burner@meduniwien.ac.at: gzip functionality
% enhanced 2013       by alexander.valentinitsch@meduniwien.ac.at: pixel dimension - for ITK Snap
% enhanced 2014-01-27 by matthew.difranco@meduniwein.ac.at: spacial units input

%% if no filename is given, generate one
if nargin==1
    filename=['volume_' datestr(now,'yyyymmddHHMMSS')];
end  

%% check if compression flag is set
if nargin<3
    opt_compress=0;
end
    
%% add .nii extension if not given
if ~strcmp(filename(end-3:end),'.nii')
    filename=[filename '.nii'];
end

% cast to uint8 if volume is logical
if islogical(volume)
    volume = uint8(volume);
end

%% generate and save nifti file
if nargin > 3
    % Save volume with voxel size (in mm)
    if strcmp(units,'um')
        pixdim = pixdim ./ 1000;
    elseif strcmp(units,'m')
        pixdim = pixdim .* 1000;
    end
    units = 'mm';
    nii = make_nii(volume,pixdim);
else
    nii = make_nii(volume);
end

% Set pixel dimenstions from input, if given
unitCode = zeros(1,'uint8');
if nargin > 4
    if any(strcmp(units,{'m','um'}))
        unitCode = bitset(unitCode,1,1);
    end
    if any(strcmp(units,{'mm','um'}))
        unitCode = bitset(unitCode,2,1);
    end
end
nii.hdr.dime.xyzt_units = unitCode;

% % Reorient volume if necessary
% view_nii(nii);
% [nii,orient] = rri_orient(nii);
% close(gcf);
% view_nii(nii);

save_nii(nii, filename);

%% compress output filename
if opt_compress
    os=computer;
    if ( strcmp( os(1:5), 'PCWIN' ) )
        system(['gzip -f9 "' filename '"']);
    else
        system(['gzip -f9 ''' filename '''']);
    end
end

