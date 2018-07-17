function [volume,header]= NIFTIimport(filename, tmpDir)
%function [volume header]= NIFTIimport(filename)
%
% filename ... The path to the file to import.
% tmpDir ... (optional) The path to the directory used for temporary files.
%            If omitted, the directory of the file to import is used.
%
% ABurner - enhanced 2011: gzip handling
%         - enhanced 2011: header
% MDiFranco - fixed 2013-10-16: temp dir for unzipping files needs to be
% unique to process to avoid contention with other processes.
% SSelendi - 2013-10-17: optional argument to set temp dir location added.
% MDiFranco - 2013-11-13: Default to using system tempdir for unzipping file, but with
%                         a unique subdirectory.

if (nargin < 2)
    destDirName = tempname;
else
    destDirName = fullfile(tmpDir,num2str(now));
end

delGZ=false;

%% to allow easy file handling, .gz must not be given as parameter
if (~exist( filename, 'file' ) && exist( [filename '.gz'], 'file' ))
    filename=[filename '.gz'];
end

%% gunzip .gz files to local temp dir
if (strcmp(filename(end-2:end),'.gz'))
    
   filenames=gunzip( filename, destDirName );
   filename=filenames{1};
   
   delGZ=true;
   
end

%% load nifti file
data = load_untouch_nii(filename);
volume = data.img;
header = data.hdr;

%% delete temp file
if delGZ
%    delete( filename );
   rmdir( destDirName , 's' );
end
