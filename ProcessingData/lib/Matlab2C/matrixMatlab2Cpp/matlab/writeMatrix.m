function writeMatrix(A, filename, typeId)
%WRITEMATRIX Dump matrix in format that my Matrix in C++ can read.
% See Matrix.h for format of file.
% Note on matrix-storage:
% (first dim. changing fastest (Matlab/Fortran style => reverse for C++))
% - 2D: X,Y (input given as A(Y,X) since this is def. in Matlab)
% - 2D with lda: LDA,X,Y (input given in same order)
% - 3D: X,Y,Z (input given in same order)
% - 3D with lda: LDA,X,Y,Z (input given in same order)
%
% USAGE: writeMatrix(A, filename, typeId)
%
% PARAMETERS:
% Input:
%   A:
%       Matrix to be dumped.
%   filename:
%       Name of resulting file.
%   typeId: (optional)
%       Type to use for dump (0 = double, 1 = single, 2 = int).
%       By default it will choose int for integer-matrices, double for
%       double-matrices and single otherwise.

% Extract info from A
if nargin < 3
    if isinteger(A)
        typeId = 2;
    elseif isa(A, 'double')
        typeId = 0;
    else
        typeId = 1;
    end
end
N = size(A);
DIM = length(N);
magic = 1234;

% choose fwrite type
switch typeId
    case 0
        fType = 'real*8';
    case 1
        fType = 'real*4';
    case 2
        fType = 'integer*4';
    otherwise
        error('Unknown type %d', typeId);
end

% write
fid = fopen(filename, 'wb');
fwrite(fid, magic, 'integer*4');
fwrite(fid, DIM, 'integer*4');
% write matrix transposed for 2D
if DIM == 2
    N = [N(2), N(1)];
    A = A';
end
fwrite(fid, N, 'integer*4');
fwrite(fid, typeId, 'integer*4');
fwrite(fid, A, fType);
fclose(fid);

end
