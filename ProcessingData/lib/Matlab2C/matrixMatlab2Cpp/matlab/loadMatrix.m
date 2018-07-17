function A = loadMatrix(filename)
%LOADMATRIX Load matrix in format that my Matrix in C++ would dump.
% See Matrix.h for format of file and writeMatrix for element-order.
%
% USAGE: A = loadMatrix(filename)
%
% PARAMETERS:
% Input:
%   filename:
%       Name of file to read
% Output:
%   A:
%       Matrix that was loaded

% get header from file
fid = fopen(filename);
magic = fread(fid, 1, 'integer*4');
assert(magic == 1234, 'Magic number does not match');
DIM = fread(fid, 1, 'integer*4');
N = fread(fid, DIM, 'integer*4');
typeId = fread(fid, 1, 'integer*4');

% choose type for fread
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

% read data and reshape
[A, count] = fread(fid, prod(N), fType);
assert(count == prod(N), 'Read %d instead of %d elements', count, prod(N));
A = reshape(A, N');
% transpose for 2D so we have A(y,x)
if DIM == 2
    A = A';
end
fclose(fid);

end
