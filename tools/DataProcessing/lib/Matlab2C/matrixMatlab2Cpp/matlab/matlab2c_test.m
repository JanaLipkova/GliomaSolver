%=========================
%
%   Test Matlab2C++
%   use Gerardo's converter
%==========================

% AIM:
% 1) define matlab matrix
% 2) store matlab matrix for input suitable for C++

function matlab2c_test

% define matrix
N = 10;
A = eye(10,10);

writeMatrix(A,'test_matrix.dat',0)

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