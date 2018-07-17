%===========================================
%
% add N voxels from border of a binary mask
%
%============================================
% INPUT:
% bnaryMask  = binary mask to be corrected
% Nvoxels    = number of voxels to be removed from the border
% OUTPUT:
% outData = input binary mask, where N voxels allong border were removed
% NOTE:
% imerode = works also in 3D but in 2D it is easier to control exact
% amount of voxels to be removed

function outData = addNvoxelsToBorder(binaryMask,Nvoxels,bShow)


% create structuring element, size 3x3 with zeros in corners
m = ones(3,3);
m(1,1) = 0;
m(1,3) = 0;
m(3,1) = 0;
m(3,3) = 0;

se = strel(m);

[Nx,Ny,Nz] = size(binaryMask);
outData = binaryMask;

for in = 1:Nvoxels
    
    for iz = 1:Nz
        outData(:,:,iz) = imdilate(outData(:,:,iz),se);
    end;
end;

if(bShow)
    vi(outData + binaryMask);
end;
