% ===================================
% 
%   Resize data from original [Ix,Iy,Iz] resolution
% to prescribed output resolution [Ox,Oy,Oz]
% =================================== 
% 
% 
function output = resize_data(input,Ox,Oy,Oz)

[Ix,Iy,Iz] = size(input);
output     = zeros(Ox,Oy,Oz);

% 1-1 mapping, just cut the padded zeros
if((Ix == 256) && (Iy == 256) && (Iz == 256))
    output = input(1:Ox, 1:Oy, 1:Oz);
else
    % 2-1 mapping, + cut the padded zeros
    if((Ix == 128) && (Iy == 128) && (Iz == 128))
        out256 = zeros(256,256,256);
        
        ho = 1/max([256-1,256-1,256-1]);
        hi = 1/max([Ix-1,Iy-1,Iz-1]);
        
        for iz = 0 : Oz-1
            for iy = 0 : Oy-1
                for ix = 0 : Ox-1
                    %compute mapped coordinates
                    mx = safeFloor(ho * ix / hi);
                    my = safeFloor(ho * iy / hi);
                    mz = safeFloor(ho * iz / hi);
                    
                    out256(ix+1,iy+1,iz+1) = input(mx+1,my+1,mz+1);
                end;
            end;
        end;

        output = out256(1:Ox, 1:Oy, 1:Oz);
        
    else
        % proper mapping based on physical spacing
        ho = 1/max([Ox-1,Oy-1,Oz-1]);
        hi = 1/max([Ix-1,Iy-1,Iz-1]);
        
        for iz = 0 : Oz-1
            for iy = 0 : Oy-1
                for ix = 0 : Ox-1
                    %compute mapped coordinates
                    mx = safeFloor(ho * ix / hi);
                    my = safeFloor(ho * iy / hi);
                    mz = safeFloor(ho * iz / hi);
                    
                    output(ix+1,iy+1,iz+1) = input(mx+1,my+1,mz+1);
                end;
            end;
        end;
        
    end;
end;



        
    
    



