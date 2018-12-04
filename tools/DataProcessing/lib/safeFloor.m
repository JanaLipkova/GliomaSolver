%===================================
%
%   Safely floor the numbers
%   - preven rounding issues of matlab
%   Example:
%   floor(0.35*2/0.1) = 6 
%   safeFloor(0.35*2/0.1) = 7
%=================================== 


function out = safeFloor(in)
eps = 1e-07;
if( (ceil(in)) - in < eps);
    out=ceil(in);
else
    out=floor(in);
end