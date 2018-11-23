%===================================
%
%   Correct empty space in path name
% -----------------------------------
%  INPUT:  path
%  OuTPUT: the same path with corrected spaces
%  UNIX ONLY!!!
% -----------------------------------
% Example:
% input=/home/My File Name/
% output=/home/My\ File\ Name
%
%===================================


function output=correctEmptySpaceInPathName(input)

spaces=find(isspace(input)==1);

if(length(spaces)>0)
    output=[];
    s=1;
    for i=1:length(spaces)
        output=[output,input(s:spaces(i)-1),'\'];
        s=spaces(i);
    end
    output=[output,input(s:end)];
else
    output=input;
end
