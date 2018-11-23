function vr = rotate90_3D(v,axis)
%function vr = rotate90_3D(v,axis)
%

    function v=switchAxis(v, axis)
        switch axis
            case 1
                v=permute(v, [2 1 3]);
            case 2
                v=permute(v, [3 2 1]);
            case 3
                v=permute(v, [1 3 2]);
        end
    end

v=switchAxis(v, axis);
vc=num2cell(v, [1,2]);
vc=cellfun(@(x) rot90(x), vc, 'uni', false);
vr=cat(3, vc{:});
vr=switchAxis(vr, axis);

end

%{
test:
    textAxis=1;
    
    v=zeros(3,3,3);
    vc=v;
    v(1,1,1)=1;

    for i=1:4
        v=rotate90_3D(v,textAxis);
        vc=vc+v;
    end

    imgHelpers.vol2movie(vc);
%}
