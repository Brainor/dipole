function a = sbc(a,dim,fix_point)
%SBC set boundary condition in specific dimension
%   a is the data
%   dim is the dimension to set periodic BC to
%   fix_point is optional array [a_left,a_right], if set, set fix BC to the x dimension
switch dim
    case 1
        switch nargin
            case 2%free-free
                a(1,:,:)=a(2,:,:);
                a(end,:,:)=a(end-1,:,:);
            case 3%fix-fix
                a(1,:,:)=fix_point(1);
                a(end,:,:)=fix_point(2);
        end
    case 2%periodic
        a(:,1,:)=a(:,end-1,:);
        a(:,end,:)=a(:,2,:);
    case 3%periodic
        a(:,:,1)=a(:,:,end-1);
        a(:,:,end)=a(:,:,2);
end
end

