%% return skew-symmetric matrix : r -> [r] or [r] - > r
function mat=skew(r)
    if (size(r) == [3,1])
        mat=zeros(3);
        mat(1,2)=-r(3);
        mat(1,3)=r(2);
        mat(2,1)=r(3);
        mat(2,3)=-r(1);
        mat(3,1)=-r(2);
        mat(3,2)=r(1);
    elseif (size(r) == [3,3])
        mat=zeros(3,1);
        mat(1,1)=r(3,2);
        mat(2,1)=r(1,3);
        mat(3,1)=r(2,1);
    else
        error('skew-symmetric matrix error: wrong input size');
    end
end