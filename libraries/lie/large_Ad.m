%% Large Adjoint Matrix of SE(3)
function AdT=large_Ad(T)
    if(size(T) == [4 4])
        R=T(1:3,1:3);
        p=T(1:3,4);
        AdT=[[R;skew(p)*R],[zeros(3,3);R]];
    else
        error('large adjoint error: wrong input size');
    end
end