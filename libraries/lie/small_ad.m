%% Small Adjoint Matrix of se(3)
function adV=small_ad(V)
    if(size(V) == [6,1])
        w=V(1:3,1);
        v=V(4:6,1);
        adV=[[skew(w);skew(v)],[zeros(3,3);skew(w)]];
    else
        error('small adjoint error: wrong input size');
    end
end